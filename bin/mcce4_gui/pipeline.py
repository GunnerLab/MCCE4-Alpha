"""
Pipeline manager for MCCE4 steps.

Runs MCCE4 computation steps as async subprocesses with real-time
log streaming, timing, and stop/resume support.
"""

from __future__ import annotations

import asyncio
import os
import time
from dataclasses import dataclass
from typing import Callable, Awaitable, Optional

from .config import MCCEConfig, StepStatus, PipelineState


@dataclass
class StepResult:
    """Tracks the status and timing of a single pipeline step."""
    step_number: int
    name: str
    status: StepStatus = StepStatus.PENDING
    start_time: float = 0.0
    end_time: float = 0.0
    return_code: int = -1
    error_message: str = ""

    @property
    def elapsed(self) -> float:
        if self.start_time == 0:
            return 0.0
        end = self.end_time if self.end_time > 0 else time.time()
        return end - self.start_time

    @property
    def elapsed_str(self) -> str:
        s = self.elapsed
        if s < 60:
            return f"{s:.0f}s"
        elif s < 3600:
            return f"{int(s // 60)}m {int(s % 60)}s"
        else:
            return f"{int(s // 3600)}h {int((s % 3600) // 60)}m"

    def to_dict(self) -> dict:
        return {
            "step_number": self.step_number,
            "name": self.name,
            "status": self.status.value,
            "elapsed": round(self.elapsed, 1),
            "elapsed_str": self.elapsed_str,
            "return_code": self.return_code,
            "error_message": self.error_message,
        }


LogCallback = Callable[[str], Awaitable[None]]


class PipelineManager:
    """
    Manages the MCCE4 4-step pipeline with async subprocess execution.

    State machine:
        IDLE ──▶ RUNNING ──▶ COMPLETED
                   │  ▲
                   ▼  │
                 PAUSED
                   │
                   ▼
                 FAILED
    """

    STEP_NAMES = {1: "Preparation", 2: "Rotamers", 3: "Energies", 4: "Monte Carlo"}

    def __init__(self, config: MCCEConfig, log_callback: Optional[LogCallback] = None):
        self.config = config
        self.log_callback = log_callback
        self.state = PipelineState.IDLE
        self.steps: list[StepResult] = []
        self.current_step_index: int = 0
        self._process: Optional[asyncio.subprocess.Process] = None
        self._stop_event = asyncio.Event()
        self._command_tuples: list[tuple[int, list[str], str]] = []

    async def _log(self, message: str) -> None:
        if self.log_callback:
            await self.log_callback(message)

    def build_commands(self) -> list[tuple[int, list[str], str]]:
        """Build (step_number, command, log_file) tuples from config."""
        cfg = self.config
        mcce_home = cfg.mcce_home
        diel = str(cfg.inner_dielectric)
        commands = []

        param_env = f"-u MCCE_HOME={mcce_home}"
        if cfg.extra_tpl:
            param_env += f",EXTRA={cfg.extra_tpl}"

        if cfg.step1.enabled:
            cmd = ["python", f"{mcce_home}/bin/step1.py", "prot.pdb", "-d", diel]
            if cfg.step1.dry:
                cmd.append("--dry")
            cmd.append(param_env)
            commands.append((1, cmd, "step1.log"))

        if cfg.step2.enabled:
            commands.append((
                2,
                ["python", f"{mcce_home}/bin/step2.py",
                 "-l", str(cfg.step2.conformer_level), "-d", diel, param_env],
                "step2.log",
            ))

        if cfg.step3.enabled:
            cmd = [
                "python", f"{mcce_home}/bin/step3.py",
                "-d", diel,
                "-do", str(cfg.outer_dielectric),
                "-p", str(cfg.step3.processes),
                "-salt", str(cfg.step3.salt),
            ]
            solver = cfg.step3.pbe_solver
            if solver != "ngpb":
                cmd.extend(["--solver", solver])
            for flag_name in ["skip_pb", "old_vdw", "fly", "refresh", "debug"]:
                if getattr(cfg.step3, flag_name, False):
                    cmd.append(f"--{flag_name.replace('_', '-')}")
            cmd.append(param_env)
            commands.append((3, cmd, "step3.log"))

        if cfg.step4.enabled:
            cmd = [
                "python", f"{mcce_home}/bin/step4.py",
                "-i", str(cfg.step4.init_ph),
                "-d", str(cfg.step4.interval),
                "-n", str(cfg.step4.titration_steps),
                "-t", cfg.step4.titration_type,
            ]
            if cfg.step4.xts:
                cmd.append("--xts")
            if cfg.step4.ms:
                cmd.append("--ms")
            cmd.append(param_env)
            commands.append((4, cmd, "step4.log"))

        return commands

    async def _prepare_pdb(self) -> bool:
        """Create prot.pdb symlink and optionally center the protein."""
        cfg = self.config
        pdb_path = os.path.abspath(cfg.pdb_file)

        if not os.path.isfile(pdb_path):
            await self._log(f"ERROR: PDB file not found: {pdb_path}")
            return False

        # Create prot.pdb symlink in the working directory
        prot_link = os.path.join(os.getcwd(), "prot.pdb")
        try:
            if os.path.exists(prot_link) or os.path.islink(prot_link):
                os.remove(prot_link)
            os.symlink(pdb_path, prot_link)
            await self._log(f">> Linked prot.pdb -> {pdb_path}")
        except OSError as e:
            await self._log(f"ERROR: Failed to create prot.pdb symlink: {e}")
            return False

        # Optionally center the protein using orientation.py
        if cfg.step1.center_protein:
            orientation_py = os.path.join(cfg.mcce_home, "bin", "orientation.py")
            if not os.path.isfile(orientation_py):
                await self._log(f"WARNING: orientation.py not found at {orientation_py}")
                return True

            await self._log(">> Centering protein with orientation.py...")
            try:
                proc = await asyncio.create_subprocess_exec(
                    "python", orientation_py, "prot.pdb",
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.STDOUT,
                )
                stdout, _ = await proc.communicate()
                if stdout:
                    for line in stdout.decode("utf-8", errors="replace").strip().split("\n"):
                        await self._log(line)

                if proc.returncode == 0:
                    # orientation.py produces prot_center.pdb — replace prot.pdb with it
                    centered = os.path.join(os.getcwd(), "prot_center.pdb")
                    if os.path.isfile(centered):
                        os.remove(prot_link)
                        os.rename(centered, prot_link)
                        await self._log(">> prot.pdb replaced with centered structure")
                    else:
                        await self._log("WARNING: prot_center.pdb not produced")
                else:
                    await self._log(f"WARNING: orientation.py exited with code {proc.returncode}")
            except Exception as e:
                await self._log(f"WARNING: orientation.py failed: {e}")

        return True

    async def run(self, from_step: int = 0) -> None:
        """Run the pipeline from the given step index."""
        if self.state == PipelineState.RUNNING:
            await self._log("Pipeline is already running.")
            return

        # Validate
        errors = self.config.validate()
        if errors:
            for err in errors:
                await self._log(f"ERROR: {err}")
            self.state = PipelineState.FAILED
            return

        # Prepare prot.pdb before building commands
        if from_step == 0:
            ok = await self._prepare_pdb()
            if not ok:
                self.state = PipelineState.FAILED
                return

        # Build commands
        self._command_tuples = self.build_commands()
        if not self._command_tuples:
            await self._log("No steps enabled.")
            return

        # Initialize step results on fresh run
        if from_step == 0:
            self.steps = [
                StepResult(step_number=num, name=self.STEP_NAMES[num])
                for num, _, _ in self._command_tuples
            ]

        self.state = PipelineState.RUNNING
        self._stop_event.clear()
        self.current_step_index = from_step

        total = len(self._command_tuples)
        await self._log(f"Starting pipeline ({from_step + 1}/{total})...")

        for i in range(from_step, total):
            if self._stop_event.is_set():
                self.state = PipelineState.PAUSED
                await self._log(f"Pipeline paused at step {i + 1}.")
                return

            step_num, cmd, log_file = self._command_tuples[i]
            self.current_step_index = i
            step = self.steps[i]
            step.status = StepStatus.RUNNING
            step.start_time = time.time()

            await self._log(f"\n{'─' * 50}")
            await self._log(f"Step {step_num}: {self.STEP_NAMES[step_num]}")
            await self._log(f"$ {' '.join(cmd)}")
            await self._log(f"{'─' * 50}")

            success = await self._run_step(cmd, log_file, step)

            step.end_time = time.time()

            if not success:
                if self._stop_event.is_set():
                    step.status = StepStatus.PENDING
                    self.state = PipelineState.PAUSED
                    await self._log(f"Step {step_num} interrupted. Resume to restart.")
                else:
                    step.status = StepStatus.FAILED
                    self.state = PipelineState.FAILED
                    await self._log(f"Step {step_num} FAILED (exit {step.return_code})")
                return

            step.status = StepStatus.COMPLETED
            await self._log(f"Step {step_num} completed in {step.elapsed_str}")

        self.state = PipelineState.COMPLETED
        total_s = sum(s.elapsed for s in self.steps if s.status == StepStatus.COMPLETED)
        await self._log(f"\nAll steps completed. Total: {total_s:.0f}s")

    async def _run_step(self, cmd: list[str], log_file: str, step: StepResult) -> bool:
        """Execute a single step. Returns True on success."""
        try:
            self._process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT,
            )

            with open(log_file, "w") as fh:
                assert self._process.stdout is not None
                async for raw in self._process.stdout:
                    if self._stop_event.is_set():
                        self._process.terminate()
                        try:
                            await asyncio.wait_for(self._process.wait(), timeout=5)
                        except asyncio.TimeoutError:
                            self._process.kill()
                        return False

                    line = raw.decode("utf-8", errors="replace").rstrip()
                    fh.write(line + "\n")
                    await self._log(line)

            await self._process.wait()
            step.return_code = self._process.returncode or 0
            return step.return_code == 0

        except FileNotFoundError:
            step.error_message = f"Command not found: {cmd[0]}"
            await self._log(f"ERROR: {step.error_message}")
            return False
        except Exception as e:
            step.error_message = str(e)
            await self._log(f"ERROR: {e}")
            return False
        finally:
            self._process = None

    async def stop(self) -> None:
        if self.state != PipelineState.RUNNING:
            return
        await self._log("Stopping...")
        self._stop_event.set()
        if self._process and self._process.returncode is None:
            self._process.terminate()
            try:
                await asyncio.wait_for(self._process.wait(), timeout=5)
            except asyncio.TimeoutError:
                self._process.kill()

    async def resume(self) -> None:
        if self.state == PipelineState.RUNNING:
            await self._log("Already running.")
            return
        await self.run(from_step=self.current_step_index)

    def get_status(self) -> dict:
        return {
            "state": self.state.value,
            "current_step_index": self.current_step_index,
            "steps": [s.to_dict() for s in self.steps],
        }
