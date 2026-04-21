"""
Configuration models for MCCE4 GUI.

Uses dataclasses for type-safe parameter handling with JSON serialization.
Minimal dependencies — no Pydantic required.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, asdict
from enum import Enum
from pathlib import Path
from typing import Optional


# ── Enums ─────────────────────────────────────────────────────────

class StepStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


class PipelineState(str, Enum):
    IDLE = "idle"
    RUNNING = "running"
    PAUSED = "paused"
    COMPLETED = "completed"
    FAILED = "failed"


# ── Step Configs ──────────────────────────────────────────────────

@dataclass
class Step1Config:
    """Step 1: Pre-run — PDB to MCCE PDB."""
    enabled: bool = True
    dry: bool = True
    center_protein: bool = True


@dataclass
class Step2Config:
    """Step 2: Rotamer conformer generation."""
    enabled: bool = True
    conformer_level: int = 1  # 1=isosteric, 2=heavy atom, 3=comprehensive


@dataclass
class Step3Config:
    """Step 3: Energy calculations (PBE solver)."""
    enabled: bool = True
    pbe_solver: str = "ngpb"
    processes: int = 1
    salt: float = 0.15
    tmp_folder: str = "/tmp"
    skip_pb: bool = False
    old_vdw: bool = False
    fly: bool = False
    refresh: bool = False
    debug: bool = False


@dataclass
class Step4Config:
    """Step 4: Monte Carlo sampling."""
    enabled: bool = True
    init_ph: float = 0.0
    interval: float = 1.0
    titration_steps: int = 15
    titration_type: str = "ph"
    xts: bool = True
    ms: bool = False


@dataclass
class SlurmConfig:
    """SLURM job submission parameters."""
    job_name: str = "mcce_job"
    cpus: int = 4
    memory_gb: int = 24
    partition: str = ""
    time_limit: str = "24:00:00"
    extra_directives: str = ""


# ── Main Config ───────────────────────────────────────────────────

@dataclass
class MCCEConfig:
    """Complete MCCE4 run configuration."""
    # Paths
    pdb_file: str = ""
    mcce_home: str = ""
    extra_tpl: str = ""
    user_param: str = "./user_param"

    # Global parameters
    inner_dielectric: float = 4.0
    outer_dielectric: float = 80.0

    # Step configs
    step1: Step1Config = field(default_factory=Step1Config)
    step2: Step2Config = field(default_factory=Step2Config)
    step3: Step3Config = field(default_factory=Step3Config)
    step4: Step4Config = field(default_factory=Step4Config)

    # Execution
    execution_mode: str = "local"
    slurm: SlurmConfig = field(default_factory=SlurmConfig)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "MCCEConfig":
        """Reconstruct config from a dictionary."""
        d = data.copy()
        # Reconstruct nested dataclasses
        for key, factory in [
            ("step1", Step1Config),
            ("step2", Step2Config),
            ("step3", Step3Config),
            ("step4", Step4Config),
            ("slurm", SlurmConfig),
        ]:
            val = d.get(key, {})
            if isinstance(val, dict):
                d[key] = factory(**val)
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})

    def save(self, path: str | Path) -> None:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, path: str | Path) -> "MCCEConfig":
        with open(path, "r") as f:
            return cls.from_dict(json.load(f))

    def auto_detect_mcce_home(self) -> str:
        """Try to find MCCE_HOME from environment or common locations."""
        # Check environment variable
        env = os.environ.get("MCCE_HOME", "")
        if env and os.path.isdir(env):
            return env

        # Check common locations
        candidates = [
            os.path.expanduser("~/MCCE4"),
            os.path.expanduser("~/MCCE4-Alpha"),
            "/opt/mcce4",
            "/usr/local/mcce4",
        ]
        for c in candidates:
            if os.path.isdir(c) and os.path.isfile(os.path.join(c, "bin", "step1.py")):
                return c
        return ""

    def validate(self) -> list[str]:
        """Return a list of validation errors (empty = valid)."""
        errors = []
        if not self.pdb_file:
            errors.append("No PDB file specified")
        elif not os.path.isfile(self.pdb_file):
            errors.append(f"PDB file not found: {self.pdb_file}")

        if not self.mcce_home:
            errors.append("MCCE_HOME not set")
        elif not os.path.isdir(self.mcce_home):
            errors.append(f"MCCE_HOME directory not found: {self.mcce_home}")
        elif not os.path.isfile(os.path.join(self.mcce_home, "bin", "step1.py")):
            errors.append(f"step1.py not found in {self.mcce_home}/bin/")

        if not any([
            self.step1.enabled, self.step2.enabled,
            self.step3.enabled, self.step4.enabled,
        ]):
            errors.append("No steps enabled")

        return errors


DEFAULT_CONFIG_PATH = Path("mcce4_gui_config.json")
