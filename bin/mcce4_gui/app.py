"""
MCCE4 GUI — FastAPI backend.

Endpoints:
  - Config:    GET/POST /api/config, save, load
  - Pipeline:  run, stop, resume, status
  - Files:     browse, upload PDB, serve PDB content
  - SLURM:     submit, status, cancel, preview script
  - Analysis:  pKa, sum_crg, conformers, available outputs
  - WebSocket: real-time log streaming
"""

from __future__ import annotations

import asyncio
import json
import os
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, WebSocket, WebSocketDisconnect, UploadFile, File, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

from .config import MCCEConfig, DEFAULT_CONFIG_PATH
from .pipeline import PipelineManager, PipelineState
from . import slurm as slurm_mod
from . import analysis as analysis_mod

# ── App ───────────────────────────────────────────────────────────

STATIC_DIR = Path(__file__).parent / "static"

app = FastAPI(title="MCCE4 GUI", version="1.0.0")


# ── State ─────────────────────────────────────────────────────────

class AppState:
    def __init__(self):
        self.config = MCCEConfig()
        self.pipeline: Optional[PipelineManager] = None
        self.ws_clients: set[WebSocket] = set()
        self.log_history: list[str] = []
        self.pipeline_task: Optional[asyncio.Task] = None
        self.slurm_job_id: str = ""

    async def broadcast_log(self, message: str) -> None:
        self.log_history.append(message)
        if len(self.log_history) > 5000:
            self.log_history = self.log_history[-5000:]
        dead = set()
        for ws in self.ws_clients:
            try:
                await ws.send_json({"type": "log", "message": message})
            except Exception:
                dead.add(ws)
        self.ws_clients -= dead

    async def broadcast_status(self) -> None:
        if not self.pipeline:
            return
        data = self.pipeline.get_status()
        dead = set()
        for ws in self.ws_clients:
            try:
                await ws.send_json({"type": "status", "data": data})
            except Exception:
                dead.add(ws)
        self.ws_clients -= dead


state = AppState()


# ── Startup ───────────────────────────────────────────────────────

@app.on_event("startup")
async def startup():
    if DEFAULT_CONFIG_PATH.exists():
        try:
            state.config = MCCEConfig.load(DEFAULT_CONFIG_PATH)
        except Exception:
            pass
    if not state.config.mcce_home:
        detected = state.config.auto_detect_mcce_home()
        if detected:
            state.config.mcce_home = detected


# ── Frontend ──────────────────────────────────────────────────────

@app.get("/", response_class=HTMLResponse)
async def index():
    p = STATIC_DIR / "index.html"
    if p.exists():
        return p.read_text()
    return HTMLResponse("<h1>mcce4_gui</h1><p>static/index.html not found</p>", status_code=500)


if STATIC_DIR.exists():
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


# ══════════════════════════════════════════════════════════════════
# CONFIG API
# ══════════════════════════════════════════════════════════════════

@app.get("/api/config")
async def get_config():
    return state.config.to_dict()


@app.post("/api/config")
async def update_config(data: dict):
    try:
        state.config = MCCEConfig.from_dict(data)
        return {"ok": True}
    except Exception as e:
        raise HTTPException(400, str(e))


@app.post("/api/config/save")
async def save_config():
    try:
        state.config.save(DEFAULT_CONFIG_PATH)
        return {"ok": True, "path": str(DEFAULT_CONFIG_PATH)}
    except Exception as e:
        raise HTTPException(500, str(e))


@app.post("/api/config/load")
async def load_config(data: dict):
    path = data.get("path", str(DEFAULT_CONFIG_PATH))
    try:
        state.config = MCCEConfig.load(path)
        return state.config.to_dict()
    except Exception as e:
        raise HTTPException(400, str(e))


# ══════════════════════════════════════════════════════════════════
# PIPELINE API
# ══════════════════════════════════════════════════════════════════

@app.get("/api/pipeline/status")
async def pipeline_status():
    if state.pipeline:
        return state.pipeline.get_status()
    return {"state": "idle", "current_step_index": 0, "steps": []}


@app.post("/api/pipeline/run")
async def pipeline_run():
    if state.pipeline and state.pipeline.state == PipelineState.RUNNING:
        raise HTTPException(409, "Pipeline already running")

    state.pipeline = PipelineManager(state.config, log_callback=_log_and_status)
    state.pipeline_task = asyncio.create_task(_run_pipeline())
    return {"ok": True}


@app.post("/api/pipeline/stop")
async def pipeline_stop():
    if not state.pipeline or state.pipeline.state != PipelineState.RUNNING:
        raise HTTPException(409, "Pipeline not running")
    await state.pipeline.stop()
    await state.broadcast_status()
    return {"ok": True}


@app.post("/api/pipeline/resume")
async def pipeline_resume():
    if not state.pipeline:
        raise HTTPException(409, "No pipeline to resume")
    if state.pipeline.state == PipelineState.RUNNING:
        raise HTTPException(409, "Already running")
    state.pipeline_task = asyncio.create_task(_run_pipeline(resume=True))
    return {"ok": True}


async def _log_and_status(msg: str):
    await state.broadcast_log(msg)
    await state.broadcast_status()


async def _run_pipeline(resume: bool = False):
    try:
        if resume:
            await state.pipeline.resume()
        else:
            await state.pipeline.run()
    except Exception as e:
        await state.broadcast_log(f"Pipeline error: {e}")
    finally:
        await state.broadcast_status()


# ══════════════════════════════════════════════════════════════════
# FILE BROWSER API
# ══════════════════════════════════════════════════════════════════

@app.get("/api/files/browse")
async def browse_files(path: str = ".", pattern: str = "*"):
    """Browse the server filesystem."""
    try:
        target = Path(path).expanduser().resolve()
        if not target.is_dir():
            target = target.parent

        dirs = []
        files = []

        # Always include parent
        parent = str(target.parent)

        for item in sorted(target.iterdir()):
            name = item.name
            if name.startswith("."):
                continue
            try:
                if item.is_dir():
                    dirs.append({"name": name, "path": str(item)})
                elif item.is_file():
                    if pattern == "*" or item.match(pattern):
                        size = item.stat().st_size
                        files.append({
                            "name": name,
                            "path": str(item),
                            "size": size,
                            "size_str": _human_size(size),
                        })
            except PermissionError:
                continue

        return {
            "current": str(target),
            "parent": parent,
            "dirs": dirs[:100],
            "files": files[:200],
        }
    except PermissionError:
        raise HTTPException(403, "Permission denied")
    except Exception as e:
        raise HTTPException(400, str(e))


@app.post("/api/upload/pdb")
async def upload_pdb(file: UploadFile = File(...)):
    if not file.filename or not file.filename.lower().endswith(".pdb"):
        raise HTTPException(400, "Must be a .pdb file")

    content = await file.read()
    save_path = Path(file.filename)
    save_path.write_bytes(content)
    state.config.pdb_file = str(save_path.resolve())
    return {"ok": True, "path": str(save_path.resolve()), "filename": file.filename}


@app.get("/api/pdb/content")
async def get_pdb_content():
    p = state.config.pdb_file
    if not p or not os.path.isfile(p):
        raise HTTPException(404, "No PDB file loaded")
    with open(p, "r") as f:
        return {"content": f.read(), "filename": os.path.basename(p)}


def _human_size(nbytes: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if nbytes < 1024:
            return f"{nbytes:.0f} {unit}" if unit == "B" else f"{nbytes:.1f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} TB"


# ══════════════════════════════════════════════════════════════════
# SLURM API
# ══════════════════════════════════════════════════════════════════

@app.get("/api/slurm/available")
async def slurm_available():
    return {"available": slurm_mod.is_slurm_available()}


@app.get("/api/slurm/preview")
async def slurm_preview():
    """Preview the sbatch script that would be submitted."""
    script = slurm_mod.generate_sbatch_script(state.config)
    return {"script": script}


@app.post("/api/slurm/submit")
async def slurm_submit():
    try:
        job_id = await slurm_mod.submit_job(state.config)
        state.slurm_job_id = job_id
        await state.broadcast_log(f"SLURM job submitted: {job_id}")
        return {"ok": True, "job_id": job_id}
    except RuntimeError as e:
        raise HTTPException(400, str(e))


@app.get("/api/slurm/status")
async def slurm_status():
    if not state.slurm_job_id:
        return {"job_id": "", "state": "NONE"}
    info = await slurm_mod.get_job_status(state.slurm_job_id)
    return info.to_dict()


@app.post("/api/slurm/cancel")
async def slurm_cancel():
    if not state.slurm_job_id:
        raise HTTPException(400, "No SLURM job to cancel")
    ok = await slurm_mod.cancel_job(state.slurm_job_id)
    return {"ok": ok}


# ══════════════════════════════════════════════════════════════════
# ANALYSIS API
# ══════════════════════════════════════════════════════════════════

@app.get("/api/analysis/available")
async def analysis_available():
    return analysis_mod.get_available_outputs()


@app.get("/api/analysis/pka")
async def analysis_pka():
    data = analysis_mod.parse_pk_out()
    if "error" in data:
        raise HTTPException(404, data["error"])
    return data


@app.get("/api/analysis/sum_crg")
async def analysis_sum_crg():
    data = analysis_mod.parse_sum_crg()
    if "error" in data:
        raise HTTPException(404, data["error"])
    return data


@app.get("/api/analysis/conformers")
async def analysis_conformers():
    data = analysis_mod.parse_head3_lst()
    if "error" in data:
        raise HTTPException(404, data["error"])
    return data


@app.get("/api/analysis/dipole")
async def analysis_dipole():
    data = analysis_mod.parse_dipole_csv()
    if "error" in data:
        raise HTTPException(404, data["error"])
    return data


@app.get("/api/analysis/dipole/status")
async def dipole_status():
    """Check if dipole CSV exists or if ms_dipole can be run."""
    return analysis_mod.can_run_dipole()


@app.post("/api/analysis/dipole/run")
async def dipole_run():
    """Run ms_dipole ensemble computation and return results."""
    status = analysis_mod.can_run_dipole()
    if status["has_csv"]:
        # Already computed — just return existing data
        return analysis_mod.parse_dipole_csv()
    if not status["can_run"]:
        raise HTTPException(400, f"Missing required files: {', '.join(status['missing'])}")

    await state.broadcast_log(">> Running dipole analysis...")

    # Run in executor to avoid blocking the event loop
    loop = asyncio.get_event_loop()

    log_messages = []
    def sync_log(msg):
        log_messages.append(msg)

    try:
        data = await loop.run_in_executor(
            None, lambda: analysis_mod.run_dipole(log_callback=sync_log)
        )
    except Exception as e:
        await state.broadcast_log(f"Dipole error: {e}")
        raise HTTPException(500, str(e))

    # Broadcast collected log messages
    for msg in log_messages:
        await state.broadcast_log(msg)

    if "error" in data:
        await state.broadcast_log(f"Dipole error: {data['error']}")
        raise HTTPException(500, data["error"])

    await state.broadcast_log(">> Dipole analysis complete.")
    return data


@app.get("/api/pdb/step2_content")
async def get_step2_pdb_content():
    """Serve step2_out.pdb content for 3D dipole visualization."""
    p = "step2_out.pdb"
    if not os.path.isfile(p):
        raise HTTPException(404, "step2_out.pdb not found")
    with open(p, "r") as f:
        return {"content": f.read(), "filename": "step2_out.pdb"}


# ══════════════════════════════════════════════════════════════════
# WEBSOCKET
# ══════════════════════════════════════════════════════════════════

@app.websocket("/ws")
async def websocket_endpoint(ws: WebSocket):
    await ws.accept()
    state.ws_clients.add(ws)

    try:
        # Send recent log history
        for line in state.log_history[-200:]:
            await ws.send_json({"type": "log", "message": line})

        # Send current status
        if state.pipeline:
            await ws.send_json({"type": "status", "data": state.pipeline.get_status()})

        # Keep alive
        while True:
            try:
                data = await ws.receive_text()
                msg = json.loads(data)
                if msg.get("type") == "ping":
                    await ws.send_json({"type": "pong"})
            except WebSocketDisconnect:
                break
    finally:
        state.ws_clients.discard(ws)
