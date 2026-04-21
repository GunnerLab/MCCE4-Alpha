"""
CLI entry point for the MCCE4 Topology Agent.
"""

import argparse
import os
import sys
import textwrap
import logging
import subprocess
from datetime import datetime

from .config import (SUPPORTED_CHARGE_METHODS, DEFAULT_CHARGE_METHOD, DEFAULT_DIELECTRICS,
                      SUPPORTED_LLM_PROVIDERS, DEFAULT_LLM_PROVIDER,
                      DIMORPHITE_PH_MIN, DIMORPHITE_PH_MAX, DIMORPHITE_PRECISION,
                      DIMORPHITE_MAX_VARIANTS, DIMORPHITE_LABEL_STATES)


def setup_logging(log_file: str, verbose: bool = False):
    """Configure logging."""
    log_level = logging.DEBUG if verbose else logging.INFO
    file_fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s",
                                  datefmt="%Y-%m-%d %H:%M:%S")
    console_fmt = logging.Formatter("%(levelname)-8s | %(message)s")
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_file, mode="w")
    fh.setLevel(logging.DEBUG); fh.setFormatter(file_fmt)
    logger.addHandler(fh)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level); ch.setFormatter(console_fmt)
    logger.addHandler(ch)


def main():
    parser = argparse.ArgumentParser(
        description="🤖 MCCE4 Topology File AI Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
        Examples:
          %(prog)s EMH.pdb                                         # Full auto (PDB input)
          %(prog)s EMH.cif                                         # CIF input (auto-converts)
          %(prog)s --lig-code EMH                                  # Ligand code only (RCSB)
          %(prog)s EMH.pdb --lig-code EMH                          # PDB + ligand code
          %(prog)s --lig-code EMH --dry-run --no-llm               # Quick test
          %(prog)s EMH.pdb --gui                                   # Web GUI (Streamlit)
          %(prog)s EMH.pdb --state-pdbs EMH_01.pdb EMH_+1.pdb     # User state PDBs
          %(prog)s EMH.pdb --no-llm                                # No LLM reasoning
          %(prog)s EMH.pdb --charge-method am1bcc                  # Override charges
          %(prog)s EMH.pdb --dry-run                               # Skip calibration
          %(prog)s EMH.pdb --llm-provider claude --api-key sk-...  # Use Claude
          %(prog)s EMH.pdb --llm-provider chatgpt --api-key sk-... # Use ChatGPT

        Install:
          pip install google-genai langgraph dimorphite-dl rdkit streamlit
          pip install anthropic openai  # optional, for Claude/ChatGPT LLM providers
          export GEMINI_API_KEY="your_free_key"   # from https://ai.google.dev
          export ANTHROPIC_API_KEY="your_key"     # for Claude provider
          export OPENAI_API_KEY="your_key"        # for ChatGPT provider
        """))

    parser.add_argument("input_file", nargs="?", default=None,
                        help="Ligand PDB or CIF file (e.g., EMH.pdb or EMH.cif). "
                             "Optional if --lig-code is provided.")
    parser.add_argument("--lig-code", default=None,
                        help="3-letter RCSB ligand code (e.g., EMH). When provided, "
                             "SMILES and metadata are fetched from RCSB. A PDB/CIF "
                             "file is optional — if omitted, the 3D structure is "
                             "built from SMILES using RDKit.")
    parser.add_argument("--state-pdbs", nargs="+", default=None,
                        help="PDB files for specific states (e.g., EMH_01.pdb EMH_+1.pdb)")
    parser.add_argument("--gui", action="store_true",
                        help="Launch web GUI (Streamlit) for interactive review")
    parser.add_argument("--ph", type=float, default=7.4, help="Target pH (default: 7.4)")
    parser.add_argument("--charge-method", default=DEFAULT_CHARGE_METHOD,
                        choices=SUPPORTED_CHARGE_METHODS, help="Charge method")
    parser.add_argument("-d", "--dielectric", nargs="+", type=int, default=DEFAULT_DIELECTRICS,
                        help="Dielectric constants (default: 2 4 8)")
    parser.add_argument("--no-llm", action="store_true", help="Disable LLM reasoning")
    parser.add_argument("--dry-run", action="store_true", help="Skip RXN calibration")

    # Dimorphite-DL options
    dimorphite_group = parser.add_argument_group("Dimorphite-DL options",
        "Control protonation state enumeration (Phase 2)")
    dimorphite_group.add_argument("--ph-min", type=float,
                                  default=DIMORPHITE_PH_MIN,
                                  help=f"Minimum pH for protonation enumeration "
                                       f"(default: {DIMORPHITE_PH_MIN})")
    dimorphite_group.add_argument("--ph-max", type=float,
                                  default=DIMORPHITE_PH_MAX,
                                  help=f"Maximum pH for protonation enumeration "
                                       f"(default: {DIMORPHITE_PH_MAX})")
    dimorphite_group.add_argument("--precision", type=float,
                                  default=DIMORPHITE_PRECISION,
                                  help=f"pKa precision factor: number of std deviations "
                                       f"from mean pKa to consider (default: {DIMORPHITE_PRECISION})")
    dimorphite_group.add_argument("--max-variants", type=int,
                                  default=DIMORPHITE_MAX_VARIANTS,
                                  help=f"Max protonation variants per compound "
                                       f"(default: {DIMORPHITE_MAX_VARIANTS})")
    dimorphite_group.add_argument("--label-states", action="store_true",
                                  default=DIMORPHITE_LABEL_STATES,
                                  help="Label output SMILES as PROTONATED/DEPROTONATED/BOTH")
    parser.add_argument("--work-dir", default=".", help="Working directory")
    parser.add_argument("-o", "--output", default=None, help="Output .ftpl filename")
    parser.add_argument("-v", "--verbose", action="store_true")

    # LLM provider options
    parser.add_argument("--llm-provider", default=DEFAULT_LLM_PROVIDER,
                        choices=SUPPORTED_LLM_PROVIDERS,
                        help=f"LLM provider (default: {DEFAULT_LLM_PROVIDER})")
    parser.add_argument("--api-key", default=None,
                        help="API key for the chosen LLM provider (overrides env vars)")

    args = parser.parse_args()

    # ── Validate inputs: need either input_file or --lig-code ──
    if args.input_file is None and args.lig_code is None:
        parser.error("Either an input PDB/CIF file or --lig-code must be provided.")

    lig_code = args.lig_code.upper() if args.lig_code else None
    pdb_path = None

    if args.input_file is not None:
        if not os.path.exists(args.input_file):
            print(f"ERROR: Input file not found: {args.input_file}")
            sys.exit(1)

        # ── CIF → PDB conversion ──
        pdb_path = args.input_file
        if args.input_file.lower().endswith(".cif"):
            pdb_path = _convert_cif_to_pdb(args.input_file)
        elif not args.input_file.lower().endswith(".pdb"):
            print(f"ERROR: Input file must be .pdb or .cif: {args.input_file}")
            sys.exit(1)

    # Determine lig_id: --lig-code takes precedence, else extract from PDB
    if lig_code:
        lig_id = lig_code
    else:
        from .tools.mcce_tools import extract_lig_id_from_pdb
        lig_id = extract_lig_id_from_pdb(pdb_path)

    # ── GUI mode: launch Streamlit ──
    if args.gui:
        if pdb_path is None:
            print("ERROR: GUI mode requires an input PDB/CIF file.")
            sys.exit(1)
        args.pdb = pdb_path  # GUI expects args.pdb
        launch_gui(args)
        return

    # ── CLI mode: run agent directly ──
    setup_logging(f"mcce4_agent_ftpl_{lig_id}.log", args.verbose)

    logging.info(f"{'='*60}")
    logging.info(f"  🤖 MCCE4 Topology Agent v5.1 (LangGraph + per-state PDBs)")
    logging.info(f"{'='*60}")
    if pdb_path:
        logging.info(f"  Input:  {os.path.abspath(pdb_path)}")
        if args.input_file and args.input_file.lower().endswith(".cif"):
            logging.info(f"  (converted from {os.path.abspath(args.input_file)})")
    else:
        logging.info(f"  Input:  --lig-code {lig_id} (no PDB file — will build from SMILES)")
    logging.info(f"  Ligand: {lig_id}   pH: {args.ph}   Method: {args.charge_method}")
    logging.info(f"  Dimorphite-DL: ph_min={args.ph_min}, ph_max={args.ph_max}, "
                 f"precision={args.precision}, max_variants={args.max_variants}"
                 f"{', label_states=True' if args.label_states else ''}")
    logging.info(f"  LLM:    {args.llm_provider}")
    logging.info(f"  Time:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"{'='*60}\n")

    # Set env to disable LLM if requested
    if args.no_llm:
        os.environ.pop("GEMINI_API_KEY", None)
        os.environ.pop("GOOGLE_API_KEY", None)
        os.environ.pop("ANTHROPIC_API_KEY", None)
        os.environ.pop("OPENAI_API_KEY", None)

    from .agent import run_agent

    final_state = run_agent(
        pdb_path=pdb_path,
        use_gui=False,
        charge_method=args.charge_method,
        dielectrics=args.dielectric,
        ph=args.ph,
        work_dir=args.work_dir,
        dry_run=args.dry_run,
        user_state_pdbs=args.state_pdbs,
        output=args.output,
        llm_provider=args.llm_provider,
        api_key=args.api_key,
        lig_code=lig_code,
        ph_min=args.ph_min,
        ph_max=args.ph_max,
        precision=args.precision,
        max_variants=args.max_variants,
        label_states=args.label_states,
    )

    # Exit code based on errors
    if final_state.get("errors"):
        sys.exit(1)


def _convert_cif_to_pdb(cif_path: str) -> str:
    """Convert a .cif file to .pdb using cif2pdb_PyMOL.

    Returns the path to the generated PDB file.
    """
    # Locate cif2pdb_PyMOL in the same bin/ directory
    bin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    converter = os.path.join(bin_dir, "cif2pdb_PyMOL")

    if not os.path.exists(converter):
        print(f"ERROR: cif2pdb_PyMOL not found at {converter}")
        print("  This script is required to convert .cif files to .pdb format.")
        sys.exit(1)

    pdb_path = os.path.splitext(cif_path)[0] + ".pdb"
    print(f"🔄 Converting {cif_path} → {pdb_path} using cif2pdb_PyMOL ...")

    try:
        result = subprocess.run(
            [sys.executable, converter, cif_path, pdb_path],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            print(f"ERROR: CIF to PDB conversion failed:")
            print(result.stderr or result.stdout)
            sys.exit(1)
        if not os.path.exists(pdb_path):
            print(f"ERROR: Conversion completed but {pdb_path} was not created.")
            sys.exit(1)
        print(f"✅ Converted successfully: {pdb_path}")
    except subprocess.TimeoutExpired:
        print(f"ERROR: CIF to PDB conversion timed out.")
        sys.exit(1)
    except FileNotFoundError:
        print(f"ERROR: Could not run cif2pdb_PyMOL. Ensure PyMOL is installed.")
        sys.exit(1)

    return pdb_path


def launch_gui(args):
    """Launch the Streamlit web GUI."""
    gui_path = os.path.join(os.path.dirname(__file__), "gui", "app.py")

    if not os.path.exists(gui_path):
        print(f"ERROR: GUI app not found at {gui_path}")
        sys.exit(1)

    port = 8501

    # Check if port is already in use and handle it
    if _is_port_in_use(port):
        print(f"⚠  Port {port} is already in use (previous Streamlit session?).")
        print()

        # Try to find the PID
        pid = _find_pid_on_port(port)
        if pid:
            print(f"   Found process PID {pid} on port {port}.")
            try:
                response = input(f"   Kill it and restart? [Y/n] ").strip().lower()
            except EOFError:
                response = "y"

            if response in ("", "y", "yes"):
                import signal
                try:
                    os.kill(pid, signal.SIGTERM)
                    print(f"   Killed PID {pid}. Waiting for port to free...")
                    import time
                    for _ in range(10):
                        time.sleep(0.5)
                        if not _is_port_in_use(port):
                            break
                    if _is_port_in_use(port):
                        os.kill(pid, signal.SIGKILL)
                        time.sleep(1)
                except ProcessLookupError:
                    pass  # Already dead
                except PermissionError:
                    print(f"   Permission denied. Run manually:")
                    print(f"     kill {pid}")
                    print(f"   Then rerun this command.")
                    sys.exit(1)
            else:
                print(f"\n   To free the port manually, run:")
                print(f"     kill {pid}")
                print(f"   Or use a different port:")
                print(f"     streamlit run {gui_path} --server.port 8502 -- {os.path.abspath(args.pdb)}")
                sys.exit(1)
        else:
            print(f"   Could not identify the process. To free port {port}, run:")
            print(f"     pkill -f 'streamlit run'")
            print(f"     # or: kill $(lsof -ti:{port})")
            print(f"   Then rerun this command.")
            sys.exit(1)

    # Pass PDB path via environment so Streamlit can access it
    os.environ["MCCE_AGENT_PDB"] = os.path.abspath(args.pdb)
    os.environ["MCCE_AGENT_PH"] = str(args.ph)
    os.environ["MCCE_AGENT_CHARGE_METHOD"] = args.charge_method

    print(f"🌐 Launching MCCE4 Topology Agent GUI...")
    print(f"   PDB: {args.pdb}")
    print(f"   Open your browser to: http://localhost:{port}")
    print(f"   (If remote: ssh -L {port}:localhost:{port} user@server)")
    print()

    try:
        subprocess.run([
            sys.executable, "-m", "streamlit", "run", gui_path,
            "--server.headless", "true",
            "--server.port", str(port),
            "--browser.gatherUsageStats", "false",
        ])
    except KeyboardInterrupt:
        print("\n  GUI stopped.")
    except FileNotFoundError:
        print("ERROR: Streamlit not installed — run: pip install streamlit")
        sys.exit(1)


def _is_port_in_use(port: int) -> bool:
    """Check if a TCP port is in use."""
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(("localhost", port)) == 0


def _find_pid_on_port(port: int):
    """Find PID of process using the given port. Returns int or None."""
    try:
        result = subprocess.run(
            ["lsof", "-ti", f":{port}"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0 and result.stdout.strip():
            # May return multiple PIDs; take the first
            return int(result.stdout.strip().split()[0])
    except (FileNotFoundError, subprocess.TimeoutExpired, ValueError):
        pass

    # Fallback: try ss
    try:
        result = subprocess.run(
            ["ss", "-tlnp", f"sport = :{port}"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            import re
            m = re.search(r'pid=(\d+)', result.stdout)
            if m:
                return int(m.group(1))
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass

    return None


if __name__ == "__main__":
    main()
