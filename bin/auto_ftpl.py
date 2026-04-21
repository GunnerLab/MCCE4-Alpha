#!/usr/bin/env python3
"""
auto_ftpl.py — AI Agent for Automated MCCE4 Topology File Creation
==================================================================

Automates the full pipeline for creating MCCE4 .ftpl topology files:
  1. Download ideal ligand PDB from Ligand Expo (or use a user-supplied PDB)
  2. Generate template .ftpl via pdb2ftpl.py
  3. Generate partial atomic charges via OpenEye QuacPac TK (or other backends)
  4. Fill 'to_be_filled' charge entries in the .ftpl template
  5. Run MCCE4 steps 1-3 for RXN (reaction field / desolvation) calibration
  6. Parse head3.lst to extract dsolv values
  7. Update CONFORMER rxn parameters and validate (dsolv ≈ 0.000)

Usage:
  auto_ftpl.py EMH -c 01 +1                          # Full auto with mmff94 (default)
  auto_ftpl.py EMH -c 01 +1 -m am1bcc                # Use AM1-BCC charges
  auto_ftpl.py EMH -c 01 +1 --pdb EMH_ideal.pdb      # Use a local PDB file
  auto_ftpl.py EMH -c 01 +1 --charges-file charges.txt  # Use pre-computed charges
  auto_ftpl.py EMH -c 01 +1 -d 4                      # Calibrate only rxn04
  auto_ftpl.py EMH -c 01 +1 --dry-run                 # Generate ftpl but skip calibration
  auto_ftpl.py EMH -c 01 +1 --skip-download           # Skip ligand download, use existing PDB

Author:  Gehan / MCCE4 Team
Version: 1.0.0
"""

import argparse
import os
import sys
import re
import subprocess
import shutil
import logging
import json
from pathlib import Path
from datetime import datetime
from typing import Optional

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────
LIGAND_EXPO_URL = "http://ligand-expo.rcsb.org/reports/{first_char}/{lig_id}/{lig_id}_ideal.pdb"

SUPPORTED_CHARGE_METHODS = [
    "mmff94",           # OpenEye MMFF94 — reliable default, works on most molecules
    "am1bcc",           # OpenEye AM1-BCC — higher quality, small molecules only
    "am1bccelf10",      # OpenEye AM1-BCC with ELF10 conformer selection
    "am1bccnosymspt",   # OpenEye AM1-BCC no symmetry, single point
    "amber",            # OpenEye Amber charges
    "amberff94",        # OpenEye Amber FF94 charges
    "antechamber",      # AmberTools antechamber AM1-BCC (free, no OpenEye needed)
    "file",             # User-supplied charge file
]

DEFAULT_CHARGE_METHOD = "mmff94"

DIELECTRIC_MAP = {2: "rxn02", 4: "rxn04", 8: "rxn08"}

# ──────────────────────────────────────────────────────────────────────────────
# Logging setup
# ──────────────────────────────────────────────────────────────────────────────
def setup_logging(log_file: str = "auto_ftpl.log", verbose: bool = False):
    """Configure logging to both file and console."""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    # Create formatters
    file_fmt = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    console_fmt = logging.Formatter("%(levelname)-8s | %(message)s")
    
    # Root logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    # File handler — always DEBUG level for full trace
    fh = logging.FileHandler(log_file, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(file_fmt)
    logger.addHandler(fh)
    
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level)
    ch.setFormatter(console_fmt)
    logger.addHandler(ch)
    
    return logger


# ──────────────────────────────────────────────────────────────────────────────
# Utility functions
# ──────────────────────────────────────────────────────────────────────────────
def run_cmd(cmd: str, description: str = "", capture: bool = False,
            cwd: Optional[str] = None) -> subprocess.CompletedProcess:
    """Run a shell command with logging and error handling.
    
    Args:
        cmd:         Shell command string to execute.
        description: Human-readable description for logging.
        capture:     If True, capture stdout/stderr for parsing.
        cwd:         Working directory for the command.
    
    Returns:
        subprocess.CompletedProcess result.
    
    Raises:
        SystemExit: If command returns non-zero exit code.
    """
    if description:
        logging.info(f"▶ {description}")
    logging.debug(f"  CMD: {cmd}")
    
    result = subprocess.run(
        cmd, shell=True,
        capture_output=capture,
        text=True,
        cwd=cwd
    )
    
    if result.returncode != 0:
        logging.error(f"  Command failed (exit code {result.returncode})")
        if capture and result.stderr:
            logging.error(f"  STDERR: {result.stderr.strip()}")
        sys.exit(1)
    
    return result


def check_dependencies(charge_method: str):
    """Verify required tools are available on PATH."""
    required = ["pdb2ftpl.py", "step1.py", "step2.py", "step3.py"]
    
    if charge_method in ("mmff94", "am1bcc", "am1bccelf10",
                         "am1bccnosymspt", "amber", "amberff94"):
        required.append("oe_assigncharges_QuacpakTK.py")
    elif charge_method == "antechamber":
        required.append("antechamber")
    
    missing = []
    for tool in required:
        if shutil.which(tool) is None:
            missing.append(tool)
    
    if missing:
        logging.error(f"Missing required tools: {', '.join(missing)}")
        logging.error("Ensure MCCE4-Alpha and relevant charge tools are on your PATH.")
        sys.exit(1)
    
    logging.info(f"✓ All dependencies found for method '{charge_method}'")


# ──────────────────────────────────────────────────────────────────────────────
# Step 1: Obtain ligand PDB
# ──────────────────────────────────────────────────────────────────────────────
def download_ligand_pdb(lig_id: str, output_path: str):
    """Download ideal PDB from RCSB Ligand Expo.
    
    Args:
        lig_id:      3-letter ligand code (e.g., 'EMH').
        output_path: Path to save the downloaded PDB file.
    """
    url = LIGAND_EXPO_URL.format(first_char=lig_id[0], lig_id=lig_id)
    logging.info(f"📥 Downloading ideal PDB for '{lig_id}' from Ligand Expo...")
    logging.debug(f"  URL: {url}")
    
    result = run_cmd(
        f'wget -q -O {output_path} "{url}"',
        description=f"Downloading {lig_id}_ideal.pdb",
        capture=True
    )
    
    # Validate the download
    if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
        logging.error(f"Download failed or file is empty: {output_path}")
        logging.error(f"Try downloading manually from: {url}")
        logging.error(f"Or supply a local PDB with: --pdb {lig_id}.pdb")
        sys.exit(1)
    
    # Quick sanity check — must contain ATOM or HETATM lines
    with open(output_path, "r") as f:
        content = f.read()
        if "ATOM" not in content and "HETATM" not in content:
            logging.error(f"Downloaded file does not appear to be a valid PDB: {output_path}")
            sys.exit(1)
    
    atom_count = sum(1 for line in content.splitlines()
                     if line.startswith(("ATOM", "HETATM")))
    logging.info(f"  ✓ Downloaded {lig_id}.pdb ({atom_count} atoms)")


def convert_cif_to_pdb(cif_path: str, pdb_path: str):
    """Convert .cif to .pdb using MCCE4's PyMOL-based converter.
    
    Args:
        cif_path: Path to input .cif file.
        pdb_path: Path to output .pdb file.
    """
    logging.info(f"🔄 Converting {cif_path} → {pdb_path} via cif2pdb_PyMOL")
    run_cmd(
        f"cif2pdb_PyMOL {cif_path}",
        description="CIF to PDB conversion"
    )
    
    # cif2pdb_PyMOL outputs to same directory with .pdb extension
    expected = cif_path.replace(".cif", ".pdb")
    if os.path.exists(expected) and expected != pdb_path:
        shutil.move(expected, pdb_path)


# ──────────────────────────────────────────────────────────────────────────────
# Step 2: Generate .ftpl template
# ──────────────────────────────────────────────────────────────────────────────
def generate_ftpl_template(lig_id: str, pdb_path: str, conformers: list,
                           ftpl_path: str):
    """Generate .ftpl template using pdb2ftpl.py.
    
    Args:
        lig_id:     3-letter ligand code.
        pdb_path:   Path to the ligand PDB file.
        conformers: List of conformer state strings (e.g., ['01', '+1']).
        ftpl_path:  Output path for the .ftpl file.
    """
    conf_str = " ".join(conformers)
    cmd = f"pdb2ftpl.py -p {pdb_path} -c {conf_str}"
    
    logging.info(f"📝 Generating .ftpl template for {lig_id} with conformers: {conf_str}")
    result = run_cmd(cmd, description="Running pdb2ftpl.py", capture=True)
    
    # pdb2ftpl.py writes to stdout
    with open(ftpl_path, "w") as f:
        f.write(result.stdout)
    
    # Count to_be_filled entries
    unfilled = result.stdout.count("to_be_filled")
    logging.info(f"  ✓ Template created: {ftpl_path} ({unfilled} charge entries to fill)")
    
    return unfilled


# ──────────────────────────────────────────────────────────────────────────────
# Step 3: Generate charges
# ──────────────────────────────────────────────────────────────────────────────
def generate_charges_openeye(pdb_path: str, method: str, lig_id: str) -> dict:
    """Generate charges using OpenEye QuacPac TK.
    
    IMPORTANT: This should be run on an isolated LIGAND PDB, not a full protein.
    AM1-BCC will fail on proteins; MMFF94 works on both but is slower on large
    molecules.
    
    Args:
        pdb_path: Path to the ligand-only PDB file.
        method:   OpenEye charge method (mmff94, am1bcc, etc.).
        lig_id:   3-letter ligand code for context.
    
    Returns:
        dict: Mapping of atom_name → charge (e.g., {' C4 ': -0.084, ...}).
    """
    logging.info(f"⚡ Generating charges via OpenEye QuacPac TK (method: {method})")
    
    result = run_cmd(
        f"oe_assigncharges_QuacpakTK.py -method {method} -in {pdb_path}",
        description=f"OpenEye {method} charge assignment",
        capture=True
    )
    
    charges = parse_openeye_output(result.stdout, lig_id)
    
    if not charges:
        logging.error(f"  No charges parsed from OpenEye output!")
        logging.error(f"  This can happen if:")
        logging.error(f"    - The input PDB is a full protein (use ligand-only PDB)")
        logging.error(f"    - The method '{method}' doesn't support this molecule")
        logging.error(f"    - OpenEye license is not configured")
        sys.exit(1)
    
    # Validate: check for all-zero charges (AM1-BCC failure mode)
    non_zero = sum(1 for v in charges.values() if abs(v) > 1e-6)
    if non_zero == 0:
        logging.warning(f"  ⚠ All charges are zero — method '{method}' likely failed!")
        logging.warning(f"  Falling back to mmff94...")
        if method != "mmff94":
            return generate_charges_openeye(pdb_path, "mmff94", lig_id)
        else:
            logging.error("  MMFF94 also returned all zeros. Check input PDB.")
            sys.exit(1)
    
    # Report charge summary
    total_charge = sum(charges.values())
    logging.info(f"  ✓ Parsed {len(charges)} atom charges (total charge: {total_charge:.3f})")
    logging.debug(f"  Charge range: [{min(charges.values()):.3f}, {max(charges.values()):.3f}]")
    
    return charges


def parse_openeye_output(output: str, lig_id: str) -> dict:
    """Parse OpenEye oe_assigncharges_QuacpakTK.py output into atom→charge dict.
    
    Expected format per line:
      Atom Name:  C4  | Symbol:  C | Charge: -0.084
    
    The challenge: OpenEye may output atoms for ALL residues (waters, etc.).
    We need to filter to only the ligand atoms and handle duplicate atom names
    from water molecules.
    
    Strategy:
      - Parse all lines matching the pattern
      - Skip obvious water atoms (H1, H2, O with typical water charges)
      - Build unique mapping using MCCE 4-char padded names
    
    Args:
        output:  Raw stdout from oe_assigncharges_QuacpakTK.py.
        lig_id:  3-letter ligand code (for logging context).
    
    Returns:
        dict: Mapping of 4-char padded atom name → charge value.
    """
    charges = {}
    pattern = re.compile(
        r"Atom Name:\s+(\S+)\s+\|\s+Symbol:\s+(\S+)\s+\|\s+Charge:\s+([-\d.]+)"
    )
    
    # Track seen atoms to detect water duplicates
    atom_counts = {}
    raw_entries = []
    
    for line in output.splitlines():
        m = pattern.search(line)
        if m:
            atom_name = m.group(1)
            symbol = m.group(2)
            charge = float(m.group(3))
            raw_entries.append((atom_name, symbol, charge))
            atom_counts[atom_name] = atom_counts.get(atom_name, 0) + 1
    
    # Filter out likely water atoms:
    # If an atom name appears many times (>5), it's probably from waters
    water_atoms = {name for name, count in atom_counts.items() if count > 5}
    if water_atoms:
        logging.debug(f"  Filtering out likely water atoms: {water_atoms}")
    
    seen_names = set()
    for atom_name, symbol, charge in raw_entries:
        if atom_name in water_atoms:
            continue
        
        # Convert to MCCE 4-character padded atom name
        mcce_name = to_mcce_atom_name(atom_name)
        
        if mcce_name not in seen_names:
            charges[mcce_name] = charge
            seen_names.add(mcce_name)
    
    return charges


def to_mcce_atom_name(name: str) -> str:
    """Convert an atom name to MCCE's 4-character padded format.
    
    MCCE uses a strict 4-character atom name format inside quotes:
      - Element symbol starts at position determined by atom count:
        * 1-char elements: ' X  ' (space-padded, left-justified from col 2)
        * 2-char elements: 'XX  ' (or with numbers like 'H21A')
      - Names from PDB typically follow:  ' C4 ', ' N9 ', 'H21A'
    
    Examples:
        'C4'   → ' C4 '
        'H21A' → 'H21A'
        'N9'   → ' N9 '
        'O20'  → ' O20'
        'H4'   → ' H4 '
    
    Args:
        name: Raw atom name string.
    
    Returns:
        4-character MCCE-formatted atom name.
    """
    name = name.strip()
    
    if len(name) >= 4:
        return name[:4]
    elif len(name) == 3:
        # Check if first char is a letter and second is a digit
        # e.g., 'C4 ' or 'O20' or 'N19'
        return f" {name}" if name[0].isalpha() and len(name) == 2 else f" {name}"
    elif len(name) == 3:
        return f" {name}"
    elif len(name) == 2:
        return f" {name} "
    elif len(name) == 1:
        return f" {name}  "
    else:
        return name.ljust(4)


def generate_charges_antechamber(pdb_path: str, lig_id: str,
                                  net_charge: int = 0) -> dict:
    """Generate AM1-BCC charges using AmberTools antechamber (free alternative).
    
    Args:
        pdb_path:   Path to the ligand PDB file.
        lig_id:     3-letter ligand code.
        net_charge: Net formal charge of the molecule.
    
    Returns:
        dict: Mapping of atom_name → charge.
    """
    logging.info(f"⚡ Generating charges via antechamber AM1-BCC (nc={net_charge})")
    
    mol2_out = f"{lig_id}_antechamber.mol2"
    
    run_cmd(
        f"antechamber -i {pdb_path} -fi pdb -o {mol2_out} -fo mol2 "
        f"-c bcc -s 2 -nc {net_charge} -at gaff",
        description="antechamber AM1-BCC charge assignment",
        capture=True
    )
    
    charges = parse_mol2_charges(mol2_out)
    
    if not charges:
        logging.error("  No charges parsed from antechamber output!")
        sys.exit(1)
    
    total = sum(charges.values())
    logging.info(f"  ✓ Parsed {len(charges)} atom charges (total: {total:.3f})")
    
    return charges


def parse_mol2_charges(mol2_path: str) -> dict:
    """Parse atom charges from a Tripos .mol2 file.
    
    The @<TRIPOS>ATOM section contains lines like:
      1 C4  7.534  0.344  17.581 c3  1 LIG -0.075900
    
    Args:
        mol2_path: Path to the .mol2 file.
    
    Returns:
        dict: Mapping of MCCE-formatted atom name → charge.
    """
    charges = {}
    in_atom_section = False
    
    with open(mol2_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue
            elif line.startswith("@<TRIPOS>"):
                in_atom_section = False
                continue
            
            if in_atom_section and line:
                parts = line.split()
                if len(parts) >= 9:
                    atom_name = parts[1]
                    charge = float(parts[8])
                    mcce_name = to_mcce_atom_name(atom_name)
                    charges[mcce_name] = charge
    
    return charges


def load_charges_from_file(charges_file: str) -> dict:
    """Load pre-computed charges from a user-supplied file.
    
    Supports multiple formats:
      - Simple: 'atom_name  charge' per line
      - CSV:    'atom_name,charge' per line
      - JSON:   {'atom_name': charge, ...}
    
    Args:
        charges_file: Path to the charges file.
    
    Returns:
        dict: Mapping of MCCE-formatted atom name → charge.
    """
    logging.info(f"📂 Loading charges from file: {charges_file}")
    
    charges = {}
    
    with open(charges_file, "r") as f:
        content = f.read().strip()
    
    # Try JSON first
    if content.startswith("{"):
        try:
            raw = json.loads(content)
            for name, charge in raw.items():
                charges[to_mcce_atom_name(name)] = float(charge)
            logging.info(f"  ✓ Loaded {len(charges)} charges from JSON")
            return charges
        except json.JSONDecodeError:
            pass
    
    # Try line-by-line (space or comma separated)
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        
        # Try comma separation, then space
        if "," in line:
            parts = line.split(",")
        else:
            parts = line.split()
        
        if len(parts) >= 2:
            atom_name = parts[0].strip().strip('"').strip("'")
            try:
                charge = float(parts[-1].strip())
                charges[to_mcce_atom_name(atom_name)] = charge
            except ValueError:
                continue
    
    logging.info(f"  ✓ Loaded {len(charges)} charges from text file")
    return charges


# ──────────────────────────────────────────────────────────────────────────────
# Step 4: Fill charges into .ftpl template
# ──────────────────────────────────────────────────────────────────────────────
def fill_ftpl_charges(ftpl_path: str, charges: dict, conformers: list,
                       lig_id: str, output_path: Optional[str] = None) -> int:
    """Replace 'to_be_filled' entries in the .ftpl with actual charge values.
    
    Handles the mapping between charge dict atom names and MCCE conformer-specific
    CHARGE lines. For molecules with multiple conformer states (e.g., 01 and +1),
    the same charge set is applied to all conformers unless state-specific charges
    are provided.
    
    Args:
        ftpl_path:   Path to the template .ftpl file.
        charges:     dict of MCCE atom name → charge.
        conformers:  List of conformer state strings.
        lig_id:      3-letter ligand code.
        output_path: Output path (defaults to overwriting ftpl_path).
    
    Returns:
        int: Number of entries that could NOT be matched.
    """
    if output_path is None:
        output_path = ftpl_path
    
    logging.info(f"✏️  Filling charges in {ftpl_path}")
    
    with open(ftpl_path, "r") as f:
        lines = f.readlines()
    
    # Pattern: CHARGE, EMH01, " C4 ": to_be_filled
    charge_pattern = re.compile(
        r'^(CHARGE,\s*' + re.escape(lig_id) + r'(\S+),\s*"(.{4})":\s*)to_be_filled(.*)$'
    )
    
    filled = 0
    unfilled = 0
    unmatched_atoms = set()
    new_lines = []
    
    for line in lines:
        m = charge_pattern.match(line)
        if m:
            prefix = m.group(1)     # 'CHARGE, EMH01, " C4 ": '
            conf_id = m.group(2)    # '01' or '+1'
            atom_name = m.group(3)  # ' C4 '
            suffix = m.group(4)     # any trailing content
            
            if atom_name in charges:
                charge_val = charges[atom_name]
                # Format charge: right-aligned, 6 chars with 3 decimal places
                charge_str = f"{charge_val:7.3f}"
                new_line = f"{prefix}{charge_str} # auto-filled{suffix}\n"
                new_lines.append(new_line)
                filled += 1
            else:
                # Try matching without padding
                stripped = atom_name.strip()
                found = False
                for cname, cval in charges.items():
                    if cname.strip() == stripped:
                        charge_str = f"{cval:7.3f}"
                        new_line = f"{prefix}{charge_str} # auto-filled{suffix}\n"
                        new_lines.append(new_line)
                        filled += 1
                        found = True
                        break
                
                if not found:
                    new_lines.append(line)  # Keep to_be_filled
                    unfilled += 1
                    unmatched_atoms.add(atom_name.strip())
        else:
            new_lines.append(line)
    
    with open(output_path, "w") as f:
        f.writelines(new_lines)
    
    logging.info(f"  ✓ Filled {filled} charge entries")
    
    if unfilled > 0:
        logging.warning(f"  ⚠ {unfilled} entries remain unfilled!")
        logging.warning(f"  Unmatched atoms: {sorted(unmatched_atoms)}")
        logging.warning(f"  Available charge keys: {sorted(charges.keys())}")
        logging.warning(f"  You may need to manually fill these in {output_path}")
    
    return unfilled


# ──────────────────────────────────────────────────────────────────────────────
# Step 5-7: RXN Calibration
# ──────────────────────────────────────────────────────────────────────────────
def run_rxn_calibration(lig_id: str, ftpl_path: str, pdb_path: str,
                         conformers: list, dielectrics: list, work_dir: str):
    """Run MCCE4 steps 1-3, parse dsolv, update rxn values, and validate.
    
    This is the calibration loop:
      1. Ensure user_param/ exists and ftpl is linked
      2. Run step1.py <pdb>, step2.py, step3.py -d <eps>
      3. Parse head3.lst for dsolv per conformer type
      4. Update CONFORMER lines in the .ftpl with rxn values
      5. Re-run step3.py to validate (dsolv → 0.000)
    
    Args:
        lig_id:      3-letter ligand code.
        ftpl_path:   Path to the .ftpl file.
        pdb_path:    Path to the ligand PDB file (passed to step1.py).
        conformers:  List of conformer state strings.
        dielectrics: List of dielectric constants to calibrate (e.g., [4]).
        work_dir:    MCCE4 working directory.
    """
    # Setup user_param directory
    user_param_dir = os.path.join(work_dir, "user_param")
    os.makedirs(user_param_dir, exist_ok=True)
    
    # Symlink the ftpl into user_param
    link_path = os.path.join(user_param_dir, os.path.basename(ftpl_path))
    if os.path.exists(link_path):
        os.remove(link_path)
    os.symlink(os.path.abspath(ftpl_path), link_path)
    logging.info(f"  🔗 Linked {ftpl_path} → user_param/")
    
    # Ensure the ligand PDB is accessible from work_dir
    pdb_basename = os.path.basename(pdb_path)
    pdb_in_workdir = os.path.join(work_dir, pdb_basename)
    if not os.path.exists(pdb_in_workdir):
        shutil.copy2(os.path.abspath(pdb_path), pdb_in_workdir)
        logging.info(f"  📄 Copied {pdb_path} → {work_dir}/")
    
    # Run steps 1-2 once (they don't depend on dielectric constant)
    run_cmd(f"step1.py {pdb_basename}",
            description=f"MCCE4 Step 1 (preparation) with {pdb_basename}",
            cwd=work_dir)
    run_cmd("step2.py", description="MCCE4 Step 2 (rotamers)", cwd=work_dir)
    
    # Calibrate rxn for each dielectric constant
    for eps in dielectrics:
        rxn_key = DIELECTRIC_MAP.get(eps)
        if not rxn_key:
            logging.warning(f"  Unsupported dielectric constant: {eps}. Skipping.")
            continue
        
        logging.info(f"\n{'='*60}")
        logging.info(f"  🔬 RXN Calibration for ε = {eps} ({rxn_key})")
        logging.info(f"{'='*60}")
        
        # Run MCCE4 step 3 with this dielectric
        run_cmd(f"step3.py -d {eps}",
                description=f"MCCE4 Step 3 (energies, ε={eps})", cwd=work_dir)
        
        # Parse head3.lst for dsolv values
        head3_path = os.path.join(work_dir, "head3.lst")
        dsolv_values = parse_head3_dsolv(head3_path, lig_id, conformers)
        
        if not dsolv_values:
            logging.error(f"  Could not parse dsolv values from {head3_path}")
            continue
        
        # Report extracted values
        for conf_type, dsolv in dsolv_values.items():
            logging.info(f"  📊 {conf_type}: dsolv = {dsolv:.3f}")
        
        # Update CONFORMER rxn values in the .ftpl
        update_ftpl_rxn(ftpl_path, dsolv_values, rxn_key)
        
        # Validation: re-run step3 to check dsolv → 0.000
        logging.info(f"\n  🔄 Validation run (expecting dsolv ≈ 0.000)...")
        run_cmd(f"step3.py -d {eps}",
                description=f"Validation Step 3 (ε={eps})", cwd=work_dir)
        
        # Check validation results
        dsolv_check = parse_head3_dsolv(head3_path, lig_id, conformers)
        all_ok = True
        for conf_type, dsolv in dsolv_check.items():
            if abs(dsolv) > 0.01:
                logging.warning(f"  ⚠ {conf_type}: dsolv = {dsolv:.3f} (expected ≈ 0.000)")
                all_ok = False
            else:
                logging.info(f"  ✓ {conf_type}: dsolv = {dsolv:.3f} ✓")
        
        if all_ok:
            logging.info(f"  🎉 RXN calibration for ε={eps} successful!")
        else:
            logging.warning(f"  ⚠ RXN calibration may need manual adjustment.")


def parse_head3_dsolv(head3_path: str, lig_id: str, conformers: list) -> dict:
    """Parse head3.lst to extract the most negative dsolv per conformer type.
    
    head3.lst format:
    iConf CONFORMER      FL occ   crg  Em0 pKa0 ne nH  vdw0  vdw1  tors  epol  dsolv  extra  history
    00001 EMH01_0000_001  f 0.00 ...                                      -8.085  ...
    
    The dsolv column is at a fixed position. We find the most negative value
    for each conformer TYPE (e.g., EMH01, EMH+1).
    
    Args:
        head3_path: Path to head3.lst file.
        conformers: List of conformer state strings (e.g., ['01', '+1']).
        lig_id:     3-letter ligand code.
    
    Returns:
        dict: Mapping of conformer type (e.g., 'EMH01') → most negative dsolv.
    """
    if not os.path.exists(head3_path):
        logging.error(f"  head3.lst not found at {head3_path}")
        return {}
    
    with open(head3_path, "r") as f:
        lines = f.readlines()
    
    # Build conformer type prefixes: EMH01, EMH+1, etc.
    conf_types = [f"{lig_id}{c}" for c in conformers]
    
    # Find the dsolv column index from the header
    header_line = None
    for line in lines:
        if "dsolv" in line.lower() and "iConf" in line:
            header_line = line
            break
    
    if not header_line:
        logging.error("  Could not find header line in head3.lst")
        return {}
    
    # Dynamically find the dsolv column index from the header
    header_parts = header_line.split()
    try:
        dsolv_col = header_parts.index("dsolv")
    except ValueError:
        logging.error("  Could not find 'dsolv' column in head3.lst header")
        logging.debug(f"  Header columns: {header_parts}")
        return {}
    
    logging.debug(f"  dsolv is at column index {dsolv_col} in head3.lst")
    
    # Parse data lines for our ligand conformers
    dsolv_values = {ct: 0.0 for ct in conf_types}
    
    for line in lines:
        for ct in conf_types:
            if ct in line and not line.strip().startswith("iConf"):
                # Parse the fixed-width format of head3.lst
                parts = line.split()
                if len(parts) > dsolv_col:
                    try:
                        dsolv = float(parts[dsolv_col])
                        # Keep the most negative value
                        if dsolv < dsolv_values[ct]:
                            dsolv_values[ct] = dsolv
                    except (ValueError, IndexError):
                        continue
    
    return dsolv_values


def update_ftpl_rxn(ftpl_path: str, dsolv_values: dict, rxn_key: str):
    """Update CONFORMER lines in .ftpl with calibrated rxn values.
    
    Finds lines like:
      CONFORMER, EMH01:  Em0=0.0, pKa0=0.00, ne=0, nH=0, rxn02= 0, rxn04= 0, rxn08= 0
    
    And updates the appropriate rxn field with the dsolv value.
    
    Args:
        ftpl_path:    Path to the .ftpl file.
        dsolv_values: dict of conformer type → dsolv value.
        rxn_key:      Which rxn to update (e.g., 'rxn04').
    """
    logging.info(f"  ✏️  Updating {rxn_key} in {ftpl_path}")
    
    with open(ftpl_path, "r") as f:
        content = f.read()
    
    for conf_type, dsolv in dsolv_values.items():
        # Pattern to find the rxn value for this conformer
        # e.g., rxn04= 0 or rxn04=  -8.085
        pattern = re.compile(
            rf"(CONFORMER,\s*{re.escape(conf_type)}:.*?{rxn_key}=\s*)([-\d.]+)"
        )
        
        match = pattern.search(content)
        if match:
            old_val = match.group(2)
            content = pattern.sub(
                rf"\g<1>{dsolv:8.3f}",
                content
            )
            logging.info(f"    {conf_type}: {rxn_key} = {old_val} → {dsolv:.3f}")
        else:
            logging.warning(f"    Could not find {rxn_key} for {conf_type} in .ftpl")
    
    with open(ftpl_path, "w") as f:
        f.write(content)


# ──────────────────────────────────────────────────────────────────────────────
# Charge summary report
# ──────────────────────────────────────────────────────────────────────────────
def print_charge_summary(charges: dict, lig_id: str, method: str):
    """Print a formatted summary table of the assigned charges.
    
    Args:
        charges: dict of atom name → charge.
        lig_id:  3-letter ligand code.
        method:  Charge method used.
    """
    total = sum(charges.values())
    
    logging.info(f"\n{'─'*60}")
    logging.info(f"  Charge Summary for {lig_id} (method: {method})")
    logging.info(f"{'─'*60}")
    logging.info(f"  {'Atom':<8} {'Charge':>8}")
    logging.info(f"  {'────':<8} {'──────':>8}")
    
    # Group by element for readability
    for name in sorted(charges.keys()):
        logging.info(f"  {name!r:<8} {charges[name]:>8.3f}")
    
    logging.info(f"  {'────':<8} {'──────':>8}")
    logging.info(f"  {'TOTAL':<8} {total:>8.3f}")
    logging.info(f"  {'Atoms':<8} {len(charges):>8d}")
    logging.info(f"{'─'*60}\n")


# ──────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ──────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="🧬 Auto-generate MCCE4 topology files (.ftpl)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s EMH -c 01 +1                           # Full auto (mmff94, calibrates rxn02/04/08)
  %(prog)s EMH -c 01 +1 -m am1bcc                 # Use AM1-BCC charges  
  %(prog)s EMH -c 01 +1 --pdb EMH.pdb             # Use local PDB file
  %(prog)s EMH -c 01 +1 --charges-file q.txt      # Use pre-computed charges
  %(prog)s EMH -c 01 +1 -d 4                      # Calibrate only rxn04
  %(prog)s EMH -c 01 +1 --dry-run                 # Generate ftpl, skip calibration
  %(prog)s EMH -c 01 -1 -m antechamber --nc -1    # antechamber with net charge -1

Charge methods:
  mmff94         OpenEye MMFF94 (default, reliable, works on most molecules)
  am1bcc         OpenEye AM1-BCC (higher quality, small molecules only)
  am1bccelf10    OpenEye AM1-BCC + ELF10 conformer selection
  am1bccnosymspt OpenEye AM1-BCC no symmetry, single point
  amber          OpenEye Amber charges
  amberff94      OpenEye Amber FF94 charges  
  antechamber    AmberTools antechamber AM1-BCC (free, no OpenEye needed)
  file           Load from a user-supplied charge file

Workflow:
  1. Download/load ideal ligand PDB
  2. Generate .ftpl template via pdb2ftpl.py
  3. Compute partial atomic charges
  4. Fill 'to_be_filled' entries in .ftpl
  5. Run MCCE4 steps 1-2 once
  6. Loop step 3 for ε=2, 4, 8 → parse dsolv → update rxn02, rxn04, rxn08
  7. Validate each rxn calibration (dsolv ≈ 0.000)
        """
    )
    
    # Positional
    parser.add_argument("ligand", type=str,
                        help="3-letter ligand code (e.g., EMH, ATP, HEM)")
    
    # Required
    parser.add_argument("-c", "--conformers", nargs="+", required=True,
                        help="Conformer states (e.g., 01 +1 -1)")
    
    # Charge options
    parser.add_argument("-m", "--method", type=str, default=DEFAULT_CHARGE_METHOD,
                        choices=SUPPORTED_CHARGE_METHODS,
                        help=f"Charge calculation method (default: {DEFAULT_CHARGE_METHOD})")
    parser.add_argument("--charges-file", type=str, default=None,
                        help="Path to pre-computed charges file (overrides -m)")
    parser.add_argument("--nc", type=int, default=0,
                        help="Net formal charge of molecule (for antechamber, default: 0)")
    
    # Input options
    parser.add_argument("--pdb", type=str, default=None,
                        help="Path to ligand PDB file (skip download)")
    parser.add_argument("--cif", type=str, default=None,
                        help="Path to ligand CIF file (will convert to PDB)")
    parser.add_argument("--skip-download", action="store_true",
                        help="Skip downloading PDB, use existing <LIG>.pdb")
    
    # Calibration options
    parser.add_argument("-d", "--dielectric", nargs="+", type=int, default=[2, 4, 8],
                        help="Dielectric constant(s) for RXN calibration (default: 2 4 8)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Generate .ftpl with charges but skip RXN calibration")
    parser.add_argument("--work-dir", type=str, default=".",
                        help="MCCE4 working directory (default: current directory)")
    
    # Output options
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output .ftpl filename (default: <LIG>.ftpl)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable debug-level logging")
    
    args = parser.parse_args()
    
    # ── Setup ──
    lig_id = args.ligand.upper()
    ftpl_path = args.output or f"{lig_id}.ftpl"
    pdb_path = args.pdb or f"{lig_id}.pdb"
    work_dir = os.path.abspath(args.work_dir)
    
    logger = setup_logging(
        log_file=f"auto_ftpl_{lig_id}.log",
        verbose=args.verbose
    )
    
    logging.info(f"{'='*60}")
    logging.info(f"  🧬 MCCE4 Topology File Agent — auto_ftpl.py")
    logging.info(f"{'='*60}")
    logging.info(f"  Ligand:      {lig_id}")
    logging.info(f"  Conformers:  {args.conformers}")
    logging.info(f"  Method:      {args.method}")
    logging.info(f"  Dielectric:  {args.dielectric}")
    logging.info(f"  Output:      {ftpl_path}")
    logging.info(f"  Work dir:    {work_dir}")
    logging.info(f"  Started:     {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"{'='*60}\n")
    
    # ── Check dependencies ──
    charge_method = "file" if args.charges_file else args.method
    if not args.dry_run:
        check_dependencies(charge_method)
    
    # ── Step 1: Obtain ligand PDB ──
    if args.cif:
        convert_cif_to_pdb(args.cif, pdb_path)
    elif args.pdb:
        pdb_path = args.pdb
        if not os.path.exists(pdb_path):
            logging.error(f"PDB file not found: {pdb_path}")
            sys.exit(1)
    elif not args.skip_download:
        download_ligand_pdb(lig_id, pdb_path)
    else:
        if not os.path.exists(pdb_path):
            logging.error(f"PDB file not found: {pdb_path}. "
                         f"Remove --skip-download to auto-download.")
            sys.exit(1)
    
    # ── Step 2: Generate .ftpl template ──
    unfilled_count = generate_ftpl_template(lig_id, pdb_path, args.conformers,
                                             ftpl_path)
    
    # ── Step 3: Generate or load charges ──
    if args.charges_file:
        charges = load_charges_from_file(args.charges_file)
    elif charge_method == "antechamber":
        charges = generate_charges_antechamber(pdb_path, lig_id, args.nc)
    else:
        # OpenEye methods
        charges = generate_charges_openeye(pdb_path, charge_method, lig_id)
    
    print_charge_summary(charges, lig_id, charge_method)
    
    # ── Step 4: Fill charges into .ftpl ──
    unfilled = fill_ftpl_charges(ftpl_path, charges, args.conformers, lig_id)
    
    if unfilled > 0:
        logging.warning(f"\n  ⚠  {unfilled} charge entries could not be auto-filled.")
        logging.warning(f"  Please fill them manually in: {ftpl_path}")
        logging.warning(f"  Then re-run with: auto_ftpl.py {lig_id} -c {' '.join(args.conformers)} "
                       f"--skip-download --dry-run")
    
    # ── Steps 5-7: RXN Calibration ──
    if args.dry_run:
        logging.info(f"\n  ⏩ Dry run — skipping RXN calibration.")
        logging.info(f"  To calibrate, run without --dry-run")
    elif unfilled > 0:
        logging.warning(f"\n  ⏩ Skipping RXN calibration (unfilled charges remain)")
    else:
        run_rxn_calibration(lig_id, ftpl_path, pdb_path, args.conformers,
                            args.dielectric, work_dir)
    
    # ── Done ──
    logging.info(f"\n{'='*60}")
    logging.info(f"  ✅ Topology file ready: {ftpl_path}")
    logging.info(f"  📋 Full log: auto_ftpl_{lig_id}.log")
    if not args.dry_run and unfilled == 0:
        logging.info(f"  📁 Copy to user_param/ for your MCCE4 simulation:")
        logging.info(f"     cp {ftpl_path} <your_mcce_run>/user_param/")
    logging.info(f"{'='*60}")


if __name__ == "__main__":
    main()
