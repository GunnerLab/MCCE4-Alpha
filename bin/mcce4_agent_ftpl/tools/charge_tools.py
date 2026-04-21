"""
Partial atomic charge generation tools.
"""

import os
import re
import logging
import shutil
import subprocess
import tempfile


def to_mcce_atom_name(name: str) -> str:
    """Convert atom name to MCCE 4-char padded format."""
    name = name.strip()
    if len(name) >= 4: return name[:4]
    elif len(name) == 3: return f" {name}"
    elif len(name) == 2: return f" {name} "
    elif len(name) == 1: return f" {name}  "
    return name.ljust(4)


def generate_charges(pdb_path: str, method: str, lig_id: str) -> dict:
    """Generate partial atomic charges. Returns MCCE atom name -> charge."""
    if method == "antechamber":
        return _charges_antechamber(pdb_path, lig_id)
    return _charges_openeye(pdb_path, method, lig_id)


def _charges_openeye(pdb_path: str, method: str, lig_id: str) -> dict:
    """Generate charges via OpenEye QuacPac TK."""
    logging.info(f"⚡ Generating charges via OpenEye (method: {method})")

    # Use absolute path so subprocess can find the file regardless of cwd
    pdb_abs = os.path.abspath(pdb_path)
    tmpdir = tempfile.mkdtemp(prefix="oe_charges_")
    try:
        result = subprocess.run(
            f"oe_assigncharges_QuacpakTK.py -method {method} -in {pdb_abs}",
            shell=True, capture_output=True, text=True, cwd=tmpdir
        )
        if result.returncode != 0:
            logging.error(f"  OpenEye failed: {result.stderr}")
            return {}
    except Exception as e:
        logging.error(f"  Could not run OpenEye: {e}")
        return {}
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    charges = _parse_openeye_output(result.stdout)

    # Validate for all-zero failure mode
    non_zero = sum(1 for v in charges.values() if abs(v) > 1e-6)
    if non_zero == 0 and method != "mmff94":
        logging.warning(f"  ⚠ All zeros from {method} — falling back to mmff94")
        return _charges_openeye(pdb_path, "mmff94", lig_id)

    if charges:
        total = sum(charges.values())
        logging.info(f"  ✓ {len(charges)} charges (total: {total:.3f})")
    return charges


def _parse_openeye_output(output: str) -> dict:
    """Parse OpenEye output, filtering water duplicates."""
    charges = {}
    pattern = re.compile(
        r"Atom Name:\s+(\S+)\s+\|\s+Symbol:\s+(\S+)\s+\|\s+Charge:\s+([-\d.]+)"
    )
    atom_counts = {}
    raw = []

    for line in output.splitlines():
        m = pattern.search(line)
        if m:
            name, sym, charge = m.group(1), m.group(2), float(m.group(3))
            raw.append((name, charge))
            atom_counts[name] = atom_counts.get(name, 0) + 1

    # Filter water atoms (appearing > 5 times)
    water_atoms = {n for n, c in atom_counts.items() if c > 5}
    seen = set()
    for name, charge in raw:
        if name in water_atoms:
            continue
        mcce = to_mcce_atom_name(name)
        if mcce not in seen:
            charges[mcce] = charge
            seen.add(mcce)
    return charges


def _charges_antechamber(pdb_path: str, lig_id: str, nc: int = 0) -> dict:
    """Generate charges via AmberTools antechamber."""
    logging.info(f"⚡ Generating charges via antechamber (nc={nc})")
    mol2 = f"{lig_id}_antechamber.mol2"

    try:
        result = subprocess.run(
            f"antechamber -i {pdb_path} -fi pdb -o {mol2} -fo mol2 "
            f"-c bcc -s 2 -nc {nc} -at gaff",
            shell=True, capture_output=True, text=True
        )
        if result.returncode != 0 or not os.path.exists(mol2):
            logging.error(f"  antechamber failed")
            return {}
    except Exception as e:
        logging.error(f"  Could not run antechamber: {e}")
        return {}

    charges = {}
    in_atoms = False
    with open(mol2) as f:
        for line in f:
            l = line.strip()
            if l.startswith("@<TRIPOS>ATOM"):
                in_atoms = True; continue
            elif l.startswith("@<TRIPOS>"):
                in_atoms = False; continue
            if in_atoms and l:
                p = l.split()
                if len(p) >= 9:
                    charges[to_mcce_atom_name(p[1])] = float(p[8])

    if charges:
        total = sum(charges.values())
        logging.info(f"  ✓ {len(charges)} charges (total: {total:.3f})")
    return charges


# ═════════════════════════════════════════════════════════════════════════════
# v4: OEChem-based ionizable site detection (fallback when RDKit fails)
# ═════════════════════════════════════════════════════════════════════════════

def get_ionizable_sites_oe(pdb_path: str) -> list:
    """Detect ionizable N/O/S sites in a PDB using OEChem.

    More robust than RDKit for unusual molecules (sp nitrile N, fused rings,
    unusual valence, etc.).  Returns the same list-of-dicts format as
    rdkit_tools.get_ionizable_sites:
      [{"idx", "name", "symbol", "type", "current_hs", "formal_charge"}]
    """
    logging.info(f"  Detecting ionizable sites via OEChem: {pdb_path}")
    try:
        from openeye import oechem
    except ImportError:
        logging.warning("  OEChem not available for ionizable site detection")
        return []

    mol = oechem.OEMol()
    ifs = oechem.oemolistream(pdb_path)
    if not oechem.OEReadMolecule(ifs, mol):
        logging.warning(f"  OEChem: cannot load {pdb_path}")
        ifs.close()
        return []
    ifs.close()

    oechem.OEAssignAromaticFlags(mol)
    oechem.OEAssignHybridization(mol)

    sites = []
    for atom in mol.GetAtoms():
        sym = oechem.OEGetAtomicSymbol(atom.GetAtomicNum())
        if sym not in ("N", "O", "S"):
            continue

        # atom.GetName() gives the PDB atom name (e.g., "N25", " N9 ")
        name = atom.GetName().strip() if atom.GetName() else f"{sym}{atom.GetIdx()}"

        n_explicit_h = sum(
            1 for bond in atom.GetBonds()
            if bond.GetNbr(atom).GetAtomicNum() == 1
        )
        hybrid = atom.GetHyb()
        fc = atom.GetFormalCharge()

        if sym == "N":
            # Skip sp N (nitrile, isonitrile, azide terminal, etc.)
            if hybrid == oechem.OEHybridization_sp:
                continue
            # Use explicit H count as ground truth for protonation state.
            # OEChem's fc can be unreliable for PDB HETATM atoms (it adds
            # implicit H based on valence rules that don't match the actual
            # neutral-PDB protonation).  A tertiary amine (no explicit H) is
            # protonatable regardless of OEChem's perceived formal charge.
            if n_explicit_h > 0 or fc == 0 or (n_explicit_h == 0 and fc == 1):
                sites.append({
                    "idx": atom.GetIdx(), "name": name, "symbol": sym,
                    "type": "protonatable_N",
                    "current_hs": n_explicit_h, "formal_charge": fc,
                })
        elif sym == "O" and n_explicit_h > 0:
            sites.append({
                "idx": atom.GetIdx(), "name": name, "symbol": sym,
                "type": "deprotonatable_OH",
                "current_hs": n_explicit_h, "formal_charge": fc,
            })
        elif sym == "S" and n_explicit_h > 0:
            sites.append({
                "idx": atom.GetIdx(), "name": name, "symbol": sym,
                "type": "deprotonatable_SH",
                "current_hs": n_explicit_h, "formal_charge": fc,
            })

    logging.info(f"  OEChem: {len(sites)} ionizable site(s): "
                 f"{[s['name'] for s in sites]}")
    return sites


# ═════════════════════════════════════════════════════════════════════════════
# v3: Per-State Charge Generation
# ═════════════════════════════════════════════════════════════════════════════

def generate_per_state_charges(state_pdbs, method, lig_id):
    """Generate charges for each protonation state PDB separately.

    Args:
        state_pdbs: {label: pdb_path}
        method: Charge method (mmff94, am1bcc, etc.)
        lig_id: Ligand identifier

    Returns:
        {label: {atom_name: charge}}
    """
    per_state = {}
    for label, pdb_path in state_pdbs.items():
        logging.info(f"  ⚡ Charges for state {label}: {pdb_path}")
        charges = generate_charges(pdb_path, method, lig_id)
        if charges:
            per_state[label] = charges
            total = sum(charges.values())
            logging.info(f"    → {len(charges)} atoms (sum={total:.3f})")
        else:
            logging.warning(f"    → No charges for state {label}")
            per_state[label] = {}
    return per_state
