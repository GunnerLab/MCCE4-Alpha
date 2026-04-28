"""
Analysis module for MCCE4 output parsing.

Parses pK.out, sum_crg.out, and other Step 4 outputs into structured
data suitable for Plotly visualization in the frontend.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Optional


def parse_pk_out(path: str = "pK.out") -> dict:
    """
    Parse pK.out which has TWO sections:

    Section 1 — pKa summary:
        pH             pKa/Em  n(slope) 1000*chi2
        NTR+A0001_        7.948     1.012     0.001
        ARG+A0001_        >14.0
        ASP-A0003_        3.008     1.004     0.005

    Section 2 — titration curves:
         ph              0.0   1.0   2.0   3.0  ...  14.0
        NTR+A0001_      1.00  1.00  1.00  ...   0.00
        ASP-A0003_      0.00 -0.01 -0.09  ...  -1.00

    Returns:
        {
            "residues": [{"name": ..., "pka": ..., "n_slope": ..., "chi2": ...}, ...],
            "ph_values": [0.0, 1.0, 2.0, ...],
            "titration_curves": {"NTR+A0001_": [1.00, 1.00, ...], ...}
        }
    """
    if not os.path.isfile(path):
        return {"error": f"File not found: {path}"}

    with open(path, "r") as f:
        content = f.read()

    if not content.strip():
        return {"error": "pK.out is empty"}

    lines = content.strip().split("\n")

    residues = []
    ph_values = []
    titration_curves = {}

    # ── Find where section 2 starts ──
    # Section 2 header: a line starting with "ph" followed by float values
    # e.g. " ph              0.0   1.0   2.0   3.0 ..."
    # Section 1 header starts with "pH" followed by text labels like "pKa/Em"
    section2_start = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        # Match: starts with "ph" (case-insensitive), then a float number
        m = re.match(r'^ph\s+([\d.-]+)', stripped, re.IGNORECASE)
        if m:
            section2_start = i
            break

    # ── Parse Section 1: pKa summary ──
    section1_end = section2_start if section2_start is not None else len(lines)
    for line in lines[:section1_end]:
        stripped = line.strip()
        if not stripped:
            continue
        # Skip header line ("pH   pKa/Em   n(slope)  1000*chi2")
        if stripped.lower().startswith("ph"):
            continue

        parts = stripped.split()
        if len(parts) < 2:
            continue

        res_name = parts[0]

        # Parse pKa — handle ">14.0", "<0.0", "nan", etc.
        pka_str = parts[1]
        pka = None
        try:
            cleaned = pka_str.lstrip("><")
            pka = float(cleaned)
        except ValueError:
            pka = None

        # Parse n(slope)
        n_slope = None
        if len(parts) > 2:
            try:
                n_slope = float(parts[2])
            except ValueError:
                pass

        # Parse 1000*chi2
        chi2 = None
        if len(parts) > 3:
            try:
                chi2 = float(parts[3])
            except ValueError:
                pass

        residues.append({
            "name": res_name,
            "pka": pka,
            "pka_display": pka_str,
            "n_slope": n_slope,
            "chi2": chi2,
        })

    # ── Parse Section 2: titration curves ──
    if section2_start is not None:
        # Parse pH values from header
        header_line = lines[section2_start].strip()
        header_parts = header_line.split()
        for part in header_parts[1:]:  # Skip "ph" label
            try:
                ph_values.append(float(part))
            except ValueError:
                break

        # Parse residue curves
        for line in lines[section2_start + 1:]:
            stripped = line.strip()
            if not stripped:
                continue

            parts = stripped.split()
            if len(parts) < 2:
                continue

            res_name = parts[0]
            curve = []
            for val_str in parts[1:]:
                try:
                    curve.append(float(val_str))
                except ValueError:
                    curve.append(None)

            if curve:
                titration_curves[res_name] = curve

    return {
        "residues": residues,
        "ph_values": ph_values,
        "titration_curves": titration_curves,
    }


def parse_sum_crg(path: str = "sum_crg.out") -> dict:
    """
    Parse sum_crg.out for net charge vs pH/Eh data.

    Format:
         ph         0.0    1.0    2.0  ...  14.0
        NTR+A0001_  1.00   1.00   1.00 ...   0.00
        ...
        Total       6.00   6.00   5.91 ...  -9.05

    Returns:
        {
            "ph_values": [...],
            "residues": {"NTR+A0001_": [crg_at_ph0, ...], ...},
            "total_charge": [total_at_ph0, ...],
        }
    """
    if not os.path.isfile(path):
        return {"error": f"File not found: {path}"}

    with open(path, "r") as f:
        lines = f.readlines()

    if not lines:
        return {"error": "sum_crg.out is empty"}

    ph_values = []
    residue_charges = {}
    total_charge = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        parts = stripped.split()

        # Detect header line: starts with "ph" followed by numbers
        if parts[0].lower() == "ph" and len(parts) > 1:
            for p in parts[1:]:
                try:
                    ph_values.append(float(p))
                except ValueError:
                    break
            continue

        if len(parts) < 2:
            continue

        res_name = parts[0]

        # Skip separator lines (e.g. "------...")
        if res_name.startswith("---"):
            continue

        charges = []
        for p in parts[1:]:
            try:
                charges.append(float(p))
            except ValueError:
                charges.append(0.0)

        # Detect summary rows
        name_lower = res_name.lower().replace("_", "")
        if name_lower in ("total", "net", "netcharge", "sum", "total:"):
            total_charge = charges
        elif name_lower in ("protons", "electrons"):
            # Store but don't add to residues
            pass
        else:
            residue_charges[res_name] = charges

    return {
        "ph_values": ph_values,
        "residues": residue_charges,
        "total_charge": total_charge,
    }


def parse_head3_lst(path: str = "head3.lst") -> dict:
    """Parse head3.lst for conformer information."""
    if not os.path.isfile(path):
        return {"error": f"File not found: {path}"}

    conformers = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("iConf"):
                continue

            parts = line.split()
            if len(parts) < 5:
                continue

            try:
                conformers.append({
                    "index": int(parts[0]) if parts[0].isdigit() else 0,
                    "name": parts[1],
                    "occ": float(parts[2]) if len(parts) > 2 else 0.0,
                    "crg": float(parts[3]) if len(parts) > 3 else 0.0,
                    "em": float(parts[4]) if len(parts) > 4 else 0.0,
                })
            except (ValueError, IndexError):
                continue

    return {"conformers": conformers}


def parse_dipole_csv(path: str = "dipole_ph_scan.csv") -> dict:
    """
    Parse dipole_ph_scan.csv produced by ms_dipole.py.

    Returns:
        {
            "ph_values": [0.0, 1.0, ...],
            "backbone_D": [34.21, ...],
            "ionizable_D": [154.79, ...],
            "full_D": [162.34, ...],
            "net_charge": [10.986, ...],
            "arrows": {
                "backbone":  [{"neg": [x,y,z], "pos": [x,y,z]}, ...],
                "ionizable": [{"neg": [x,y,z], "pos": [x,y,z]}, ...],
                "full":      [{"neg": [x,y,z], "pos": [x,y,z]}, ...],
            },
            "vectors": {
                "backbone":  [[x,y,z], ...],
                "ionizable": [[x,y,z], ...],
                "full":      [[x,y,z], ...],
            },
            "quadrupole_eigenvalues": [[e1,e2,e3], ...],
        }
    """
    import csv

    if not os.path.isfile(path):
        return {"error": f"File not found: {path}"}

    with open(path, "r") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        return {"error": "dipole_ph_scan.csv is empty"}

    ph_values = []
    backbone_D = []
    ionizable_D = []
    full_D = []
    net_charge = []
    arrows_bb = []
    arrows_ion = []
    arrows_full = []
    vectors_bb = []
    vectors_ion = []
    vectors_full = []
    quad_eig = []

    for row in rows:
        ph_values.append(float(row["pH"]))
        backbone_D.append(float(row["backbone_D"]))
        ionizable_D.append(float(row["ionizable_D"]))
        full_D.append(float(row["full_D"]))
        net_charge.append(float(row["net_charge"]))

        arrows_bb.append({
            "neg": [float(row["bb_neg_x"]), float(row["bb_neg_y"]), float(row["bb_neg_z"])],
            "pos": [float(row["bb_pos_x"]), float(row["bb_pos_y"]), float(row["bb_pos_z"])],
        })
        arrows_ion.append({
            "neg": [float(row["ion_neg_x"]), float(row["ion_neg_y"]), float(row["ion_neg_z"])],
            "pos": [float(row["ion_pos_x"]), float(row["ion_pos_y"]), float(row["ion_pos_z"])],
        })
        arrows_full.append({
            "neg": [float(row["full_neg_x"]), float(row["full_neg_y"]), float(row["full_neg_z"])],
            "pos": [float(row["full_pos_x"]), float(row["full_pos_y"]), float(row["full_pos_z"])],
        })
        vectors_bb.append([float(row["bb_x"]), float(row["bb_y"]), float(row["bb_z"])])
        vectors_ion.append([float(row["ion_x"]), float(row["ion_y"]), float(row["ion_z"])])
        vectors_full.append([float(row["full_x"]), float(row["full_y"]), float(row["full_z"])])
        quad_eig.append([float(row["Q_eig1"]), float(row["Q_eig2"]), float(row["Q_eig3"])])

    return {
        "ph_values": ph_values,
        "backbone_D": backbone_D,
        "ionizable_D": ionizable_D,
        "full_D": full_D,
        "net_charge": net_charge,
        "arrows": {
            "backbone": arrows_bb,
            "ionizable": arrows_ion,
            "full": arrows_full,
        },
        "vectors": {
            "backbone": vectors_bb,
            "ionizable": vectors_ion,
            "full": vectors_full,
        },
        "quadrupole_eigenvalues": quad_eig,
    }


def can_run_dipole() -> dict:
    """Check whether ms_dipole can be run (required files exist, CSV does not)."""
    has_csv = os.path.isfile("dipole_ph_scan.csv")
    has_step2 = os.path.isfile("step2_out.pdb")
    has_head3 = os.path.isfile("head3.lst")
    has_fort38 = os.path.isfile("fort.38")
    ready = has_step2 and has_head3 and has_fort38
    missing = []
    if not has_step2:
        missing.append("step2_out.pdb")
    if not has_head3:
        missing.append("head3.lst")
    if not has_fort38:
        missing.append("fort.38")
    return {
        "has_csv": has_csv,
        "can_run": ready,
        "missing": missing,
    }


def run_dipole(log_callback=None) -> dict:
    """
    Run the ms_dipole ensemble computation in the current directory.

    Imports from mcce_dipole and produces dipole_ph_scan.csv.
    Returns parsed dipole data on success.
    """
    import sys
    import numpy as np

    # Ensure MCCE_bin is importable
    mcce_bin_dir = str(Path(__file__).resolve().parent.parent.parent / "MCCE_bin")
    if mcce_bin_dir not in sys.path:
        sys.path.insert(0, mcce_bin_dir)

    from mcce_dipole.parsers import parse_step2_pdb, parse_head3lst, parse_fort38
    from mcce_dipole.compute import compute_from_ensemble
    from mcce_dipole.visualize import generate_ph_scan_csv

    def log(msg):
        if log_callback:
            log_callback(msg)

    # Check files
    for name in ("step2_out.pdb", "head3.lst", "fort.38"):
        if not os.path.isfile(name):
            return {"error": f"Required file not found: {name}"}

    log("Dipole: parsing step2_out.pdb...")
    conformers, all_atoms = parse_step2_pdb("step2_out.pdb")
    log(f"Dipole: {len(conformers)} conformers, {len(all_atoms)} atoms")

    log("Dipole: parsing head3.lst...")
    head3_data = parse_head3lst("head3.lst")
    log(f"Dipole: {len(head3_data)} entries")

    log("Dipole: parsing fort.38...")
    ph_values, conf_ids, occupancies = parse_fort38("fort.38")
    log(f"Dipole: {len(ph_values)} pH values, {len(conf_ids)} conformers")

    log("Dipole: computing pH-dependent dipole moments...")
    results = compute_from_ensemble(
        conformers, head3_data,
        ph_values, conf_ids, occupancies
    )

    csv_path = "dipole_ph_scan.csv"
    generate_ph_scan_csv(results, csv_path)
    log(f"Dipole: wrote {csv_path}")

    return parse_dipole_csv(csv_path)


def can_run_pdb_dipole() -> dict:
    """Check whether PDB standard protonation dipole can be computed."""
    has_step1 = os.path.isfile("step1_out.pdb")
    has_param = os.path.isdir("param") and os.path.isfile(os.path.join("param", "mcce.tpl"))
    missing = []
    if not has_step1:
        missing.append("step1_out.pdb")
    if not has_param:
        missing.append("param/mcce.tpl")
    return {
        "can_run": has_step1 and has_param,
        "missing": missing,
    }


def run_pdb_dipole(log_callback=None) -> dict:
    """
    Compute PDB standard protonation dipole from step1_out.pdb + param/mcce.tpl.

    Returns a dict with single-state dipole data formatted for the GUI,
    including arrow endpoints and quadrupole data.
    """
    import sys
    import numpy as np

    mcce_bin_dir = str(Path(__file__).resolve().parent.parent.parent / "MCCE_bin")
    if mcce_bin_dir not in sys.path:
        sys.path.insert(0, mcce_bin_dir)

    from mcce_dipole.parsers import parse_step2_pdb, parse_tpl_charges
    from mcce_dipole.compute import compute_from_step1

    def log(msg):
        if log_callback:
            log_callback(msg)

    if not os.path.isfile("step1_out.pdb"):
        return {"error": "step1_out.pdb not found"}
    param_dir = "param"
    if not os.path.isfile(os.path.join(param_dir, "mcce.tpl")):
        return {"error": "param/mcce.tpl not found"}

    log("PDB dipole: parsing step1_out.pdb...")
    conformers, all_atoms = parse_step2_pdb("step1_out.pdb")
    log(f"PDB dipole: {len(conformers)} conformers, {len(all_atoms)} atoms")

    log("PDB dipole: loading topology charges...")
    tpl_charges = parse_tpl_charges(param_dir)
    log(f"PDB dipole: {len(tpl_charges)} charge entries")

    log("PDB dipole: computing standard protonation dipole...")
    results = compute_from_step1(conformers, all_atoms, tpl_charges)

    center = results["center"].tolist()
    arrow_scale = 0.1

    data = {
        "mode": "pdb",
        "center": center,
        "net_charge": float(results["net_charge"]),
        "arrows": {},
        "vectors": {},
        "quadrupole_eigenvalues": results["quadrupole_eigenvalues"].tolist(),
        "quadrupole_eigenvectors": results["quadrupole_eigenvectors"].T.tolist(),
    }

    for key, label in [("backbone_dipole", "backbone"),
                       ("ionizable_dipole", "ionizable"),
                       ("full_dipole", "full")]:
        mu = results[key]
        mag = float(np.linalg.norm(mu))
        data[f"{label}_D"] = mag

        if mag > 0.1:
            direction = mu / mag
            neg = [center[i] - direction[i] * mag * arrow_scale for i in range(3)]
            pos = [center[i] + direction[i] * mag * arrow_scale for i in range(3)]
        else:
            neg = center[:]
            pos = center[:]
        data["arrows"][label] = {"neg": neg, "pos": pos}
        data["vectors"][label] = mu.tolist()

    log("PDB dipole: done.")
    return data


def get_available_outputs() -> dict:
    """Check which MCCE4 output files exist in the current directory."""
    files = {
        "pK.out": os.path.isfile("pK.out"),
        "sum_crg.out": os.path.isfile("sum_crg.out"),
        "head3.lst": os.path.isfile("head3.lst"),
        "step1_out.pdb": os.path.isfile("step1_out.pdb"),
        "step2_out.pdb": os.path.isfile("step2_out.pdb"),
        "energies": os.path.isdir("energies"),
        "run.prm": os.path.isfile("run.prm"),
        "run.log": os.path.isfile("run.log"),
        "dipole_ph_scan.csv": os.path.isfile("dipole_ph_scan.csv"),
        "fort.38": os.path.isfile("fort.38"),
    }
    return {"files": files, "workdir": os.getcwd()}
