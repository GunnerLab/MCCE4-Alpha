#!/usr/bin/env python
"""
ms_dipole.py - MCCE4 Dipole & Quadrupole Moment Analysis

Computes four electrostatic moments from the MCCE4 microstate ensemble,
Boltzmann-weighted over conformer occupancies at each pH:

  1. Backbone Dipole Moment
     Sum of q*r over backbone atoms (N, CA, C, O, H, HA) only.
     Reflects the intrinsic polarity from aligned peptide bonds — largely
     pH-independent since backbone charges are fixed.

  2. Ionizable Protein Dipole Moment
     Sum of q*r over ionizable residue atoms (ASP, GLU, LYS, ARG, HIS,
     CYS, TYR, NTR, CTR). This is the most pH-sensitive component —
     it changes as titratable groups gain/lose protons.

  3. Full Protein Dipole Moment
     Sum of q*r over ALL atoms (backbone + sidechain + ions). The total
     charge asymmetry of the protein at each pH.

  4. Quadrupole Moment Tensor
     The traceless symmetric 3x3 tensor Q_ab = Sum qi(3*ri_a*ri_b - |ri|^2*d_ab).
     Describes the shape of the charge distribution beyond the dipole.
     Reported as eigenvalues (principal moments) in e*A^2.

All quantities use partial charges embedded in step2_out.pdb and conformer
occupancies from fort.38 Monte Carlo sampling.

Units: Dipole in Debye (1 e*A = 4.803 D), Quadrupole in e*A^2.

Usage:
  ms_dipole.py                        # ensemble from cwd
  ms_dipole.py --pdb protein.pdb      # PDB standard protonation
  ms_dipole.py --pqr state_0001.pqr   # single microstate PQR
  ms_dipole.py --dir /path/to/run     # ensemble from another dir
  ms_dipole.py --ph 4.0               # PyMOL snapshot at pH 4
"""

import argparse
import os
import subprocess
import sys
import numpy as np

# Resolve imports — works from MCCE_bin/ or anywhere via PATH/symlink
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from mcce_dipole.parsers import (parse_step2_pdb, parse_head3lst, parse_fort38,
                                 parse_pqr, parse_tpl_charges)
from mcce_dipole.compute import (compute_from_pqr, compute_from_ensemble,
                                 compute_from_step1, dipole_magnitude)
from mcce_dipole.visualize import generate_pymol_script, generate_ph_scan_csv


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def find_file(directory, *candidates):
    """Return the first existing file path from candidate names, or None."""
    for name in candidates:
        path = os.path.join(directory, name)
        if os.path.exists(path):
            return path
    return None


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------
def print_header(title):
    print(f"\n{'=' * 60}")
    print(f"  {title}")
    print(f"{'=' * 60}")


def print_single_state(results, label=""):
    """Print formatted single-state dipole results with coordinates."""
    if label:
        print_header(label)

    center = results["center"]
    print(f"\n  Geometric center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")

    for name, key in [("Backbone", "backbone_dipole"),
                      ("Ionizable", "ionizable_dipole"),
                      ("Full Protein", "full_dipole")]:
        mu = results[key]
        mag = np.linalg.norm(mu)
        direction = mu / mag if mag > 0.1 else np.zeros(3)
        print(f"\n  {name} Dipole Moment:")
        print(f"    |mu| = {mag:.2f} D")
        print(f"    mu = ({mu[0]:8.2f}, {mu[1]:8.2f}, {mu[2]:8.2f}) D")
        if mag > 0.1:
            print(f"    dir  = ({direction[0]:.4f}, {direction[1]:.4f}, {direction[2]:.4f})")

    print(f"\n  Quadrupole Tensor (e*A^2):")
    Q = results["quadrupole_tensor"]
    for i in range(3):
        print(f"    [{Q[i, 0]:10.2f}  {Q[i, 1]:10.2f}  {Q[i, 2]:10.2f}]")
    eigvals = results["quadrupole_eigenvalues"]
    print(f"    Eigenvalues: ({eigvals[0]:.2f}, {eigvals[1]:.2f}, {eigvals[2]:.2f})")

    print(f"\n  Net Charge: {results['net_charge']:.3f} e")
    print(f"  Arrow direction: (-) to (+)")
    print()


def print_ph_table(results):
    """Print compact pH scan table with magnitudes and vectors."""
    ph = results["ph_values"]
    print(f"\n{'=' * 72}")
    print(f"  pH-Dependent Dipole Moments (Debye)")
    print(f"{'=' * 72}")
    print(f"  {'pH':>5s}  {'Backbone':>10s}  {'Ionizable':>10s}  {'Full':>10s}  {'Net Crg':>10s}")
    print(f"  {'-' * 5}  {'-' * 10}  {'-' * 10}  {'-' * 10}  {'-' * 10}")

    for i in range(len(ph)):
        mag_bb = np.linalg.norm(results["backbone_dipole"][i])
        mag_ion = np.linalg.norm(results["ionizable_dipole"][i])
        mag_full = np.linalg.norm(results["full_dipole"][i])
        nq = results["net_charge"][i]
        print(f"  {ph[i]:5.1f}  {mag_bb:10.2f}  {mag_ion:10.2f}  {mag_full:10.2f}  {nq:10.3f}")
    print()

    # Dipole vector components: negative end <-> positive end
    center = results["center"]
    print(f"  Geometric center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    print(f"\n{'=' * 100}")
    print(f"  pH-Dependent Dipole Vectors: (-) end  <-->  (+) end   (Debye)")
    print(f"  Convention: mu = Sum(q*r) points from (-) to (+)")
    print(f"{'=' * 100}")

    for i in range(len(ph)):
        bb = results["backbone_dipole"][i]
        ion = results["ionizable_dipole"][i]
        full = results["full_dipole"][i]

        print(f"\n  pH {ph[i]:.1f}:")
        for label, mu in [("Backbone ", bb),
                          ("Ionizable", ion),
                          ("Full     ", full)]:
            mag = np.linalg.norm(mu)
            if mag > 0.1:
                direction = mu / mag
                neg_end = center - direction * mag * 0.1
                pos_end = center + direction * mag * 0.1
                print(f"    {label}  ({neg_end[0]:8.2f},{neg_end[1]:8.2f},{neg_end[2]:8.2f})"
                      f"  <-->  ({pos_end[0]:8.2f},{pos_end[1]:8.2f},{pos_end[2]:8.2f})"
                      f"   |mu|={mag:.1f} D")
            else:
                print(f"    {label}  (magnitude too small: {mag:.2f} D)")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        prog="ms_dipole.py",
        description="MCCE4 Dipole & Quadrupole Moment Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Computes from the MCCE4 conformer ensemble (Boltzmann-weighted at each pH):

  1. Backbone Dipole       Peptide bond polarity (N,CA,C,O,H,HA)
  2. Ionizable Dipole      pH-dependent charge redistribution
                           (ASP,GLU,LYS,ARG,HIS,CYS,TYR,NTR,CTR)
  3. Full Protein Dipole   Total charge asymmetry (all atoms)
  4. Quadrupole Tensor     Charge distribution shape (traceless 3x3,
                           eigenvalues in e*A^2)

Required files (auto-detected in run directory):
  step2_out.pdb   atom coordinates + partial charges
  head3.lst       conformer metadata
  fort.38         Monte Carlo occupancies at each pH

Examples:
  ms_dipole.py                        # ensemble from cwd
  ms_dipole.py --dir /path/to/run     # ensemble from another dir
  ms_dipole.py --pqr state_0001.pqr   # single microstate PQR
  ms_dipole.py --ph 4.0               # PyMOL snapshot at pH 4
  ms_dipole.py -o lysozyme            # custom output prefix

PyMOL controls (type 'dipole_help' after loading the .pml script):
  show_backbone / hide_backbone       solo_backbone
  show_ionizable / hide_ionizable     solo_ionizable
  show_full / hide_full               solo_full
  show_quadrupole / hide_quadrupole   solo_quadrupole
  show_all_dipoles / hide_all_dipoles
        """,
    )

    parser.add_argument("--pdb", type=str, default=None,
                        help="Input PDB file for standard protonation dipole. "
                             "Uses step1_out.pdb + param/mcce.tpl charges.")
    parser.add_argument("--pqr", type=str, default=None,
                        help="PQR file for single microstate calculation. "
                             "If omitted, runs ensemble mode.")
    parser.add_argument("--dir", type=str, default=".",
                        help="MCCE4 run directory (default: current directory)")
    parser.add_argument("--pdb_display", type=str, default=None,
                        help="PDB/PQR for PyMOL protein display "
                             "(default: step2_out.pdb or the PQR file)")
    parser.add_argument("--ph", type=float, default=7.0,
                        help="pH for PyMOL snapshot (default: 7.0)")
    parser.add_argument("--arrow_scale", type=float, default=0.1,
                        help="Arrow length: Debye to Angstrom (default: 0.1)")
    parser.add_argument("-o", "--output_prefix", type=str, default="dipole",
                        help="Output file prefix (default: 'dipole')")

    args = parser.parse_args()

    # ===================================================================
    # PQR MODE — single microstate
    # ===================================================================
    if args.pqr:
        if not os.path.exists(args.pqr):
            sys.exit(f"Error: PQR file not found: {args.pqr}")

        print(f"\n  Reading PQR: {args.pqr}")
        pqr_atoms = parse_pqr(args.pqr)
        print(f"  Atoms: {len(pqr_atoms)}")

        results = compute_from_pqr(pqr_atoms)
        print_single_state(results, label=f"Microstate: {os.path.basename(args.pqr)}")

        # Generate PyMOL script
        pdb_file = args.pdb_display or args.pqr
        pml_path = f"{args.output_prefix}_pymol.pml"
        generate_pymol_script(pdb_file, results, pml_path,
                              ph_index=None, arrow_scale=args.arrow_scale)

        pse_path = pml_path.replace(".pml", ".pse")
        print(f"  Output:")
        print(f"    PyMOL script  {pml_path}")
        print(f"    PyMOL session {pse_path}  (open directly in PyMOL)")
        print(f"\n  -> pymol {pml_path}")
        print(f"  -> Type 'dipole_help' in PyMOL for toggle commands\n")
        return

    # ===================================================================
    # PDB MODE — standard protonation from step1_out.pdb + topology
    # ===================================================================
    if args.pdb:
        if not os.path.exists(args.pdb):
            sys.exit(f"Error: PDB file not found: {args.pdb}")

        mcce_dir = args.dir

        # Check for step1_out.pdb; run step1.py if missing
        step1_path = os.path.join(mcce_dir, "step1_out.pdb")
        if not os.path.exists(step1_path):
            print(f"\n  step1_out.pdb not found — running step1.py...")
            step1_script = os.path.join(os.path.dirname(SCRIPT_DIR), "bin", "step1.py")
            if not os.path.exists(step1_script):
                step1_script = "step1.py"
            try:
                result = subprocess.run(
                    [sys.executable, step1_script, args.pdb],
                    cwd=mcce_dir, capture_output=True, text=True, timeout=300
                )
                if result.returncode != 0:
                    print(f"\n  step1.py failed (exit code {result.returncode}):")
                    if result.stderr:
                        print(result.stderr[:500])
                    sys.exit("  Fix these issues and run ms_dipole.py --pdb again.")
            except FileNotFoundError:
                sys.exit("  Error: step1.py not found. Run step1 manually first.")

        if not os.path.exists(step1_path):
            sys.exit(f"  Error: step1_out.pdb still not found after running step1.py")

        # Load topology charges
        param_dir = os.path.join(mcce_dir, "param")
        tpl_path = os.path.join(param_dir, "mcce.tpl")
        if not os.path.exists(tpl_path):
            sys.exit(f"  Error: {tpl_path} not found. Run step1.py first to generate topology.")

        print_header("PDB Standard Protonation Dipole Analysis")
        print(f"  Input PDB: {args.pdb}")

        print(f"\n  Loading topology charges from {tpl_path}...")
        tpl_charges = parse_tpl_charges(param_dir)
        print(f"  Topology entries: {len(tpl_charges)}")

        print(f"  Parsing step1_out.pdb...", end=" ", flush=True)
        conformers, all_atoms = parse_step2_pdb(step1_path)
        print(f"{len(conformers)} conformers, {len(all_atoms)} atoms")

        print(f"  Computing dipole with standard protonation charges...")
        results = compute_from_step1(conformers, all_atoms, tpl_charges)
        print_single_state(results, label="Standard Protonation State")

        # Generate outputs with _pdb suffix
        pdb_file = args.pdb_display or args.pdb
        pml_path = f"{args.output_prefix}_pdb_pymol.pml"
        generate_pymol_script(pdb_file, results, pml_path,
                              ph_index=None, arrow_scale=args.arrow_scale)

        pse_path = pml_path.replace(".pml", ".pse")
        print(f"  Output:")
        print(f"    PyMOL script  {pml_path}")
        print(f"    PyMOL session {pse_path}  (open directly in PyMOL)")
        print(f"\n  -> pymol {pml_path}")
        print(f"  -> Type 'dipole_help' in PyMOL for toggle commands\n")
        return

    # ===================================================================
    # ENSEMBLE MODE — default
    # ===================================================================
    mcce_dir = args.dir

    # Locate required files
    step2_path = find_file(mcce_dir, "step2_out.pdb")
    head3_path = find_file(mcce_dir, "head3.lst")
    fort38_path = find_file(mcce_dir, "fort.38")

    missing = []
    if not step2_path: missing.append("step2_out.pdb")
    if not head3_path: missing.append("head3.lst")
    if not fort38_path: missing.append("fort.38")

    if missing:
        print(f"\n  Error: Required files not found in {os.path.abspath(mcce_dir)}:")
        for f in missing:
            print(f"    - {f}")
        print(f"\n  Run from a completed MCCE4 directory, or use --dir.")
        print(f"  For single microstate: ms_dipole.py --pqr <file.pqr>")
        sys.exit(1)

    # --- Parse ---
    print_header("MCCE4 Ensemble Dipole Analysis")
    print(f"  Run dir: {os.path.abspath(mcce_dir)}")

    print(f"\n  Parsing step2_out.pdb...", end=" ", flush=True)
    conformers, all_atoms = parse_step2_pdb(step2_path)
    print(f"{len(conformers)} conformers, {len(all_atoms)} atoms")

    print(f"  Parsing head3.lst...", end=" ", flush=True)
    head3_data = parse_head3lst(head3_path)
    print(f"{len(head3_data)} entries")

    print(f"  Parsing fort.38...", end=" ", flush=True)
    ph_values, conf_ids, occupancies = parse_fort38(fort38_path)
    print(f"{len(ph_values)} pH values, {len(conf_ids)} conformers")

    # Check conformer overlap
    fort38_set = set(conf_ids)
    pdb_set = set(conformers.keys())
    matched = fort38_set & pdb_set
    print(f"  Conformer match: {len(matched)}/{len(fort38_set)} from fort.38 found in step2_out.pdb")

    # --- Compute ---
    print(f"\n  Computing pH-dependent dipole moments...")
    results = compute_from_ensemble(
        conformers, head3_data,
        ph_values, conf_ids, occupancies
    )

    # --- Print results ---
    print_ph_table(results)

    # --- Write outputs ---
    csv_path = f"{args.output_prefix}_ph_scan.csv"
    generate_ph_scan_csv(results, csv_path)

    ph_idx = np.argmin(np.abs(ph_values - args.ph))
    actual_ph = ph_values[ph_idx]
    pdb_file = args.pdb_display or find_file(mcce_dir, "prot_center.pdb", "prot.pdb")
    if not pdb_file:
        print(f"\n  Error: No protein PDB found for visualization.")
        print(f"  Looked for prot_center.pdb and prot.pdb in {os.path.abspath(mcce_dir)}")
        print(f"  Please provide one with: ms_dipole.py --pdb <your_protein.pdb>")
        sys.exit(1)
    print(f"  Protein PDB for visualization: {os.path.basename(pdb_file)}")
    pml_path = f"{args.output_prefix}_pH{actual_ph:.0f}_pymol.pml"
    generate_pymol_script(pdb_file, results, pml_path,
                          ph_index=ph_idx, arrow_scale=args.arrow_scale)

    pse_path = pml_path.replace(".pml", ".pse")
    print(f"  Output:")
    print(f"    pH scan CSV   {csv_path}")
    print(f"    PyMOL script  {pml_path}  (pH {actual_ph:.1f})")
    print(f"    PyMOL session {pse_path}  (open directly in PyMOL)")
    print(f"\n  -> pymol {pml_path}")
    print(f"  -> Type 'dipole_help' in PyMOL for toggle commands\n")


if __name__ == "__main__":
    main()
