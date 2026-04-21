"""
__main__.py - Run as: python -m mcce_dipole

Primary entry point is ms_dipole.py. This module provides the same
functionality when running the package directly.
"""

import argparse
import os
import subprocess
import sys
import numpy as np

from .parsers import parse_step2_pdb, parse_head3lst, parse_fort38, parse_pqr, parse_tpl_charges
from .compute import compute_from_pqr, compute_from_ensemble, compute_from_step1, dipole_magnitude
from .visualize import generate_pymol_script, generate_ph_scan_csv


def find_file(directory, *candidates):
    for name in candidates:
        path = os.path.join(directory, name)
        if os.path.exists(path):
            return path
    return None


def main():
    parser = argparse.ArgumentParser(
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
  python -m mcce_dipole                        # ensemble from cwd
  python -m mcce_dipole --dir /path/to/run     # from another dir
  python -m mcce_dipole --pqr state_0001.pqr   # single microstate
  python -m mcce_dipole --ph 4.0               # snapshot at pH 4
        """,
    )

    parser.add_argument("--pqr", type=str, default=None,
                        help="PQR file for single microstate calculation")
    parser.add_argument("--pdb", type=str, default=None, metavar="FILE",
                        help="Input PDB for single-structure dipole using "
                             "standard protonation charges from step1/topology")
    parser.add_argument("--dir", type=str, default=".",
                        help="MCCE4 run directory (default: cwd)")
    parser.add_argument("--pdb_display", type=str, default=None,
                        help="PDB for PyMOL display (default: input file)")
    parser.add_argument("--ph", type=float, default=7.0,
                        help="pH for PyMOL snapshot (default: 7.0)")
    parser.add_argument("--arrow_scale", type=float, default=0.1,
                        help="Arrow length: Debye to Angstrom (default: 0.1)")
    parser.add_argument("-o", "--output_prefix", type=str, default="dipole",
                        help="Output file prefix (default: 'dipole')")

    args = parser.parse_args()

    if args.pqr:
        if not os.path.exists(args.pqr):
            sys.exit(f"Error: PQR file not found: {args.pqr}")
        print(f"\n  Reading PQR: {args.pqr}")
        pqr_atoms = parse_pqr(args.pqr)
        print(f"  Atoms: {len(pqr_atoms)}")
        results = compute_from_pqr(pqr_atoms)

        pdb_file = args.pdb_display or args.pqr
        pml_path = f"{args.output_prefix}_pymol.pml"
        generate_pymol_script(pdb_file, results, pml_path,
                              ph_index=None, arrow_scale=args.arrow_scale)
        print(f"  PyMOL script: {pml_path}\n")
        return

    if args.pdb:
        pdb_input = args.pdb
        if not os.path.exists(pdb_input):
            sys.exit(f"Error: PDB file not found: {pdb_input}")

        mcce_dir = args.dir
        step1_path = os.path.join(mcce_dir, "step1_out.pdb")
        param_path = os.path.join(mcce_dir, "param")

        if not os.path.exists(step1_path):
            print(f"\n  step1_out.pdb not found, running step1.py on {pdb_input}...")
            result = subprocess.run(
                ["step1.py", pdb_input],
                cwd=os.path.abspath(mcce_dir),
                capture_output=True, text=True
            )
            if result.returncode != 0 or not os.path.exists(step1_path):
                print(f"\n  Error: step1.py failed.")
                if result.stdout:
                    print(f"  stdout: {result.stdout.strip()}")
                if result.stderr:
                    print(f"  stderr: {result.stderr.strip()}")
                print(f"\n  Please fix these issues and run again.")
                sys.exit(1)

        if not os.path.exists(param_path):
            sys.exit(f"Error: param/ directory not found in {os.path.abspath(mcce_dir)}.")

        conformers, all_atoms = parse_step2_pdb(step1_path)
        tpl_charges = parse_tpl_charges(param_path)
        results = compute_from_step1(conformers, all_atoms, tpl_charges)

        pdb_file = args.pdb_display or pdb_input
        prefix = f"{args.output_prefix}_pdb"
        pml_path = f"{prefix}_pymol.pml"
        generate_pymol_script(pdb_file, results, pml_path,
                              ph_index=None, arrow_scale=args.arrow_scale)
        print(f"  PyMOL script: {pml_path}\n")
        return

    mcce_dir = args.dir
    step2_path = find_file(mcce_dir, "step2_out.pdb")
    head3_path = find_file(mcce_dir, "head3.lst")
    fort38_path = find_file(mcce_dir, "fort.38")

    missing = [f for f, p in [("step2_out.pdb", step2_path),
                                ("head3.lst", head3_path),
                                ("fort.38", fort38_path)] if not p]
    if missing:
        print(f"\n  Error: Missing in {os.path.abspath(mcce_dir)}: {', '.join(missing)}")
        print(f"  Use --dir or --pqr.")
        sys.exit(1)

    conformers, all_atoms = parse_step2_pdb(step2_path)
    head3_data = parse_head3lst(head3_path)
    ph_values, conf_ids, occupancies = parse_fort38(fort38_path)

    print(f"\n  {len(conformers)} conformers, {len(all_atoms)} atoms, "
          f"{len(ph_values)} pH values")

    results = compute_from_ensemble(conformers, head3_data,
                                    ph_values, conf_ids, occupancies)

    csv_path = f"{args.output_prefix}_ph_scan.csv"
    generate_ph_scan_csv(results, csv_path)

    ph_idx = np.argmin(np.abs(ph_values - args.ph))
    actual_ph = ph_values[ph_idx]
    pdb_file = args.pdb_display or step2_path
    pml_path = f"{args.output_prefix}_pH{actual_ph:.0f}_pymol.pml"
    generate_pymol_script(pdb_file, results, pml_path,
                          ph_index=ph_idx, arrow_scale=args.arrow_scale)

    print(f"  CSV: {csv_path}")
    print(f"  PyMOL: {pml_path} (pH {actual_ph:.1f})\n")


if __name__ == "__main__":
    main()
