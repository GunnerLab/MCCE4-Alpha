"""
compute.py - Dipole and quadrupole moment calculations.

Physics:
  Dipole moment:     mu = Sum qi * ri           (vector, 3 components)
  Quadrupole tensor: Q_ab = Sum qi (3*ria*rib - |ri|^2 * d_ab)   (traceless symmetric 3x3)

  For MCCE4 ensemble at a given pH:
    mu(pH) = Sum_conf [ occ(conf, pH) * Sum_atoms_in_conf ( q_atom * r_atom ) ]

  where occupancies come from fort.38 Monte Carlo output, and charges are
  read directly from step2_out.pdb (embedded partial charges from topology).

Units:
  - Input coordinates: Angstrom (from PDB)
  - Input charges: elementary charge (e)
  - Output dipole: Debye  (1 e*A = 4.80321 D)
  - Output quadrupole: e*A^2 (Buckingham)
"""

import numpy as np
from . import E_ANG_TO_DEBYE, IONIZABLE_RES, BACKBONE_ATOMS


def _dipole_raw(coords, charges):
    """
    Compute raw dipole moment vector in e*A.

    Parameters:
        coords:  np.array of shape (N, 3)
        charges: np.array of shape (N,)

    Returns:
        mu: np.array of shape (3,) in e*A
    """
    return np.sum(charges[:, None] * coords, axis=0)


def _quadrupole_raw(coords, charges, center):
    """
    Compute traceless quadrupole tensor in e*A^2.

    Parameters:
        coords:  np.array of shape (N, 3)
        charges: np.array of shape (N,)
        center:  np.array of shape (3,) — reference point

    Returns:
        Q: np.array of shape (3, 3) — traceless quadrupole tensor
    """
    r = coords - center
    r2 = np.sum(r ** 2, axis=1)  # |r|^2 for each atom
    Q = np.zeros((3, 3))
    for alpha in range(3):
        for beta in range(3):
            delta = 1.0 if alpha == beta else 0.0
            Q[alpha, beta] = np.sum(
                charges * (3.0 * r[:, alpha] * r[:, beta] - r2 * delta)
            )
    return Q


def compute_center_of_geometry(coords):
    """Geometric center of a set of coordinates."""
    return np.mean(coords, axis=0)


# ---------------------------------------------------------------------------
# PQR-based (single state) calculations
# ---------------------------------------------------------------------------
def compute_from_pqr(pqr_atoms, backbone_only_names=None):
    """
    Compute dipole and quadrupole from a parsed PQR file (single state).

    Parameters:
        pqr_atoms: list of dicts from parsers.parse_pqr()
        backbone_only_names: set of atom names for backbone selection

    Returns:
        dict with keys: backbone_dipole, full_dipole, ionizable_dipole, quadrupole, etc.
    """
    if backbone_only_names is None:
        backbone_only_names = BACKBONE_ATOMS

    all_coords = np.array([a["coords"] for a in pqr_atoms])
    all_charges = np.array([a["charge"] for a in pqr_atoms])
    center = compute_center_of_geometry(all_coords)

    coords_centered = all_coords - center

    # --- Full protein dipole ---
    mu_full = _dipole_raw(coords_centered, all_charges)

    # --- Backbone dipole ---
    bb_mask = np.array([a["name"] in backbone_only_names for a in pqr_atoms])
    mu_backbone = _dipole_raw(coords_centered[bb_mask], all_charges[bb_mask])

    # --- Ionizable dipole ---
    ion_mask = np.array([a["res_name"] in IONIZABLE_RES for a in pqr_atoms])
    mu_ionizable = _dipole_raw(coords_centered[ion_mask], all_charges[ion_mask])

    # --- Quadrupole (full protein) ---
    Q = _quadrupole_raw(all_coords, all_charges, center)
    eigenvalues, eigenvectors = np.linalg.eigh(Q)

    return {
        "center": center,
        "full_dipole": mu_full * E_ANG_TO_DEBYE,
        "backbone_dipole": mu_backbone * E_ANG_TO_DEBYE,
        "ionizable_dipole": mu_ionizable * E_ANG_TO_DEBYE,
        "quadrupole_tensor": Q,
        "quadrupole_eigenvalues": eigenvalues,
        "quadrupole_eigenvectors": eigenvectors,
        "net_charge": float(np.sum(all_charges)),
    }


# ---------------------------------------------------------------------------
# MCCE4 ensemble calculations (pH-dependent)
# ---------------------------------------------------------------------------
def compute_from_ensemble(conformers, head3_data,
                          ph_values, conf_ids, occupancies):
    """
    Compute pH-dependent dipole moments from the MCCE4 conformer ensemble.

    Charges are read directly from the Atom objects (embedded in step2_out.pdb),
    so no external topology file is needed for the calculation.

    Parameters:
        conformers:  dict {conf_id: [Atom, ...]} from parsers.parse_step2_pdb()
        head3_data:  dict {conf_id: {...}} from parsers.parse_head3lst()
        ph_values:   np.array of pH values from fort.38
        conf_ids:    list of conformer IDs from fort.38
        occupancies: np.array (n_ph, n_conf) from fort.38

    Returns:
        results: dict with pH-resolved dipole data
    """
    # Pre-compute per-conformer data: coordinates, charges, dipole contribution
    conf_dipoles = {}      # {conf_id: dipole_vector in e*A}
    conf_atom_data = {}    # {conf_id: (coords_array, charges_array)}

    # Collect ALL atom coordinates for geometric center
    all_coords_list = []
    for conf_id, atoms in conformers.items():
        for atom in atoms:
            all_coords_list.append(atom.coords)

    if not all_coords_list:
        raise ValueError("No atoms found in step2_out.pdb")

    global_center = compute_center_of_geometry(np.array(all_coords_list))

    # Build per-conformer dipole vectors using embedded charges
    for conf_id, atoms in conformers.items():
        coords = np.array([a.coords for a in atoms])
        charges = np.array([a.charge for a in atoms])

        coords_centered = coords - global_center
        mu = _dipole_raw(coords_centered, charges)
        conf_dipoles[conf_id] = mu
        conf_atom_data[conf_id] = (coords, charges)

    # Mapping from fort.38 conformer index → conf_id
    conf_id_to_idx = {cid: i for i, cid in enumerate(conf_ids)}

    # Classify conformers
    def is_backbone(conf_id):
        return len(conf_id) >= 5 and conf_id[3:5] == "BK"

    def is_ionizable(conf_id):
        return conf_id[:3] in IONIZABLE_RES

    # Compute at each pH
    n_ph = len(ph_values)
    results = {
        "ph_values": ph_values,
        "center": global_center,
        "full_dipole": np.zeros((n_ph, 3)),
        "backbone_dipole": np.zeros((n_ph, 3)),
        "ionizable_dipole": np.zeros((n_ph, 3)),
        "quadrupole_tensor": np.zeros((n_ph, 3, 3)),
        "quadrupole_eigenvalues": np.zeros((n_ph, 3)),
        "net_charge": np.zeros(n_ph),
    }

    for ph_idx in range(n_ph):
        mu_full = np.zeros(3)
        mu_bb = np.zeros(3)
        mu_ion = np.zeros(3)
        total_charge = 0.0

        # For quadrupole: accumulate weighted atoms
        all_weighted_coords = []
        all_weighted_charges = []

        for conf_id, mu_conf in conf_dipoles.items():
            # Get occupancy for this conformer at this pH
            if conf_id in conf_id_to_idx:
                occ = occupancies[ph_idx, conf_id_to_idx[conf_id]]
            elif conf_id in head3_data:
                occ = head3_data[conf_id]["occ"]  # fallback: fixed occ
            else:
                # Backbone conformers have occ=1 by default
                if is_backbone(conf_id):
                    occ = 1.0
                else:
                    occ = 0.0

            if occ < 1e-8:
                continue

            # Weighted dipole contribution
            weighted_mu = occ * mu_conf
            mu_full += weighted_mu

            if is_backbone(conf_id):
                mu_bb += weighted_mu
            if is_ionizable(conf_id):
                mu_ion += weighted_mu

            # Net charge
            if conf_id in head3_data:
                total_charge += occ * head3_data[conf_id]["crg"]

            # Accumulate for quadrupole
            coords, charges = conf_atom_data[conf_id]
            all_weighted_coords.append(coords)
            all_weighted_charges.append(occ * charges)

        results["full_dipole"][ph_idx] = mu_full * E_ANG_TO_DEBYE
        results["backbone_dipole"][ph_idx] = mu_bb * E_ANG_TO_DEBYE
        results["ionizable_dipole"][ph_idx] = mu_ion * E_ANG_TO_DEBYE
        results["net_charge"][ph_idx] = total_charge

        # Quadrupole from weighted atoms
        if all_weighted_coords:
            qp_coords = np.concatenate(all_weighted_coords, axis=0)
            qp_charges = np.concatenate(all_weighted_charges, axis=0)
            Q = _quadrupole_raw(qp_coords, qp_charges, global_center)
            results["quadrupole_tensor"][ph_idx] = Q
            eigvals, _ = np.linalg.eigh(Q)
            results["quadrupole_eigenvalues"][ph_idx] = eigvals

    return results


# ---------------------------------------------------------------------------
# PDB standard-protonation (single state from step1_out.pdb + topology)
# ---------------------------------------------------------------------------
STANDARD_PROTONATION = {
    "ARG": "+1",   # pKa ~12.5, protonated
    "LYS": "+1",   # pKa ~10.5, protonated
    "HIS": "01",   # pKa ~6.0, neutral
    "ASP": "-1",   # pKa ~3.9, deprotonated
    "GLU": "-1",   # pKa ~4.1, deprotonated
    "NTR": "+1",   # pKa ~8.0, protonated
    "CTR": "-1",   # pKa ~3.2, deprotonated
    "TYR": "01",   # pKa ~10.1, neutral
    "CYS": "01",   # pKa ~8.3, neutral
}


def compute_from_step1(conformers, all_atoms, tpl_charges):
    """
    Compute dipole/quadrupole from step1_out.pdb using standard protonation
    charges from the topology (param/mcce.tpl).

    For ionizable sidechain atoms, charges are looked up under the standard
    protonation conf_key (e.g. ARG+1) first, falling back to the original
    conf_key (e.g. ARG01) for atoms that exist in both forms.

    Note: step1_out.pdb only contains atoms from the 01 conformer geometry.
    Protons that exist only in +1/-1 forms (e.g. HE on ARG+1) are absent,
    so the net charge will be systematically underestimated for ionizable
    residues whose standard state differs from 01.

    Parameters:
        conformers: dict from parsers.parse_step2_pdb() (works on step1_out.pdb)
        all_atoms:  list of Atom objects from parse_step2_pdb()
        tpl_charges: dict {(conf_key, atom_name): charge} from parse_tpl_charges()

    Returns:
        dict — same format as compute_from_pqr()
    """
    pqr_atoms = []
    unmatched = 0

    for atom in all_atoms:
        if atom.res_name in STANDARD_PROTONATION and atom.conf_type != "BK":
            std_type = STANDARD_PROTONATION[atom.res_name]
            std_key = atom.res_name + std_type
            orig_key = atom.conf_key

            if (std_key, atom.name) in tpl_charges:
                charge = tpl_charges[(std_key, atom.name)]
            elif (orig_key, atom.name) in tpl_charges:
                charge = tpl_charges[(orig_key, atom.name)]
            else:
                charge = 0.0
                unmatched += 1
        else:
            conf_key = atom.conf_key
            charge = tpl_charges.get((conf_key, atom.name), 0.0)
            if (conf_key, atom.name) not in tpl_charges:
                unmatched += 1

        pqr_atoms.append({
            "name": atom.name,
            "res_name": atom.res_name,
            "chain": atom.chain,
            "res_seq": atom.res_seq,
            "coords": atom.coords,
            "charge": charge,
            "radius": atom.radius,
        })

    if unmatched > 0:
        import sys
        print(f"  Note: {unmatched}/{len(all_atoms)} atoms not in topology charge table "
              f"(non-polar atoms, charge=0.0)", file=sys.stderr)

    return compute_from_pqr(pqr_atoms)


def dipole_magnitude(mu_vector):
    """Magnitude of a dipole vector (or array of vectors)."""
    if mu_vector.ndim == 1:
        return np.linalg.norm(mu_vector)
    return np.linalg.norm(mu_vector, axis=1)
