"""
Protonation state enumeration using Dimorphite-DL.

Naming convention:
  - Single state per charge: 01 (neutral), +1 (protonated), -1 (deprotonated)
  - Multiple states at same charge: +a/+b, -a/-b, 0a/0b, etc.
"""

import logging
from typing import List
from collections import Counter
from ..models import ConformerState, make_labels_unique
from ..config import CHARGE_TO_CONF


def enumerate_protonation_states(smiles: str, ph: float = 7.4,
                                  pdb_path: str = None,
                                  ph_min: float = 6.4,
                                  ph_max: float = 8.4,
                                  precision: float = 1.0,
                                  max_variants: int = 128,
                                  label_states: bool = False) -> List[ConformerState]:
    """Enumerate protonation states at target pH using Dimorphite-DL.

    Dimorphite-DL is an open-source tool that enumerates small-molecule
    ionization states (protonation/deprotonation) based on a pH range.
    It accepts SMILES input and returns protonated SMILES variants.

    This is distinct from Phase 1 ionizable site detection, which uses
    RDKit heuristics on the 3D PDB structure to identify *potential*
    protonation sites (N, O, S atoms). Dimorphite-DL actually *enumerates*
    the chemically relevant states at the target pH.

    Keeps ALL unique protonation states (multiple per charge).
    Applies naming: single charge → 01/+1/-1, multiple → +a/+b.
    Does NOT auto-add states — trusts Dimorphite-DL output, with
    LLM reasoning phase handling validation and trimming.

    Args:
        smiles: Input SMILES string.
        ph: Target pH (used for state rationale text).
        pdb_path: Optional path to PDB (unused, kept for API compat).
        ph_min: Minimum pH for protonation enumeration (default: 6.4).
        ph_max: Maximum pH for protonation enumeration (default: 8.4).
        precision: pKa precision factor (std devs from mean pKa).
        max_variants: Max protonation variants per compound.
        label_states: Label output SMILES as PROTONATED/DEPROTONATED/BOTH.

    Returns:
        List of ConformerState objects for all unique protonation states.
    """
    logging.info(f"  Dimorphite-DL: enumerating protonation states from SMILES")
    logging.info(f"    ph_min={ph_min}, ph_max={ph_max}")
    logging.info(f"    precision={precision}, max_variants={max_variants}"
                 f"{', label_states=True' if label_states else ''}")
    logging.info(f"    Input SMILES: {smiles}")

    state_smiles = [smiles]  # fallback

    try:
        from dimorphite_dl import protonate_smiles
        results = protonate_smiles(
            smiles,
            ph_min=ph_min,
            ph_max=ph_max,
            precision=precision,
            max_variants=max_variants,
            label_states=label_states,
        )
        if results:
            state_smiles = list(set(results))
        logging.info(f"  ✓ Dimorphite-DL returned {len(state_smiles)} unique state(s)")
        for i, smi in enumerate(state_smiles, 1):
            logging.info(f"    [{i}] {smi}")
    except ImportError:
        logging.warning("  ⚠ Dimorphite-DL not installed (pip install dimorphite-dl)")
        logging.warning("    Falling back to input SMILES as single neutral state")
    except Exception as e:
        logging.warning(f"  ⚠ Dimorphite-DL error: {e}")
        logging.warning("    Falling back to input SMILES as single neutral state")

    # Compute formal charges and build ConformerState objects
    states = _smiles_to_states(state_smiles, ph)

    # Trust Dimorphite-DL output — do NOT auto-add a -1 state.
    # The LLM reasoning phase will validate and trim states as needed.

    # Ensure neutral (01) is always present
    has_neutral = any(s.charge == 0 for s in states)
    if not has_neutral:
        states.insert(0, ConformerState(
            label="01", charge=0, nH=0, smiles=smiles,
            source="dimorphite", rationale="Neutral reference (auto-added)"
        ))

    # Apply naming convention via make_labels_unique
    states = make_labels_unique(states)

    for s in states:
        logging.info(f"  📋 {s.label}: charge={s.charge:+d}, nH={s.nH:+d}, "
                     f"smiles={s.smiles[:50] if s.smiles else 'N/A'}...")

    return states


def _auto_add_deprotonation(pdb_path: str, ph: float) -> ConformerState:
    """Check for deprotonatable sites (O-H, S-H, N-H) and create a -1 state.

    Returns a ConformerState or None.
    """
    try:
        from .charge_tools import get_ionizable_sites_oe
        sites = get_ionizable_sites_oe(pdb_path)
        deprot_candidates = [
            s for s in sites
            if s.get("current_hs", 0) > 0 and s.get("type", "").startswith("deprotonatable")
        ]
        if not deprot_candidates:
            # Also check N-H sites (secondary amines can be deprotonated)
            deprot_candidates = [
                s for s in sites
                if s.get("current_hs", 0) > 0
            ]
        if deprot_candidates:
            # Prefer O-H (most commonly deprotonatable), then S-H, then N-H
            for preferred_sym in ("O", "S", "N"):
                for c in deprot_candidates:
                    if c.get("symbol") == preferred_sym:
                        return ConformerState(
                            label="-1", charge=-1, nH=-1,
                            smiles="", source="auto-detected",
                            site_atom=c["name"],
                            rationale=f"Deprotonation at {c['name']} ({c['type']})",
                            proton_exchange=f"-H from {c['name']}",
                        )
            # Fallback: first candidate
            c = deprot_candidates[0]
            return ConformerState(
                label="-1", charge=-1, nH=-1,
                smiles="", source="auto-detected",
                site_atom=c["name"],
                rationale=f"Deprotonation at {c['name']} ({c['type']})",
                proton_exchange=f"-H from {c['name']}",
            )
    except Exception as e:
        logging.warning(f"  Auto-deprotonation detection failed: {e}")
    return None


def _smiles_to_states(smiles_list: list, ph: float) -> List[ConformerState]:
    """Convert SMILES list to ConformerState objects.

    Keeps ALL unique SMILES (multiple per charge level).
    Labels are assigned using MCCE convention:
      - Single state at charge 0  → "01"
      - Single state at charge +1 → "+1"
      - Single state at charge -1 → "-1"
      - Multiple states at charge +1 → "+1", "+1" (make_labels_unique → "+a", "+b")
    """
    states = []
    seen_smiles: set = set()

    try:
        from rdkit import Chem

        for smi in smiles_list:
            if smi in seen_smiles:
                continue
            seen_smiles.add(smi)

            mol = Chem.MolFromSmiles(smi)
            charge = Chem.GetFormalCharge(mol) if mol else 0

            label = CHARGE_TO_CONF.get(charge, f"{charge:+d}"[-2:])
            nH = charge  # relative to neutral

            states.append(ConformerState(
                label=label, charge=charge,
                nH=nH if charge != 0 else 0,
                smiles=smi, source="dimorphite",
                rationale=f"Dimorphite-DL at pH {ph}"
            ))

    except ImportError:
        # No RDKit fallback
        states.append(ConformerState(
            label="01", charge=0, nH=0, smiles=smiles_list[0],
            source="fallback", rationale="No RDKit — assumed neutral"
        ))

    return make_labels_unique(states)
