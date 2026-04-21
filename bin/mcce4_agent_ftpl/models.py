"""
Data models for the MCCE4 Topology Agent v3.

v3 additions to ConformerState:
  - h_added / h_removed: track which H atoms differ from neutral
  - per_state_charges: charges computed on this state's PDB
  - rxn_values / dsolv_values: per-state calibration results
"""

from dataclasses import dataclass, field, asdict
from typing import Optional, TypedDict


@dataclass
class ConformerState:
    """Represents a single protonation/redox conformer state.

    v3: Each state has its own PDB with correct H atoms.
    """
    label: str              # e.g., "01", "+1", "-1"
    charge: int = 0         # formal charge
    nH: int = 0             # protons relative to neutral
    pka: Optional[float] = None  # solution pKa
    smiles: str = ""        # SMILES for this state
    pdb_path: Optional[str] = None  # PDB file for this state (v3: per-state)
    source: str = "auto"    # "user", "dimorphite", "llm", "auto"
    site_atom: Optional[str] = None  # atom where protonation occurs
    rationale: str = ""     # why this state was included

    # ── v3: per-state PDB fields ──
    h_added: list = field(default_factory=list)     # H atom names added vs neutral
    h_removed: list = field(default_factory=list)   # H atom names removed vs neutral
    per_state_charges: dict = field(default_factory=dict)  # atom_name → charge
    rxn_values: dict = field(default_factory=dict)  # {rxn02: v, rxn04: v, rxn08: v}
    dsolv_values: dict = field(default_factory=dict)  # {2: dsolv, 4: dsolv, 8: dsolv}

    # ── v3: provenance & proton exchange info ──
    proton_exchange: str = ""   # Which protons differ vs neutral (e.g., "+H on N1 piperidine")
    llm_model: str = ""         # LLM that produced this state (e.g., "gemini-2.5-flash")
    references: list = field(default_factory=list)  # Literature refs (DOIs, URLs, titles)

    def to_dict(self) -> dict:
        return asdict(self)

    def __repr__(self):
        return (f"ConformerState({self.label}, charge={self.charge:+d}, "
                f"nH={self.nH:+d}, source={self.source})")


def make_labels_unique(states: list) -> list:
    """Ensure all conformer state labels are distinct and ≤ 2 characters.

    5-Character Conformer Naming Convention
    ========================================
    MCCE conformer type names are lig_id (3 chars) + label (≤ 2 chars) = 5 chars.

    Naming rules:
      - Neutral states:  indexed as 01, 02, ..., 0n
      - Unique ionic states (single conformer at that charge):
            identified by polarity + magnitude, e.g., +1, -1, +2, -2
      - Redundant ionic states (multiple conformers with identical net charge):
            transition to alpha-numeric code to preserve the index.
            Uppercase letters encode positive magnitudes:
                A = +1, B = +2, C = +3, ...
            Lowercase letters encode negative magnitudes:
                a = -1, b = -2, c = -3, ...
            The final digit is the conformer index (1-based):
                +1, +1  →  A1, A2   (two conformers at charge +1)
                -1, -1  →  a1, a2   (two conformers at charge -1)
                +2, +2  →  B1, B2   (two conformers at charge +2)
                -2, -2  →  b1, b2   (two conformers at charge -2)
      - A lone state at a given charge keeps its original label.

    Works with both ConformerState objects (mutated in place) and dicts
    (replaced with a shallow copy).  Returns the modified list.
    """
    import logging
    from collections import Counter

    # ── Magnitude → letter mappings ──
    # Positive: +1→A, +2→B, +3→C, ...
    # Negative: -1→a, -2→b, -3→c, ...
    def _charge_letter(charge: int) -> str:
        mag = abs(charge)
        if charge > 0:
            return chr(ord('A') + mag - 1)  # +1→A, +2→B, +3→C
        else:
            return chr(ord('a') + mag - 1)  # -1→a, -2→b, -3→c

    def _get_charge(s) -> int:
        if isinstance(s, dict):
            return int(s.get('charge', 0) or 0)
        return int(getattr(s, 'charge', 0) or 0)

    def _get_label(s) -> str:
        return s['label'] if isinstance(s, dict) else s.label

    # Count how many states share the same label
    label_counts = Counter(_get_label(s) for s in states)

    # Also count how many states share the same charge (for disambiguation)
    charge_counts = Counter(_get_charge(s) for s in states)

    label_index: dict = {}     # base_label → next index
    charge_index: dict = {}    # charge → next index (for disambiguation)
    neutral_index = 0          # for multiple neutral states

    result = []
    for s in states:
        label = _get_label(s)
        charge = _get_charge(s)

        if label_counts[label] > 1 or charge_counts[charge] > 1:
            # Need disambiguation
            if charge == 0:
                # Neutral: 01, 02, 03, ...
                neutral_index += 1
                new_label = f"0{neutral_index}"
            else:
                # Ionic with duplicates: use alpha-numeric code
                idx = charge_index.get(charge, 0)
                charge_index[charge] = idx + 1
                letter = _charge_letter(charge)
                new_label = f"{letter}{idx + 1}"
            if new_label != label:
                logging.info(f"  Label disambiguated: '{label}' → '{new_label}'")
        else:
            new_label = label
            # Single neutral state keeps "01"
            if charge == 0 and new_label in ("01",):
                pass  # keep as-is

        if isinstance(s, dict):
            s = {**s, 'label': new_label}
        else:
            s.label = new_label
        result.append(s)
    return result


def sort_conformer_labels(labels):
    """Sort conformer labels: neutral first, then positive, then negative.

    Label categories (5-character naming convention):
      Neutral:  starts with '0'  (01, 02, 03, ...)
      Positive: starts with '+' or uppercase A-Z (e.g., +1, +2, A1, A2, B1)
      Negative: starts with '-' or lowercase a-z (e.g., -1, -2, a1, a2, b1)

    Within each category, labels are sorted alphabetically.
    Example: ['01', '02', '+1', 'A1', 'A2', '-1', 'a1', 'a2']
    """
    def _sort_key(label):
        if label.startswith('0'):
            return (0, label)  # neutral
        elif label.startswith('+') or (label[0].isupper() and label[0] != 'M'):
            return (1, label)  # positive
        elif label.startswith('-') or label[0].islower():
            return (2, label)  # negative
        else:
            return (3, label)  # unknown — put last
    return sorted(labels, key=_sort_key)


class AgentState(TypedDict, total=False):
    """State object for the LangGraph agent.

    This is passed between nodes in the agent graph.
    Each node reads/writes fields as needed.
    """
    # ── Inputs ──
    pdb_path: str                   # Input ligand PDB file
    lig_id: str                     # 3-letter ligand code
    ph: float                       # Target pH
    charge_method: str              # Charge calculation method
    dielectrics: list               # Dielectric constants for RXN
    work_dir: str                   # Working directory

    # ── Molecule info ──
    smiles: str                     # Canonical SMILES
    name: str                       # Molecule name
    formula: str                    # Molecular formula
    formal_charge: int              # Net formal charge from RCSB

    # ── Conformer states ──
    states: list                    # List of ConformerState dicts
    conformer_labels: list          # List of label strings (e.g., ["01", "+1"])
    user_state_pdbs: list           # User-supplied state PDB paths

    # ── v3: Per-state PDBs ──
    state_pdbs: dict                # label → pdb_path (each with correct H atoms)
    h_diffs: dict                   # label → {"added": [...], "removed": [...]}
    per_state_connects: dict        # label → ftpl_path from pdb2ftpl

    # ── Generated files ──
    ftpl_path: str                  # Output .ftpl file path
    charges: dict                   # atom_name -> charge mapping
    per_state_charges: dict         # label → {atom_name: charge}

    # ── Validation ──
    hybridization_issues: list      # List of hybridization mismatches
    charge_issues: list             # List of charge validation issues
    rxn_values: dict                # conformer_type -> {rxn02: v, rxn04: v, rxn08: v}

    # ── Agent control ──
    phase: str                      # Current phase name
    errors: list                    # Accumulated errors
    warnings: list                  # Accumulated warnings
    messages: list                  # LLM conversation messages
    needs_user_review: bool         # Whether to pause for GUI review
    user_approved: bool             # Whether user approved states
    complete: bool                  # Whether agent is done
    dry_run: bool                   # Skip RXN calibration
