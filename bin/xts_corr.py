#!/usr/bin/env python3
"""
Created on Oct 15 2025
@author: Gehan Ranepura

Name: xts_corr.py
Applies conformational entropy corrections to non-amino acid residues for file fort.38

ENTROPY CORRECTION METHOD:
==========================
This script corrects for the bias that arises when different charge states have
different numbers of conformers. Charge states with more conformers would otherwise
appear more probable simply due to having more "counts" in the ensemble.

The correction penalizes conformers belonging to charge states with multiple conformers.

EXAMPLE with GLU conformers:
----------------------------
Original probabilities:
  GLU01A0292_001: 0.005  }
  GLU01A0292_002: 0.007  } charge state 01, N=2
  GLU02A0292_003: 0.016  }
  GLU02A0292_004: 0.012  } charge state 02, N=2
  GLU-1A0292_005: 0.960  } charge state -1, N=1

Step 1: Convert probabilities to relative free energies
    ΔG_i/RT = -ln(P_i/P_ref)

    Using the most probable conformer as reference (P_ref = 0.960):

    GLU01A0292_001: -RT * ln(0.005/0.960) = 3.112
    GLU01A0292_002: -RT * ln(0.007/0.960) = 2.913
    GLU02A0292_003: -RT * ln(0.016/0.960) = 2.423
    GLU02A0292_004: -RT * ln(0.012/0.960) = 2.594
    GLU-1A0292_005: -RT * ln(0.960/0.960) = 0.000

Step 2: Conformers within a residue are grouped according to charge value. For example, GLU01 and GLU02 are both neutral (charge=0) so they form one group, while GLU-1 (charge=-1) forms another.
    For each group, entropy is computed as:
    E_TS = -RT * Σ(P_i * ln(P_i))

    Normalize each group's probabilities:
    Neutral group (charge=0): 0.005 + 0.007 + 0.016 + 0.012 = 0.04
    Charged group (charge=-1): 0.960

    GLU01A0292_001: 0.005/0.04 = 0.125
    GLU01A0292_002: 0.007/0.04 = 0.175
    GLU02A0292_003: 0.016/0.04 = 0.400
    GLU02A0292_004: 0.012/0.04 = 0.300
    GLU-1A0292_005: 0.960/0.960 = 1.000

    Shannon entropy for each group:
    Neutral (charge=0): -RT * [ (0.125 * ln(0.125)) + (0.175 * ln(0.175)) + (0.400 * ln(0.400)) + (0.300 * ln(0.300)) ] = 0.765
    Charged (charge=-1): -RT * [ (1.000 * ln(1.000) ] = 0.000

Step 3: Each conformer's energy is corrected using its charge group's entropy:
    E_i_corr = E_i + E_TS

    Corrected energies:
      GLU01A0292_001: 3.112 + 0.765 = 3.877
      GLU01A0292_002: 2.913 + 0.765 = 3.678
      GLU02A0292_003: 2.423 + 0.765 = 3.188
      GLU02A0292_004: 2.594 + 0.765 = 3.359
      GLU-1A0292_005: 0.000 + 0.000 = 0.000

Step 4: Convert back to probabilities using Boltzmann distribution
    P_i = exp(-ΔG_corrected/RT) / Σ_j exp(-ΔG_corrected/RT)

    Boltzmann factors (exp(-ΔG/RT)):
      GLU01A0292_001: exp(-3.877/RT) = exp(-6.548) = 0.00143
      GLU01A0292_002: exp(-3.678/RT) = exp(-6.212) = 0.00200
      GLU02A0292_003: exp(-3.188/RT) = exp(-5.385) = 0.00458
      GLU02A0292_004: exp(-3.359/RT) = exp(-5.674) = 0.00343
      GLU-1A0292_005: exp(-0.000/RT) = exp(-0.000) = 1.00000

    Partition function Z = 1.01144

    Final corrected probabilities:
      GLU01A0292_001: 0.00143/1.01144 = 0.00141 (0.14%)  [was 0.50%]
      GLU01A0292_002: 0.00200/1.01144 = 0.00198 (0.20%)  [was 0.70%]
      GLU02A0292_003: 0.00458/1.01144 = 0.00453 (0.45%)  [was 1.60%]
      GLU02A0292_004: 0.00343/1.01144 = 0.00339 (0.34%)  [was 1.20%]
      GLU-1A0292_005: 1.00000/1.01144 = 0.98868 (98.87%) [was 96.00%]

RESULT: All neutral conformers (GLU01 and GLU02) share the same entropy penalty
(0.765) because they belong to the same charge state (charge=0). The charged state
(GLU-1) with only 1 conformer receives no penalty and becomes more dominant.

EFFECT:
- Conformers are grouped by their numeric charge value (not exact charge state string)
- Charge states with more conformers (larger N) get penalized more
- The penalty depends on the Shannon entropy: -Σ(P_i × ln(P_i)) within each charge group
- For equal probabilities: N=2 → 0.41 kcal/mol, N=3 → 0.65 kcal/mol, N=4 → 0.82 kcal/mol
- This lowers the probability of conformers in charge states with many conformers
- Charge states with single conformers (N=1) get no penalty (entropy = 0)

CONSTANTS:
- R = 1.987×10⁻³ kcal·K⁻¹·mol⁻¹
- T = 298 K
- RT = 0.592 kcal/mol
- Energy penalty = RT × Shannon_entropy kcal/mol
  where Shannon_entropy = -Σ(P_i × ln(P_i)) for normalized probabilities within each charge group

NOTE:
- Actual charge values are read from head3.lst (column 5) for accurate sum_crg.out generation
- If head3.lst is not available, charges are parsed from conformer names
"""

import re
import math
import argparse
from collections import defaultdict
from typing import Dict, List, Tuple

# Constants
R = 1.987e-3  # kcal·K⁻¹·mol⁻¹
T = 298.0     # K
RT = R * T    # = 0.592 kcal/mol

# Amino acids to exclude from entropy correction (unless --doAAs flag is used)
AMINO_ACIDS_LIST = {
    'ACE', 'NME', 'CTR', 'NTR', 'CTG', 'NTG',
    'ALA', 'ARG', 'ASN', 'ASP', 'CYD', 'GLN', 'GLU',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL', 'CYS'
}


class Conformer:
    """Represents a single conformer"""
    def __init__(self, name: str, probabilities: List[float], actual_charge: float = None):
        self.full_name = name
        self.probabilities = probabilities
        self.original_probs = probabilities.copy()
        self.actual_charge = actual_charge  # Actual charge from head3.lst

        # Parse the conformer name
        parsed = self._parse_name(name)
        if parsed:
            self.res_type = parsed['res_type']
            self.charge_state = parsed['charge_state']
            self.chain = parsed['chain']
            self.res_num = parsed['res_num']
            self.conf_num = parsed['conf_num']
            self.full_res_id = f"{self.res_type}{self.chain}{self.res_num}"
            # Group by numeric charge value (not exact charge state string)
            self.charge_value = parse_charge_state(self.charge_state)
            self.charge_group = f"{self.res_type}_charge{self.charge_value}"
        else:
            raise ValueError(f"Cannot parse conformer name: {name}")

    @staticmethod
    def _parse_name(name: str) -> Dict[str, str]:
        """Parse conformer name: e.g., '0WN+1A1101_001'"""
        parts = name.split('_')
        if len(parts) != 2:
            return None

        conf_num = parts[1]
        main_part = parts[0]

        # Pattern: RESTYPE(3 chars) + CHARGE_STATE + CHAIN(1 letter) + RESNUM(4 digits)
        match = re.match(r'^([A-Z0-9]{3})([\+\-]?\w+)([A-Z])(\d{4})$', main_part)
        if not match:
            return None

        return {
            'res_type': match.group(1),
            'charge_state': match.group(2),
            'chain': match.group(3),
            'res_num': match.group(4),
            'conf_num': conf_num
        }


def parse_charge_state(charge_state: str) -> float:
    """Parse charge state string to numeric value

    Examples:
        '+1' -> 1.0
        '+2' -> 2.0
        '-'  -> -1.0
        '-2' -> -2.0
        '01' -> 0.0
        '02' -> 0.0
        '+a' -> 2.0 (doubly protonated)
        '+b' -> 3.0 (triply protonated)
    """
    if not charge_state:
        return 0.0

    # Handle numeric formats
    if charge_state[0] in '+-':
        # '+1', '-1', '+2', etc.
        if len(charge_state) > 1 and charge_state[1:].isdigit():
            return float(charge_state)
        # Just '+' means +1, just '-' means -1
        elif charge_state == '+':
            return 1.0
        elif charge_state == '-':
            return -1.0
        # Handle letter codes: '+a' = +2, '+b' = +3, etc.
        elif len(charge_state) > 1 and charge_state[1].isalpha():
            letter_charge = ord(charge_state[1].lower()) - ord('a') + 2
            return float(letter_charge) if charge_state[0] == '+' else -float(letter_charge)

    # Handle '01', '02' formats (neutral)
    if charge_state.startswith('0'):
        return 0.0

    # Default to 0
    return 0.0


def read_head3_charges(head3_file: str) -> Dict[str, float]:
    """Read actual charge values from head3.lst file
    
    Returns:
        Dictionary mapping conformer names to their actual charges
    """
    charges = {}
    
    import os
    if not os.path.exists(head3_file):
        print(f"Warning: {head3_file} not found. Will use parsed charges from conformer names.")
        return charges
    
    with open(head3_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('iConf'):
                continue
            
            parts = line.split()
            if len(parts) < 5:
                continue
            
            conformer_name = parts[1]
            try:
                charge = float(parts[4])
                charges[conformer_name] = charge
            except (ValueError, IndexError):
                continue
    
    print(f"Read charges for {len(charges)} conformers from {head3_file}")
    return charges


def apply_entropy_correction(conformers: List[Conformer], ph_index: int) -> Tuple[List[float], Dict]:
    """Apply Shannon entropy correction to penalize charge states with more conformers

    This corrects for the bias where charge states with more conformers appear more
    probable simply because they have more "counts" in the ensemble.

    Method (all energies in dimensionless units of RT):
    1. Group conformers by numeric charge value
    2. Convert probabilities to relative free energies: ΔG/RT = -ln(P_i/P_ref)
    3. Calculate Shannon entropy for each group: E_TS/RT = -Σ(P_i * ln(P_i))
    4. Add entropy penalty: ΔG_corrected/RT = ΔG_original/RT + E_TS/RT
    5. Convert back to probabilities: P_i = exp(-ΔG_corrected/RT) / Z
    """

    # Group conformers by charge state
    charge_groups = defaultdict(list)
    for conf in conformers:
        charge_groups[conf.charge_group].append(conf)

    # Step 1: Find reference (most probable conformer)
    probs = [c.probabilities[ph_index] for c in conformers]
    p_ref = max(probs)

    # Step 2: Calculate relative free energies using reference (dimensionless: ΔG/RT)
    relative_energies = []
    for conf in conformers:
        prob = conf.probabilities[ph_index]

        # ΔG/RT = -ln(P_i/P_ref)
        if prob > 1e-10:
            relative_energy = -math.log(prob / p_ref)
        else:
            # For very small probabilities, use a large but finite energy
            relative_energy = 30.0

        relative_energies.append(relative_energy)

    # Step 3: Calculate Shannon entropy for each charge group (dimensionless: E_TS/RT)
    group_info = {}
    
    for charge_group, confs in charge_groups.items():
        num_conformers = len(confs)
        
        # Get probabilities for this group
        group_probs = [c.probabilities[ph_index] for c in confs]
        total_prob = sum(group_probs)
        
        # Calculate Shannon entropy: E_TS/RT = -Σ(P_i * ln(P_i))
        # where P_i are normalized probabilities within the group
        entropy_penalty_dimensionless = 0.0
        
        if total_prob > 1e-10 and num_conformers > 1:
            # Normalize probabilities within the group
            normalized_probs = [p / total_prob for p in group_probs]
            
            # Calculate Shannon entropy (dimensionless, in units of RT)
            for p_norm in normalized_probs:
                if p_norm > 1e-10:
                    entropy_penalty_dimensionless -= p_norm * math.log(p_norm)
        
        group_info[charge_group] = {
            'total_prob': total_prob,
            'num_conformers': num_conformers,
            'entropy_penalty_dimensionless': entropy_penalty_dimensionless,
            'entropy_penalty_kcal': entropy_penalty_dimensionless * RT
        }

    # Step 4: Apply entropy penalty to each conformer (dimensionless)
    # ΔG_corrected/RT = ΔG_original/RT + E_TS/RT
    corrected_energies = []
    for i, conf in enumerate(conformers):
        penalty = group_info[conf.charge_group]['entropy_penalty_dimensionless']
        corrected_energy = relative_energies[i] + penalty  # ADD penalty
        corrected_energies.append(corrected_energy)

    # Step 5: Calculate Boltzmann factors: exp(-ΔG/RT)
    boltzmann_factors = []
    for e in corrected_energies:
        if e > 100:  # Prevent overflow
            boltzmann_factors.append(0.0)
        else:
            # e is already in units of RT (dimensionless), so we compute exp(-e)
            # which is equivalent to exp(-ΔG/RT)
            boltzmann_factors.append(math.exp(-e))

    # Step 6: Calculate partition function and normalize to get probabilities
    partition = sum(boltzmann_factors)
    if partition < 1e-100:
        # All probabilities are essentially zero, return original
        return [c.probabilities[ph_index] for c in conformers], group_info

    # P_i = exp(-ΔG_i/RT) / Z
    new_probs = [bf / partition for bf in boltzmann_factors]

    return new_probs, group_info


def process_fort38(input_file: str, output_file: str, log_file: str, amino_acids: set, head3_file: str = 'head3.lst'):
    """Process fort.38 file and apply entropy corrections"""

    # Read actual charges from head3.lst
    conformer_charges = read_head3_charges(head3_file)

    # Read the file
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Parse header
    header = lines[0].strip().split()
    ph_values = [float(x) for x in header[1:]]
    num_ph = len(ph_values)

    print(f"Processing {len(lines)-1} conformers at {num_ph} pH values")

    # Parse all conformers
    all_conformers = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        name = parts[0]
        probs = [float(x) for x in parts[1:]]

        try:
            # Get actual charge from head3.lst if available
            actual_charge = conformer_charges.get(name, None)
            conf = Conformer(name, probs, actual_charge)
            all_conformers.append(conf)
        except ValueError as e:
            print(f"Warning: {e}")
            continue

    # Group by residue
    residue_groups = defaultdict(list)
    for conf in all_conformers:
        residue_groups[conf.full_res_id].append(conf)

    # Process each residue
    log_entries = []
    num_corrected = 0

    for res_id, res_conformers in residue_groups.items():
        res_type = res_conformers[0].res_type
        is_amino_acid = res_type in amino_acids

        if not is_amino_acid:
            num_corrected += 1
            print(f"Applying entropy correction to {res_id} (type: {res_type})")

            # Apply correction at each pH
            for ph_idx in range(num_ph):
                new_probs, group_info = apply_entropy_correction(res_conformers, ph_idx)

                # Update probabilities
                for i, conf in enumerate(res_conformers):
                    conf.probabilities[ph_idx] = new_probs[i]

                # Log significant corrections
                for charge_group, info in group_info.items():
                    if info['total_prob'] > 1e-10 and info['num_conformers'] > 1:
                        log_entries.append({
                            'residue': res_id,
                            'res_type': res_type,
                            'ph': ph_values[ph_idx],
                            'charge_group': charge_group,
                            'num_conformers': info['num_conformers'],
                            'entropy_penalty_dimensionless': info['entropy_penalty_dimensionless'],
                            'entropy_penalty_kcal': info['entropy_penalty_kcal']
                        })

    if amino_acids:
        print(f"\nCorrected {num_corrected} non-amino acid residues")
    else:
        print(f"\nCorrected {num_corrected} residues (including amino acids)")

    # Write output file - match exact fort.38 format
    with open(output_file, 'w') as f:
        # Write header
        f.write(' ph            ')
        f.write(''.join(f'{ph:>5.1f} ' for ph in ph_values))
        f.write('\n')

        # Write conformer data
        for conf in all_conformers:
            f.write(conf.full_name)
            f.write(' ')
            f.write(''.join(f'{p:5.3f} ' for p in conf.probabilities))
            f.write('\n')

    print(f"Wrote corrected file to: {output_file}")

    # Write log file
    with open(log_file, 'w') as f:
        f.write('SHANNON ENTROPY CORRECTION LOG\n')
        f.write('=' * 80 + '\n')
        if amino_acids:
            f.write('Applied to non-amino acid residues only\n\n')
        else:
            f.write('Applied to ALL residues (including amino acids)\n\n')

        f.write('METHOD:\n')
        f.write('-------\n')
        f.write('Corrects for bias when charge states have different numbers of conformers.\n')
        f.write('Conformers are grouped by their NUMERIC CHARGE VALUE (e.g., all neutral\n')
        f.write('states with charge=0 are in one group, regardless of state designation).\n')
        f.write('Charge states with more conformers are penalized using Shannon entropy\n')
        f.write('to account for the entropic cost of the conformational distribution.\n\n')

        f.write('PROCEDURE:\n')
        f.write('1. Group conformers by NUMERIC charge value (e.g., charge=0, charge=-1)\n')
        f.write('2. Convert probabilities to free energies: ΔG/RT = -ln(P_i/P_ref)\n')
        f.write('3. Normalize probabilities within each charge group\n')
        f.write('4. Calculate Shannon entropy for each group: E_TS/RT = -Σ(P_i * ln(P_i))\n')
        f.write('5. Add entropy penalty: ΔG_corrected/RT = ΔG_original/RT + E_TS/RT\n')
        f.write('6. Convert back to probabilities: P_i = exp(-ΔG_corrected/RT) / Z\n\n')

        f.write('CONSTANTS:\n')
        f.write(f'R = {R} kcal·K⁻¹·mol⁻¹\n')
        f.write(f'T = {T} K\n')
        f.write(f'RT = {RT:.3f} kcal/mol\n\n')

        f.write('PENALTIES (Shannon entropy):\n')
        f.write('N=1: Entropy = 0.000 → No penalty (0.000 kcal/mol)\n')
        f.write('N=2 (equal probs): Entropy = 0.693 → Penalty = 0.410 kcal/mol\n')
        f.write('N=3 (equal probs): Entropy = 1.099 → Penalty = 0.651 kcal/mol\n')
        f.write('N=4 (equal probs): Entropy = 1.386 → Penalty = 0.821 kcal/mol\n')
        f.write('Note: Actual penalties depend on the probability distribution within each group\n\n')

        f.write('Amino acids excluded from correction:\n')
        if amino_acids:
            f.write(', '.join(sorted(amino_acids)) + '\n\n')
        else:
            f.write('(None - all residues corrected)\n\n')

        f.write('=' * 80 + '\n')
        if amino_acids:
            f.write('NON-AMINO ACID RESIDUES PROCESSED:\n')
        else:
            f.write('ALL RESIDUES PROCESSED:\n')
        f.write('=' * 80 + '\n\n')

        # List all non-amino acid residues
        for res_id, res_conformers in sorted(residue_groups.items()):
            res_type = res_conformers[0].res_type
            if res_type not in amino_acids:
                f.write(f'Residue: {res_id} (Type: {res_type})\n')
                f.write(f'  Conformers:\n')
                for conf in res_conformers:
                    f.write(f'    - {conf.full_name} (charge group: {conf.charge_group})\n')
                f.write('\n')

        f.write('=' * 80 + '\n')
        f.write('ENTROPY PENALTIES APPLIED:\n')
        f.write('=' * 80 + '\n\n')

        if not log_entries:
            f.write('No penalties applied.\n')
            f.write('(All charge state groups had single conformers)\n\n')
        else:
            f.write(f'Total penalties applied: {len(log_entries)}\n\n')

            # Group log entries by residue
            by_residue = defaultdict(list)
            for entry in log_entries:
                by_residue[entry['residue']].append(entry)

            for res_id, entries in sorted(by_residue.items()):
                f.write(f"Residue: {res_id} ({entries[0]['res_type']})\n")
                for entry in entries:
                    f.write(f"  pH {entry['ph']:.1f}, Charge Group: {entry['charge_group']}\n")
                    f.write(f"    Number of conformers (N): {entry['num_conformers']}\n")
                    f.write(f"    Shannon entropy: {entry['entropy_penalty_dimensionless']:.4f}\n")
                    f.write(f"    Energy penalty: {entry['entropy_penalty_kcal']:.4f} kcal/mol\n\n")

        f.write('=' * 80 + '\n')
        f.write('SUMMARY:\n')
        f.write('=' * 80 + '\n\n')
        f.write(f"Total residues processed: {len(residue_groups)}\n")
        if amino_acids:
            f.write(f"Non-amino acid residues corrected: {num_corrected}\n")
            f.write(f"Amino acid residues (unchanged): {len(residue_groups) - num_corrected}\n")
        else:
            f.write(f"All residues corrected (including amino acids): {num_corrected}\n")
        f.write(f"Charge groups with penalties (N>1): {len(set(e['charge_group'] for e in log_entries))}\n")

    print(f"Wrote log file to: {log_file}")

    return all_conformers, ph_values, residue_groups


def generate_sum_crg(conformers: List[Conformer], ph_values: List[float],
                     residue_groups: Dict, output_file: str, original_sum_crg: str = 'sum_crg.out'):
    """Generate sum_crg.out file from corrected conformer probabilities"""

    print(f"\nGenerating {output_file}...")

    # Read original Protons and Electrons rows if file exists
    original_protons = None
    original_electrons = None
    original_entries = set()  # Track which entries were in the original file

    import os
    if os.path.exists(original_sum_crg):
        with open(original_sum_crg, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('Protons'):
                    original_protons = line.strip()
                elif line.startswith('Electrons'):
                    original_electrons = line.strip()
                elif line.strip() and not line.startswith('ph') and not line.startswith('-') and not line.startswith('Net_Charge'):
                    # Extract the residue identifier (first column)
                    parts = line.split()
                    if parts:
                        original_entries.add(parts[0])

    # Calculate net charges for each residue grouped by charge sign
    residue_charges = {}

    for res_id, res_conformers in residue_groups.items():
        res_type = res_conformers[0].res_type
        chain = res_conformers[0].chain
        res_num = res_conformers[0].res_num

        # Group conformers by charge sign (positive or negative only)
        charge_sign_groups = {'+': [], '-': []}

        for conf in res_conformers:
            # Use actual charge from head3.lst if available, otherwise parse from name
            if conf.actual_charge is not None:
                charge_val = conf.actual_charge
            else:
                charge_val = parse_charge_state(conf.charge_state)
            
            # Only consider charges with absolute value > 0.001 (skip neutrals)
            if charge_val > 0.001:
                charge_sign_groups['+'].append((conf, charge_val))
            elif charge_val < -0.001:
                charge_sign_groups['-'].append((conf, charge_val))
            # Skip neutral (close to 0) charge states

        # For each charge sign group that has conformers, create a line
        for sign in ['+', '-']:
            conf_list = charge_sign_groups[sign]
            if not conf_list:
                continue

            # Create identifier: RESTYPE + SIGN + CHAIN + RESNUM + '_'
            identifier = f"{res_type}{sign}{chain}{res_num}_"

            # Calculate charge contributions at each pH
            charges = []
            for ph_idx in range(len(ph_values)):
                total_charge = sum(charge_val * conf.probabilities[ph_idx]
                                  for conf, charge_val in conf_list)
                charges.append(total_charge)

            # Only include this entry if it was in the original sum_crg.out
            # OR if no original file exists, include all with significant charge
            if original_entries:
                if identifier in original_entries:
                    residue_charges[identifier] = charges
            else:
                # No original file - include entries with significant charge
                max_abs_charge = max(abs(c) for c in charges)
                if max_abs_charge > 0.005:
                    residue_charges[identifier] = charges

    # Sort residue identifiers by chain and residue number
    def sort_key(identifier):
        chain = identifier[-6]
        res_num = identifier[-5:-1]
        return (chain, int(res_num))

    sorted_residues = sorted(residue_charges.keys(), key=sort_key)

    # Calculate net charges
    net_charges = [0.0] * len(ph_values)

    for res_id in sorted_residues:
        for ph_idx in range(len(ph_values)):
            net_charges[ph_idx] += residue_charges[res_id][ph_idx]

    # Write output file
    with open(output_file, 'w') as f:
        # Write header
        f.write(' ph            ')
        f.write(''.join(f'{ph:>5.1f} ' for ph in ph_values))
        f.write('\n')

        # Write residue data
        for res_id in sorted_residues:
            f.write(f'{res_id:<15}')
            for charge in residue_charges[res_id]:
                f.write(f'{charge:5.2f} ')
            f.write('\n')

        # Write separator
        f.write('-' * (15 + len(ph_values) * 6) + '\n')

        # Write Net_Charge row
        f.write(f'{"Net_Charge":<15}')
        for charge in net_charges:
            f.write(f'{charge:5.2f} ')
        f.write('\n')

        # Write original Protons and Electrons rows
        if original_protons:
            f.write(original_protons + '\n')
        else:
            print("Warning: Could not find Protons row in original sum_crg.out")

        if original_electrons:
            f.write(original_electrons + '\n')
        else:
            print("Warning: Could not find Electrons row in original sum_crg.out")

    print(f"Wrote sum_crg file to: {output_file}")


if __name__ == '__main__':
    import sys
    import os

    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Apply Shannon entropy corrections to fort.38 file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
This script will:
  - Read the fort.38 file
  - Read actual charge values from head3.lst
  - Apply Shannon entropy corrections to non-amino acids (by default)
  - Penalize charge states with more conformers using Shannon entropy
  - Generate xts_fort.38 (corrected file)
  - Generate xts_sum_crg.out (corrected charge summary)
  - Generate entropy_correction.log (detailed log)

Examples:
  python xts_corr.py                      # Process only non-amino acids
  python xts_corr.py --all                # Process all residues including amino acids
  python xts_corr.py input.txt --all
  python xts_corr.py --head3 myhead3.lst  # Specify custom head3.lst path
        '''
    )

    parser.add_argument(
        'input_file',
        nargs='?',
        default='fort.38',
        help='Path to fort.38 file (default: fort.38)'
    )

    parser.add_argument(
        '--head3',
        default='head3.lst',
        help='Path to head3.lst file (default: head3.lst)'
    )

    parser.add_argument(
        '--all',
        action='store_true',
        help='Apply corrections to all residues including amino acids (default: only non-amino acids)'
    )

    args = parser.parse_args()
    input_file = args.input_file
    head3_file = args.head3

    # Set AMINO_ACIDS based on flag
    if args.all:
        AMINO_ACIDS = {}  # Empty set means no residues are excluded
        print("Mode: Correcting ALL residues (including amino acids)")
    else:
        AMINO_ACIDS = AMINO_ACIDS_LIST  # Use the full list to exclude amino acids
        print("Mode: Correcting non-amino acids only")

    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found!")
        print("\nPlease ensure fort.38 exists in the current directory,")
        print("or specify the correct path as an argument.")
        sys.exit(1)

    output_file = 'xts_fort.38'
    sum_crg_file = 'xts_sum_crg.out'
    log_file = 'entropy_correction.log'

    print("=" * 60)
    print("SHANNON ENTROPY CORRECTION FOR fort.38")
    print("=" * 60)
    print(f"Input file: {input_file}")
    print(f"Charge file: {head3_file}")
    print("\nMethod: Penalizing charge states using Shannon entropy")
    print("Grouping: By numeric charge value (e.g., all charge=0 together)")
    print("Penalty: ΔG += -RT × Σ(P_i × ln(P_i)) within each charge group")
    print()

    try:
        conformers, ph_values, residue_groups = process_fort38(input_file, output_file, log_file, AMINO_ACIDS, head3_file)
        generate_sum_crg(conformers, ph_values, residue_groups, sum_crg_file)

        print("\n" + "=" * 60)
        print("COMPLETE!")
        print("=" * 60)
        print(f"\nGenerated files:")
        print(f"  - {output_file} (corrected probabilities)")
        print(f"  - {sum_crg_file} (corrected charge summary)")
        print(f"  - {log_file} (detailed correction log)")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
