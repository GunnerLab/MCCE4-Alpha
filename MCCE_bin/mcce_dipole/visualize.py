"""
visualize.py - PyMOL visualization for MCCE4 dipole moments.

Chemistry convention dipole arrows:
  (+) end = cross/plus mark (tail) -- colored in + color
  (-) end = arrowhead (tip)        -- colored in - color

Each dipole arrow is two-toned: + color on one half, - color on the other.

Color scheme:
  backbone_dipole:   yellow (+) / dark gold (-)
  ionizable_dipole:  green (+) / magenta (-)
  full_dipole:       blue (+) / red (-)
  quadrupole:        cyan axes
"""

import os
import subprocess
import numpy as np
from . import E_ANG_TO_DEBYE

CGO_CYLINDER = 9.0

# Two-tone color scheme: (positive_color, negative_color)
DIPOLE_COLORS = {
    "backbone":   ((1.0, 0.85, 0.0),  (0.8, 0.55, 0.0)),   # yellow / dark gold
    "ionizable":  ((0.0, 0.8, 0.2),   (0.85, 0.2, 0.85)),  # green / magenta
    "full":       ((0.3, 0.5, 1.0),   (1.0, 0.2, 0.2)),    # blue / red
}

QUAD_COLOR = (0.0, 0.85, 1.0)  # cyan


def _normalize(v):
    mag = np.linalg.norm(v)
    return v / mag if mag > 1e-10 else np.zeros(3)


def _perp_vectors(d):
    """Two unit vectors perpendicular to d and to each other."""
    ref = np.array([1.0, 0.0, 0.0]) if abs(d[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(d, ref)
    u /= np.linalg.norm(u)
    v = np.cross(d, u)
    return u, v


def _cgo_dipole(center, direction, magnitude, color_pos, color_neg,
                shaft_radius=0.5, head_radius=1.4, head_fraction=0.12,
                cross_size=2.0, cross_radius=0.3,
                scale=0.1, min_length=5.0, n_taper=6, n_grad=8):
    """
    Gradient two-color dipole arrow with cross at (+) and arrowhead at (-).

    Color smoothly transitions from (+) color at the positive end to
    (-) color at the negative end. Each CYLINDER segment uses PyMOL's
    built-in two-color interpolation for smooth blending.

    Chemistry convention: arrow from (+) toward (-).
    mu vector points from (-) to (+), so:
      (+) end = center + direction * half   (cross mark here)
      (-) end = center - direction * half   (arrowhead here)
    """
    half = max(magnitude * scale, min_length) * 0.5
    head_len = half * head_fraction * 2

    pos_tip = center + direction * half            # (+) end
    neg_tip = center - direction * half            # (-) end
    head_start = center - direction * (half - head_len)  # arrowhead starts

    rp, gp, bp = color_pos
    rn, gn, bn = color_neg
    full_len = half * 2

    vals = []

    # Gradient shaft: n_grad segments from pos_tip to head_start
    # Color goes from color_pos at pos_tip to color_neg at head_start
    shaft_len = full_len - head_len
    for i in range(n_grad):
        # Fraction along the shaft (0 = pos_tip, 1 = head_start)
        t0 = i / n_grad
        t1 = (i + 1) / n_grad
        p0 = pos_tip - direction * (shaft_len * t0)
        p1 = pos_tip - direction * (shaft_len * t1)
        # Interpolate colors
        r0 = rp + (rn - rp) * t0
        g0 = gp + (gn - gp) * t0
        b0 = bp + (bn - bp) * t0
        r1 = rp + (rn - rp) * t1
        g1 = gp + (gn - gp) * t1
        b1 = bp + (bn - bp) * t1
        vals.extend([CGO_CYLINDER,
                     p0[0], p0[1], p0[2],
                     p1[0], p1[1], p1[2],
                     shaft_radius,
                     r0, g0, b0,
                     r1, g1, b1])

    # Arrowhead: tapered cylinders from head_start to neg_tip (full neg color)
    for i in range(n_taper):
        f0 = i / n_taper
        f1 = (i + 1) / n_taper
        s0 = head_start + (neg_tip - head_start) * f0
        s1 = head_start + (neg_tip - head_start) * f1
        r_seg = head_radius * (1.0 - (f0 + f1) * 0.5)
        if r_seg < 0.05:
            r_seg = 0.05
        vals.extend([CGO_CYLINDER,
                     s0[0], s0[1], s0[2],
                     s1[0], s1[1], s1[2],
                     r_seg, rn, gn, bn, rn, gn, bn])

    # (+) cross mark: two short perpendicular bars at pos_tip
    u, v = _perp_vectors(direction)
    for perp in [u, v]:
        bar_start = pos_tip - perp * cross_size * 0.5
        bar_end = pos_tip + perp * cross_size * 0.5
        vals.extend([CGO_CYLINDER,
                     bar_start[0], bar_start[1], bar_start[2],
                     bar_end[0], bar_end[1], bar_end[2],
                     cross_radius, rp, gp, bp, rp, gp, bp])

    return vals


def _format_cgo(vals):
    return "[" + ", ".join(f"{v:.3f}" for v in vals) + "]"


def generate_pymol_script(protein_pdb, results, output_path, ph_index=None,
                          arrow_scale=0.1):
    """
    Generate .pml with two-color dipole arrows.
    (+) end has cross mark, (-) end has arrowhead.
    """
    center = results["center"]

    if ph_index is not None:
        mu_bb = results["backbone_dipole"][ph_index]
        mu_ion = results["ionizable_dipole"][ph_index]
        mu_full = results["full_dipole"][ph_index]
        Q_eigvals = results["quadrupole_eigenvalues"][ph_index]
        Q_eigvecs = results.get("quadrupole_eigenvectors", None)
        if Q_eigvecs is None:
            Q = results["quadrupole_tensor"][ph_index]
            Q_eigvals, Q_eigvecs = np.linalg.eigh(Q)
        ph_label = f"pH {results['ph_values'][ph_index]:.1f}"
        net_q = results["net_charge"][ph_index]
    else:
        mu_bb = results["backbone_dipole"]
        mu_ion = results["ionizable_dipole"]
        mu_full = results["full_dipole"]
        Q_eigvals = results["quadrupole_eigenvalues"]
        Q_eigvecs = results["quadrupole_eigenvectors"]
        ph_label = "single state"
        net_q = results["net_charge"]

    mag_bb = np.linalg.norm(mu_bb)
    mag_ion = np.linalg.norm(mu_ion)
    mag_full = np.linalg.norm(mu_full)

    dir_bb = _normalize(mu_bb) if mag_bb > 0.1 else np.zeros(3)
    dir_ion = _normalize(mu_ion) if mag_ion > 0.1 else np.zeros(3)
    dir_full = _normalize(mu_full) if mag_full > 0.1 else np.zeros(3)

    L = []
    L.append(f"# MCCE4 Dipole Moment Visualization -- {ph_label}")
    L.append(f"# Generated by mcce_dipole v1.0")
    L.append(f"# Convention: cross (+) --|--======> arrowhead (-)")
    L.append("")
    L.append("from pymol import cmd")
    L.append("")

    # 1. Protein
    L.append(f'cmd.load("{protein_pdb}", "protein")')
    L.append('cmd.dss("protein")')
    L.append('cmd.hide("everything", "protein")')
    L.append('cmd.show("cartoon", "protein")')
    L.append('cmd.set("cartoon_fancy_helices", 1)')
    L.append('cmd.color("palegreen", "protein")')
    L.append('cmd.color("paleyellow", "protein and ss s")')
    L.append('cmd.set("cartoon_transparency", 0.15)')
    L.append('cmd.bg_color("black")')
    L.append("")

    # Center marker
    L.append(f'cmd.pseudoatom("center_marker", '
             f'pos=[{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}])')
    L.append('cmd.show("spheres", "center_marker")')
    L.append('cmd.set("sphere_scale", 0.5, "center_marker")')
    L.append('cmd.color("white", "center_marker")')
    L.append("")

    # Dipole arrows
    dipoles = [
        ("backbone_dipole",  mu_bb,  mag_bb,  dir_bb,
         DIPOLE_COLORS["backbone"],  "yellow(+) / gold(-)"),
        ("ionizable_dipole", mu_ion, mag_ion, dir_ion,
         DIPOLE_COLORS["ionizable"], "green(+) / magenta(-)"),
        ("full_dipole",      mu_full, mag_full, dir_full,
         DIPOLE_COLORS["full"],      "blue(+) / red(-)"),
    ]

    for obj_name, mu_vec, mag, direction, (c_pos, c_neg), color_desc in dipoles:
        if mag < 0.1:
            L.append(f"# {obj_name}: too small ({mag:.2f} D), skipped")
            continue

        cgo = _cgo_dipole(center, direction, mag, c_pos, c_neg,
                          scale=arrow_scale)

        L.append(f"# {obj_name}: {mag:.1f} D  {color_desc}")
        L.append(f"cmd.load_cgo({_format_cgo(cgo)}, '{obj_name}')")
        L.append("")

    # Quadrupole (cyan, 3 axes as single object)
    if Q_eigvecs is not None:
        r, g, b = QUAD_COLOR
        quad_cgo = []
        for i in range(3):
            eigval = Q_eigvals[i]
            eigvec = Q_eigvecs[:, i]
            axis_len = max(abs(eigval) * 0.02, 3.0)
            dv = _normalize(eigvec)
            start = center - dv * axis_len
            end = center + dv * axis_len
            quad_cgo.extend([CGO_CYLINDER,
                             start[0], start[1], start[2],
                             end[0], end[1], end[2],
                             0.3, r, g, b, r, g, b])
        L.append(f"# quadrupole: 3 principal axes (cyan)")
        L.append(f"cmd.load_cgo({_format_cgo(quad_cgo)}, 'quadrupole')")
        L.append("")

    # Rendering
    L.append('cmd.set("ray_shadow", "on")')
    L.append('cmd.set("antialias", 2)')
    L.append('cmd.set("specular", 0.3)')
    L.append('cmd.center("protein")')
    L.append('cmd.zoom("protein", buffer=10)')
    L.append("")

    # Summary with coordinates
    L.append(f'print("")')
    L.append(f'print("  MCCE4 Dipole Visualization -- {ph_label}")')
    L.append(f'print("  Convention: cross (+) --|--======> arrow (-)")')
    L.append(f'print("  -----------------------------------------")')
    L.append(f'print("  Geometric center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")')
    L.append(f'print("")')
    L.append(f'print("  backbone_dipole:  |mu| = {mag_bb:7.1f} D")')
    L.append(f'print("    yellow (+) / gold (-)")')
    L.append(f'print("    mu  = ({mu_bb[0]:.2f}, {mu_bb[1]:.2f}, {mu_bb[2]:.2f}) D")')
    L.append(f'print("    dir = ({dir_bb[0]:.4f}, {dir_bb[1]:.4f}, {dir_bb[2]:.4f})")')
    L.append(f'print("")')
    L.append(f'print("  ionizable_dipole: |mu| = {mag_ion:7.1f} D")')
    L.append(f'print("    green (+) / magenta (-)")')
    L.append(f'print("    mu  = ({mu_ion[0]:.2f}, {mu_ion[1]:.2f}, {mu_ion[2]:.2f}) D")')
    L.append(f'print("    dir = ({dir_ion[0]:.4f}, {dir_ion[1]:.4f}, {dir_ion[2]:.4f})")')
    L.append(f'print("")')
    L.append(f'print("  full_dipole:      |mu| = {mag_full:7.1f} D")')
    L.append(f'print("    blue (+) / red (-)")')
    L.append(f'print("    mu  = ({mu_full[0]:.2f}, {mu_full[1]:.2f}, {mu_full[2]:.2f}) D")')
    L.append(f'print("    dir = ({dir_full[0]:.4f}, {dir_full[1]:.4f}, {dir_full[2]:.4f})")')
    L.append(f'print("")')
    L.append(f'print("  quadrupole:       3 axes (cyan)")')
    L.append(f'print("    eigenvalues = ({Q_eigvals[0]:.1f}, {Q_eigvals[1]:.1f}, {Q_eigvals[2]:.1f}) e*A^2")')
    L.append(f'print("")')
    L.append(f'print("  Net charge: {net_q:.2f} e")')
    L.append(f'print("")')

    # Toggle commands
    L.append("try:")
    L.append("    def _solo(keep):")
    L.append('        for o in ["backbone_dipole","ionizable_dipole","full_dipole","quadrupole"]:')
    L.append("            cmd.disable(o) if o != keep else cmd.enable(o)")
    L.append('    def show_backbone(): cmd.enable("backbone_dipole")')
    L.append('    def hide_backbone(): cmd.disable("backbone_dipole")')
    L.append('    def show_ionizable(): cmd.enable("ionizable_dipole")')
    L.append('    def hide_ionizable(): cmd.disable("ionizable_dipole")')
    L.append('    def show_full(): cmd.enable("full_dipole")')
    L.append('    def hide_full(): cmd.disable("full_dipole")')
    L.append('    def show_quadrupole(): cmd.enable("quadrupole")')
    L.append('    def hide_quadrupole(): cmd.disable("quadrupole")')
    L.append('    def show_all_dipoles():')
    L.append('        for o in ["backbone_dipole","ionizable_dipole","full_dipole","quadrupole"]: cmd.enable(o)')
    L.append('    def hide_all_dipoles():')
    L.append('        for o in ["backbone_dipole","ionizable_dipole","full_dipole","quadrupole"]: cmd.disable(o)')
    L.append('    def solo_backbone(): _solo("backbone_dipole")')
    L.append('    def solo_ionizable(): _solo("ionizable_dipole")')
    L.append('    def solo_full(): _solo("full_dipole")')
    L.append('    def solo_quadrupole(): _solo("quadrupole")')
    L.append('    for fn in [show_backbone, hide_backbone, show_ionizable, hide_ionizable,')
    L.append('               show_full, hide_full, show_quadrupole, hide_quadrupole,')
    L.append('               show_all_dipoles, hide_all_dipoles,')
    L.append('               solo_backbone, solo_ionizable, solo_full, solo_quadrupole]:')
    L.append('        cmd.extend(fn.__name__, fn)')
    L.append("except:")
    L.append("    pass")
    L.append("")

    with open(output_path, "w") as fh:
        fh.write("\n".join(L) + "\n")

    # Generate .pse by running PyMOL headlessly
    pse_path = output_path.replace(".pml", ".pse")
    abs_pml = os.path.abspath(output_path)
    abs_pse = os.path.abspath(pse_path)
    try:
        subprocess.run(
            ["pymol", "-cq", abs_pml, "-d", f'cmd.save("{abs_pse}")'],
            capture_output=True, timeout=60,
        )
    except FileNotFoundError:
        pass  # pymol not installed; .pml still usable
    except subprocess.TimeoutExpired:
        pass

    return output_path


def generate_ph_scan_csv(results, output_path):
    """Write pH scan CSV with magnitudes, vector components, and endpoints."""
    ph = results["ph_values"]
    center = results["center"]
    with open(output_path, "w") as fh:
        fh.write("pH,backbone_D,ionizable_D,full_D,net_charge,"
                 "bb_x,bb_y,bb_z,ion_x,ion_y,ion_z,full_x,full_y,full_z,"
                 "bb_neg_x,bb_neg_y,bb_neg_z,bb_pos_x,bb_pos_y,bb_pos_z,"
                 "ion_neg_x,ion_neg_y,ion_neg_z,ion_pos_x,ion_pos_y,ion_pos_z,"
                 "full_neg_x,full_neg_y,full_neg_z,full_pos_x,full_pos_y,full_pos_z,"
                 "Q_eig1,Q_eig2,Q_eig3\n")
        for i in range(len(ph)):
            mag_bb = np.linalg.norm(results["backbone_dipole"][i])
            mag_ion = np.linalg.norm(results["ionizable_dipole"][i])
            mag_full = np.linalg.norm(results["full_dipole"][i])
            bb = results["backbone_dipole"][i]
            ion = results["ionizable_dipole"][i]
            full = results["full_dipole"][i]
            nq = results["net_charge"][i]
            qe = results["quadrupole_eigenvalues"][i]

            # Compute negative and positive endpoints for each dipole
            endpoints = []
            for mu, mag in [(bb, mag_bb), (ion, mag_ion), (full, mag_full)]:
                if mag > 0.1:
                    d = mu / mag
                    neg = center - d * mag * 0.1
                    pos = center + d * mag * 0.1
                else:
                    neg = center.copy()
                    pos = center.copy()
                endpoints.extend([neg, pos])

            fh.write(f"{ph[i]:.1f},{mag_bb:.2f},{mag_ion:.2f},{mag_full:.2f},"
                     f"{nq:.3f},"
                     f"{bb[0]:.2f},{bb[1]:.2f},{bb[2]:.2f},"
                     f"{ion[0]:.2f},{ion[1]:.2f},{ion[2]:.2f},"
                     f"{full[0]:.2f},{full[1]:.2f},{full[2]:.2f},"
                     f"{endpoints[0][0]:.2f},{endpoints[0][1]:.2f},{endpoints[0][2]:.2f},"
                     f"{endpoints[1][0]:.2f},{endpoints[1][1]:.2f},{endpoints[1][2]:.2f},"
                     f"{endpoints[2][0]:.2f},{endpoints[2][1]:.2f},{endpoints[2][2]:.2f},"
                     f"{endpoints[3][0]:.2f},{endpoints[3][1]:.2f},{endpoints[3][2]:.2f},"
                     f"{endpoints[4][0]:.2f},{endpoints[4][1]:.2f},{endpoints[4][2]:.2f},"
                     f"{endpoints[5][0]:.2f},{endpoints[5][1]:.2f},{endpoints[5][2]:.2f},"
                     f"{qe[0]:.2f},{qe[1]:.2f},{qe[2]:.2f}\n")
    return output_path
