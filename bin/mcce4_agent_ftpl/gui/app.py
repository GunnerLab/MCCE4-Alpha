"""
Streamlit Web GUI for the MCCE4 Topology Agent.

Features:
  - 7-phase progress tracker in sidebar
  - Persistent scrollable log panel with protonation state highlighting
  - 2D molecule depiction with protonation site highlighting
  - Conformer state table with add/edit/remove
  - Molecule editor via Ketcher (if streamlit-ketcher installed)
  - Runs in browser — works over SSH tunnel, no X11 needed

Launch:
  streamlit run mcce4_agent_ftpl/gui/app.py -- EMH.pdb
  # Or via the CLI:
  mcce4_agent_ftpl.py EMH.pdb --gui
"""

import os
import sys
import json
import base64
import logging
import tempfile
import io
from pathlib import Path

# Ensure package is importable
PACKAGE_DIR = str(Path(__file__).resolve().parent.parent.parent)
if PACKAGE_DIR not in sys.path:
    sys.path.insert(0, PACKAGE_DIR)

import streamlit as st

from mcce4_agent_ftpl.models import ConformerState, AgentState
from mcce4_agent_ftpl.config import (
    GUI_TITLE, SUPPORTED_CHARGE_METHODS, DEFAULT_CHARGE_METHOD, DEFAULT_DIELECTRICS,
    AGENT_PHASES, PHASE_NAME_MAP,
)
from mcce4_agent_ftpl.tools.rdkit_tools import mol_to_svg, get_mol_from_pdb, get_mol_from_smiles
from mcce4_agent_ftpl.tools.mcce_tools import extract_lig_id_from_pdb
import html as _html


# ──────────────────────────────────────────────────────────────────────────────
# Log highlighting
# ──────────────────────────────────────────────────────────────────────────────

def _highlight_log(log_text: str) -> str:
    """Add HTML color highlighting to log text for protonation changes.

    Green: added H, protonation events, success markers
    Red: removed H, deprotonation events, errors
    Yellow: warnings
    Cyan: phase headers, calibration info
    Magenta: dsolv/rxn values
    Bold: state labels like EMH01, EMH+1, EMH+a
    """
    lines = []
    for line in _html.escape(log_text).splitlines():
        # Phase headers — cyan bold
        if any(kw in line for kw in ["PHASE 1", "PHASE 2", "PHASE 3", "PHASE 4",
                                      "PHASE 5", "PHASE 6", "PHASE 7"]):
            line = f'<span style="color:#4fc1ff;font-weight:bold;">{line}</span>'
        elif any(kw in line for kw in ["PHASE", "═", "─", "🔬", "🎉", "✅"]):
            line = f'<span style="color:#4fc1ff;">{line}</span>'
        # Protonation additions — green
        elif any(kw in line for kw in ["+H ", "Added H", "added on", "protonate",
                                        "🟢", "h_added", "Add H:", "Added "]):
            line = f'<span style="color:#4ec94e;font-weight:bold;">▸ {line}</span>'
        # Deprotonation removals — red
        elif any(kw in line for kw in ["-H ", "Removed H", "removed from", "deprotonate",
                                        "🔴", "h_removed", "Remove H:", "Removed "]):
            line = f'<span style="color:#f44747;font-weight:bold;">▸ {line}</span>'
        # State label lines — highlight the label
        elif any(kw in line for kw in ["State +", "State -", "State 0",
                                        "state +", "state -", "state 0"]):
            line = f'<span style="color:#dcdcaa;">{line}</span>'
        # Naming/label lines
        elif any(kw in line for kw in ["Label disambiguated", "→", "label=",
                                        "+a", "+b", "-a", "-b", "0a", "0b"]):
            line = f'<span style="color:#ce9178;">{line}</span>'
        # Warnings
        elif any(kw in line for kw in ["⚠", "WARNING", "warning"]):
            line = f'<span style="color:#cca700;">{line}</span>'
        # Errors
        elif any(kw in line for kw in ["❌", "ERROR", "failed", "FAILED"]):
            line = f'<span style="color:#f44747;font-weight:bold;">{line}</span>'
        # dsolv/rxn calibration values — magenta
        elif any(kw in line for kw in ["📊", "dsolv", "rxn0", "rxn_"]):
            line = f'<span style="color:#c586c0;">{line}</span>'
        # Charge info
        elif any(kw in line for kw in ["⚡", "Charges", "charges", "charge="]):
            line = f'<span style="color:#9cdcfe;">{line}</span>'
        # Conformer labels in output (EMH01, EMH+1, etc.)
        elif any(kw in line for kw in ["Conformer", "conformer", "CONFLIST"]):
            line = f'<span style="color:#dcdcaa;">{line}</span>'
        lines.append(line)
    return "\n".join(lines)


def _render_log_panel(log_text: str, title: str = "Agent Log",
                      max_height: int = 500, key_suffix: str = ""):
    """Render a scrollable, syntax-highlighted log panel."""
    if not log_text or not log_text.strip():
        return
    st.markdown(f"**📋 {title}** _(scroll to review all output)_")
    st.markdown(
        f'<div id="log-panel-{key_suffix}" style="max-height:{max_height}px; '
        f'overflow-y:auto; background:#1e1e1e; color:#d4d4d4; padding:12px; '
        f'border-radius:6px; font-family:\'Consolas\',\'Courier New\',monospace; '
        f'font-size:13px; line-height:1.5; white-space:pre-wrap; '
        f'border:1px solid #333;">{_highlight_log(log_text)}</div>',
        unsafe_allow_html=True,
    )


# ──────────────────────────────────────────────────────────────────────────────
# Log capture helper
# ──────────────────────────────────────────────────────────────────────────────

def _start_log_capture():
    """Start capturing log output. Returns (StringIO, handler)."""
    log_capture = io.StringIO()
    log_handler = logging.StreamHandler(log_capture)
    log_handler.setLevel(logging.INFO)
    log_handler.setFormatter(logging.Formatter("%(message)s"))
    logging.getLogger().addHandler(log_handler)
    return log_capture, log_handler


def _stop_log_capture(log_capture, log_handler):
    """Stop capturing and return text."""
    logging.getLogger().removeHandler(log_handler)
    return log_capture.getvalue()


# ──────────────────────────────────────────────────────────────────────────────
# Page config
# ──────────────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title=GUI_TITLE,
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ──────────────────────────────────────────────────────────────────────────────
# Session state init
# ──────────────────────────────────────────────────────────────────────────────
def init_session():
    """Initialize session state variables."""
    defaults = {
        "agent_state": None,
        "pdb_path": None,
        "lig_id": None,
        "states": [],
        "phase": "upload",
        "approved": False,
        "running": False,
        "analysis_done": False,
        "log_output": "",
        "phase_status": {},  # {phase_id: "pending"|"running"|"done"|"error"}
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


init_session()


# ──────────────────────────────────────────────────────────────────────────────
# Sidebar: Settings + 7-Phase Progress Tracker
# ──────────────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.title("🤖 MCCE4 Agent")
    st.markdown("---")

    # Settings
    st.subheader("⚙️ Settings")
    ph = st.slider("Target pH", 0.0, 14.0, 7.4, 0.1)
    charge_method = st.selectbox("Charge Method", SUPPORTED_CHARGE_METHODS,
                                  index=SUPPORTED_CHARGE_METHODS.index(DEFAULT_CHARGE_METHOD))

    dielectric_opts = st.multiselect("Dielectric Constants", [2, 4, 8], default=[2, 4, 8])
    dry_run = st.checkbox("Dry Run (skip RXN calibration)", value=False)

    st.markdown("---")

    # ── 7-Phase Progress Tracker ──
    st.subheader("📊 Pipeline Progress")

    phase_status = st.session_state.get("phase_status", {})
    agent_state = st.session_state.get("agent_state")
    current_internal_phase = ""
    if agent_state:
        current_internal_phase = agent_state.get("phase", "")

    # Determine which phase ID is current
    current_phase_id = PHASE_NAME_MAP.get(current_internal_phase)

    for p in AGENT_PHASES:
        pid = p["id"]
        status = phase_status.get(pid, "pending")

        if status == "done":
            icon = "✅"
            color = "#4ec94e"
        elif status == "running":
            icon = "🔄"
            color = "#4fc1ff"
        elif status == "error":
            icon = "❌"
            color = "#f44747"
        else:
            icon = "⬜"
            color = "#666"

        st.markdown(
            f'<div style="padding:3px 0; color:{color}; font-size:13px;">'
            f'{icon} <b>Phase {p["num"]}</b>: {p["name"]}'
            f'</div>',
            unsafe_allow_html=True,
        )

    st.markdown("---")

    # Naming convention legend
    st.subheader("🏷️ Naming Convention")
    st.markdown("""
    <div style="font-size:12px; line-height:1.6;">
    <b>Single state per charge:</b><br>
    &nbsp;&nbsp;<code>01</code> = neutral<br>
    &nbsp;&nbsp;<code>+1</code> = protonated (+1)<br>
    &nbsp;&nbsp;<code>-1</code> = deprotonated (-1)<br>
    <br>
    <b>Multiple states at same charge:</b><br>
    &nbsp;&nbsp;<code>+a, +b</code> = two +1 states<br>
    &nbsp;&nbsp;<code>-a, -b</code> = two -1 states<br>
    &nbsp;&nbsp;<code>0a, 0b</code> = two neutral states<br>
    <br>
    Names are disambiguated based on<br>
    which atom site is protonated.
    </div>
    """, unsafe_allow_html=True)

    st.markdown("---")

    # Ligand info
    lig_id = st.session_state.get("lig_id")
    if lig_id:
        st.caption(f"Ligand: **{lig_id}**")
    states = st.session_state.get("states", [])
    if states:
        st.caption(f"{len(states)} conformer state(s)")


# ──────────────────────────────────────────────────────────────────────────────
# Main content
# ──────────────────────────────────────────────────────────────────────────────
st.title("🧬 MCCE4 Topology File Agent")
st.markdown("Create MCCE4 topology files (.ftpl) with AI-powered protonation state analysis")

# ── Tab layout: Input | Conformer States | Output | Log ──
tab_input, tab_states, tab_output, tab_log = st.tabs([
    "📂 Input", "🧪 Conformer States & Editor", "📄 Output", "📋 Log"
])


# ──────────────────────────────────────────────────────────────────────────────
# TAB 1: Input
# ──────────────────────────────────────────────────────────────────────────────
with tab_input:
    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader("Upload Ligand PDB")

        uploaded_pdb = st.file_uploader("Ligand PDB file", type=["pdb"],
                                         key="pdb_upload")
        if uploaded_pdb:
            # Save to temp file
            tmp_dir = tempfile.mkdtemp()
            pdb_path = os.path.join(tmp_dir, uploaded_pdb.name)
            with open(pdb_path, "wb") as f:
                f.write(uploaded_pdb.getvalue())
            st.session_state["pdb_path"] = pdb_path
            st.session_state["lig_id"] = extract_lig_id_from_pdb(pdb_path)
            st.success(f"✓ Loaded: {uploaded_pdb.name} (Ligand: {st.session_state['lig_id']})")

        # Or specify path directly
        pdb_path_input = st.text_input("Or enter PDB file path on server:",
                                        placeholder="/path/to/EMH.pdb")
        if pdb_path_input and os.path.exists(pdb_path_input):
            st.session_state["pdb_path"] = pdb_path_input
            st.session_state["lig_id"] = extract_lig_id_from_pdb(pdb_path_input)
            st.success(f"✓ Found: {pdb_path_input} (Ligand: {st.session_state['lig_id']})")

        st.markdown("---")
        st.subheader("Optional: State-specific PDBs")
        st.caption("Upload PDB files for specific protonation states (e.g., EMH_01.pdb, EMH_+1.pdb)")
        state_pdbs = st.file_uploader("State PDB files", type=["pdb"],
                                       accept_multiple_files=True, key="state_pdbs")

    with col2:
        st.subheader("Molecule Preview")
        if st.session_state.get("pdb_path"):
            svg = mol_to_svg(st.session_state["pdb_path"], size=(500, 400))
            if svg:
                st.image(svg, use_container_width=True)
            else:
                st.info("RDKit not available for 2D preview. Install: `pip install rdkit`")
        else:
            st.info("Upload or specify a PDB file to see the molecule.")

    # ── Analyze button ──
    if st.session_state.get("pdb_path"):
        if st.button("🧠 Analyze Protonation States", type="primary", use_container_width=True):
            st.session_state["phase"] = "analyzing"
            # Reset phase status
            st.session_state["phase_status"] = {}

            from mcce4_agent_ftpl.agent import run_agent

            # Show real-time progress below the button
            status = st.status("🤖 Agent is analyzing the molecule...", expanded=True)

            log_capture, log_handler = _start_log_capture()

            with status:
                # Phase 1
                st.session_state["phase_status"]["phase1"] = "running"
                st.write("**Phase 1:** Molecule Intelligence — fetching info from RCSB...")

                agent_state = run_agent(
                    st.session_state["pdb_path"],
                    use_gui=True,
                    charge_method=charge_method,
                    dielectrics=dielectric_opts,
                    ph=ph,
                    dry_run=dry_run,
                )

                # Update phase statuses based on what completed
                for pid in ["phase1", "phase2", "phase3", "phase4"]:
                    st.session_state["phase_status"][pid] = "done"

                # Show what happened
                if agent_state.get("smiles"):
                    smi = agent_state.get('smiles', '')
                    display_smi = f"{smi[:80]}..." if len(smi) > 80 else smi
                    st.write(f"✅ **SMILES:** `{display_smi}`")
                if agent_state.get("name"):
                    st.write(f"✅ **Name:** {agent_state['name']}")

                st.write("**Phase 2:** Protonation State Enumeration")
                states = agent_state.get("states", [])
                if states:
                    st.write(f"✅ Found **{len(states)}** conformer state(s):")
                    for s in states:
                        label = s.get('label', '?')
                        charge = int(s.get('charge', 0) or 0)
                        pka = s.get('pka')
                        pka_str = ""
                        if pka is not None:
                            try:
                                pka_clean = str(pka).lstrip("~≈><≥≤ ")
                                pka_str = f", pKa≈{float(pka_clean):.1f}"
                            except (ValueError, TypeError):
                                pka_str = f", pKa={pka}"
                        source = s.get('source', '?')

                        # Color-code by charge
                        if charge > 0:
                            badge = f'<span style="color:#4ec94e;font-weight:bold;">⊕ {label}</span>'
                        elif charge < 0:
                            badge = f'<span style="color:#f44747;font-weight:bold;">⊖ {label}</span>'
                        else:
                            badge = f'<span style="color:#4fc1ff;font-weight:bold;">◯ {label}</span>'

                        proton_ex = s.get('proton_exchange', '')
                        ex_str = f" — _{proton_ex}_" if proton_ex else ""
                        st.markdown(
                            f"&nbsp;&nbsp; {badge} (charge={charge:+d}{pka_str}) "
                            f"— {source}{ex_str}",
                            unsafe_allow_html=True,
                        )

                # Phase 4: show per-state PDB generation results
                state_pdbs_result = agent_state.get("state_pdbs", {})
                if state_pdbs_result:
                    st.write(f"**Phase 4:** Generated **{len(state_pdbs_result)}** per-state PDB(s)")
                    h_diffs = agent_state.get("h_diffs", {})
                    for label, pdb_path in state_pdbs_result.items():
                        diff = h_diffs.get(label, {})
                        added = diff.get("added", [])
                        removed = diff.get("removed", [])
                        diff_parts = []
                        if added:
                            added_on = diff.get("added_on", {})
                            for h in added:
                                parent = added_on.get(h, "?")
                                diff_parts.append(
                                    f'<span style="color:#4ec94e;">+{h} on {parent}</span>'
                                )
                        if removed:
                            removed_from = diff.get("removed_from", {})
                            for h in removed:
                                parent = removed_from.get(h, "?")
                                diff_parts.append(
                                    f'<span style="color:#f44747;">-{h} from {parent}</span>'
                                )
                        diff_html = ", ".join(diff_parts) if diff_parts else ""
                        st.markdown(
                            f"&nbsp;&nbsp; **{label}**: `{os.path.basename(pdb_path)}` "
                            f"{diff_html}",
                            unsafe_allow_html=True,
                        )

                if agent_state.get("warnings"):
                    st.write("**⚠️ Warnings:**")
                    for w in agent_state["warnings"]:
                        st.warning(w)

                # Show captured log in scrollable container
                log_text = _stop_log_capture(log_capture, log_handler)
                if log_text.strip():
                    _render_log_panel(log_text, "Analysis Log", 400, "analysis")
                    st.session_state["log_output"] = log_text

            status.update(label="✅ Analysis complete — switch to **Conformer States** tab to review",
                          state="complete", expanded=True)

            st.session_state["agent_state"] = agent_state
            st.session_state["states"] = agent_state.get("states", [])
            st.session_state["phase"] = "review"
            st.session_state["analysis_done"] = True

            # Rerun so other tabs pick up the new state
            st.rerun()

    # Show success banner (persists after rerun)
    if st.session_state.get("analysis_done"):
        st.success("🎉 Analysis complete! Click the **Conformer States** tab above to review and approve.")


# ──────────────────────────────────────────────────────────────────────────────
# TAB 2: Conformer States & Editor (v3)
# ──────────────────────────────────────────────────────────────────────────────
with tab_states:
    states = st.session_state.get("states", [])
    agent_state = st.session_state.get("agent_state")

    if not states:
        st.info("Run analysis first to see proposed conformer states.")
    else:
        lig_id = st.session_state.get("lig_id", "?")
        h_diffs = agent_state.get("h_diffs", {}) if agent_state else {}
        state_pdbs = agent_state.get("state_pdbs", {}) if agent_state else {}

        st.subheader(f"Protonation States for {lig_id}")

        # ══════════════════════════════════════════════════════════════════
        # SECTION 1: Hydrogen comparison table + per-state focused views
        # ══════════════════════════════════════════════════════════════════

        # 1a. Hydrogen comparison table
        if h_diffs:
            st.markdown("#### Hydrogen Differences vs Neutral")
            st.caption(
                "🟢 Green = H added (protonation) | "
                "🔴 Red = H removed (deprotonation) | "
                "— = no change (neutral reference)"
            )
            try:
                from mcce4_agent_ftpl.tools.rdkit_tools import generate_state_comparison_table
                comp_rows = generate_state_comparison_table(states, h_diffs, lig_id)
                if comp_rows:
                    import pandas as pd
                    comp_df = pd.DataFrame(comp_rows)
                    st.dataframe(comp_df, use_container_width=True, hide_index=True)
            except Exception as e:
                logging.warning(f"Comparison table failed: {e}")

        # 1b. Per-state views side by side
        if state_pdbs and len(state_pdbs) > 0:
            st.markdown("---")
            st.markdown(
                "#### Per-State Structures\n"
                "All atoms labeled with PDB names. "
                "🟢 **Green** = H atom added vs neutral. "
                "🔴 **Red** = H atom removed vs neutral."
            )

            cols_per_row = min(3, len(states))
            for row_start in range(0, len(states), cols_per_row):
                row_cols = st.columns(cols_per_row)
                for i, col in enumerate(row_cols):
                    idx = row_start + i
                    if idx >= len(states):
                        break
                    s = states[idx]
                    sd = s if isinstance(s, dict) else s.to_dict()
                    label = sd.get("label", "?")
                    charge = int(sd.get("charge", 0) or 0)
                    diff = h_diffs.get(label, {})
                    h_added = diff.get("added", []) or sd.get("h_added", [])
                    h_removed = diff.get("removed", []) or sd.get("h_removed", [])

                    with col:
                        # Color-coded header
                        if charge > 0:
                            hdr_color = "#4ec94e"
                            charge_icon = "⊕"
                        elif charge < 0:
                            hdr_color = "#f44747"
                            charge_icon = "⊖"
                        else:
                            hdr_color = "#4fc1ff"
                            charge_icon = "◯"

                        st.markdown(
                            f'<div style="background:{hdr_color}22; border-left:4px solid {hdr_color}; '
                            f'padding:6px 10px; border-radius:4px; margin-bottom:8px;">'
                            f'<b>{charge_icon} {lig_id}{label}</b> (charge: {charge:+d})'
                            f'</div>',
                            unsafe_allow_html=True,
                        )

                        # Determine which PDB to render
                        pdb_for_state = state_pdbs.get(label) or sd.get("pdb_path")
                        if not pdb_for_state or not os.path.exists(str(pdb_for_state)):
                            pdb_for_state = st.session_state.get("pdb_path")

                        rendered = False
                        if pdb_for_state and os.path.exists(str(pdb_for_state)):
                            try:
                                from mcce4_agent_ftpl.tools.rdkit_tools import (
                                    mol_to_svg_with_h_diff,
                                )
                                svg = mol_to_svg_with_h_diff(
                                    pdb_for_state, h_added, h_removed,
                                    size=(550, 450),
                                )
                                if svg:
                                    st.image(svg, use_container_width=True)
                                    rendered = True
                            except Exception as e:
                                st.warning(f"Render: {e}")

                        if not rendered:
                            st.info(f"No structure available for {label}")

                        # H-diff annotations below each molecule
                        if h_added:
                            added_on = diff.get("added_on", {})
                            for h in h_added:
                                parent = added_on.get(h, "?")
                                st.markdown(f"🟢 **{h}** added on **{parent}**")
                        if h_removed:
                            removed_from = diff.get("removed_from", {})
                            for h in h_removed:
                                parent = removed_from.get(h, "?")
                                st.markdown(f"🔴 **{h}** removed from **{parent}**")
                        if label in ("01", "00"):
                            st.caption("_(neutral reference)_")
                        proton_ex = sd.get("proton_exchange", "")
                        if proton_ex:
                            st.caption(f"Exchange: {proton_ex}")

        else:
            # No state_pdbs — render neutral PDB for all states
            neutral_pdb = st.session_state.get("pdb_path")
            if neutral_pdb:
                svg = mol_to_svg(neutral_pdb, size=(450, 350))
                if svg:
                    st.image(svg, use_container_width=True)
                if len(states) > 1:
                    state_cols = st.columns(min(len(states), 3))
                    for i, s in enumerate(states):
                        with state_cols[i % 3]:
                            smi = s.get("smiles", "") if isinstance(s, dict) else s.smiles
                            if smi:
                                mol = get_mol_from_smiles(smi)
                                if mol:
                                    svg_s = mol_to_svg(mol, size=(250, 200))
                                    if svg_s:
                                        lbl = s.get("label", "?") if isinstance(s, dict) else s.label
                                        chg = s.get("charge", 0) if isinstance(s, dict) else s.charge
                                        st.image(svg_s, caption=f"{lbl} ({chg:+d})",
                                                 use_container_width=True)

        # ══════════════════════════════════════════════════════════════════
        # SECTION 1b: PyMOL visualization (full 3D with all H atoms)
        # ══════════════════════════════════════════════════════════════════
        if state_pdbs and len(state_pdbs) > 1:
            st.markdown("---")
            st.subheader("🔬 3D Visualization (PyMOL)")
            st.markdown(
                "For **full control** over hydrogen visualization, use PyMOL "
                "with the generated per-state PDB files. Each PDB has the "
                "correct H atoms for that protonation state."
            )

            pymol_script = agent_state.get("pymol_script") if agent_state else None

            col_pymol, col_files = st.columns([1, 1])

            with col_pymol:
                if pymol_script and os.path.exists(pymol_script):
                    with open(pymol_script) as f:
                        pml_content = f.read()
                    st.download_button(
                        "📥 Download PyMOL script (.pml)",
                        data=pml_content,
                        file_name=os.path.basename(pymol_script),
                        mime="text/plain",
                    )
                    st.caption("Run: `pymol " + os.path.basename(pymol_script) + "`")

                    with st.expander("Preview PyMOL script"):
                        st.code(pml_content[:2000], language="python")
                else:
                    st.info("PyMOL script will be generated after analysis completes.")

            with col_files:
                st.markdown("**Per-state PDB files:**")
                for label, pdb_path in sorted(state_pdbs.items()):
                    if os.path.exists(str(pdb_path)):
                        with open(pdb_path) as f:
                            pdb_content = f.read()
                        n_atoms = pdb_content.count("\nHETATM") + pdb_content.count("\nATOM")
                        n_h = sum(1 for line in pdb_content.splitlines()
                                  if line.startswith(("ATOM", "HETATM")) and
                                  line[76:78].strip() == "H")
                        st.download_button(
                            f"📥 {os.path.basename(pdb_path)} ({n_atoms} atoms, {n_h} H)",
                            data=pdb_content,
                            file_name=os.path.basename(pdb_path),
                            mime="chemical/x-pdb",
                            key=f"dl_{label}",
                        )

        # ══════════════════════════════════════════════════════════════════
        # SECTION 2: Enriched conformer state table
        # ══════════════════════════════════════════════════════════════════
        st.markdown("---")
        st.subheader("Conformer State Table")

        import pandas as pd
        table_data = []
        for s in states:
            sd = s if isinstance(s, dict) else s.to_dict()
            row = {
                "Label": sd.get("label", "?"),
                "Charge": sd.get("charge", 0),
                "nH": sd.get("nH", 0),
                "pKa": sd.get("pka"),
                "Proton Exchange": sd.get("proton_exchange", ""),
                "H Added": ", ".join(sd.get("h_added", [])) or "—",
                "H Removed": ", ".join(sd.get("h_removed", [])) or "—",
                "Source": sd.get("source", ""),
                "LLM Model": sd.get("llm_model", ""),
                "Rationale": sd.get("rationale", ""),
                "References": "; ".join(sd.get("references", [])) if sd.get("references") else "—",
            }
            table_data.append(row)

        df = pd.DataFrame(table_data)
        edited_df = st.data_editor(
            df,
            num_rows="dynamic",
            use_container_width=True,
            column_config={
                "Label": st.column_config.TextColumn("Label", help="e.g., 01, +1, -1, +a, +b"),
                "Charge": st.column_config.NumberColumn("Charge", format="%+d"),
                "nH": st.column_config.NumberColumn("nH", help="Protons relative to neutral"),
                "pKa": st.column_config.NumberColumn("pKa", format="%.1f"),
                "Proton Exchange": st.column_config.TextColumn(
                    "Proton Exchange", help="Which protons differ from neutral",
                    width="medium"),
                "H Added": st.column_config.TextColumn("H Added", disabled=True),
                "H Removed": st.column_config.TextColumn("H Removed", disabled=True),
                "Source": st.column_config.TextColumn("Source"),
                "LLM Model": st.column_config.TextColumn("LLM Model", disabled=True),
                "Rationale": st.column_config.TextColumn("Rationale", width="medium"),
                "References": st.column_config.TextColumn(
                    "References", help="Literature sources",
                    width="large"),
            },
            key="state_editor"
        )

        # Merge edits back into state dicts
        if edited_df is not None:
            old_labels = [
                (s.get("label", "?") if isinstance(s, dict) else s.label)
                for s in states
            ]
            label_remap = {}

            for i, row in edited_df.iterrows():
                if i < len(states):
                    orig = states[i] if isinstance(states[i], dict) else states[i].to_dict()
                    new_label = str(row.get("Label", orig.get("label", "?")))
                    old_label = orig.get("label", "?")

                    if new_label != old_label:
                        label_remap[old_label] = new_label
                        logging.info(f"  Label change: {old_label} → {new_label}")

                    orig["label"] = new_label
                    try:
                        orig["charge"] = int(row.get("Charge", orig.get("charge", 0)) or 0)
                    except (ValueError, TypeError):
                        pass
                    try:
                        orig["nH"] = int(row.get("nH", orig.get("nH", 0)) or 0)
                    except (ValueError, TypeError):
                        pass
                    orig["pka"] = row.get("pKa")
                    orig["rationale"] = row.get("Rationale", "")
                    orig["proton_exchange"] = row.get("Proton Exchange", "")
                    states[i] = orig

            st.session_state["states"] = states

            # Propagate label renames to agent_state keys
            if label_remap and agent_state:
                if agent_state.get("state_pdbs"):
                    new_pdbs = {}
                    for old_lbl, pdb_path in agent_state["state_pdbs"].items():
                        new_lbl = label_remap.get(old_lbl, old_lbl)
                        new_pdbs[new_lbl] = pdb_path
                    agent_state["state_pdbs"] = new_pdbs

                if agent_state.get("h_diffs"):
                    new_diffs = {}
                    for old_lbl, diff in agent_state["h_diffs"].items():
                        new_lbl = label_remap.get(old_lbl, old_lbl)
                        new_diffs[new_lbl] = diff
                    agent_state["h_diffs"] = new_diffs

                agent_state["conformer_labels"] = [
                    (s.get("label") if isinstance(s, dict) else s.label)
                    for s in states
                ]

                st.session_state["agent_state"] = agent_state

        # Apply changes button
        if st.button("🔄 Apply Table Changes", help="Refresh page after editing labels or parameters"):
            st.rerun()

        # ══════════════════════════════════════════════════════════════════
        # SECTION 3: Add custom state via ionizable site selector
        # ══════════════════════════════════════════════════════════════════
        st.markdown("---")
        st.subheader("🧪 Add / Remove Protonation State")
        st.caption(
            "Select an ionizable site on the molecule to create a new "
            "protonation state. RDKit handles the chemistry automatically."
        )

        # Get ionizable sites from agent state
        ionizable = []
        if agent_state:
            ionizable = agent_state.get("_ionizable_sites", [])
            if not ionizable:
                try:
                    from mcce4_agent_ftpl.tools.rdkit_tools import (
                        get_ionizable_sites, get_mol_from_pdb as _gm
                    )
                    pdb = st.session_state.get("pdb_path")
                    if pdb:
                        mol = _gm(pdb, remove_hs=False)
                        if mol:
                            ionizable = get_ionizable_sites(mol)
                except Exception:
                    pass

        col_proto, col_depro = st.columns(2)

        with col_proto:
            st.markdown("**➕ Protonate a site** (add H)")
            proto_sites = [s for s in ionizable if "protonatable" in s.get("type", "")]
            if proto_sites:
                options = [
                    f"{s['name']} ({s['symbol']}) — {s['type'].replace('_', ' ')}"
                    for s in proto_sites
                ]
                selected = st.selectbox("Protonate at:", options, key="proto_site")
                new_label = st.text_input("Label:", value="+1", key="proto_label")
                if st.button("➕ Add protonated state", key="btn_proto"):
                    site_idx = options.index(selected)
                    site = proto_sites[site_idx]
                    new_state = {
                        "label": new_label, "charge": 1,
                        "nH": 1, "pka": None,
                        "smiles": "", "source": "user",
                        "site_atom": site["name"],
                        "rationale": f"User-added: protonate {site['name']}",
                        "proton_exchange": f"+H on {site['name']} ({site['type'].replace('_', ' ')})",
                        "pdb_path": None, "h_added": [], "h_removed": [],
                        "llm_model": "", "references": [],
                    }
                    states.append(new_state)
                    st.session_state["states"] = states
                    st.success(f"Added state '{new_label}': +H on {site['name']}")
                    st.rerun()
            else:
                st.info("No protonatable sites detected (run analysis first)")

        with col_depro:
            st.markdown("**➖ Deprotonate a site** (remove H)")
            depro_sites = [s for s in ionizable if "deprotonatable" in s.get("type", "")]
            if depro_sites:
                options = [
                    f"{s['name']} ({s['symbol']}) — {s['type'].replace('_', ' ')}"
                    for s in depro_sites
                ]
                selected = st.selectbox("Deprotonate at:", options, key="depro_site")
                new_label = st.text_input("Label:", value="-1", key="depro_label")
                if st.button("➖ Add deprotonated state", key="btn_depro"):
                    site_idx = options.index(selected)
                    site = depro_sites[site_idx]
                    new_state = {
                        "label": new_label, "charge": -1,
                        "nH": -1, "pka": None,
                        "smiles": "", "source": "user",
                        "site_atom": site["name"],
                        "rationale": f"User-added: deprotonate {site['name']}",
                        "proton_exchange": f"-H from {site['name']} ({site['type'].replace('_', ' ')})",
                        "pdb_path": None, "h_added": [], "h_removed": [],
                        "llm_model": "", "references": [],
                    }
                    states.append(new_state)
                    st.session_state["states"] = states
                    st.success(f"Added state '{new_label}': -H from {site['name']}")
                    st.rerun()
            else:
                st.info("No deprotonatable sites detected (run analysis first)")

        # Optional SMILES fallback for advanced users
        with st.expander("Advanced: add state from SMILES", expanded=False):
            manual_smiles = st.text_input(
                "SMILES for custom state:",
                placeholder="e.g., [NH2+]1CCCCC1",
                key="manual_smiles"
            )
            if manual_smiles:
                try:
                    from rdkit import Chem
                    mol = Chem.MolFromSmiles(manual_smiles)
                    if mol:
                        charge = Chem.GetFormalCharge(mol)
                        svg_custom = mol_to_svg(mol, size=(300, 200))
                        if svg_custom:
                            st.image(svg_custom, use_container_width=False)
                        st.info(f"Formal charge: {charge:+d}")
                        smiles_label = st.text_input("Label:",
                                                      value=f"{charge:+d}" if charge else "01",
                                                      key="smiles_label")
                        if st.button("➕ Add", key="add_smiles"):
                            new_state = {
                                "label": smiles_label, "charge": charge,
                                "nH": charge, "pka": None,
                                "smiles": manual_smiles,
                                "source": "user-smiles",
                                "rationale": "Custom SMILES entry",
                                "proton_exchange": "",
                                "pdb_path": None, "site_atom": None,
                                "h_added": [], "h_removed": [],
                                "llm_model": "", "references": [],
                            }
                            states.append(new_state)
                            st.session_state["states"] = states
                            st.success(f"Added '{smiles_label}' (charge={charge:+d})")
                            st.rerun()
                    else:
                        st.error("Invalid SMILES")
                except ImportError:
                    st.error("RDKit required — pip install rdkit")

        # ══════════════════════════════════════════════════════════════════
        # Approve / Run
        # ══════════════════════════════════════════════════════════════════
        st.markdown("---")
        col_cancel, col_spacer, col_approve = st.columns([1, 2, 1])

        with col_cancel:
            if st.button("✖ Cancel", type="secondary", use_container_width=True):
                st.session_state["phase"] = "upload"
                st.session_state["states"] = []
                st.session_state["phase_status"] = {}
                st.rerun()

        with col_approve:
            if st.button("✔ Approve & Generate .ftpl", type="primary", use_container_width=True):
                st.session_state["phase"] = "running"
                st.session_state["approved"] = True

                from mcce4_agent_ftpl.agent import (
                    node_generate_template, node_assign_charges,
                    node_rxn_calibration, node_done,
                )

                agent_state = st.session_state["agent_state"]
                agent_state["states"] = st.session_state["states"]
                agent_state["conformer_labels"] = [
                    (s["label"] if isinstance(s, dict) else s.label)
                    for s in st.session_state["states"]
                ]
                agent_state["user_approved"] = True
                agent_state["needs_user_review"] = False

                # Set up log capture for generation phases
                gen_log_capture, gen_log_handler = _start_log_capture()

                status = st.status("🔧 Generating topology file...", expanded=True)
                generation_ok = True

                with status:
                    # Phase 5: Template generation
                    st.session_state["phase_status"]["phase5"] = "running"
                    st.write("**Phase 5:** Generating per-state ftpl templates (pdb2ftpl.py)...")
                    try:
                        agent_state = node_generate_template(agent_state)
                        st.session_state["phase_status"]["phase5"] = "done"
                        if agent_state.get("errors"):
                            st.session_state["phase_status"]["phase5"] = "error"
                            for e in agent_state["errors"]:
                                st.error(f"❌ {e}")
                            generation_ok = False
                        else:
                            ftpl = agent_state.get("ftpl_path", "?")
                            st.write(f"✅ Template: `{ftpl}`")
                    except Exception as e:
                        st.session_state["phase_status"]["phase5"] = "error"
                        st.error(f"❌ Template generation failed: {e}")
                        generation_ok = False

                    # Phase 6: Charge assignment
                    if generation_ok:
                        st.session_state["phase_status"]["phase6"] = "running"
                        st.write("**Phase 6:** Computing per-state charges (OpenEye QuacPac TK)...")
                        try:
                            agent_state = node_assign_charges(agent_state)
                            st.session_state["phase_status"]["phase6"] = "done"
                            if agent_state.get("complete"):
                                st.session_state["phase_status"]["phase6"] = "error"
                                for e in agent_state.get("errors", []):
                                    st.error(f"❌ {e}")
                                generation_ok = False
                            else:
                                n_states = len(agent_state.get("per_state_charges", {}))
                                st.write(f"✅ Charges computed for {n_states} state(s)")
                        except Exception as e:
                            st.session_state["phase_status"]["phase6"] = "error"
                            st.error(f"❌ Charge assignment failed: {e}")
                            generation_ok = False

                    # Phase 7: RXN calibration
                    if generation_ok and not agent_state.get("dry_run", False):
                        st.session_state["phase_status"]["phase7"] = "running"
                        st.write("**Phase 7:** RXN calibration — MCCE step1 + step2, "
                                 "then step3.py for each dielectric (2, 4, 8)...")
                        st.caption(
                            "For each dielectric, the lowest dsolv across all conformers "
                            "of each state is used to calibrate rxn."
                        )
                        try:
                            agent_state = node_rxn_calibration(agent_state)
                            rxn = agent_state.get("rxn_values", {})
                            if isinstance(rxn, dict) and "error" not in str(rxn):
                                st.session_state["phase_status"]["phase7"] = "done"
                                st.write(f"✅ RXN calibration complete")
                                # Show rxn summary
                                for rxn_key, vals in rxn.items():
                                    if isinstance(vals, dict):
                                        parts = [f"{ct}: {v:.3f}" for ct, v in vals.items()]
                                        st.caption(f"  {rxn_key}: {', '.join(parts)}")
                            else:
                                st.session_state["phase_status"]["phase7"] = "error"
                                st.warning(f"⚠ RXN calibration: {rxn}")
                        except Exception as e:
                            st.session_state["phase_status"]["phase7"] = "error"
                            st.warning(f"⚠ RXN calibration failed: {e}")
                    elif agent_state.get("dry_run", False):
                        st.write("⏩ Dry run — RXN calibration skipped")

                    # Done
                    agent_state = node_done(agent_state)

                    if agent_state.get("warnings"):
                        for w in agent_state["warnings"]:
                            st.warning(f"⚠ {w}")

                    # Show generation log
                    gen_log_text = _stop_log_capture(gen_log_capture, gen_log_handler)
                    if gen_log_text.strip():
                        _render_log_panel(gen_log_text, "Generation Log", 500, "generation")
                        # Append to session log
                        prev_log = st.session_state.get("log_output", "")
                        st.session_state["log_output"] = (
                            prev_log + "\n\n" +
                            "=" * 60 + "\n" +
                            "  GENERATION PHASES (5-7)\n" +
                            "=" * 60 + "\n\n" +
                            gen_log_text
                        )

                if generation_ok:
                    status.update(label="✅ Topology file generated!", state="complete",
                                  expanded=True)
                else:
                    status.update(label="❌ Generation failed — check errors above",
                                  state="error", expanded=True)

                st.session_state["agent_state"] = agent_state
                st.session_state["phase"] = "complete"
                st.rerun()


# ──────────────────────────────────────────────────────────────────────────────
# TAB 3: Output
# ──────────────────────────────────────────────────────────────────────────────
with tab_output:
    st.subheader("📄 Generated Topology File")

    agent_state = st.session_state.get("agent_state")
    if agent_state and agent_state.get("ftpl_path"):
        ftpl_path = agent_state["ftpl_path"]

        if os.path.exists(ftpl_path):
            with open(ftpl_path) as f:
                ftpl_content = f.read()

            # Download button
            st.download_button(
                "⬇️ Download .ftpl",
                data=ftpl_content,
                file_name=os.path.basename(ftpl_path),
                mime="text/plain",
                type="primary"
            )

            # Show content
            st.code(ftpl_content, language="text", line_numbers=True)

            # Show summary
            if agent_state.get("rxn_values"):
                st.subheader("RXN Calibration Results")
                st.json(agent_state["rxn_values"])

            if agent_state.get("warnings"):
                st.subheader("⚠️ Warnings")
                for w in agent_state["warnings"]:
                    st.warning(w)
        else:
            st.warning(f"File not found: {ftpl_path}")
    else:
        st.info("Complete the analysis and approve states to generate the topology file.")


# ──────────────────────────────────────────────────────────────────────────────
# TAB 4: Full Session Log (persistent, always available)
# ──────────────────────────────────────────────────────────────────────────────
with tab_log:
    st.subheader("📋 Full Session Log")
    st.caption(
        "Complete log of all agent phases. Scroll to review output at any time. "
        "Protonation state changes are highlighted: "
        "🟢 green = H added, 🔴 red = H removed, "
        "🔵 cyan = phase headers, 🟣 magenta = dsolv/rxn values."
    )

    full_log = st.session_state.get("log_output", "")
    if full_log.strip():
        # Download log button
        col_dl, col_clear = st.columns([1, 3])
        with col_dl:
            lig_id = st.session_state.get("lig_id", "LIG")
            st.download_button(
                "⬇️ Download Log",
                data=full_log,
                file_name=f"mcce4_agent_{lig_id}.log",
                mime="text/plain",
            )
        with col_clear:
            if st.button("🗑 Clear Log"):
                st.session_state["log_output"] = ""
                st.rerun()

        # Render the full log with max height and scrolling
        _render_log_panel(full_log, "Session Log", 700, "full_session")

        # Show state naming summary if available
        states = st.session_state.get("states", [])
        if states:
            st.markdown("---")
            st.subheader("🏷️ State Naming Summary")
            for s in states:
                sd = s if isinstance(s, dict) else s.to_dict()
                label = sd.get("label", "?")
                charge = int(sd.get("charge", 0) or 0)
                rationale = sd.get("rationale", "")
                proton_ex = sd.get("proton_exchange", "")

                if charge > 0:
                    badge_color = "#4ec94e"
                    badge_icon = "⊕"
                elif charge < 0:
                    badge_color = "#f44747"
                    badge_icon = "⊖"
                else:
                    badge_color = "#4fc1ff"
                    badge_icon = "◯"

                st.markdown(
                    f'<div style="padding:4px 0;">'
                    f'<span style="background:{badge_color}33; color:{badge_color}; '
                    f'padding:2px 8px; border-radius:3px; font-weight:bold;">'
                    f'{badge_icon} {lig_id}{label}</span> '
                    f'(charge={charge:+d}) '
                    f'{"— " + proton_ex if proton_ex else ""} '
                    f'<span style="color:#888;">| {rationale[:80]}</span>'
                    f'</div>',
                    unsafe_allow_html=True,
                )
    else:
        st.info(
            "No log output yet. Run the analysis from the **Input** tab to "
            "start generating logs. All phase output will appear here."
        )


# ──────────────────────────────────────────────────────────────────────────────
# Footer
# ──────────────────────────────────────────────────────────────────────────────
st.markdown("---")
st.caption("MCCE4 Topology Agent v5.0 | GunnerLab | "
           "[Tutorial](https://gunnerlab.github.io/mcce4_tutorial/docs/topology/)")
