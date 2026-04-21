"""
mcce4_agent_ftpl — MCCE4 Topology File AI Agent Package
=========================================================

An AI agent for creating MCCE4 topology files (.ftpl) from ligand PDB files.
Uses Dimorphite-DL for protonation state enumeration, RDKit for chemistry
validation, Google Gemini (free tier) for reasoning, and LangGraph for
agentic orchestration.

Install:
    pip install google-genai langgraph dimorphite-dl rdkit streamlit

Usage:
    # CLI
    mcce4_agent_ftpl.py EMH.pdb --gui

    # Python
    from mcce4_agent_ftpl.agent import run_agent
    run_agent("EMH.pdb")
"""

__version__ = "5.0.0"
__author__ = "Gehan / MCCE4 Team (GunnerLab)"
