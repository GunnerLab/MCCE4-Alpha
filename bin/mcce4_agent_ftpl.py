#!/usr/bin/env python3
"""
mcce4_agent_ftpl.py — MCCE4 Topology File AI Agent
====================================================

Launcher script. Place this in MCCE4-Alpha/bin/ alongside the
mcce4_agent_ftpl/ package directory.

Usage:
  mcce4_agent_ftpl.py EMH.pdb                                     # Full auto (PDB input)
  mcce4_agent_ftpl.py EMH.cif                                     # CIF input (auto-converts)
  mcce4_agent_ftpl.py --lig-code EMH                               # Ligand code only (fetches from RCSB)
  mcce4_agent_ftpl.py EMH.pdb --lig-code EMH                      # PDB + ligand code
  mcce4_agent_ftpl.py --lig-code EMH --dry-run --no-llm            # Quick test from RCSB
  mcce4_agent_ftpl.py EMH.pdb --gui                               # Web GUI
  mcce4_agent_ftpl.py EMH.pdb --state-pdbs EMH_01.pdb EMH_+1.pdb  # User states
  mcce4_agent_ftpl.py EMH.pdb --no-llm --dry-run                  # Minimal run
  mcce4_agent_ftpl.py EMH.pdb --llm-provider claude               # Use Claude LLM
  mcce4_agent_ftpl.py EMH.pdb --llm-provider chatgpt              # Use ChatGPT LLM
  mcce4_agent_ftpl.py EMH.pdb --llm-provider claude --api-key KEY  # Custom API key

Install dependencies:
  pip install google-genai langgraph dimorphite-dl rdkit streamlit
  pip install anthropic openai  # optional, for Claude/ChatGPT providers
  export GEMINI_API_KEY="your_free_key"   # from https://ai.google.dev
  export ANTHROPIC_API_KEY="your_key"     # for Claude provider
  export OPENAI_API_KEY="your_key"        # for ChatGPT provider
"""

import sys
import os

# Ensure the package is findable from bin/
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from mcce4_agent_ftpl.cli import main

if __name__ == "__main__":
    main()
