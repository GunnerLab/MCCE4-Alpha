# First Time Using MCCE?

Consider viewing the [MCCE4 Tutorial](README_MCCE4_Tutorial.md).

# MCCE Command Line Tools
<!--
 Keep listing in alphabetical order;
 If a README file is provided, format the 1st part of the link (the one in square brackets)
 so that the tool name is enclosed in backticks, followed by a very short description, e.g.:
 [`tcms`: Tautomeric, topN Charge Microstates](README_tcms.md);
-->

## [Benchmarking](README_benchmark.md)
This tool comes with multiple commands to: prepare the run folders for a dataset of pdbs, launch the calculations, 
analyze the results, and compare two datasets.
Main Commands: `bench_setup`, `bench_launch`, `bench_analyze`, `bench_compare`.

### [Graphical workflows](README_benchmark_workflows.md)
### [Tutorial with use cases](README_benchmark_tutorial.md)
### [Cheat Sheet] WIP

### Help in a notebook:
This notebook: [analysis.ipynb](analyze.ipynb) shows how a user can access the code base. It demonstrates a specific  
use-case: how to selectively redo the figures of an analysis or comparison. Note that the data paths given will have to be
amended to your own and that the setting of the 'MC4' path will depend on the location were the notebook is ultimately used.

## `getpdb`
Download one or more (bioassembly) pdb files from the RSCB Protein Data Bank using `mcce4.protinfo` download function.
__Usage__: getpdb 1ots 4LZT

## [`glossary`: Glossary of MCCE parameters](README_glossary.md)

## [`ProtInfo`](README_ProtInfo.md)
This is the app that creates the 'prerun report' in the benchmarking tool.
Does the same at the command line for a single protein (with option to download it).

## [`tcms`: Tautomeric, topN Charge Microstates](README_tcms.md)
This tool outputs:
  - A file listing the topN tautomeric charge states, along with
    their related properties: energy (E), net charge (sum_crg), count, and occupancy (occ).
  - The <topN> pdb files of each charge state.

---
[Home](README.md) | [Feature Requests](Features.md) | [Developer Manual](DevManual.md) | [User Manual](UserManual.md)
