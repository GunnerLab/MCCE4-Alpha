
#!/usr/bin/env python

"""
Module: pbs_interfaces.py

PB Solver Interface Module

"""

import logging
import math
import time
import numpy as np
import os
from pathlib import Path
from pyexpat import features
import struct
import subprocess
import sys


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class PBS_TEMPLATE:
    """
    Template interface. This serves as a template for writing your own interface.
    
    To test out this template, run this command in a working directory with step2_out.pdb file
        step3.py -c 101 101 --fly --debug -s template
    This command will
    * calculate conformer 101, you can choose other number.
    * calculate reference reaction field energy on the fly
    * keep pbe solver working directory for debugging purpose
    * call a dummy pbe solver named template

    When writing interface, add a interface connection to step3.py around line 392.

    The PBE solver working directory is saved and the location is reported at the end of step3.py.

    The energy look up table is in energies/*.raw.
    """

    def __init__(self):
        """Initialize default values and allow to be overwritten from run_options."""
        self.exe = "pbs_template.py"   # PBE solver executable, should be made available by the execution environment
        self.epsilon_prot = 4.0         # default dielectric constant for protein
        self.epsilon_solv = 80.0        # default dielectric constant for solvent
        return

    def write_pqr(self, bound):
        """
        This is a help function. It writes 3 pqr files in PBE solver working directory
            * floating side chain: float.pqr
            * single conformation boundary: single.pqr
            * multi conformation boundary: multi.pqr
        """
        lines = self.xyzrcp2prq(bound.float_bnd_xyzrcp)
        open("float.pqr", "w").writelines(lines)

        lines = self.xyzrcp2prq(bound.single_bnd_xyzrcp)
        open("single.pqr", "w").writelines(lines)

        lines = self.xyzrcp2prq(bound.multi_bnd_xyzrcp)
        open("multi.pqr", "w").writelines(lines)

    def xyzrcp2prq(self, xyzrcp_records):
        """
        Convert xyzrcp records to pqr lines
        """
        lines = []
        serial = 1
        for record in xyzrcp_records:
            line = "ATOM  %5d ATOM RES %5d    %8.3f%8.3f%8.3f%8.3f%12.3f\n" % (serial, serial, record.x, record.y,
                                                                             record.z, record.r, record.c)
            serial += 1
            lines.append(line)
        return lines

    def write_run_options(self, run_options):
        """
        Write command options passed to this object into file run_options.txt.
        """
        with open("run_options.txt", "w") as fh:
            for key, value in vars(run_options).items():
                fh.write("%s: %s %s\n" % (key, str(value), str(type(value))))

    def run(self, bound, run_options):
        """
        This is the PBE solver's main program. The solver will
        Input
            * bound - 3 boundary in xyzrcp format
            * run_options - dictionary of run options
        Output
            * update xyzrcp records p value
            * return (rxn0, rxn)
        """
        # Convert xyzrcp records in bound to pqr files
        self.write_pqr(bound)
        # Show what run options are passed in
        self.write_run_options(run_options)

        # Calculate reference reaction field energy, using float boundary condition
        # There is no need to update site potential
        rxn0 = 0.0
        if run_options.fly:
            # Set up PBE solver at float bnd condition
            # Use bound.float_bnd_xyzrcp
            # ...

            # Run PBE solver and save the log to result
            result = subprocess.run([self.exe], capture_output=True, text=True)

            # Obtain rxn0 from the log file as an example of extracting information from stdout
            lines = result.stdout.split("\n")
            for line in lines:
                if "RXN =" in line:  # find unique pattern to extract the value
                    fields = line.split("=")
                    rxn0 = float(fields[1])
                    break

        # Calculate site potential at single bnd condition
        # Need to calculate rxn and update site potential
        rxn = 0.0
        # Set up PBE solver at single bnd condition
        # Use bound.single_bnd_xyzrcp

        # Run PBE solver and save the log to result
        result = subprocess.run([self.exe], capture_output=True, text=True)
        # update p in bound.single_bnd_xyzrcp
        # Obtain rxn
        rxn = -4.500

        # Calculate site potential at multi bnd condition
        # Need to update site potential, No need to calculate rxn
        # Set up PBE solver at multu bnd condition
        # Use bound.multi_bnd_xyzrcp

        # Run PBE solver and save the log to result
        result = subprocess.run([self.exe], capture_output=True, text=True)
        # update p in bound.multi_bnd_xyzrcp

        return (rxn0, rxn)


class PBS_DELPHI:
    """
    Legacy delphi interface
    """
    def __init__(self):
        # consider loading these parameters from a file
        self.exe = "delphi"   # This has to be made available by the execution environment
        self.radius_probe = 1.4
        self.grids_per_ang = 2.0
        self.epsilon_prot = 4.0  # default value, be overwritten by run_options
        self.epsilon_solv = 80.0
        self.ionrad = 2.0
        self.salt = 0.15
        self.grids_delphi = 65
        self.KCAL2KT = 1.688
        time.sleep(2)
        return
    
    def depth(self, bound):
        # determine delphi focusing depth
        x_min = x_max = bound.single_bnd_xyzrcp[0].x
        y_min = y_max = bound.single_bnd_xyzrcp[0].y
        z_min = z_max = bound.single_bnd_xyzrcp[0].z
        for p in bound.single_bnd_xyzrcp[1:]:
            if x_min > p.x: x_min = p.x
            if x_max < p.x: x_max = p.x
            if y_min > p.y: y_min = p.y
            if y_max < p.y: y_max = p.y
            if z_min > p.z: z_min = p.z
            if z_max < p.z: z_max = p.z

        dx = x_max - x_min
        dy = y_max - y_min
        dz = z_max - z_min
        dm = max(dx, dy, dz)
        # expand the largest dimension by the probe radius and safety:
        dm += self.radius_probe * 2 + 3.4

        # scale is a multiplier on grid_per_ang required to reach the target resolution
        scale = self.grids_per_ang/(self.grids_delphi/(2*dm))
        if scale <= 1.0:
            depth = 1
        else:
            depth = math.ceil(math.log(scale)/math.log(2.0)) + 1

        #print(dm, depth, scale)
        return depth

    def write_fort15(self, xyzrcp):
        i = 1
        with open("fort.15", "w") as fh:
            for p in xyzrcp:
                header = "ATOM      0  O   LYS %5d    " % i
                fh.write("%-30s%8.3f%8.3f%8.3f\n" % (header, p.x, p.y, p.z))
                i += 1
        return

    def write_fort13(self, xyzrcp):
        struct_fmt = '=ifffffi'
        with open("fort.13", "wb") as fh:
            for p in xyzrcp:
                record_unf = struct.pack(struct_fmt, 20, p.x, p.y, p.z, p.r, p.c, 20)
                fh.write(record_unf)
        return

    def collect_phi(self, depth, xyzrcp):
        # collect results from the log
        try:
            lines = open("run01.frc", "r").readlines()
        except OSError:
            logger.critical("Could not open Delphi output file run01.frc.")
            sys.exit(1)

        for counter in range(len(xyzrcp)):
            #line_segment = lines[12 + counter][19:]
            #print(f"line_segment: {line_segment}")
            phi = float(lines[12+counter][19:30])
            #print(f"phi: {phi}")
            xyzrcp[counter].p = phi

        # If the potential is non 0 in focusing runs, update
        for i in range(1, depth):
            frc_name = "run%02d.frc" % (i+1)
            try:
                lines = open(frc_name, "r").readlines()
            except OSError:
                logger.critical("Could not open Delphi output file %s." % frc_name)
                sys.exit(1)

            for counter in range(len(xyzrcp)):
                phi = float(lines[12 + counter][19:30])
                if abs(phi) > 0.0001:
                    xyzrcp[counter].p = phi

        return


    def collect_rxn(self, log):
        log_lines = log.split("\n")
        found = False
        for line in log_lines:
            if "corrected reaction field energy:" in line:
                rxn = float(line[34:].split()[0]) / self.KCAL2KT
                found = True

        if not found:
            logger.critical("Did not detect corrected reaction field energy line. Delphi failed!")
            sys.exit(1)

        return rxn

    def run(self, bound, run_options):
        """PBE solver interface for delphi.
        Generate site p in both boundary conditions and
        return reference rxn0, and rxn in single boundary condition.
        Input:
            bound - dielectric boundary object
            run_options - command options in dictionary
        Returns:
         -  A 2-tuple: rxn0, rxn: the most negative value of the focusing runs 
        """
        # snippets to check the input and environment
        # Current working directory
        # cwd = os.getcwd()
        # print(cwd)
        # What are in bound
        #print(vars(bound))

        # Caclulate rxn0 using float boundary condition
        depth = self.depth(bound)
        logger.info("Delphi focusing depth: %d" % depth)
        rxns = []
        # floating boundary condition
        self.write_fort13(bound.float_bnd_xyzrcp)
        self.write_fort15(bound.float_bnd_xyzrcp)

        # fort.27
        center = [0.0, 0.0, 0.0]
        weight = 0.0
        for p in bound.float_bnd_xyzrcp:
            w = abs(p.c)
            if w > 0.00001:
                center[0] += p.x * w
                center[1] += p.y * w
                center[2] += p.z * w
                weight += w

        if weight > 0.000001:
            center = [c/(weight+0.000001) for c in center]
        else:
            logger.error("PB solver shouldn't run a conformer with no charged atom.")

        with open("fort.27", "w") as fh:
            fh.write("ATOM  %5d  C   CEN  %04d    %8.3f%8.3f%8.3f\n" % (1, 1, center[0], center[1], center[2]))

        # fort.10
        self.epsilon_prot = run_options.d
        self.epsilon_solv = run_options.do
        self.salt=run_options.salt
        with open("fort.10", "w") as fh:
            fh.write("gsize=%d\n" % self.grids_delphi)
            fh.write("scale=%.2f\n" % (self.grids_per_ang/2**(depth-1)))
            fh.write("in(unpdb,file=\"fort.13\")\n")
            fh.write("indi=%.1f\n" % self.epsilon_prot)
            fh.write("exdi=%.1f\n" % self.epsilon_solv)
            fh.write("ionrad=%.1f\n" % self.ionrad)
            fh.write("salt=%.2f\n" % self.salt)
            fh.write("bndcon=2\n")
            fh.write("center(777, 777, 0)\n")
            fh.write("out(frc,file=\"run01.frc\")\n")
            fh.write("out(phi,file=\"run01.phi\")\n")
            fh.write("site(a,c,p)\n")
            fh.write("energy(g,an,sol)\n")   # g for grid energy, sol for corrected rxn

        # 1st delphi run
        result = subprocess.run([self.exe], capture_output=True, text=True)
        if result.returncode != 0:
            logger.critical(f"Delphi failed with error:\n{result.stderr}")
            sys.exit(1)

        rxns.append(self.collect_rxn(result.stdout))

        # subsequent delphi runs
        for i in range(1, depth):
            with open("fort.10", "w") as fh:
                fh.write("gsize=%d\n" % self.grids_delphi)
                fh.write("scale=%.2f\n" % (self.grids_per_ang/2**(depth-1-i)))
                fh.write("in(unpdb,file=\"fort.13\")\n")
                fh.write("in(phi,file=\"run%02d.phi\")\n" % i)
                fh.write("indi=%.1f\n" % self.epsilon_prot)
                fh.write("exdi=%.1f\n" % self.epsilon_solv)
                fh.write("ionrad=%.1f\n" % self.ionrad)
                fh.write("salt=%.2f\n" % self.salt)
                fh.write("bndcon=3\n")
                fh.write("center(777, 777, 0)\n")
                fh.write("out(frc,file=\"run%02d.frc\")\n" % (i+1))
                fh.write("out(phi,file=\"run%02d.phi\")\n" % (i+1))
                fh.write("site(a,c,p)\n")
                fh.write("energy(g,an,sol)\n")  # g for grid energy, sol for corrected rxn

            result = subprocess.run([self.exe], capture_output=True, text=True)
            if result.returncode != 0:
                logger.critical(f"Delphi failed with error:\n{result.stderr}")
                sys.exit(1)

            rxns.append(self.collect_rxn(result.stdout))
        rxn0 = min(rxns)

        # single side chain boundary condition
        depth = self.depth(bound)
        logger.info("Delphi focusing depth: %d" % depth)
        rxns = []
        # single side chain boundary condition
        # fort.13
        # The first run starts with fort.13 as dielectric boundary, the following
        # runs will be focusing runs, using the phi as input
        #
        self.write_fort13(bound.single_bnd_xyzrcp)
        self.write_fort15(bound.single_bnd_xyzrcp)

        # fort.27
        center = [0.0, 0.0, 0.0]
        weight = 0.0
        for p in bound.single_bnd_xyzrcp:
            w = abs(p.c)
            if w > 0.00001:
                center[0] += p.x * w
                center[1] += p.y * w
                center[2] += p.z * w
                weight += w

        if weight > 0.000001:
            center = [c/(weight+0.000001) for c in center]
        else:
            logger.error("PB solver shouldn't run a conformer with no charged atom.")

        with open("fort.27", "w") as fh:
            fh.write("ATOM  %5d  C   CEN  %04d    %8.3f%8.3f%8.3f\n" % (1, 1, center[0], center[1], center[2]))

        # fort.10
        self.epsilon_prot = run_options.d
        self.epsilon_solv = run_options.do
        with open("fort.10", "w") as fh:
            fh.write("gsize=%d\n" % self.grids_delphi)
            fh.write("scale=%.2f\n" % (self.grids_per_ang/2**(depth-1)))
            fh.write("in(unpdb,file=\"fort.13\")\n")
            fh.write("indi=%.1f\n" % self.epsilon_prot)
            fh.write("exdi=%.1f\n" % self.epsilon_solv)
            fh.write("ionrad=%.1f\n" % self.ionrad)
            fh.write("salt=%.2f\n" % self.salt)
            fh.write("bndcon=2\n")
            fh.write("center(777, 777, 0)\n")
            fh.write("out(frc,file=\"run01.frc\")\n")
            fh.write("out(phi,file=\"run01.phi\")\n")
            fh.write("site(a,c,p)\n")
            fh.write("energy(g,an,sol)\n")   # g for grid energy, sol for corrected rxn

        # 1st delphi run
        result = subprocess.run([self.exe], capture_output=True, text=True)
        if result.returncode != 0:
            logger.critical(f"Delphi failed with error:\n{result.stderr}")
            sys.exit(1)

        rxns.append(self.collect_rxn(result.stdout))

        # subsequent delphi runs
        for i in range(1, depth):
            with open("fort.10", "w") as fh:
                fh.write("gsize=%d\n" % self.grids_delphi)
                fh.write("scale=%.2f\n" % (self.grids_per_ang/2**(depth-1-i)))
                fh.write("in(unpdb,file=\"fort.13\")\n")
                fh.write("in(phi,file=\"run%02d.phi\")\n" % i)
                fh.write("indi=%.1f\n" % self.epsilon_prot)
                fh.write("exdi=%.1f\n" % self.epsilon_solv)
                fh.write("ionrad=%.1f\n" % self.ionrad)
                fh.write("salt=%.2f\n" % self.salt)
                fh.write("bndcon=3\n")
                fh.write("center(777, 777, 0)\n")
                fh.write("out(frc,file=\"run%02d.frc\")\n" % (i+1))
                fh.write("out(phi,file=\"run%02d.phi\")\n" % (i+1))
                fh.write("site(a,c,p)\n")
                fh.write("energy(g,an,sol)\n")  # g for grid energy, sol for corrected rxn

            result = subprocess.run([self.exe], capture_output=True, text=True)
            if result.returncode != 0:
                logger.critical(f"Delphi failed with error:\n{result.stderr}")
                sys.exit(1)

            rxns.append(self.collect_rxn(result.stdout))

        # collect results from frc files
        self.collect_phi(depth, bound.single_bnd_xyzrcp)

        # multi side chain boundary condition
        # fort.13 for first run
        self.write_fort13(bound.multi_bnd_xyzrcp)
        self.write_fort15(bound.multi_bnd_xyzrcp)

        # fort.27
        center = [0.0, 0.0, 0.0]
        weight = 0.0
        for p in bound.multi_bnd_xyzrcp:
            w = abs(p.c)
            if w > 0.00001:
                center[0] += p.x * w
                center[1] += p.y * w
                center[2] += p.z * w
                weight += w

        if weight > 0.000001:
            center = [c/(weight+0.000001) for c in center]
        else:
            logger.error("PB solver shouldn't run a conformer with no charged atom.")

        with open("fort.27", "w") as fh:
            fh.write("ATOM  %5d  C   CEN  %04d    %8.3f%8.3f%8.3f\n" % (1, 1, center[0], center[1], center[2]))

        # fort.10
        self.epsilon_prot = run_options.d
        self.epsilon_solv = run_options.do
        with open("fort.10", "w") as fh:
            fh.write("gsize=%d\n" % self.grids_delphi)
            fh.write("scale=%.2f\n" % (self.grids_per_ang/2**(depth-1)))
            fh.write("in(unpdb,file=\"fort.13\")\n")
            fh.write("indi=%.1f\n" % self.epsilon_prot)
            fh.write("exdi=%.1f\n" % self.epsilon_solv)
            fh.write("ionrad=%.1f\n" % self.ionrad)
            fh.write("salt=%.2f\n" % self.salt)
            fh.write("bndcon=2\n")
            fh.write("center(777, 777, 0)\n")
            fh.write("out(frc,file=\"run01.frc\")\n")
            fh.write("out(phi,file=\"run01.phi\")\n")
            fh.write("site(a,c,p)\n")
            fh.write("energy(g,an,sol)\n")   # g for grid energy, sol for corrected rxn

        # 1st delphi run
        result = subprocess.run([self.exe], capture_output=True, text=True)
        if result.returncode != 0:
            logger.critical(f"Delphi failed with error:\n{result.stderr}")
            sys.exit(1)

        # subsequent delphi runs
        for i in range(1, depth):
            with open("fort.10", "w") as fh:
                fh.write("gsize=%d\n" % self.grids_delphi)
                fh.write("scale=%.2f\n" % (self.grids_per_ang/2**(depth-1-i)))
                fh.write("in(unpdb,file=\"fort.13\")\n")
                fh.write("in(phi,file=\"run%02d.phi\")\n" % i)
                fh.write("indi=%.1f\n" % self.epsilon_prot)
                fh.write("exdi=%.1f\n" % self.epsilon_solv)
                fh.write("ionrad=%.1f\n" % self.ionrad)
                fh.write("salt=%.2f\n" % self.salt)
                fh.write("bndcon=3\n")
                fh.write("center(777, 777, 0)\n")
                fh.write("out(frc,file=\"run%02d.frc\")\n" % (i+1))
                fh.write("out(phi,file=\"run%02d.phi\")\n" % (i+1))
                fh.write("site(a,c,p)\n")
                fh.write("energy(g,an,sol)\n")  # g for grid energy, sol for corrected rxn
            result = subprocess.run([self.exe], capture_output=True, text=True)
            if result.returncode != 0:
                logger.critical(f"Delphi failed with error:\n{result.stderr}")
                sys.exit(1)

        # collect results from frc files
        self.collect_phi(depth, bound.multi_bnd_xyzrcp)
        #print(f"depth: {depth}, rxn0: {rxn0}, rxns: {rxns}")
        return (rxn0, min(rxns))


class PBS_NGPB:
    """
    NGPB interface
    """
    def __init__(self,instance_name):
        # consider loading these parameters from a file
        #self.exe = ['mpirun', '-np', '1', 'poisson_boltzmann', '--potfile', 'options.pot', '--pqrfile']
        self.exe = ['apptainer', 'exec', f'instance://{instance_name}', 'ngpb', '--potfile', 'options.pot', '--pqrfile']
        
        self.mesh_shape = 3
        self.scale = 1.5
        self.outlevel =1
        self.cx_foc = 0
        self.cy_foc = 0
        self.cz_foc = 0
        self.n_grid = 50
        self.perfil1 = 0.98
        self.perfil2 = 0.95
        self.bc_type = 2
        self.epsilon_prot = 4.0  
        self.epsilon_solv = 80.0
        self.salt = 0.150
        self.KCAL2KT = 1.688
        self.number_core = 1
        self.T = 298.15
        self.calc_energy = 1
        self.atoms_write = 1
        self.stern_layer_surf = 1
        self.surface_type = 0 
        self.stern_layer_thickness = 2.0
        self.radius_probe = 1.4
        # get environment variables
        self.my_env = os.environ.copy()

        return
    
    def write_option_file(self, center, number_grid, energy, atoms_w):
        with open("options.pot", "w") as fh:
            fh.write("#Parameters file\n")

            fh.write("[mesh]\n")
            fh.write("mesh_shape=%d\n" % self.mesh_shape)
            fh.write("scale=%.2f\n" % self.scale)
            fh.write("perfil1=%.2f\n" % self.perfil1)
            fh.write("perfil2=%.2f\n" % self.perfil2)
            fh.write("cx_foc=%.2f\n" % center[0])
            fh.write("cy_foc=%.2f\n" % center[1])
            fh.write("cz_foc=%.2f\n" % center[2])
            fh.write("n_grid=%.2f\n" % number_grid)
            fh.write("outlevel=%d\n" % self.outlevel)
            fh.write("[../]\n")

            fh.write("[model]\n")
            fh.write("linearized=1\n")
            fh.write("bc_type=2\n")
            fh.write("molecular_dielectric_constant=%.1f\n" % self.epsilon_prot)
            fh.write("solvent_dielectric_constant=%.1f\n" % self.epsilon_solv)
            fh.write("ionic_strength=%.2f\n" % self.salt)
            fh.write("T=%.2f\n" % self.T)
            fh.write("calc_energy=%d\n" % energy)
            fh.write("atoms_write=%d\n" %atoms_w)
            fh.write("[../]\n")

            fh.write("[surface]\n")
            fh.write("surface_type=%d\n" % self.surface_type)
            fh.write("surface_parameter=%.3f\n" % self.radius_probe)
            fh.write("stern_layer_surf=%d\n" % self.stern_layer_surf)
            fh.write("stern_layer_thickness=%.3f\n" % self.stern_layer_thickness)
            fh.write("[../]\n")

            fh.write("[algorithm]\n")
            fh.write("linear_solver = lis\n")
            fh.write("solver_options = -p\ ssor\ -ssor_omega\ 0.51\ -i\ cgs\ -tol\ 1.e-4\ -print\ 2\ -conv_cond\ 2\ -tol_w\ 0 \n")
            fh.write("[../]\n")

    def write_pqr(self, bound):
        i = 1
        with open("float.pqr", "w") as fh:
            for p in bound.float_bnd_xyzrcp:
                header = "ATOM      %5d  O   LYS 0    " % i
                fh.write("%-30s%8.3f%8.3f%8.3f%8.4f%8.4f\n" % (header, p.x, p.y, p.z,p.c,p.r))
                i += 1
        i = 1
        with open("single.pqr", "w") as fh:
            for p in bound.single_bnd_xyzrcp:
                header = "ATOM      %5d  O   LYS 0    " % i
                fh.write("%-30s%8.3f%8.3f%8.3f%8.4f%8.4f\n" % (header, p.x, p.y, p.z,p.c,p.r))
                i += 1
        i = 1
        with open("multi.pqr", "w") as fh:
            for p in bound.multi_bnd_xyzrcp:
                header = "ATOM      %5d  O   LYS 0    " % i
                fh.write("%-30s%8.3f%8.3f%8.3f%8.4f%8.4f\n" % (header, p.x, p.y, p.z,p.c,p.r))
                i += 1
        return
    
    def collect_energy(self, log: str):
        #log_lines = log.split("\n")
        #logger.info("Number of log_lines %d" % len(log_lines))

        found_pol  = False
        pol = None
        for line in log.splitlines():
            if "Polarization energy" in line:
                pol = float(line.split(":")[1].strip().split()[0]) / self.KCAL2KT
                found_pol = True
                break

        if not found_pol:
            logger.critical(f"Did not detect Polarization energy line. NGPB failed!\n{log}\n")
            sys.exit(1)

        return pol

    def collect_phi(self, xyzrcp):
        """Collect results from the log"""
        try:
            with open("phi_on_atoms.txt") as fh:
                lines = fh.readlines()
        except OSError:
            logger.critical("Could not open NGPB output file phi_on_atoms.txt.")
            sys.exit(1)

        expected_lines = len(xyzrcp)
        if len(lines) < expected_lines:
            logger.critical(
                "NGPB output phi_on_atoms.txt has fewer lines than expected: "
                "%d lines for %d atoms." % (len(lines), expected_lines)
            )
            sys.exit(1)

        for counter in range(expected_lines):
            parts = lines[counter].split()
            if not parts:
                logger.critical(
                    "NGPB output phi_on_atoms.txt has an empty line at line %d."
                    % (counter + 1)
                )
                sys.exit(1)

            try:
                phi = float(parts[-1])
            except ValueError:
                logger.critical(
                    "Could not parse NGPB atom potential from phi_on_atoms.txt "
                    "line %d: %r" % (counter + 1, lines[counter].rstrip("\n"))
                )
                sys.exit(1)
            xyzrcp[counter].p = phi

        return
    
    def run(self, bound, run_options: dict):
        """PBE solver interface for NGPB. 
        Returns rxn0 and rxn as the most negative value of the focusing runs.
        It will generate site p in both boundary conditions and
        return reference rxn0, and rxn in single boundary condition.
        Args:
            bound (Exchange?): dielectric boundary object
            run_options (dict) - command options in dictionary
        Returns:
          A 2-tuple: rxn0, rxn
        """
        # snippets to check the input and environment
        # Current working directory
        # cwd = os.getcwd()
        # print(cwd)
        # What are in bound
        #print(vars(bound))

        # Caclulate rxn0 using float boundary condition
        rxn0 = 0.0
        self.write_pqr(bound)
        self.epsilon_prot = run_options.d
        self.epsilon_solv = run_options.do
        self.salt=run_options.salt
        num_spacing = 20
        if run_options.fly:
            weight = 0.0
            center = [0.0, 0.0, 0.0]
            for p in bound.float_bnd_xyzrcp:
                w = abs(p.c)
                if w > 0.00001:
                    center[0] += p.x * w
                    center[1] += p.y * w
                    center[2] += p.z * w
                    weight += w

            if weight > 0.000001:
                center = [c/(weight+0.000001) for c in center]
            else:
                logger.error("PB solver shouldn't run a conformer with no charged atom.")
            
            self.write_option_file(center,num_spacing,1,0)

            # ngpb run
            command = []
            command = self.exe.copy()
            command.append('float.pqr')
            result = subprocess.run(command, capture_output=True, text=True, env=self.my_env)
            rxn0 = self.collect_energy(result.stdout)

        rxn = 0.0
        # single side chain boundary condition
        weight = 0.0
        center = [0.0, 0.0, 0.0]
        for p in bound.single_bnd_xyzrcp:
            w = abs(p.c)
            if w > 0.00001:
                center[0] += p.x * w
                center[1] += p.y * w
                center[2] += p.z * w
                weight += w

        if weight > 0.000001:
            center = [c/(weight+0.000001) for c in center]
        else:
            logger.error("PB solver shouldn't run a conformer with no charged atom.")
        
        self.write_option_file(center, num_spacing, 1, 1)
        command = []
        command = self.exe.copy()
        command.append('single.pqr')
        result = subprocess.run(command, capture_output=True, text=True, env=self.my_env)
        rxn = self.collect_energy(result.stdout)
        self.collect_phi(bound.single_bnd_xyzrcp)

        # multi side chain boundary condition
        weight = 0.0
        center = [0.0, 0.0, 0.0]
        for p in bound.multi_bnd_xyzrcp:
            w = abs(p.c)
            if w > 0.00001:
                center[0] += p.x * w
                center[1] += p.y * w
                center[2] += p.z * w
                weight += w

        if weight > 0.000001:
            center = [c/(weight+0.000001) for c in center]
        else:
            logger.error("PB solver shouldn't run a conformer with no charged atom.")

        self.write_option_file(center,num_spacing,0,1);
        command = []
        command = self.exe.copy()
        command.append('multi.pqr')
        result = subprocess.run(command, capture_output=True, text=True, env=self.my_env)
        # Test result, if result is error, exit
        if result.returncode != 0:
            logger.critical(f"NGPB failed with error:\n{result.stderr}")
            sys.exit(1)

        self.collect_phi(bound.multi_bnd_xyzrcp)

        return rxn0, rxn


class PBS_ZAP:
    """
    Zap interface
    """
    def __init__(self):
        try:
            from openeye import oechem, oezap
            self.oechem = oechem
            self.oezap = oezap
        except ModuleNotFoundError:
            logger.critical(("Error: Module openeye.oechem & openeye.zap not found\n"
                          "***** NOTE: If calling ZAP, please ensure a dedicated conda enviroment "
                          "for openeye and an openeye license path is setup prior to running step3 ****\n"
                          )
            )
            sys.exit(1)

        # consider loading these parameters from a file
        self.KCalsPerKT = 0.59
        self.KCalsPerSqAngstrom = 0.025
        self.scale = 2
        self.epsilon_prot = 4.0  
        self.epsilon_solv = 80.0
        self.salt = 0.050
        self.KCAL2KT = 1.688
        self.T = 298.15
        self.radius_probe = 1.4
        self.InterfaceData = """ """
        #self.pdbfile = ""
        return

    def InterfaceData_writer(self, pdbfile):
        self.InterfaceData = f"""
        !PARAMETER -in
          !TYPE string
          !BRIEF Input molecule file.
          !REQUIRED false
          !DEFAULT {pdbfile}
          !KEYLESS 1
        !END
        
        !PARAMETER -file_charges
          !TYPE bool
          !DEFAULT false
          !BRIEF Use partial charges from input file rather than calculating with MMFF.
        !END

        !PARAMETER -calc_type
          !TYPE string
          !DEFAULT solvent_only
          !LEGAL_VALUE solvent_only
          !LEGAL_VALUE remove_self
          !LEGAL_VALUE coulombic
          !LEGAL_VALUE breakdown
          !BRIEF Choose type of atom potentials to calculate
        !END

        !PARAMETER -atomtable
          !TYPE bool
          !DEFAULT true
          !BRIEF Output a table of atom potentials
        !END

        !PARAMETER -epsin
          !TYPE float
          !BRIEF Inner dielectric
          !DEFAULT {self.epsilon_prot}
          !LEGAL_RANGE 0.0 100.0
        !END

        !PARAMETER -grid_spacing
          !TYPE float
          !DEFAULT 0.5
          !BRIEF Spacing between grid points (Angstroms)
          !LEGAL_RANGE 0.1 2.0
        !END

        !PARAMETER -boundary
          !ALIAS -buffer
          !TYPE float
          !DEFAULT 2.0
          !BRIEF Extra buffer outside extents of molecule.
          !LEGAL_RANGE 0.1 10.0
        !END
        """
        return True

    def SetupInterface(self, itf, InterfaceData):
        self.oechem.OEConfigure(itf, self.InterfaceData)
        
        return True      

    def xyzrc_to_pdb(self, bound):
        i = 1
        element = "X"
        with open("float.pdb", "w") as fh:
            for p in bound.float_bnd_xyzrcp:
                fh.write("ATOM  {:5d}  {:2s}  MOL     1    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.3f}          {:>2s}\n".format(
                i, element, p.x, p.y, p.z,p.r, p.c, element))
                i += 1
        i = 1
        with open("single.pdb", "w") as fh:
            for p in bound.single_bnd_xyzrcp:
                fh.write("ATOM  {:5d}  {:2s}  MOL     1    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.3f}          {:>2s}\n".format(
                i, element, p.x, p.y, p.z,p.r, p.c, element))
                i += 1
        i = 1
        with open("multi.pdb", "w") as fh:
            for p in bound.multi_bnd_xyzrcp:
                fh.write("ATOM  {:5d}  {:2s}  MOL     1    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.3f}          {:>2s}\n".format(
                i, element, p.x, p.y, p.z,p.r, p.c, element))
                i += 1
        return
    
    def collect_phi(self, phi_zap, xyzrcp):
        # collect results from the log
        for counter in range(len(xyzrcp)):
            xyzrcp[counter].p = phi_zap[counter]

    def run(self, bound, run_options):        
        """PBE solver interface for ZAP. 
        Input:
            bound - dielectric boundary object
            run_options - command options in dictionary
        Returns a 2-tuple: rxn0, rxn
            It will generate site p in both boundary conditions and
            return reference rxn0, and rxn in single boundary condition.
        """
        # Calculate rxn0 using float boundary condition
        rxn0 = 0.0
        self.epsilon_prot = run_options.d
        self.epsilon_solv = run_options.do
        salt_concentration = run_options.salt
        if salt_concentration > 0.05:
           logger.warning(("\n!!!!! WARNING ZAP !!!!!\nZAP is unstable at salt_concentrations "
                           "over 0.05\nFor better results please set flag when running STEP3: "
                           "$ step3.py -s ZAP -salt 0.05 \n!!!!!!!!!!\n")
           )
        self.xyzrc_to_pdb(bound)

        #### Here we set target focusing of the conformer ---> Results are worse with focusing
        #itf = oechem.OEInterface()
        #self.InterfaceData_writer("float.pdb")
        #if not self.SetupInterface(itf, self.InterfaceData):
        #       logger.critical("PB solver Interface not recognized.")
        #mol_focus = oechem.OEGraphMol()
        #ifs = oechem.oemolistream()
        #ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_DELPHI)
        #if not ifs.open(itf.GetString("-in")):
        #       oechem.OEThrow.Fatal("Unable to open %s for reading" % itf.GetString("-in"))
        #oechem.OEReadMolecule(ifs, mol_focus)

        if run_options.fly:
            weight = 0.0
            center = [0.0, 0.0, 0.0]
            for p in bound.float_bnd_xyzrcp:
                w = abs(p.c)
                if w > 0.00001:
                    center[0] += p.x * w
                    center[1] += p.y * w
                    center[2] += p.z * w
                    weight += w

            if weight > 0.000001:
                center = [c/(weight+0.000001) for c in center]
            else:
                logger.error("PB solver shouldn't run a conformer with no charged atom.")
            
            # Zap run: float
            itf = self.oechem.OEInterface()
            self.InterfaceData_writer("float.pdb")
            if not self.SetupInterface(itf, self.InterfaceData):
               logger.critical("PB solver Interface not recognized.")
               sys.exit(1)
            
            mol = self.oechem.OEGraphMol()
            ifs = self.oechem.oemolistream()
            ifs.SetFlavor(self.oechem.OEFormat_PDB, self.oechem.OEIFlavor_PDB_DELPHI)
            if not ifs.open(itf.GetString("-in")):
                self.oechem.OEThrow.Fatal("Unable to open %s for reading" % itf.GetString("-in"))

            self.oechem.OEReadMolecule(ifs, mol)
            zap = self.oezap.OEZap()
            zap.SetMolecule(mol)
            zap.SetInnerDielectric(itf.GetFloat("-epsin"))
            zap.SetBoundarySpacing(itf.GetFloat("-boundary"))
            zap.SetGridSpacing(itf.GetFloat("-grid_spacing"))
            zap.SetSaltConcentration(salt_concentration)
           
            solv = zap.CalcSolvationEnergy()
            print(f"Solv(kcal) = {self.KCalsPerKT * solv:.3f}")

            rxn0 = solv * self.KCalsPerKT
 
        ########################################################
        # single side chain boundary condition
        rxn = 0.0
        weight = 0.0
        center = [0.0, 0.0, 0.0]
        for p in bound.single_bnd_xyzrcp:
            w = abs(p.c)
            if w > 0.00001:
                center[0] += p.x * w
                center[1] += p.y * w
                center[2] += p.z * w
                weight += w

        if weight > 0.000001:
            center = [c/(weight+0.000001) for c in center]
        else:
            logger.error("PB solver shouldn't run a conformer with no charged atom.")
        
        # Zap run: Single
        itf = self.oechem.OEInterface() 
        self.InterfaceData_writer("single.pdb")
        if not self.SetupInterface(itf, self.InterfaceData):
            logger.critical("PB solver Interface not recognized.")
            sys.exit(1)

        mol = self.oechem.OEGraphMol()
        ifs = self.oechem.oemolistream()
        ifs.SetFlavor(self.oechem.OEFormat_PDB, self.oechem.OEIFlavor_PDB_DELPHI)
        if not ifs.open(itf.GetString("-in")):
            self.oechem.OEThrow.Fatal("Unable to open %s for reading" % itf.GetString("-in"))

        self.oechem.OEReadMolecule(ifs, mol)
        zap = self.oezap.OEZap()
        zap.SetMolecule(mol)

        # Uncomment below line for focusing ---> Results are worse
        #zap.SetFocusTarget(mol_focus)     
        zap.SetInnerDielectric(itf.GetFloat("-epsin"))
        zap.SetBoundarySpacing(itf.GetFloat("-boundary"))
        zap.SetGridSpacing(itf.GetFloat("-grid_spacing"))
        zap.SetSaltConcentration(salt_concentration)

        solv = zap.CalcSolvationEnergy()
        print(f"Solv(kcal) = {self.KCalsPerKT * solv:.3f}")
        rxn = solv * self.KCalsPerKT

        apot = self.oechem.OEFloatArray(mol.GetMaxAtomIdx())
        zap.CalcAtomPotentials(apot)
        self.collect_phi(apot, bound.single_bnd_xyzrcp)

        ########################################################
        # multi side chain boundary condition
        weight = 0.0
        center = [0.0, 0.0, 0.0]
        for p in bound.multi_bnd_xyzrcp:
            w = abs(p.c)
            if w > 0.00001:
                center[0] += p.x * w
                center[1] += p.y * w
                center[2] += p.z * w
                weight += w

        if weight > 0.000001:
            center = [c/(weight+0.000001) for c in center]
        else:
            logger.error("PB solver shouldn't run a conformer with no charged atom.")
        
        # Zap run: Multi
        itf = self.oechem.OEInterface()
        self.InterfaceData_writer("multi.pdb")
        if not self.SetupInterface(itf, self.InterfaceData):
            logger.critical("PB solver Interface not recognized.")
            sys.exit(1)

        mol = self.oechem.OEGraphMol()
        ifs = self.oechem.oemolistream()
        ifs.SetFlavor(self.oechem.OEFormat_PDB, self.oechem.OEIFlavor_PDB_DELPHI)
        if not ifs.open(itf.GetString("-in")):
               self.oechem.OEThrow.Fatal("Unable to open %s for reading" % itf.GetString("-in"))

        self.oechem.OEReadMolecule(ifs, mol)
        zap = self.oezap.OEZap()
        zap.SetMolecule(mol)
        # Uncomment below line for focusing ----> Results are worse
        #zap.SetFocusTarget(mol_focus)
        zap.SetInnerDielectric(itf.GetFloat("-epsin"))
        zap.SetBoundarySpacing(itf.GetFloat("-boundary"))
        zap.SetGridSpacing(itf.GetFloat("-grid_spacing"))
        zap.SetSaltConcentration(salt_concentration)

        apot = self.oechem.OEFloatArray(mol.GetMaxAtomIdx())
        zap.CalcAtomPotentials(apot)
        self.collect_phi(apot, bound.multi_bnd_xyzrcp)

        return (rxn0, rxn)

class PBS_ML:
    """
    Machine Learning PBS interface
    """ 
    def __init__(self):
        self.epsilon_prot = 4.0  # default value, be overwritten by run_options
        self.epsilon_solv = 80.0

         
    def run(self, bound, run_options):
        """PBE solver interface for Machine Learning.
        Generate site p in both boundary conditions and
        return reference rxn0, and rxn in single boundary condition.
        Input:
            bound - dielectric boundary object
            run_options - command options in dictionary
        Returns:
            -  A 2-tuple: rxn0, rxn: the most negative value of the focusing runs 
        """
        from ml_pbs import pqr2rxn, pqr2pot, K_rxn


        # Update dielectric constants for ML
        e_prot = run_options.d
        e_solv = run_options.do

        # Use float_bnd_xyzrcp for rxn0
        pqrs = np.array([[p.x, p.y, p.z, p.c, p.r] for p in bound.float_bnd_xyzrcp])
        rxn0 =  pqr2rxn(pqrs, e_prot, e_solv)

        # Use single_bnd_xyzrcp for rxn
        pqrs = np.array([[p.x, p.y, p.z, p.c, p.r] for p in bound.single_bnd_xyzrcp])
        rxn = pqr2rxn(pqrs, e_prot, e_solv)
        site_potentials = pqr2pot(pqrs, e_prot, e_solv)
        for atom, pot in zip(bound.single_bnd_xyzrcp, site_potentials):
            atom.p = pot

        # Use multi_bnd_xyzrcp for site potential
        pqrs = np.array([[p.x, p.y, p.z, p.c, p.r] for p in bound.multi_bnd_xyzrcp])
        site_potentials = pqr2pot(pqrs, e_prot, e_solv)
        for atom, pot in zip(bound.multi_bnd_xyzrcp, site_potentials):
            atom.p = pot

        return rxn0, rxn


class PBS_APBS:
    """
    APBS interface
    """
    def __init__(self):
        # consider loading these parameters from a file
        self.KCalsPerKT = 0.59
        self.KCalsPerSqAngstrom = 0.025
        self.scale = 2
        self.epsilon_prot = 4.0
        self.epsilon_solv = 80.0
        self.salt = 0.05
        self.KJOL2KCAL = 4.184
        self.T = 297.334198  # seriously ??
        self.radius_probe = 1.4
        #self.pdbfile = ""
        return

    def write_pqr(self, bound):
        i = 1
        with open("float.pqr", "w") as fh:
            for p in bound.float_bnd_xyzrcp:
                header = "ATOM      %5d  O   LYS 0    " % i
                fh.write("%-30s%8.3f%8.3f%8.3f%8.4f%8.4f\n" % (header, p.x, p.y, p.z,p.c,p.r))
                i += 1
        i = 1
        with open("single.pqr", "w") as fh:
            for p in bound.single_bnd_xyzrcp:
                header = "ATOM      %5d  O   LYS 0    " % i
                fh.write("%-30s%8.3f%8.3f%8.3f%8.4f%8.4f\n" % (header, p.x, p.y, p.z,p.c,p.r))
                i += 1
        i = 1
        with open("multi.pqr", "w") as fh:
            for p in bound.multi_bnd_xyzrcp:
                header = "ATOM      %5d  O   LYS 0    " % i
                fh.write("%-30s%8.3f%8.3f%8.3f%8.4f%8.4f\n" % (header, p.x, p.y, p.z,p.c,p.r))
                i += 1
        return

    def collect_energy(self, log):
        log_lines = log.split("\n")
        found_pol  = False

        for line in log_lines:
            if "Global net ELEC energy =" in line:
                pol = float(line[27:].split()[0])
                found_pol = True

        if not found_pol:
           logger.critical("Did not detect Polarization energy line. APBS failed!")
           print(log)
           sys.exit(1)

        return pol

    def collect_phi(self, xyzrcp):
        # collect results from the log
        try:
            lines = open("pot.txt", "r").readlines()
        except OSError:
            logger.critical("Could not open APBS output file pot.txt.")
            sys.exit(1)

        for counter in range(len(xyzrcp)):
            phi = float(lines[counter+4][:])
            xyzrcp[counter].p = phi

        return

    def get_molecule_dimensions_from_pqr(self, pqr, padding=5.0, grid_spacing=0.5):
        # Read the PQR file and extract atom coordinates
        x,y,z,r = [], [], [], []
        with open(pqr, 'r') as f:
             for line in f:
                 if line.startswith('ATOM'):
                    x.append(float(line[30:38].strip()))
                    y.append(float(line[38:46].strip()))
                    z.append(float(line[46:54].strip()))
                    r.append(float(line[62:70].strip()))

        min_x, min_y, min_z = min(x), min(y), min(z)
        max_x, max_y, max_z, max_r = max(x), max(y), max(z), max(r)
        #print(f"min_xyz = {min_x, min_y, min_z}, max_xyz = {max_x, max_y, max_z, max_r}")

        length_x, length_y, length_z = (max_x - min_x) + 2*max_r,  (max_y - min_y) + 2*max_r, (max_z - min_z) + 2*max_r
        cent_x, cent_y, cent_z = (max_x + min_x) / 2,  (max_y + min_y) / 2, (max_z + min_z) / 2
        print(f"geo_center = {cent_x}, {cent_y}, {cent_z}")
        cglen = max([length_x, length_y, length_z]) + 2*padding
        fglen = max([length_x, length_y, length_z]) * 0.25
        dime = math.ceil(cglen/grid_spacing + 1)
        dime = 100

        return dime, cglen, fglen

    #def calculate_grid_parameters(self, min_coords, max_coords, padding=5.0, coarseness=2.0):
    #    lengths = max_coords - min_coords + 2 * padding
    #    dime = np.ceil(lengths / coarseness).astype(int) + 1  # Grid dimensions
    #    print(dime)
    #    return dime

    def APBS_input_writer(self, pqr, dime, cglen, fglen, cent):
         input_string = f"""
         # READ IN MOLECULES
         read
           mol pqr {pqr}
         end
         elec name solv                        # Electrostatics calculation on the solvated state
           mg-manual                           # Specify the mode for APBS to run
           dime {dime} {dime} {dime}           # The grid dimensions
           nlev 4                              # Multigrid level parameter
           grid 0.5 0.5 0.5                    # Grid spacing
           cglen  {cglen} {cglen} {cglen}
           fglen  {fglen} {fglen} {fglen}
           #gcent  {cent[0]} {cent[1]} {cent[2]} # Center manually on molecule/residue
           #cgcent {cent[0]} {cent[1]} {cent[2]}
           #fgcent {cent[0]} {cent[1]} {cent[2]}
           gcent  mol 1                        # Center the grid on molecule 1
           cgcent mol 1
           fgcent mol 1
           mol 1                               # Perform the calculation on molecule 1
           lpbe                                # Solve the linearized Poisson-Boltzmann equation
           bcfl mdh                            # Use all multipole moments when calculating the potential
           pdie {self.epsilon_prot}            # Solute dielectric
           sdie {self.epsilon_solv}            # Solvent dielectric
           chgm spl2                           # Spline-based discretization of the delta functions
           srfm spl2                           # Molecular surface definition
           srad 1.4                            # Solvent probe radius (for molecular surface)
           swin 0.3                            # Solvent surface spline window (not used here)
           sdens 10.0                          # Sphere density for accessibility object
           temp 298.15                         # Temperature
           calcenergy total                    # Calculate energies
           calcforce no                        # Do not calculate forces
           ion charge +1 conc {self.salt} radius 2.0 # Positive ion (e.g., Na+)
           ion charge -1 conc {self.salt} radius 2.0 # Negative ion (e.g., Cl-)
           write atompot flat pot              # write potential on atoms
         end

         elec name ref # Calculate potential for reference (vacuum) state
            mg-manual
            dime {dime} {dime} {dime}
            nlev 4
            grid 0.5 0.5 0.5
            cglen  {cglen} {cglen} {cglen}
            fglen  {fglen} {fglen} {fglen}
            #gcent  {cent[0]} {cent[1]} {cent[2]} # Center manually on molecule/residue
            #cgcent {cent[0]} {cent[1]} {cent[2]}
            #fgcent {cent[0]} {cent[1]} {cent[2]}
            gcent  mol 1                        # Center the grid on molecule 1
            cgcent mol 1
            fgcent mol 1
            mol 1
            lpbe
            bcfl mdh
            pdie {self.epsilon_prot}
            sdie {self.epsilon_prot}
            chgm spl2
            srfm spl2
            srad 1.4
            swin 0.3
            sdens 10.0
            temp 298.15
            calcenergy total
            calcforce no
            ion charge +1 conc {self.salt} radius 2.0
            ion charge -1 conc {self.salt} radius 2.0
            #write atompot flat ref_pot
         end


         # Calculate solvation energy
         print ElecEnergy solv - ref end
         quit
         """
         return input_string

    def run_apbs(self, pqr, dime, cglen, fglen, cent):
        # print("Welcome to the APBS Nightmare for MCCE4")

        # Create the APBS input and convert to a file
        apbs_input = self.APBS_input_writer(pqr, dime, cglen, fglen, cent)
        with open("apbs_input.in", "w") as f:
             f.write(apbs_input)

        # Run the APBS command
        result = subprocess.run(["apbs", "apbs_input.in"], capture_output=True, text=True)
        print(result.stdout)
        print(result.stderr)

        solv = self.collect_energy(result.stdout)

        return solv

    #def run0(self):
    #    dime = self.get_molecule_dimensions_from_pqr()
    #    self.run_apbs(dime)
    #    self.collect_phi()

    def run(self, bound, run_options):
        """PBE solver interface for ZAP.
        Input:
            bound - dielectric boundary object
            run_options - command options in dictionary
        Return value:
            It will generate site p in both boundary conditions and
            return reference rxn0, and rxn in single boundary condition.
        """
        self.epsilon_prot = run_options.d
        self.epsilon_solv = run_options.do
        self.salt = run_options.d.salt

        # Calculate rxn0 using float boundary condition
        rxn0 = 0.0
        salt_concentration = 0.05  # Example: 0.15 M (150 mM)
        self.write_pqr(bound)

        if run_options.fly:
            weight = 0.0
            center = [0.0, 0.0, 0.0]
            for p in bound.float_bnd_xyzrcp:
                w = abs(p.c)
                if w > 0.00001:
                    center[0] += p.x * w
                    center[1] += p.y * w
                    center[2] += p.z * w
                    weight += w

            if weight > 0.000001:
                center = [c/(weight+0.000001) for c in center]
            else:
                logger.error("PB solver shouldn't run a conformer with no charged atom.")
            print(f"charged_center = {center}")

            # APBS run: float
            pqr = "float.pqr"
            dime, cglen, fglen = self.get_molecule_dimensions_from_pqr(pqr)
            solv = self.run_apbs(pqr, dime, cglen, fglen, center)
            #print(f"Solv(kcal) = {solv / self.KJOL2KCAL}")

            rxn0 = solv / self.KJOL2KCAL

        ########################################################
        # single side chain boundary condition
        rxn = 0.0
        weight = 0.0
        center = [0.0, 0.0, 0.0]
        for p in bound.single_bnd_xyzrcp:
            w = abs(p.c)
            if w > 0.00001:
                center[0] += p.x * w
                center[1] += p.y * w
                center[2] += p.z * w
                weight += w

        if weight > 0.000001:
            center = [c/(weight+0.000001) for c in center]
        else:
            logger.error("PB solver shouldn't run a conformer with no charged atom.")
        print(f"charged_center = {center}")

        # APBS run: Single
        pqr = "single.pqr"
        dime, cglen, fglen = self.get_molecule_dimensions_from_pqr(pqr)
        solv = self.run_apbs(pqr, dime, cglen, fglen, center)
        self.collect_phi(bound.single_bnd_xyzrcp)

        #print(f"Solv(kcal) = {solv / self.KJOL2KCAL}")

        rxn = solv / self.KJOL2KCAL

        ########################################################
        # multi side chain boundary condition
        weight = 0.0
        center = [0.0, 0.0, 0.0]
        for p in bound.multi_bnd_xyzrcp:
            w = abs(p.c)
            if w > 0.00001:
                center[0] += p.x * w
                center[1] += p.y * w
                center[2] += p.z * w
                weight += w

        if weight > 0.000001:
            center = [c/(weight+0.000001) for c in center]
        else:
            logger.error("PB solver shouldn't run a conformer with no charged atom.")
        print(f"charged_center = {center}")

        # APBS run: Multi
        pqr = "multi.pqr"
        dime, cglen, fglen = self.get_molecule_dimensions_from_pqr(pqr)
        solv = self.run_apbs(pqr, dime, cglen, fglen, center)
        self.collect_phi(bound.multi_bnd_xyzrcp)

        return rxn0, rxn
