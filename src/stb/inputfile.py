#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
try:
    from importlib.metadata import version as _pkg_version
    VERSION = _pkg_version("stb_suite")
except Exception:
    VERSION = "1.9.5"    
from stb.cli import color_text, show_intro, print_info, print_ok, print_warn, print_error

import numpy as np
import math
import sys
import os
import shutil
import argparse



CALC_RELAX_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include structure.fdf

#%block Geometry.Constraints
  #stress 1  # Fixes XX 
  #stress 2  # Fixes YY
  #stress 3  # Fixes ZZ 
  #stress 4  # Fixes YZ 
  #stress 5  # Fixes XZ 
  #stress 6  # Fixes XY 
#%endblock Geometry.Constraints

## Optional Parameters
  # %block Supercell
  # %endblock Supercell
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       DZP
PAO.EnergyShift     0.01 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH 
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [6  6  6]
Mesh.CutOff           320  Ry  
FilterCutoff          150  Ry

## Optional Parameters
  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-5  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.

## Optional Parameters
  # SCF.MustConverge true 
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   false

## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## STRUCTURE RELAXATION (MOLECULAR DYNAMICS - MD)
## ===================================================================
##Essential Parameters
MD.TypeOfRun            CG
MD.VariableCell         true
MD.MaxForceTol          0.01 eV/Ang
MD.Steps                150

## Optional Parameters
  # MD.RelaxCellOnly true
  # Constant.Volume false
  # MD.MaxStressTol 1 GPa
  # MD.MaxDispl 0.2 Bohr
  # MD.PreconditionVariableCell 5 Ang
  # MD.Broyden.History.Steps 5
  # MD.Broyden.Cycle.On.Mait true
  # Target.Pressure 0 GPa
  # MD.RemoveIntramáscularPressure false

## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates 
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize 
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.  #(Use previous Density Matrix if available)

## Optional Parameters (Set to false or inactive)
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.  
  # ON.UseSaveLWF  .false.  

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteMDHistory          .true.
WriteCoorStep           .true.
WriteForces             .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.

## Optional Parameters (Interface flags)
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""
# --- End of Relax Template ---


# --- Internal Total Energy Template ---
CALC_TOTAL_ENERGY_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include struct.fdf

## Optional Parameters
  # %block Supercell
  # %endblock Supercell
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       DZP
PAO.EnergyShift     0.02 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [4  13  1]
Mesh.CutOff           320  Ry   
FilterCutoff          150  Ry

## Optional Parameters
  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-5  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.

## Optional Parameters
  # SCF.MustConverge true
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   .false.
## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.
## Optional Parameters (Set to false or inactive)
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.
  # ON.UseSaveLWF  .false.

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteMDHistory          .true.
WriteCoorStep           .true.
WriteForces             .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.
## Optional Parameters (Interface flags)
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""
# --- End of Total Energy Template ---


# --- Internal AIMD Template ---
CALC_AIMD_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include struct.fdf

%block Supercell
   2   0   0
   0   2   0
   0   0   2
%endblock Supercell

#%block Geometry.Constraints
  #stress 1  # Fixes XX 
  #stress 2  # Fixes YY
  #stress 3  # Fixes ZZ 
  #stress 4  # Fixes YZ 
  #stress 5  # Fixes XZ 
  #stress 6  # Fixes XY 
#%endblock Geometry.Constraints


## Optional Parameters
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       SZ
PAO.EnergyShift     0.02 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [1  1  1]
Mesh.CutOff           320  Ry    
FilterCutoff          150  Ry

## Optional Parameters
  
  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-4  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.

## Optional Parameters
  # SCF.MustConverge true
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   .false.
## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## STRUCTURE RELAXATION (MOLECULAR DYNAMICS - MD)
## ===================================================================
##Essential Parameters
MD.TypeOfRun  Nose
MD.InitialTimeStep  1
MD.FinalTimeStep  5000
MD.LengthTimeStep  1.0  fs
MD.InitialTemperature  300.0  K
MD.TargetTemperature   300.0  K
MD.NoseMass  30.0  Ry*fs**2

## Optional Parameters
  # MD.RelaxCellOnly true
  # Constant.Volume false
  # MD.MaxStressTol 1 GPa
  # MD.MaxDispl 0.2 Bohr
  # MD.PreconditionVariableCell 5 Ang
  # MD.Broyden.History.Steps 5
  # MD.Broyden.Cycle.On.Mait true
  # Target.Pressure 0 GPa
  # MD.RemoveIntramáscularPressure false

## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.
MD.UseSaveXV            .true.
## Optional Parameters (Set to false or inactive)
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.
  # ON.UseSaveLWF  .false.

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteMDHistory          .true.
WriteCoorStep           .true.
WriteForces             .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.
## Optional Parameters (Interface flags)
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""
# --- End of AIMD Template ---


# --- Internal Band Structure Template ---
CALC_BANDS_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include struct.fdf

## Optional Parameters
  #%block Supercell  %endblock Supercell
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       DZP
PAO.EnergyShift     0.02 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [1  1  1]
Mesh.CutOff           320  Ry 
FilterCutoff          150  Ry

## Optional Parameters

  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-5  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.

## Optional Parameters
  # SCF.MustConverge true
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   .false.
## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## BAND STRUCTURE AND DENSITY OF STATES
## ===================================================================
##Essential Parameters
#___
#IMPORTANT: k-path file created by STB-SUITE.
#For manual k-path use optional
#Parameters "%block BandLines" and comment the bellow line.
%include kpath_bs.fdf
#___

%PDOS.kgrid_Monkhorst_Pack
 20     0    0  0.0
  0    20    0  0.0
  0     0   20  0.0
%end PDOS.kgrid_Monkhorst_Pack

%block  ProjectedDensityOfStates
-20.0  20.0  0.01  8000  eV
%endblock  ProjectedDensityOfStates

## Optional Parameters
#BandLinesScale  ReciprocalLatticeVectors
# %block BandLines
#  Important: Include here the k-path
# %endblock BandLines


## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.
## Optional Parameters (Set to false or inactive)
  # MD.UseSaveXV  .true.
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.
  # ON.UseSaveLWF .false.

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.
## Optional Parameters (Interface flags)
  #WriteMDHistory          .true.
  #WriteCoorStep           .true.
  #WriteForces             .true.
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""


def parse_structure_fdf(filename):
    """
    Parses a .fdf (Siesta) file and returns the lattice vectors and
    the list of chemical species symbols.
    Modified to 'raise Exception' instead of 'sys.exit'.
    """
    in_lattice_block = False
    in_species_block = False
    all_values = []
    species_list = []
    lattice_constant = 1.0 # Default value

    try:
        with open(filename, 'r') as f:
            for line in f:
                cleaned_line = line.split('#', 1)[0].strip()
                if not cleaned_line: continue
                parts = cleaned_line.split()
                lower_line = cleaned_line.lower()

                if lower_line.startswith('latticeconstant'):
                    try:
                        lattice_constant = float(parts[1])
                    except (IndexError, ValueError):
                        raise ValueError(f"'LatticeConstant' line malformed: {line}")
                    continue 

                if lower_line == '%block latticevectors':
                    in_lattice_block = True; continue 
                if lower_line == '%endblock latticevectors':
                    in_lattice_block = False; continue 
                if lower_line == '%block chemicalspecieslabel':
                    in_species_block = True; continue
                if lower_line == '%endblock chemicalspecieslabel':
                    in_species_block = False; continue

                if in_lattice_block:
                    for part in parts:
                        try: all_values.append(float(part))
                        except ValueError: pass
                
                if in_species_block:
                    try:
                        symbol = parts[2] # ex: 'C', 'B', 'N'
                        if symbol not in species_list:
                            species_list.append(symbol)
                    except IndexError: pass
        
        if len(all_values) != 9:
            raise ValueError(f"Expected 9 values in LatticeVectors block, found {len(all_values)}.")
            
        lattice = np.array(all_values).reshape(3, 3) * lattice_constant
        
        if not species_list:
            raise ValueError("No chemical species found in ChemicalSpeciesLabel block.")

        return lattice, species_list

    except FileNotFoundError:
        raise FileNotFoundError(f"Structure file '{filename}' not found.")
    except Exception as e:
        # Re-raise other errors
        raise e

def compute_monkhorts(cella, cellb, cellc, k_density):
    """
    Calculates the reciprocal vectors and the number of Monkhorst-Pack divisions.
    Modified to 'raise' errors.
    """
    volume = np.dot(cella, np.cross(cellb, cellc))
    
    if abs(volume) < 1e-9:
        raise ValueError("Cell volume is zero. Check lattice vectors.")
        
    b1 = 2 * np.pi * np.cross(cellb, cellc) / volume
    b2 = 2 * np.pi * np.cross(cellc, cella) / volume
    b3 = 2 * np.pi * np.cross(cella, cellb) / volume

    lengths = [np.linalg.norm(b) for b in (b1, b2, b3)]
    divisions = [max(1, math.ceil(length / k_density)) for length in lengths]
    return divisions

def _section(title: str) -> None:
    """Print a styled section separator."""
    bar = color_text("─" * 60, 'white')
    print(f"\n{bar}")
    print(color_text(f"  {title}", 'bold'))
    print(bar)


def copy_pseudopotentials(species_list, pp_path):
    """
    Copies required .psml or .psf files from pp_path to the current directory.
    Prioritises .psml. Returns a list of (ok, message) tuples.
    """
    _section("Pseudopotentials")

    if not os.path.isdir(pp_path):
        print_warn(f"PP path '{pp_path}' is not a valid directory — skipping copy.")
        return

    copied, failed = 0, 0
    for symbol in species_list:
        psml = f"{symbol}.psml"
        psf  = f"{symbol}.psf"
        src_psml = os.path.join(pp_path, psml)
        src_psf  = os.path.join(pp_path, psf)

        if os.path.exists(src_psml):
            try:
                shutil.copy2(src_psml, psml)
                print_ok(f"Copied  {psml}")
                copied += 1
            except Exception as e:
                print_error(f"Failed to copy {psml}: {e}")
                failed += 1
        elif os.path.exists(src_psf):
            try:
                shutil.copy2(src_psf, psf)
                print_ok(f"Copied  {psf}")
                copied += 1
            except Exception as e:
                print_error(f"Failed to copy {psf}: {e}")
                failed += 1
        else:
            print_warn(f"Not found: {psml} or {psf} in '{pp_path}'")
            failed += 1

    print()
    if failed == 0:
        print_ok(f"All pseudopotentials copied ({copied}/{len(species_list)})")
    else:
        print_warn(f"Copied {copied}/{len(species_list)} — {failed} missing")


def generate_calculation(struct_file, chosen_mode, pp_path):
    """Main generation logic: parse structure, fill template, write calc.fdf."""
    output_file = "calc.fdf"

    # ── Header ────────────────────────────────────────────────────
    bar = color_text("═" * 60, 'white')
    print(f"\n{bar}")
    print(color_text("  SIESTA Input File Generator", 'bold'))
    print(bar)

    def _row(label, value):
        print(f"  {color_text(label, 'cyan')}  {value}")

    _row("Structure :", struct_file)
    _row("Mode      :", color_text(chosen_mode, 'yellow'))
    _row("PP path   :", pp_path if pp_path else color_text("(none)", 'yellow'))
    _row("Output    :", output_file)
    print()

    try:
        # ── 1. Validate input ──────────────────────────────────────
        if not os.path.exists(struct_file):
            print_error(f"Structure file '{struct_file}' not found.")
            return

        # ── 2. Select template ─────────────────────────────────────
        _section("Template")
        templates = {
            'relax':        ('Structural Relaxation',  CALC_RELAX_TEMPLATE),
            'relax+d3':     ('Structural Relaxation + DFT-D3', CALC_RELAX_TEMPLATE),
            'total_energy': ('Total Energy',           CALC_TOTAL_ENERGY_TEMPLATE),
            'total_energy+d3': ('Total Energy + DFT-D3', CALC_TOTAL_ENERGY_TEMPLATE),
            'aimd':         ('Ab Initio MD',           CALC_AIMD_TEMPLATE),
            'aimd+d3':      ('Ab Initio MD + DFT-D3',  CALC_AIMD_TEMPLATE),
            'bands':        ('Band Structure',         CALC_BANDS_TEMPLATE),
            'bands+d3':     ('Band Structure + DFT-D3', CALC_BANDS_TEMPLATE),
        }
        label, template_string = templates[chosen_mode]
        print_info(f"Template  : {label}")

        use_d3 = chosen_mode.endswith('+d3')
        d3_status = color_text('ENABLED', 'green') if use_d3 else color_text('disabled', 'yellow')
        print_info(f"DFT-D3    : {d3_status}")
        d3_line_new = f"DFTD3                   {'.true.' if use_d3 else '.false.'}"

        # ── 3. Parse structure ─────────────────────────────────────
        _section("Structure")
        lattice, species = parse_structure_fdf(struct_file)
        vol = abs(np.dot(lattice[0], np.cross(lattice[1], lattice[2])))
        a, b, c = [np.linalg.norm(v) for v in lattice]
        print_info(f"Species   : {color_text(', '.join(species), 'yellow')}")
        print_info(f"Cell a/b/c: {a:.4f}  {b:.4f}  {c:.4f}  Å")
        print_info(f"Volume    : {vol:.4f}  Å³")

        # ── 4. K-grid ──────────────────────────────────────────────
        _section("K-Grid")
        replace_kgrid = chosen_mode not in ('aimd', 'aimd+d3')
        kgrid_line_new = ""
        if replace_kgrid:
            k_density = 0.2
            kd = compute_monkhorts(lattice[0], lattice[1], lattice[2], k_density)
            kgrid_line_new = f"kgrid.MonkhorstPack   [{kd[0]}  {kd[1]}  {kd[2]}]"
            print_info(f"Density   : {k_density} Å⁻¹  →  grid {color_text(f'{kd[0]} × {kd[1]} × {kd[2]}', 'green')}")
        else:
            print_info("AIMD mode — k-grid kept from template (not modified)")

        # ── 5. Write output ────────────────────────────────────────
        _section("Writing Output")
        include_line_new = f"%include {os.path.basename(struct_file)}"
        output_lines = []
        for line in template_string.splitlines(keepends=True):
            low = line.strip().lower()
            if low.startswith('%include') and 'struct' in low:
                output_lines.append(include_line_new + '\n')
            elif low.startswith('kgrid.monkhorstpack') and replace_kgrid:
                output_lines.append(kgrid_line_new + '\n')
            elif low.startswith('dftd3'):
                output_lines.append(d3_line_new + '\n')
            else:
                output_lines.append(line)

        with open(output_file, 'w') as f:
            f.writelines(output_lines)

        print_info(f"Writing   : {output_file}  ({len(output_lines)} lines)")
        print()
        print_ok(f"'{output_file}' generated successfully.")

        # ── 6. Pseudopotentials ────────────────────────────────────
        if pp_path:
            copy_pseudopotentials(species, pp_path)

        # ── Footer ─────────────────────────────────────────────────
        print()
        print(color_text("═" * 60, 'white'))
        print(color_text("  Done. Edit calc.fdf to adjust calculation parameters.", 'bold'))
        print(color_text("═" * 60, 'white') + "\n")

    except Exception as e:
        print()
        print_error(f"An error occurred: {e}")
        print_error("Operation aborted.")

# --- Main Function (Argument Parser) ---
def main():
    """
    Processes the command-line arguments and calls the generation logic.
    """
    # --- Definition of valid modes ---
    mode_list = [
        'total_energy', 'total_energy+d3',
        'relax', 'relax+d3',
        'aimd', 'aimd+d3',
        'bands', 'bands+d3'
    ]

    # --- Configure argparse ---
    parser = argparse.ArgumentParser(
        description="Script to prepare a Siesta input file (.fdf).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s struct.fdf -t relax+d3 -p /path/to/pps\n"
               "  %(prog)s system.fdf -t bands\n"
               "  %(prog)s system.fdf --type aimd\n"
    )
    
    # Argument 1: Structure File (Required, Positional)
    parser.add_argument(
        "structure_file",
        type=str,
        help="Path to the structure file (e.g., struct.fdf)"
    )
    
    # Argument 2: Calculation Mode (Required, with Flag)
    parser.add_argument(
        "-t", "--type",
        type=str,
        required=True, # <-- The -t/--type flag is now mandatory
        choices=mode_list,
        dest="calc_type", # The value will be stored in 'args.calc_type'
        help=f"Calculation mode. Valid choices: {', '.join(mode_list)}"
    )
    
    # Argument 3: PPs Path (Optional, with flag)
    parser.add_argument(
        "-p", "--pp-path",
        type=str,
        default="",
        dest="pp_path",
        help="Path to the pseudopotentials folder - Only psml format (optional)"
    )

    
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-inputfile {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro()

    generate_calculation(args.structure_file, args.calc_type, args.pp_path)

# --- Script Entry Point ---
if __name__ == "__main__":
    main()
