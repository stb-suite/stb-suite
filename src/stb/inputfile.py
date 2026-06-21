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
from stb.cli import color_text, show_intro

import numpy as np
import math
import sys
import os
import shutil
import argparse



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

def copy_pseudopotentials(species_list, pp_path):
    """
    Copies the required .psml or .psf files from the 'pp_path' folder to
    the current working directory. Prioritizes .psml if both exist.
    Returns a list of warnings/errors.
    """
    warnings = []
    if not os.path.isdir(pp_path):
        warnings.append(f"Warning: PP path '{pp_path}' is not a valid directory. PPs not copied.")
        return warnings

    print("\n--- Copying Pseudopotentials ---")
    copied_files = 0
    for symbol in species_list:
        psml_filename = f"{symbol}.psml"
        psf_filename = f"{symbol}.psf"
        
        source_psml = os.path.join(pp_path, psml_filename)
        source_psf = os.path.join(pp_path, psf_filename)
        
        dest_psml = os.path.join(os.getcwd(), psml_filename)
        dest_psf = os.path.join(os.getcwd(), psf_filename)
        
        if os.path.exists(source_psml):
            try:
                shutil.copy2(source_psml, dest_psml)
                print(color_text(f"  [OK] Copied: {psml_filename}", 'green'))
                copied_files += 1
            except Exception as e:
                warnings.append(f"  [ERROR] Failed to copy {psml_filename}: {e}")
                
        elif os.path.exists(source_psf):
            try:
                shutil.copy2(source_psf, dest_psf)
                print(color_text(f"  [OK] Copied: {psf_filename}", 'green'))
                copied_files += 1
            except Exception as e:
                warnings.append(f"  [ERROR] Failed to copy {psf_filename}: {e}")
                
        else:
            warnings.append(f"  [WARNING] Neither {psml_filename} nor {psf_filename} found in {pp_path}")

    if copied_files == len(species_list):
        print("All pseudopotentials copied successfully.")
    else:
        warnings.append("Warning: Some PP files were not found or could not be copied.")
    print("----------------------------------")
    
    return warnings

# --- Main Generation Logic Function ---
def generate_calculation(struct_file, chosen_mode, pp_path):
    """
    Main function that executes the file generation logic.
    """
    try:
        # --- 1. Input Validation ---
        if not os.path.exists(struct_file):
            print(f"Error: Structure file '{struct_file}' not found.")
            return # Abort the function

        if not pp_path:
            print(color_text("[INFO] PP path is blank. Skipping pseudopotential copy.", 'cyan'))
        else:
            print(color_text(f"[INFO] PP path set to: {pp_path}", 'cyan'))

        print(color_text(f"[INFO] Mode selected: {chosen_mode}", 'cyan'))

        # --- 2. Start Generation Logic ---
        print(f"\nProcessing '{struct_file}'...")
        
        output_file = "calc.fdf" # Fixed output name
        
        # --- Select Template ---
        template_string = ""
        if chosen_mode in ['relax', 'relax+d3']:
            template_string = CALC_RELAX_TEMPLATE
            print(color_text("[INFO] Using 'Relax' template.", 'cyan'))
        elif chosen_mode in ['total_energy', 'total_energy+d3']:
            template_string = CALC_TOTAL_ENERGY_TEMPLATE
            print(color_text("[INFO] Using 'Total Energy' template.", 'cyan'))
        elif chosen_mode in ['aimd', 'aimd+d3']:
            template_string = CALC_AIMD_TEMPLATE
            print(color_text("[INFO] Using 'AIMD' template.", 'cyan'))
        elif chosen_mode in ['bands', 'bands+d3']:
            template_string = CALC_BANDS_TEMPLATE
            print(color_text("[INFO] Using 'Bands' template.", 'cyan'))
        
        # --- Set D3 flag ---
        d3_flag = ".false."
        if chosen_mode in ['relax+d3', 'total_energy+d3', 'aimd+d3', 'bands+d3']:
            d3_flag = ".true."
            print(color_text("[INFO] DFT-D3 (van der Waals) correction will be ENABLED.", 'cyan'))
        else:
            # 'relax', 'total_energy', 'aimd', or 'bands'
            print(color_text("[INFO] DFT-D3 (van der Waals) correction will be DISABLED.", 'cyan'))
        
        d3_line_new = f"DFTD3                   {d3_flag}"
        template_lines = template_string.splitlines(keepends=True)

        # Parse structure and get species
        lattice, species = parse_structure_fdf(struct_file)
        print(f"  Species found: {', '.join(species)}")

        # --- K-Grid Conditional Logic ---
        replace_kgrid = False
        kgrid_line_new = "" # Initialize
        
        if chosen_mode in ['aimd', 'aimd+d3']:
            # For AIMD, do nothing, keep the K-grid from the template
            replace_kgrid = False
            print(color_text("[INFO] AIMD mode selected. K-grid will NOT be modified.", 'cyan'))
        else:
            # For Total Energy, Relax, and Bands, calculate the K-grid
            replace_kgrid = True
            print(color_text("[INFO] Calculating K-grid (density = 0.2 1/Å)...", 'cyan'))
            k_density = 0.2
            kgrid_divs = compute_monkhorts(lattice[0], lattice[1], lattice[2], k_density)
            kgrid_line_new = f"kgrid.MonkhorstPack   [{kgrid_divs[0]}  {kgrid_divs[1]}  {kgrid_divs[2]}]"
            print(f"  Suggested K-grid: {kgrid_divs[0]} {kgrid_divs[1]} {kgrid_divs[2]}")

        
        # Use only the base file name for the include
        include_line_new = f"%include {os.path.basename(struct_file)}"

        # Process the template and replace lines
        print(f"Writing '{output_file}'...")
        output_lines = []
        for line in template_lines:
            line_stripped_lower = line.strip().lower()
            
            # Replace the structure include line
            if line_stripped_lower.startswith('%include') and 'struct' in line_stripped_lower:
                output_lines.append(include_line_new + '\n')
            
            # Replace the k-grid (if not AIMD)
            elif line_stripped_lower.startswith('kgrid.monkhorstpack'):
                if replace_kgrid:
                    output_lines.append(kgrid_line_new + '\n')
                else:
                    # Keep the original template line (for AIMD)
                    output_lines.append(line)
            
            # Replace the D3 flag
            elif line_stripped_lower.startswith('dftd3'):
                output_lines.append(d3_line_new + '\n')
            
            else:
                # Keep the original template line
                output_lines.append(line)
        
        # Write the output file
        with open(output_file, 'w') as f_output:
            f_output.writelines(output_lines)
        
        print(color_text(f"\n[OK] File '{output_file}' generated successfully.", 'green'))

        # --- 5. Copy Pseudopotentials ---
        pp_warnings = []
        if pp_path:
            pp_warnings = copy_pseudopotentials(species, pp_path)
        
        if pp_warnings:
            print("\n--- Pseudopotential Log ---")
            for warning in pp_warnings:
                print(warning)

    except Exception as e:
        # Catch any errors (e.g., File not found, Zero cell volume)
        print(f"\n--- ERROR ---")
        print(f"An error occurred: {e}")
        print("Operation aborted.")

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

    if args.intro == True:
        show_intro()

    print("\n" + color_text("Create a input file from structure in fdf format (converted by STB):", 'bold'))
    print("-"*60)

    # Process the arguments provided on the command line
    args = parser.parse_args()

    # Call the main logic function with the processed arguments
    print("--- Siesta Calculation Preparer ---")
    generate_calculation(args.structure_file, args.calc_type, args.pp_path)

# --- Script Entry Point ---
if __name__ == "__main__":
    main()
