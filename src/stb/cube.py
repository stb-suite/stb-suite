#!/usr/bin/env python

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
from stb.cli import COLORS, color_text, show_intro, print_info, print_ok, print_warn, print_error

import os
import sys
import argparse
import numpy as np

# Try to import sisl
try:
    import sisl
except ImportError:
    print()
    print_error("sisl library not found.")
    print("Please install it using: pip install sisl")
    sys.exit(1)

def get_cell_matrix(lattice_obj):
    """
    Robustly retrieves the unit cell matrix from sisl object.
    Supports: .cell, .matrix, .vectors, .get_cell()
    """
    attributes = ['cell', 'matrix', 'vectors']
    
    for attr in attributes:
        if hasattr(lattice_obj, attr):
            return getattr(lattice_obj, attr)
            
    if hasattr(lattice_obj, 'get_cell') and callable(lattice_obj.get_cell):
        return lattice_obj.get_cell()
        
    try:
        return np.array(lattice_obj)
    except:
        return None

def main():
    parser = argparse.ArgumentParser(
        description="Convert SIESTA Grid files (VT, VH, RHO, BADER) to Gaussian Cube format.",
        epilog="Example usage:\n"
               "  stb_convert --label graphene --type RHO\n"
               "  stb_convert -l system -t VT -o potential.cube",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-l", "--label", help="System Label (required).", required=True)
    parser.add_argument("-t", "--type", help="File type to convert.", choices=['VT', 'VH', 'RHO', 'BADER'], required=True)
    parser.add_argument("-o", "--output", help="Custom output filename. Default: label_type.cube", default=None)
    
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    args = parser.parse_args()

    # Filenames
    input_file = f"{args.label}.{args.type}"
    xv_file = f"{args.label}.XV"
    output_file = args.output if args.output else f"{args.label}_{args.type}.cube"

    if args.intro:
        show_intro()

    print("\n" + color_text("GRID CONVERTER:", 'bold'))
    print("-" * 60)

    # --- 1. Check Input ---
    if not os.path.exists(input_file):
        print_error(f"Input file '{input_file}' not found.")
        sys.exit(1)

    # --- 2. Load Grid ---
    print()
    print_info(f"Reading SIESTA file: {color_text(input_file, 'cyan')}...")
    try:
        sile = sisl.get_sile(input_file)
        grid = sile.read_grid()
    except Exception as e:
        print(f"{COLORS['red']}[FATAL ERROR] Could not read grid: {e}{COLORS['reset']}")
        sys.exit(1)

    # --- 3. Load Geometry (Optional) ---
    if os.path.exists(xv_file):
        print_info(f"Found geometry file: {color_text(xv_file, 'cyan')}. Mapping atoms...")
        try:
            geom = sisl.get_sile(xv_file).read_geometry()
            grid.geometry = geom
        except Exception as e:
            print_warn(f"Failed to read geometry from .XV: {e}")
    else:
        print_warn(f"No .XV file found. Atoms will be missing in the .cube output.")

    # --- 4. Metadata Display ---
    print_info(f"Grid Metadata:")
    print(f"       Shape:   {grid.shape}")
    
    cell = get_cell_matrix(grid.lattice)
    if cell is not None:
        print(f"       Lattice: [{cell[0,0]:.2f}, {cell[1,1]:.2f}, {cell[2,2]:.2f}] Ang")
    
    if args.type.upper() == "RHO":
        # Approximate integral check
        total_val = grid.grid.sum() * grid.volume / (grid.shape[0] * grid.shape[1] * grid.shape[2])
        print(f"       Integral: {total_val:.4f} e (Check if close to total electrons)")

    # --- 5. Write Output ---
    print_info(f"Converting to Cube format...")
    try:
        grid.write(output_file)
        print_ok(f"File saved as: {color_text(output_file, 'green')}")
    except Exception as e:
        print_error(f"Failed to write output file: {e}")
        sys.exit(1)

    print()
    print_info("Complete job!")
    print("\n"+"-"*60)
    # Funny closing phrase
    print(color_text("Conversion complete. It's hip to be square (cube).\n\n", 'bold'))

if __name__ == "__main__":
    main()
