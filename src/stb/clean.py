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
from stb.cli import color_text, show_intro, print_info

import os
import argparse
from typing import List
import numpy as np
import re

def should_delete(file, allowed_exts):
    return os.path.isfile(file) and os.path.splitext(file)[1] not in allowed_exts

def main():
    parser = argparse.ArgumentParser(description="Remove all files except those with specified extensions.")
    
    parser.add_argument(
        '--keep', nargs='+', default=['.psml', '.psf', '.fdf', '.sh'],
        help="List of extensions to keep (e.g. .fdf .sh)"
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help="Only print the files that would be removed"
    )
    parser.add_argument(
        '--no-confirm', '-n', action='store_true',
        help="Do not ask for confirmation before deleting files"
    )
    parser.add_argument(
        '--path', default='.',
        help="Directory to clean (default: current directory)"
    )
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    allowed_exts = set(args.keep)
    
    if args.intro == True:
        show_intro()
    print("\n" + color_text("Clean:", 'bold'))
    print("-"*60)
    print("Remove all files except those with specified extensions.\n")

    for file in os.listdir(args.path):
        full_path = os.path.join(args.path, file)
        if should_delete(full_path, allowed_exts):
            if args.dry_run:
                print(f"[Dry-run] Would remove: {file}")
            elif args.no_confirm:
                os.remove(full_path)
                print_info(f"Removed: {file}")
            else:
                answer = input(f"[CONFIRM] Delete {file}? [y/N] ").strip().lower()
                if answer == 'y':
                    os.remove(full_path)
                    print_info(f"Removed: {file}")
                else:
                    print_info(f"Skipped: {file}")
    
    print()
    print_info("Complete job!")
    print("\n"+"-"*60)
    print(color_text("Cleaned up! May the deleted files rest in pieces.\n\n", 'bold'))

if __name__ == "__main__":
    main()

