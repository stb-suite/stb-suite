#!/usr/bin/env python3

import os
from time import sleep

try:
    from importlib.metadata import version as _pkg_version
    VERSION = _pkg_version("stb_suite")
except Exception:
    VERSION = "1.9.5"

COLORS = {
    'reset':     '\033[0m',
    'cyan':      '\033[96m',
    'blue':      '\033[94m',
    'green':     '\033[92m',
    'yellow':    '\033[93m',
    'red':       '\033[91m',
    'bold':      '\033[1m',
    'underline': '\033[4m',
    'bg_red':    '\033[41m',
    'white':     '\033[97m',
    'magenta':   '\033[95m',
}


def color_text(text: str, color: str) -> str:
    return f"{COLORS[color]}{text}{COLORS['reset']}"


def print_info(msg: str) -> None:
    print(color_text(f"[INFO]  {msg}", 'cyan'))

def print_ok(msg: str) -> None:
    print(color_text(f"[OK]    {msg}", 'green'))

def print_warn(msg: str) -> None:
    print(color_text(f"[WARN]  {msg}", 'yellow'))

def print_error(msg: str) -> None:
    print(color_text(f"[ERROR] {msg}", 'red'))


def show_intro() -> None:
    os.system('cls' if os.name == 'nt' else 'clear')

    logo = color_text(r"""
.----------------.  .----------------.  .----------------.
| .--------------. || .--------------. || .--------------. |
| |    _______   | || |  _________   | || |   ______     | |
| |   /  ___  |  | || | |  _   _  |  | || |  |_   _ \    | |
| |  |  (__ \_|  | || | |_/ | | \_|  | || |    | |_) |   | |
| |   '.___`-.   | || |     | |      | || |    |  __'.   | |
| |  |`\____) |  | || |    _| |_     | || |   _| |__) |  | |
| |  |_______.'  | || |   |_____|    | || |  |_______/   | |
| |              | || |              | || |              | |
| '--------------' || '--------------' || '--------------' |
 '----------------'  '----------------'  '----------------'
 """, 'cyan')

    description = [
        "Siesta ToolBox Suite",
        "A comprehensive toolkit for SIESTA DFT simulations",
        f"Version {VERSION} | University of Brasilia - 2025",
        "Developed by Dr. Carlos M. O. Bastos"
    ]

    print(logo)
    print("\n" + "=" * 60)
    for line in description:
        print(line.center(60))
        sleep(0.2)
    print("=" * 60 + "\n")
