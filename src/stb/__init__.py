try:
    from importlib.metadata import version
    __version__ = version("stb_suite")
except Exception:
    __version__ = "1.9.5"

from stb.cli import COLORS, color_text, show_intro


def __getattr__(name):
    """Lazy-load heavy submodules on first access."""
    _lazy = {
        "read_bands":       ("stb.bands",          "read_data"),
        "cbm_vbm":          ("stb.bands",          "cbm_vbm"),
        "process_pdos_xml": ("stb.dos",            "process_pdos_xml"),
        "compute_monkhorts":("stb.kgrid",          "compute_monkhorts"),
        "analyze_mechanics":("stb.strain_analysis","analyze_mechanics"),
    }
    if name in _lazy:
        module_path, attr = _lazy[name]
        import importlib
        mod = importlib.import_module(module_path)
        return getattr(mod, attr)
    raise AttributeError(f"module 'stb' has no attribute {name!r}")


__all__ = [
    "__version__",
    "COLORS",
    "color_text",
    "show_intro",
    "read_bands",
    "cbm_vbm",
    "process_pdos_xml",
    "compute_monkhorts",
    "analyze_mechanics",
]
