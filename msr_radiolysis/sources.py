
from typing import Dict
from .utils import g_to_source

def build_sources(g_values: Dict[str, float], dose_rate_J_m3_s: float) -> Dict[str,float]:
    """Return dict of species -> S_i [mol/m^3/s]"""
    sources = {}
    for s, G in g_values.items():
        sources[s] = g_to_source(G, dose_rate_J_m3_s)
    return sources
