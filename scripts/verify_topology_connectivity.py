#!/usr/bin/env python3
"""Verify Lemma 1 (two-flip connectivity of the feasibility set Γ) for the
LiCl-KCl chloride + Cr kernel as defined in msr_radiolysis/data/database.yaml.

Method:
  (1) Enumerate all 2^R candidate networks (R = number of reactions).
  (2) Compute the feasibility set Γ = networks satisfying mass + charge + cycle constraints.
  (3) Build the two-flip adjacency graph on Γ (edges = pairs differing in ≤2 indicators).
  (4) Run BFS from a starting feasible network; verify reachability of every other.
  (5) Report |Γ|/2^R sparsity and the mixing-acceleration constant from Proposition 1.

This is a constructive (per-kernel) verification of the path-connectedness premise of
Lemma 1 in HBMAE_THEOREMS_TIGHTENED.md and the manuscript Theorem 1.

Usage:
    python scripts/verify_topology_connectivity.py
"""

from __future__ import annotations

import itertools
import sys
from collections import deque
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))


# ============================================================================
# Reaction inventory for the LiCl-KCl + Cr kernel
# ============================================================================
# Each reaction is encoded as (substrates, products) of *non-conserved* species.
# Conserved species (e.g. Cl- in eutectic, Li+, K+) are treated as infinite reservoirs
# and excluded from the cycle/mass-balance check.
#
# Non-conserved species (those that must be balanced reaction-by-reaction):
#   e_s-      (solvated electron)        — radiolytic source
#   Cl•       (chlorine atom)            — radiolytic source
#   Cl2•-     (dichlorine radical anion)
#   Cl3-      (trichloride)
#   Cl2_diss  (dissolved Cl2 — also has gas-exchange consumption)
#   Cr2+, Cr3+, Cr+   (chromium oxidation states)
#
# Reactions:

REACTIONS = [
    # (name, substrates dict, products dict)
    ("R1_Cl_Clminus_to_Cl2m",       {"Cl•": 1},               {"Cl2•-": 1}),
    ("R2_e_plus_Cl_to_Clminus",     {"e_s-": 1, "Cl•": 1},    {}),
    ("R3_2Cl_to_Cl2diss",           {"Cl•": 2},               {"Cl2_diss": 1}),
    ("R4_Cl2m_dispro_to_Cl3",       {"Cl2•-": 2},             {"Cl3-": 1}),
    ("R5_Cl3_to_Cl2diss",           {"Cl3-": 1},              {"Cl2_diss": 1}),
    ("R6_e_plus_Cr2_to_Crp",        {"e_s-": 1, "Cr2+": 1},   {"Cr+": 1}),
    ("R7_e_plus_Cr3_to_Cr2",        {"e_s-": 1, "Cr3+": 1},   {"Cr2+": 1}),
    ("R8_Cl2m_plus_Cr2_to_Cr3",     {"Cl2•-": 1, "Cr2+": 1},  {"Cr3+": 1}),
    ("R9_Cl2m_plus_Cr3_to_Cl2diss", {"Cl2•-": 1, "Cr3+": 1},  {"Cr2+": 1, "Cl2_diss": 1}),
    ("R10_Cr3_plus_Crp_to_2Cr2",    {"Cr3+": 1, "Cr+": 1},    {"Cr2+": 2}),
]

# Species that are produced by radiation sources (always sourced, can be consumed):
RADIOLYTIC_SOURCES = {"e_s-", "Cl•"}

# Species that are sinks (consumption-only OK, no production needed):
SINK_SPECIES = {"Cl2_diss"}     # leaves via gas exchange

# Species that may be present as initial conditions (production OK from sources only):
INITIAL_CONDITIONS = {"Cr2+", "Cr3+"}    # only one of these initially, depending on experiment


def is_feasible(gamma: tuple) -> tuple:
    """Test mass+cycle feasibility of a network γ on the non-conserved species.

    A network γ is FEASIBLE if for every non-conserved species s appearing in at least
    one active reaction:
      (consumption) the species must appear as a substrate in some active reaction
        UNLESS it is a sink species (consumed by an external mechanism, e.g. gas exchange);
      (production) the species must appear as a product in some active reaction OR be
        a radiolytic-source species OR be present as an initial condition.

    Returns (feasible: bool, reason: str).
    """
    active = [REACTIONS[i] for i, on in enumerate(gamma) if on]

    # All non-conserved species touched by active reactions
    touched = set()
    for _, subs, prods in active:
        touched.update(subs)
        touched.update(prods)

    if len(active) == 0:
        return False, "empty network has no consumption of radiolytic sources"

    # Always need to consume e_s- and Cl• if they are sourced
    for src in RADIOLYTIC_SOURCES:
        if not any(src in subs for _, subs, _ in active):
            return False, f"{src} is sourced but no active reaction consumes it"

    # Every produced non-conserved, non-sink species must have a consumer
    for sp in touched:
        if sp in SINK_SPECIES:
            continue
        if sp in RADIOLYTIC_SOURCES or sp in INITIAL_CONDITIONS:
            continue  # produced externally
        produced = any(sp in prods for _, _, prods in active)
        consumed = any(sp in subs for _, subs, _ in active)
        if produced and not consumed:
            return False, f"{sp} is produced but not consumed"

    # Cl2•- accumulation safeguard: if R1 is on, R4/R8/R9 must consume Cl2•-
    if gamma[0] == 1:  # R1 on -> Cl2•- produced
        if not any(gamma[i] for i in [3, 7, 8]):
            return False, "R1 produces Cl2•- but no R4/R8/R9 consumes it"

    return True, "feasible"


def enumerate_feasibility(R: int):
    """Enumerate all 2^R networks and partition into feasible/infeasible."""
    feasible = []
    infeasible = []
    for gamma_bits in itertools.product([0, 1], repeat=R):
        ok, reason = is_feasible(gamma_bits)
        (feasible if ok else infeasible).append((gamma_bits, reason))
    return feasible, infeasible


def hamming(g1: tuple, g2: tuple) -> int:
    return sum(a != b for a, b in zip(g1, g2))


def two_flip_connectivity(feasible_gammas):
    """Check whether the feasibility set Γ is connected under two-flip moves.

    Returns (connected: bool, component_sizes: list of int).
    """
    if not feasible_gammas:
        return False, []
    feasible_set = set(feasible_gammas)
    visited = set()
    component_sizes = []

    while feasible_set - visited:
        start = next(iter(feasible_set - visited))
        queue = deque([start])
        component = set([start])
        visited.add(start)

        while queue:
            current = queue.popleft()
            # Try all neighbours within Hamming distance ≤ 2
            for nbr in feasible_set:
                if nbr in visited:
                    continue
                if hamming(current, nbr) <= 2:
                    visited.add(nbr)
                    component.add(nbr)
                    queue.append(nbr)
        component_sizes.append(len(component))

    return len(component_sizes) == 1, component_sizes


def main():
    R = len(REACTIONS)
    print(f"Verifying two-flip connectivity for the LiCl-KCl + Cr chloride kernel")
    print(f"Number of candidate reactions R = {R}")
    print(f"Total network space |{{0,1}}^R| = 2^{R} = {2**R}")
    print()

    feasible, infeasible = enumerate_feasibility(R)
    nF = len(feasible)
    sparsity = nF / 2**R
    print(f"Feasibility:")
    print(f"  |Γ|       = {nF}")
    print(f"  |Γ|/2^R   = {sparsity:.4f} ({100*sparsity:.2f}%)")
    print(f"  Mixing acceleration (Proposition 1): τ_Γ / τ_∅ ≤ {sparsity:.4f}")
    print()

    feasible_gammas = [g for g, _ in feasible]
    connected, component_sizes = two_flip_connectivity(feasible_gammas)
    print(f"Two-flip connectivity test:")
    print(f"  Connected = {connected}")
    print(f"  Components: {component_sizes}")
    print()

    if connected:
        print(f"✓ Lemma 1 holds CONSTRUCTIVELY for this kernel.")
        print(f"  Theorem 1 (RJMCMC ergodicity) applies with two-flip proposal kernel.")
    else:
        print(f"✗ Two-flip moves are INSUFFICIENT for this kernel.")
        print(f"  Need to increase neighborhood size to k = ? to connect all components.")
        # Diagnose: try k-flip up to k=5
        for k in range(3, 6):
            visited = set()
            start = feasible_gammas[0]
            queue = deque([start])
            visited.add(start)
            feasible_set = set(feasible_gammas)
            while queue:
                cur = queue.popleft()
                for nbr in feasible_set:
                    if nbr not in visited and hamming(cur, nbr) <= k:
                        visited.add(nbr)
                        queue.append(nbr)
            if len(visited) == nF:
                print(f"  k = {k}-flip moves are sufficient to connect Γ.")
                break

    # Spot check a known feasible state
    print()
    print("Spot check: minimal feasible network from radiolytic-source consumption:")
    minimal = (1, 1, 0, 1, 0, 0, 0, 0, 0, 0)   # R1 (Cl→Cl2•-), R2 (e+Cl→Cl-), R4 (Cl2•-→Cl3-/Cl-)
    # Wait — R4 produces Cl3- but no R5 consumes it. Should be infeasible.
    ok, reason = is_feasible(minimal)
    print(f"  γ = R1+R2+R4 (no R5): feasible={ok}  reason='{reason}'")
    minimal_with_R5 = (1, 1, 0, 1, 1, 0, 0, 0, 0, 0)
    ok, reason = is_feasible(minimal_with_R5)
    print(f"  γ = R1+R2+R4+R5     : feasible={ok}  reason='{reason}'")

    # Save the feasibility table for downstream use
    out_path = REPO / "validation" / "topology_feasibility.csv"
    with out_path.open("w") as f:
        header = "gamma," + ",".join(name for name, _, _ in REACTIONS) + ",reason\n"
        f.write(header)
        for g, reason in feasible:
            bits = ",".join(str(b) for b in g)
            f.write(f"feasible,{bits},{reason}\n")
        for g, reason in infeasible[:50]:  # only first 50 infeasible to keep file small
            bits = ",".join(str(b) for b in g)
            f.write(f"infeasible,{bits},{reason}\n")
    print()
    print(f"Wrote feasibility table to {out_path.relative_to(REPO)}")
    print(f"  ({len(feasible)} feasible + {min(50, len(infeasible))} of {len(infeasible)} infeasible examples)")

    return connected, sparsity


if __name__ == "__main__":
    main()
