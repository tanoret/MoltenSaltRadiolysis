#!/usr/bin/env python3
"""Read selected finite-volume fields from native STAR-CCM+ ``.sim`` files.

The files used for the SABRE capsule calculations were saved by STAR-CCM+
20.06.010 in parallel. They contain a compact ASCII object state following a set
of binary arrays. This module implements a deliberately narrow reader for the
saved objects needed by the radiolysis post-processing:

* region and finite-volume-region metadata;
* partitioned cell scalar fields;
* partitioned vertex coordinates;
* interior and boundary face connectivity; and
* STAR's ghost-to-owner maps, used to reconstruct global cell polygons.

It is not intended as a general replacement for STAR-CCM+. The reader checks
all structural assumptions and raises a descriptive error if a file does not
match the serialization used by the supplied simulations.
"""
from __future__ import annotations

import argparse
import ast
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

_HEADER_RE = re.compile(
    br"\{'Binary': .*?'StatePosition': (\d+)L, 'Version': (\d+)\}\s*\n",
    re.S,
)
_LONG_RE = re.compile(r"(?<![A-Za-z0-9_\.])(-?\d+)L\b")


def _clean_literal(text: str) -> str:
    text = _LONG_RE.sub(r"\1", text)
    text = re.sub(r"\btrue\b", "True", text)
    text = re.sub(r"\bfalse\b", "False", text)
    text = re.sub(r"\bnull\b", "None", text)
    return text


@dataclass(frozen=True)
class PartitionMap:
    """Map STAR partition-local indices onto a compact global ordering."""

    n_element: np.ndarray
    n_interior: np.ndarray
    ghosts: tuple[tuple[tuple[int, int], ...], ...]

    @property
    def offsets(self) -> np.ndarray:
        return np.concatenate(([0], np.cumsum(self.n_interior[:-1], dtype=np.int64)))

    @property
    def global_size(self) -> int:
        return int(np.sum(self.n_interior))

    def global_index(self, partition: int, local_index: int) -> int:
        if partition < 0 or partition >= len(self.n_interior):
            raise IndexError(f"partition {partition} is outside the partition map")
        n_int = int(self.n_interior[partition])
        if local_index < 0 or local_index >= int(self.n_element[partition]):
            raise IndexError(
                f"local index {local_index} is outside partition {partition} "
                f"with {int(self.n_element[partition])} elements"
            )
        if local_index < n_int:
            return int(self.offsets[partition] + local_index)
        owner_partition, owner_local = self.ghosts[partition][local_index - n_int]
        if owner_local >= int(self.n_interior[owner_partition]):
            raise ValueError(
                "ghost owner points outside the owner's interior range: "
                f"({owner_partition}, {owner_local})"
            )
        return int(self.offsets[owner_partition] + owner_local)


@dataclass
class SaltMesh:
    """Extracted salt mesh and cell-centred STAR solution fields."""

    source_path: Path
    simulation_name: str
    version: dict[str, Any]
    region_name: str
    fields: dict[str, np.ndarray]
    vertices_xyz_m: np.ndarray
    cell_vertex_ids: tuple[np.ndarray, ...]
    cell_centroids_xr_m: np.ndarray
    cell_area_m2: np.ndarray
    cell_volume_m3: np.ndarray

    @property
    def n_cells(self) -> int:
        return len(self.cell_vertex_ids)

    @property
    def total_volume_m3(self) -> float:
        return float(np.sum(self.cell_volume_m3))

    def volume_average(self, values: np.ndarray) -> float:
        values = np.asarray(values, dtype=float)
        if values.shape != (self.n_cells,):
            raise ValueError(f"expected {self.n_cells} values, got {values.shape}")
        denominator = float(np.sum(self.cell_volume_m3))
        if denominator <= 0.0:
            raise ValueError("nonpositive mesh volume")
        return float(np.sum(values * self.cell_volume_m3) / denominator)

    def summary(self) -> dict[str, Any]:
        temperature = self.fields.get("Temperature")
        liquid = 1.0 - self.fields.get("Solidity", np.zeros(self.n_cells))
        energy = self.fields.get("UserSpecifiedEnergySource")
        density = self.fields.get("Density")
        out: dict[str, Any] = {
            "file": self.source_path.name,
            "simulation": self.simulation_name,
            "star_release_number": self.version.get("ReleaseNumber"),
            "star_presentation_version": self.version.get("PresentationVersion"),
            "star_release_date": self.version.get("ReleaseDate"),
            "star_machine_config": self.version.get("MachineConfig"),
            "region": self.region_name,
            "cells": self.n_cells,
            "vertices": len(self.vertices_xyz_m),
            "volume_m3": self.total_volume_m3,
            "x_min_m": float(np.min(self.cell_centroids_xr_m[:, 0])),
            "x_max_m": float(np.max(self.cell_centroids_xr_m[:, 0])),
            "r_min_m": float(np.min(self.cell_centroids_xr_m[:, 1])),
            "r_max_m": float(np.max(self.cell_centroids_xr_m[:, 1])),
            "liquid_fraction_min": float(np.min(liquid)),
            "liquid_fraction_max": float(np.max(liquid)),
            "liquid_fraction_volume_average": self.volume_average(liquid),
        }
        if temperature is not None:
            out.update(
                {
                    "temperature_min_K": float(np.min(temperature)),
                    "temperature_max_K": float(np.max(temperature)),
                    "temperature_volume_average_K": self.volume_average(temperature),
                }
            )
        if energy is not None:
            out.update(
                {
                    "energy_source_min_W_m3": float(np.min(energy)),
                    "energy_source_max_W_m3": float(np.max(energy)),
                    "energy_source_volume_average_W_m3": self.volume_average(energy),
                }
            )
        if density is not None:
            out.update(
                {
                    "density_min_kg_m3": float(np.min(density)),
                    "density_max_kg_m3": float(np.max(density)),
                    "density_volume_average_kg_m3": self.volume_average(density),
                    "salt_mass_kg": float(np.sum(density * self.cell_volume_m3)),
                }
            )
        return out


__all__ = [name for name in globals() if not name.startswith("__")]
