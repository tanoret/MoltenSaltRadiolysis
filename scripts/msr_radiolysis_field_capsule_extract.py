"""Couple supplied STAR-CCM+ fields to the radiolysis closures."""
from __future__ import annotations
from msr_radiolysis_field_core import *
from starccm_sim_extract import StarSimReader

def _cell_polygons(mesh: SaltMesh) -> list[np.ndarray]:
    return [mesh.vertices_xyz_m[ids, :2] * 1000.0 for ids in mesh.cell_vertex_ids]


def _field_collection(
    ax: plt.Axes,
    mesh: SaltMesh,
    values: np.ndarray,
    *,
    cmap: str,
    norm: Normalize,
) -> PolyCollection:
    collection = PolyCollection(
        _cell_polygons(mesh),
        array=np.asarray(values, dtype=float),
        cmap=cmap,
        norm=norm,
        edgecolors="none",
        linewidths=0.0,
        rasterized=True,
    )
    ax.add_collection(collection)
    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Axial coordinate [mm]")
    return collection


def extract_capsule_case(case: CapsuleCase, sim_dir: Path) -> dict[str, Any]:
    sim_path = sim_dir / case.file_name
    reader = StarSimReader(sim_path)
    requested_fields = (
        "Temperature",
        "Solidity",
        "VolumeFraction",
        "VoidFraction",
        "Density",
        "UserSpecifiedEnergySource",
        "U_Velocity",
        "V_Velocity",
    )
    salt = reader.extract_salt_mesh("Salt", requested_fields)
    plenum = reader.extract_salt_mesh("Plenum", requested_fields)
    temperature = salt.fields["Temperature"]
    density = salt.fields.get("Density", np.full(salt.n_cells, case.composition.density_user_kg_m3))
    liquid = 1.0 - salt.fields.get("Solidity", np.zeros(salt.n_cells))
    q_user = np.full(salt.n_cells, case.user_power_W_cc * 1e6) * liquid
    composition_state = case.composition.composition_state(density)
    intermediate = steady_intermediates(
        case.composition,
        temperature,
        q_user,
        composition_state["halide_mol_m3"],
    )

    # Temperature-dependent kinetic redistribution of the empirical stable
    # source.  The total G-based source is preserved exactly by volume-weight
    # normalization; only its spatial location changes.
    network_rate = np.maximum(intermediate["network_stable_halogen_mol_m3_s"], 0.0)
    mean_network = volume_average(salt, network_rate)
    if mean_network <= 0.0:
        stable_weight = np.ones(salt.n_cells)
    else:
        stable_weight = network_rate / mean_network
    empirical_stable_source = g_source(case.composition.stable_g_value, q_user) * stable_weight

    salt_mass_kg = float(np.sum(density * salt.cell_volume_m3))
    uranium_mass_density = composition_state["u_mass_kg_m3"]
    uranium_mass_kg = float(np.sum(uranium_mass_density * salt.cell_volume_m3))
    thermal_power_W = float(np.sum(q_user * salt.cell_volume_m3))
    t3 = irradiation_time_s(case.burnup_3_GWd_tU, uranium_mass_kg, thermal_power_W)
    t4 = irradiation_time_s(case.burnup_4_GWd_tU, uranium_mass_kg, thermal_power_W)
    chi = gas_capacitance_mol_Pa(salt, plenum, case.composition.henry)

    # beta_ref is redistributed with the same modeled stable-branch weight so
    # its volume average remains BETA_REDOX_CAPSULE.
    beta_local = BETA_REDOX_CAPSULE * stable_weight
    net_redox_rate = beta_local * intermediate["source_pair_mol_m3_s"]
    redox3 = uranium_redox_state(
        composition_state["total_u_mol_m3"],
        composition_state["nominal_u3_fraction"],
        temperature,
        net_redox_rate,
        t3,
    )
    redox4 = uranium_redox_state(
        composition_state["total_u_mol_m3"],
        composition_state["nominal_u3_fraction"],
        temperature,
        net_redox_rate,
        t4,
    )

    stable_source_total_mol_s = float(np.sum(empirical_stable_source * salt.cell_volume_m3))
    gross3 = stable_source_total_mol_s * t3
    gross4 = stable_source_total_mol_s * t4
    source_pair_total_mol_s = float(
        np.sum(intermediate["source_pair_mol_m3_s"] * salt.cell_volume_m3)
    )
    max_pair3 = 0.5 * source_pair_total_mol_s * t3
    max_pair4 = 0.5 * source_pair_total_mol_s * t4

    return {
        "case": case,
        "salt": salt,
        "plenum": plenum,
        "temperature_K": temperature,
        "density_kg_m3": density,
        "liquid_fraction": liquid,
        "q_user_W_m3": q_user,
        "composition_state": composition_state,
        "intermediate": intermediate,
        "stable_weight": stable_weight,
        "empirical_stable_source_mol_m3_s": empirical_stable_source,
        "beta_local": beta_local,
        "net_redox_rate_mol_m3_s": net_redox_rate,
        "redox3": redox3,
        "redox4": redox4,
        "salt_mass_kg": salt_mass_kg,
        "uranium_mass_kg": uranium_mass_kg,
        "thermal_power_W": thermal_power_W,
        "time3_s": t3,
        "time4_s": t4,
        "gas_capacitance_mol_Pa": chi,
        "stable_source_total_mol_s": stable_source_total_mol_s,
        "gross_stable3_mol": gross3,
        "gross_stable4_mol": gross4,
        "maximum_pair3_mol": max_pair3,
        "maximum_pair4_mol": max_pair4,
    }


def write_capsule_cell_data(data: dict[str, Any], results: Path) -> Path:
    case: CapsuleCase = data["case"]
    salt: SaltMesh = data["salt"]
    path = results / "capsule_cell_fields" / f"{case.key}.npz"
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        x_m=salt.cell_centroids_xr_m[:, 0],
        r_m=salt.cell_centroids_xr_m[:, 1],
        volume_m3=salt.cell_volume_m3,
        temperature_K=data["temperature_K"],
        density_kg_m3=data["density_kg_m3"],
        liquid_fraction=data["liquid_fraction"],
        power_density_W_m3=data["q_user_W_m3"],
        electron_mol_m3=data["intermediate"]["electron_mol_m3"],
        primary_radical_mol_m3=data["intermediate"]["primary_radical_mol_m3"],
        oxidizing_intermediate_mol_m3=data["intermediate"]["oxidizing_intermediate_mol_m3"],
        polyhalide_mol_m3=data["intermediate"]["polyhalide_mol_m3"],
        network_stable_halogen_mol_m3_s=data["intermediate"]["network_stable_halogen_mol_m3_s"],
        empirical_stable_halogen_mol_m3_s=data["empirical_stable_source_mol_m3_s"],
        beta_local=data["beta_local"],
        net_redox_rate_mol_m3_s=data["net_redox_rate_mol_m3_s"],
        delta_E_3cycle_V=data["redox3"]["delta_E_V"],
        delta_E_4cycle_V=data["redox4"]["delta_E_V"],
        U3_fraction_3cycle=data["redox3"]["u3_fraction_final"],
        U3_fraction_4cycle=data["redox4"]["u3_fraction_final"],
    )
    return path

