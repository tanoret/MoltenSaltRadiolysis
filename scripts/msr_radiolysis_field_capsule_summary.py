"""Machine-readable summaries and orchestration for capsule studies."""
from __future__ import annotations
from msr_radiolysis_field_core import *
from msr_radiolysis_field_capsule_extract import extract_capsule_case, write_capsule_cell_data
from msr_radiolysis_field_capsule_plots import *

def write_capsule_summaries(cases: list[dict[str, Any]], results: Path) -> dict[str, Any]:
    case_rows: list[dict[str, Any]] = []
    cycle_rows: list[dict[str, Any]] = []
    redox_rows: list[dict[str, Any]] = []
    for data in cases:
        case: CapsuleCase = data["case"]
        salt: SaltMesh = data["salt"]
        plenum: SaltMesh = data["plenum"]
        T = data["temperature_K"]
        density = data["density_kg_m3"]
        stored_source = salt.fields.get("UserSpecifiedEnergySource", np.full(salt.n_cells, np.nan))
        velocity_u = salt.fields.get("U_Velocity")
        velocity_v = salt.fields.get("V_Velocity")
        velocity_rms = math.nan
        if velocity_u is not None and velocity_v is not None:
            velocity_rms = math.sqrt(
                volume_average(salt, velocity_u**2 + velocity_v**2)
            )
        case_rows.append(
            {
                "case": case.key,
                "capsule": case.capsule,
                "salt": case.composition.key,
                "sim_file": case.file_name,
                "star_release_number": data["salt"].version.get("ReleaseNumber", "20.06.010"),
                "star_presentation_version": data["salt"].version.get("PresentationVersion", "2510.0001"),
                "salt_cells": salt.n_cells,
                "plenum_cells": plenum.n_cells,
                "salt_volume_cm3": salt.total_volume_m3 * 1e6,
                "plenum_volume_cm3": plenum.total_volume_m3 * 1e6,
                "salt_mass_g": data["salt_mass_kg"] * 1000.0,
                "uranium_mass_g": data["uranium_mass_kg"] * 1000.0,
                "user_power_W_cc": case.user_power_W_cc,
                "stored_power_W_cc": volume_average(salt, stored_source) / 1e6,
                "thermal_power_W": data["thermal_power_W"],
                "temperature_min_K": float(np.min(T)),
                "temperature_volume_average_K": volume_average(salt, T),
                "temperature_max_K": float(np.max(T)),
                "density_min_kg_m3": float(np.min(density)),
                "density_volume_average_kg_m3": volume_average(salt, density),
                "density_max_kg_m3": float(np.max(density)),
                "liquid_fraction_min": float(np.min(data["liquid_fraction"])),
                "liquid_fraction_volume_average": volume_average(salt, data["liquid_fraction"]),
                "velocity_rms_m_s": velocity_rms,
                "plenum_temperature_capacity_average_K": plenum.total_volume_m3
                / float(np.sum(plenum.cell_volume_m3 / plenum.fields["Temperature"])),
                "gas_capacitance_mol_Pa": data["gas_capacitance_mol_Pa"],
                "electron_volume_average_mol_m3": volume_average(salt, data["intermediate"]["electron_mol_m3"]),
                "oxidant_volume_average_mol_m3": volume_average(salt, data["intermediate"]["oxidizing_intermediate_mol_m3"]),
                "oxidant_p05_mol_m3": volume_quantile(salt, data["intermediate"]["oxidizing_intermediate_mol_m3"], 0.05),
                "oxidant_p95_mol_m3": volume_quantile(salt, data["intermediate"]["oxidizing_intermediate_mol_m3"], 0.95),
                "electron_effective_lifetime_volume_average_s": volume_average(salt, data["intermediate"]["tau_e_s"]),
                "oxidant_effective_lifetime_volume_average_s": volume_average(salt, data["intermediate"]["tau_oxidant_s"]),
                "temperature_oxidant_volume_correlation": volume_correlation(
                    salt, T, data["intermediate"]["oxidizing_intermediate_mol_m3"]
                ),
                "temperature_stable_source_volume_correlation": volume_correlation(
                    salt, T, data["empirical_stable_source_mol_m3_s"]
                ),
                "temperature_redox_cycle4_volume_correlation": volume_correlation(
                    salt, T, data["redox4"]["delta_E_V"]
                ),
                "stable_source_p95_over_p05": volume_quantile(
                    salt, data["empirical_stable_source_mol_m3_s"], 0.95
                ) / max(
                    volume_quantile(salt, data["empirical_stable_source_mol_m3_s"], 0.05),
                    1.0e-300,
                ),
                "stable_source_total_mol_s": data["stable_source_total_mol_s"],
            }
        )
        for cycle, burnup, gas, time_s, gross, pair_max, redox in (
            (3, case.burnup_3_GWd_tU, case.gas_3_mol, data["time3_s"], data["gross_stable3_mol"], data["maximum_pair3_mol"], data["redox3"]),
            (4, case.burnup_4_GWd_tU, case.gas_4_mol, data["time4_s"], data["gross_stable4_mol"], data["maximum_pair4_mol"], data["redox4"]),
        ):
            chi = data["gas_capacitance_mol_Pa"]
            cycle_rows.append(
                {
                    "case": case.key,
                    "cycle": cycle,
                    "burnup_GWd_tU": burnup,
                    "equivalent_full_power_days": time_s / SECONDS_PER_DAY,
                    "reported_total_gas_mol": gas,
                    "reported_total_gas_pressure_equivalent_Pa": gas / chi,
                    "reported_total_gas_equivalent_G_molecules_100eV": gas
                    / max(data["thermal_power_W"] * time_s * G_TO_MOL_PER_J, 1.0e-300),
                    "gross_empirical_halogen_mol": gross,
                    "empirical_stable_halogen_G_molecules_100eV": case.composition.stable_g_value,
                    "gross_empirical_halogen_pressure_Pa": gross / chi,
                    "maximum_primary_pair_halogen_mol": pair_max,
                    "reported_gas_over_empirical_halogen": gas / max(gross, 1e-300),
                    "reported_gas_over_maximum_pair_halogen": gas / max(pair_max, 1e-300),
                    "conditional_halogen_survival_upper_bound": min(1.0, gas / max(gross, 1e-300)),
                    "redox_beta_volume_average": BETA_REDOX_CAPSULE,
                    "delta_E_volume_average_mV": 1000.0 * volume_average(salt, redox["delta_E_V"]),
                    "delta_E_p05_mV": 1000.0 * volume_quantile(salt, redox["delta_E_V"], 0.05),
                    "delta_E_p95_mV": 1000.0 * volume_quantile(salt, redox["delta_E_V"], 0.95),
                    "U3_fraction_consumed_volume_average": volume_average(salt, redox["fraction_u3_consumed"]),
                }
            )
        for beta in (-1e-4, -1e-5, -1e-6, -1e-7, 0.0, 1e-7, 1e-6, 1e-5, 1e-4):
            beta_local = beta * data["stable_weight"]
            rate = beta_local * data["intermediate"]["source_pair_mol_m3_s"]
            redox = uranium_redox_state(
                data["composition_state"]["total_u_mol_m3"],
                data["composition_state"]["nominal_u3_fraction"],
                T,
                rate,
                data["time4_s"],
            )
            redox_rows.append(
                {
                    "case": case.key,
                    "cycle": 4,
                    "net_redox_beta_volume_average": beta,
                    "delta_E_volume_average_mV": 1000.0 * volume_average(salt, redox["delta_E_V"]),
                    "delta_E_p05_mV": 1000.0 * volume_quantile(salt, redox["delta_E_V"], 0.05),
                    "delta_E_p95_mV": 1000.0 * volume_quantile(salt, redox["delta_E_V"], 0.95),
                    "fraction_U3_consumed_volume_average": volume_average(salt, redox["fraction_u3_consumed"]),
                    "fraction_U4_reduced_volume_average": volume_average(salt, redox["fraction_u4_reduced"]),
                }
            )
    _write_csv(results / "capsule_star_field_summary.csv", case_rows)
    _write_csv(results / "capsule_cycle_summary.csv", cycle_rows)
    _write_csv(results / "capsule_redox_sensitivity.csv", redox_rows)
    return {"case_rows": case_rows, "cycle_rows": cycle_rows, "redox_rows": redox_rows}


def run_capsule_study(sim_dir: Path, figures: Path, results: Path) -> dict[str, Any]:
    cases = [extract_capsule_case(case, sim_dir) for case in CAPSULE_CASES]
    for data in cases:
        write_capsule_cell_data(data, results)
    plot_capsule_temperature(cases, figures)
    plot_capsule_intermediates(cases, figures)
    plot_capsule_stable_sources(cases, figures)
    plot_capsule_redox(cases, figures)
    plot_capsule_fast_evolution(cases, figures)
    plot_capsule_long_evolution(cases, figures)
    summaries = write_capsule_summaries(cases, results)
    return {"cases": cases, **summaries}

