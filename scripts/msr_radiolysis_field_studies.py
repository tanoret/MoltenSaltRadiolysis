"""Generate reactor-scale and STAR-CCM+ field-resolved radiolysis studies."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import argparse
import json
from msr_radiolysis_field_core import *
from msr_radiolysis_field_reactor import run_msr_offgas_study, run_msr_redox_study, run_msr_intermediate_study
from msr_radiolysis_field_capsules import extract_capsule_case, run_capsule_study

def run_all(sim_dir: Path, output_dir: Path) -> dict[str, Any]:
    figures = output_dir / 'figures'
    results = output_dir / 'results'
    figures.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    manifest = {'model_boundaries': {'star_state': 'frozen saved fields; no resolved irradiation-time sequence', 'stable_halogen': 'empirical G-value source redistributed spatially by repository kinetic branching', 'redox': 'relative U(IV)/U(III) Nernst shift from a net-equivalent branching envelope', 'intermediates': 'repository homogeneous kernels plus explicit pseudo-first-order scavenging closures', 'reported_gas': 'total gas composition unspecified; used only as an inventory ceiling/context'}, 'constants': {'G_Cl2_effective_molecules_100eV': G_CL2_EFFECTIVE, 'G_F2_empirical_molecules_100eV': G_F2_EMPIRICAL, 'tau_e_reference_s': TAU_E_REF_S, 'tau_oxidant_reference_s': TAU_OX_REF_S, 'beta_redox_capsule': BETA_REDOX_CAPSULE, 'KH_Cl2_mol_m3_Pa': KH_CL2, 'KH_F2_mol_m3_Pa': KH_F2}}
    offgas = run_msr_offgas_study(figures, results)
    redox = run_msr_redox_study(figures, results)
    intermediates = run_msr_intermediate_study(figures, results)
    capsule = run_capsule_study(sim_dir, figures, results)
    manifest['capsules'] = capsule['case_rows']
    (results / 'study_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    return {'manifest': manifest, 'offgas': offgas, 'redox': redox, 'intermediates': intermediates, 'capsule': capsule}

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sim-dir', type=Path, default=Path('/mnt/data'), help='Directory containing the four supplied STAR-CCM+ .sim files')
    parser.add_argument('--output-dir', type=Path, default=HERE.parent, help='Report package root; figures/ and results/ are written below it')
    args = parser.parse_args()
    missing = [case.file_name for case in CAPSULE_CASES if not (args.sim_dir / case.file_name).is_file()]
    if missing:
        parser.error('missing STAR-CCM+ files: ' + ', '.join(missing))
    run_all(args.sim_dir, args.output_dir)
    print(args.output_dir / 'results' / 'study_manifest.json')
if __name__ == '__main__':
    main()
