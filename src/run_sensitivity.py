"""Reproduce OLS diagnostics and matched-estimator candidate comparisons."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys
import traceback
import pandas as pd

from lb_sensitivity import (first_file, read_primary_inputs, regenerate_ols,
                            reproduce_power, candidate_score_effects, sha256)


def prepare_support(project: Path, out: Path, stage: str, reference: Path | None) -> Path:
    """Collect local inputs without redistributing them in the source distribution.

    An explicit reference directory supplies independent expected tables. With
    no such directory, metadata are read from the input project and no expected
    OLS result is assumed.
    """
    if reference is not None:
        if not reference.is_dir():
            raise FileNotFoundError('Reference directory not found: ' + str(reference))
        return reference
    support = out / 'supporting_inputs'
    support.mkdir()
    inter = project / '04_intermediate'
    tables = project / '06_tables'
    sources = {}
    if stage in ('all','ols','effects'):
        sources['Population_Definition.csv'] = first_file(inter, ('population_definition.csv',))
        sources['Accession_Phenotypes.csv'] = first_file(inter, (
            'phenotype_all_accession_traits.csv','phenotype_all_accession_traits_v2.csv',
            'phenotype_all_accession_traits_v02.csv'))
    if stage in ('all','ols'):
        sources['Structure_PCs.csv'] = first_file(inter, ('structure_cache/PRIMARY_N100/pcs.csv',))
        sources['Structure_Markers.csv'] = first_file(inter, ('structure_cache/common_ld_pruned_structure_markers.csv',))
    if stage in ('all','effects'):
        sources['Top_Ranked_Candidates.csv'] = first_file(tables, ('Top_Ranked_Candidates.csv',))
        sources['LOO_Candidates.csv'] = first_file(tables, (
            'Leave_One_Location_Out_Candidates.csv', 'Final_Leave_One_Location_Out_Candidates.csv'))
    manifest = {}
    for name, source in sources.items():
        shutil.copyfile(source, support/name)
        manifest[name] = {'relative_source':str(source.relative_to(project)), 'sha256':sha256(source)}
    if stage in ('all','power'):
        source = first_file(tables, ('Power.csv',))
        power = pd.read_csv(source)
        old = 'minimum_detectable_standardized_per_allele_beta'
        new = 'minimum_detectable_beta_residual_SD_per_allele'
        if new not in power and old in power:
            power = power.rename(columns={old:new})
        if new not in power:
            raise ValueError('Power input lacks the detectable-effect column')
        power = power.drop(columns=[c for c in power if 'variance_explained' in c])
        power['method_note'] = ('Independent observations; residual SD per allele; additive genotype '
                                'variance 2f(1-f); no genomic covariance or phenotype-estimation uncertainty.')
        power.to_csv(support/'Power.csv', index=False)
        manifest['Power.csv'] = {'relative_source':str(source.relative_to(project)), 'sha256':sha256(source),
                                 'unit':'residual SD per allele'}
    (out/'Supporting_Input_Sources.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    return support


def run(project: Path, output: Path, stage: str = 'all', reference: Path | None = None,
        gemma_executable: str | None = None) -> dict:
    project, output = project.resolve(), output.resolve()
    if not project.is_dir():
        raise FileNotFoundError('Input project does not exist')
    if output == project or project in output.parents or output in project.parents:
        raise ValueError('Use a separate output directory outside the input project')
    if output.exists():
        raise FileExistsError('Use a new output directory: ' + str(output))
    output.mkdir(parents=True)
    summary = {'software_version':'2.0.2', 'stage':stage,
               'independent_reference_directory_supplied':reference is not None,
               'input_project_modified':False}
    try:
        supporting = prepare_support(project,output,stage,reference)
        if stage in ('all','power'):
            print('BEGIN power', flush=True)
            summary['power'] = reproduce_power(supporting,output)
            print('END power',flush=True)
        if stage in ('all','ols'):
            print('BEGIN OLS diagnostics and PC reconstruction',flush=True)
            inputs = read_primary_inputs(project,supporting,output)
            try:
                summary['OLS'] = regenerate_ols(inputs,supporting,output)
                summary['input_sha256'] = inputs['hashes']
            finally:
                inputs['geno']._mmap.close()
            print('END OLS',flush=True)
        if stage in ('all','effects'):
            print('BEGIN candidate-only GEMMA mode 3',flush=True)
            summary['candidate_effects'] = candidate_score_effects(project,supporting,output,gemma_executable)
            print('END candidate effects',flush=True)
        outcomes = [v for v in summary.values() if isinstance(v,dict) and 'passed' in v]
        summary['requested_calculations_passed'] = bool(outcomes and all(v['passed'] for v in outcomes))
        return summary
    except Exception as exc:
        summary['requested_calculations_passed'] = False
        summary['error_type'] = type(exc).__name__
        summary['error'] = str(exc)
        raise
    finally:
        (output/'Calculation_Summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--reference-dir',type=Path)
    parser.add_argument('--gemma-executable')
    parser.add_argument('--stage',choices=['all','ols','effects','power'],default='all')
    args = parser.parse_args()
    summary = run(args.project_root,args.output_dir,args.stage,args.reference_dir,args.gemma_executable)
    print(json.dumps(summary,indent=2),flush=True)
    return 0 if summary['requested_calculations_passed'] else 2


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception:
        traceback.print_exc(file=sys.stdout)
        raise SystemExit(1)
