"""Evaluate filtering choices on labelled synthetic expression profiles.

This measures detection/retention tradeoffs, not runtime or real-data accuracy.
Labels are assigned by the data-generating intervention, never by a margin.
Run from the repository with PYTHONPATH=. and a supported Python version.
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from amalgkit.cross_species_computation import calculate_correlation_within_group
from amalgkit.cross_species_filter import _apply_csfilter_outlier_flags
from amalgkit.per_species_python import _apply_within_group_filter, _should_stop_within_group_filter


def make_expression(seed, scenario, species_count=5, replicates=5, genes=200):
    """Generate log-scale expression; a technical swap and a biological shift
    intentionally have the same observed expression in the paired scenarios.
    """
    rng = np.random.default_rng(seed)
    common = rng.normal(12, 2, genes)
    tissue = rng.normal(0, 2, (2, genes))
    columns, rows = {}, []
    for species in range(species_count):
        species_effect = rng.normal(0, 0.4, genes)
        for group in range(2):
            for rep in range(replicates):
                run = f's{species}g{group}r{rep}'
                shifted = species == 0 and group == 0
                bad = False
                biological = False
                source_group = group
                if scenario in {'individual_swap', 'sparse_swap'} and shifted and rep == 0:
                    source_group, bad = 1, True
                if scenario in {'lineage_shift', 'species_artifact'} and shifted:
                    source_group = 1
                    bad = scenario == 'species_artifact'
                    biological = scenario == 'lineage_shift'
                values = common + tissue[source_group] + species_effect + rng.normal(0, 0.5, genes)
                if scenario in {'low_support', 'sparse_swap'} and shifted and rep == 0:
                    # Sparse observation without a technical corruption of the
                    # observed values. Labels distinguish support from quality.
                    values[2:] = np.nan
                columns[f'S{species}_{run}'] = values
                rows.append(dict(run=run, species_tag=f'S{species}', sample_group=f'G{group}',
                                 bioproject=f'P{species}_{rep % 3}', exclusion='no',
                                 bad=bad, biological=biological))
    return pd.DataFrame(columns), pd.DataFrame(rows)


def cs_flags(matrix, metadata, options):
    result = calculate_correlation_within_group(
        metadata, matrix, 'corrected',
        reference_exclusion=options.get('reference_exclusion', 'run'),
        min_common_genes=options.get('min_common_genes', 0),
    )
    for metric in ('within_group_cor', 'max_nongroup_cor'):
        result[metric + '_uncorrected'] = result[metric + '_corrected']
    result = _apply_csfilter_outlier_flags(
        result, outlier_method='robust_margin',
        small_group_policy=options.get('small_group_policy', 'margin_fallback'),
        robust_z_scope=options.get('robust_z_scope', 'sample_group'),
    )
    return result['exclusion'].ne('no').to_numpy()


def ws_flags(matrix, metadata, options):
    matrix = matrix.copy()
    matrix.columns = metadata['run'].tolist()
    args = SimpleNamespace(**options)
    for round_index in range(len(metadata) + 1):
        next_matrix, next_metadata, excluded = _apply_within_group_filter(
            matrix, metadata, args, ['G0', 'G1'],
        )
        stop = _should_stop_within_group_filter(
            matrix, next_matrix, excluded, round_index + 1, options.get('max_filter_iterations'),
        )
        matrix, metadata = next_matrix, next_metadata
        if stop:
            return metadata['exclusion'].ne('no').to_numpy()
    raise AssertionError('Finite monotone exclusion must terminate')


def evaluate(seeds=12):
    variants = {
        'csfilter': {
            'baseline': {}, 'retain_small': {'small_group_policy': 'retain'},
            'minimum_50': {'min_common_genes': 50},
            'exclude_species': {'reference_exclusion': 'species'},
            'within_species_z': {'robust_z_scope': 'species_group'},
        },
        'wsfilter': {
            'baseline': {}, 'retain_small': {'small_group_policy': 'retain'},
            'minimum_50': {'min_common_genes': 50},
            'exclude_project': {'reference_exclusion': 'bioproject'},
            'single_pass': {'max_filter_iterations': 1},
        },
    }
    records = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        warnings.simplefilter('ignore', UserWarning)
        for command, choices in variants.items():
            for replicates in (2, 5):
                for scenario in ('clean', 'individual_swap', 'lineage_shift', 'species_artifact', 'low_support', 'sparse_swap'):
                    totals = {name: dict(tp=0, fp=0, fn=0, tn=0, biological_removed=0, biological_total=0)
                              for name in choices}
                    for seed in range(seeds):
                        matrix, metadata = make_expression(
                            seed, scenario, species_count=5 if command == 'csfilter' else 1,
                            replicates=replicates,
                        )
                        truth = metadata['bad'].to_numpy()
                        bio = metadata['biological'].to_numpy()
                        for name, options in choices.items():
                            flags = (cs_flags if command == 'csfilter' else ws_flags)(matrix, metadata, options)
                            result = totals[name]
                            for key, count in dict(tp=(flags & truth).sum(), fp=(flags & ~truth).sum(),
                                                   fn=(~flags & truth).sum(), tn=(~flags & ~truth).sum(),
                                                   biological_removed=(flags & bio).sum(),
                                                   biological_total=bio.sum()).items():
                                result[key] += int(count)
                    for name, counts in totals.items():
                        records.append(dict(command=command, replicates=replicates, scenario=scenario,
                                            variant=name, **counts))
    return dict(seeds=list(range(seeds)), genes=200, scale='synthetic log expression',
                limitations='Synthetic interventions only; no independent real quality labels. '
                            'Paired lineage_shift and species_artifact have identical inputs and different truths.',
                variants=variants, results=records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seeds', type=int, default=12)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.seeds < 1:
        parser.error('--seeds must be positive')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(evaluate(args.seeds), indent=2) + '\n')
