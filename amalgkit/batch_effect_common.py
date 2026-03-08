from dataclasses import dataclass, field
import os
from typing import Any, Dict, List, Optional

import numpy
import pandas


def normalize_run_ids(values):
    normalized = []
    seen = set()
    for value in values:
        if value is None:
            continue
        run_id = str(value).strip()
        if run_id == '':
            continue
        if run_id in seen:
            continue
        seen.add(run_id)
        normalized.append(run_id)
    return normalized


@dataclass
class BatchEffectResult:
    backend: str
    method: str
    skip_reason: str = ''
    stable: Optional[bool] = None
    corrected_run_ids: List[str] = field(default_factory=list)
    uncorrected_run_ids: List[str] = field(default_factory=list)
    resolved_sva_nsv: Optional[int] = None
    resolved_sva_B: Optional[int] = None
    resolved_ruv_k: Optional[int] = None
    resolved_ruv_controls: Optional[int] = None
    negative_values_before_clip: Optional[int] = None
    negative_values_after_clip: Optional[int] = None
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_jsonable(self):
        payload = {
            'backend': self.backend,
            'method': self.method,
            'skip_reason': self.skip_reason,
            'stable': self.stable,
            'corrected_run_ids': normalize_run_ids(self.corrected_run_ids),
            'uncorrected_run_ids': normalize_run_ids(self.uncorrected_run_ids),
            'resolved_sva_nsv': self.resolved_sva_nsv,
            'resolved_sva_B': self.resolved_sva_B,
            'resolved_ruv_k': self.resolved_ruv_k,
            'resolved_ruv_controls': self.resolved_ruv_controls,
            'negative_values_before_clip': self.negative_values_before_clip,
            'negative_values_after_clip': self.negative_values_after_clip,
        }
        payload.update(self.extra)
        return payload


def initialize_batch_info(run_ids=(), batch_effect_alg='no'):
    normalized_runs = normalize_run_ids(run_ids)
    return {
        'batch_effect_alg_requested': str(batch_effect_alg),
        'batch_effect_alg_applied': str(batch_effect_alg),
        'corrected_runs': [],
        'uncorrected_runs': normalized_runs,
        'resolved_sva_nsv': None,
        'resolved_sva_B': None,
        'sva_estimation_method': None,
        'sva_stable': None,
        'resolved_ruv_k': None,
        'resolved_ruv_controls': None,
        'ruv_baseline_score': None,
        'ruv_selected_score': None,
        'ruv_selected_penalized_score': None,
        'ruv_penalty': None,
        'skip_reason': 'not_run',
    }


def _batch_info_to_mapping(batch_info):
    if isinstance(batch_info, BatchEffectResult):
        return batch_info.to_jsonable()
    return dict(batch_info)


def _batch_info_value(batch_info, *keys, default=None):
    mapping = _batch_info_to_mapping(batch_info)
    for key in keys:
        if key in mapping:
            value = mapping[key]
            if value is not None:
                return value
    return default


def _summary_scalar(value):
    return numpy.nan if value is None else value


def annotate_metadata_with_batch_info(metadata_df, batch_info, run_column='run'):
    if run_column not in metadata_df.columns:
        raise ValueError('Missing required metadata column: {}'.format(run_column))
    annotated = metadata_df.copy()
    if annotated.shape[0] == 0:
        annotated.loc[:, 'batch_corrected'] = pandas.Series(dtype=object)
        annotated.loc[:, 'batch_alg_used'] = pandas.Series(dtype=object)
        return annotated
    corrected_runs = normalize_run_ids(
        _batch_info_value(batch_info, 'corrected_runs', 'corrected_run_ids', default=[])
    )
    used_alg = str(
        _batch_info_value(
            batch_info,
            'batch_effect_alg_applied',
            'method',
            default='no',
        )
    ).strip()
    if used_alg == '':
        used_alg = 'no'
    run_values = annotated.loc[:, run_column].astype(str)
    is_corrected = run_values.isin(corrected_runs)
    annotated.loc[:, 'batch_corrected'] = numpy.where(is_corrected, 'yes', 'no')
    annotated.loc[:, 'batch_alg_used'] = numpy.where(is_corrected, used_alg, 'no')
    return annotated


def build_batch_effect_summary_dataframe(
    batch_info,
    scientific_name,
    random_seed_value=None,
):
    corrected_runs = normalize_run_ids(
        _batch_info_value(batch_info, 'corrected_runs', 'corrected_run_ids', default=[])
    )
    uncorrected_runs = normalize_run_ids(
        _batch_info_value(batch_info, 'uncorrected_runs', 'uncorrected_run_ids', default=[])
    )
    requested_alg = str(
        _batch_info_value(
            batch_info,
            'batch_effect_alg_requested',
            'backend',
            'method',
            default='no',
        )
    )
    applied_alg = str(
        _batch_info_value(
            batch_info,
            'batch_effect_alg_applied',
            'method',
            default=requested_alg,
        )
    )
    if random_seed_value is None:
        random_seed = 'auto'
    else:
        try:
            random_seed = 'auto' if bool(pandas.isna(random_seed_value)) else str(random_seed_value)
        except TypeError:
            random_seed = str(random_seed_value)
    return pandas.DataFrame(
        [
            {
                'scientific_name': str(scientific_name),
                'batch_effect_alg_requested': requested_alg,
                'batch_effect_alg_applied': applied_alg,
                'random_seed': random_seed,
                'resolved_sva_nsv': _summary_scalar(_batch_info_value(batch_info, 'resolved_sva_nsv')),
                'resolved_sva_B': _summary_scalar(_batch_info_value(batch_info, 'resolved_sva_B')),
                'sva_estimation_method': _summary_scalar(_batch_info_value(batch_info, 'sva_estimation_method')),
                'sva_stable': _summary_scalar(_batch_info_value(batch_info, 'sva_stable', 'stable')),
                'resolved_ruv_k': _summary_scalar(_batch_info_value(batch_info, 'resolved_ruv_k')),
                'resolved_ruv_controls': _summary_scalar(_batch_info_value(batch_info, 'resolved_ruv_controls')),
                'ruv_baseline_score': _summary_scalar(_batch_info_value(batch_info, 'ruv_baseline_score')),
                'ruv_selected_score': _summary_scalar(_batch_info_value(batch_info, 'ruv_selected_score')),
                'ruv_selected_penalized_score': _summary_scalar(_batch_info_value(batch_info, 'ruv_selected_penalized_score')),
                'ruv_penalty': _summary_scalar(_batch_info_value(batch_info, 'ruv_penalty')),
                'corrected_run_count': int(len(corrected_runs)),
                'corrected_runs': '|'.join(corrected_runs),
                'uncorrected_run_count': int(len(uncorrected_runs)),
                'uncorrected_runs': '|'.join(uncorrected_runs),
                'skip_reason': _summary_scalar(_batch_info_value(batch_info, 'skip_reason', default='')),
            }
        ]
    )


def write_batch_effect_summary_tsv(
    batch_info,
    scientific_name,
    species_tag,
    dir_tsv,
    random_seed_value=None,
):
    os.makedirs(dir_tsv, exist_ok=True)
    summary_df = build_batch_effect_summary_dataframe(
        batch_info=batch_info,
        scientific_name=scientific_name,
        random_seed_value=random_seed_value,
    )
    requested_alg = str(summary_df.loc[0, 'batch_effect_alg_requested'])
    out_path = os.path.join(
        dir_tsv,
        '{}.{}.batch_effect_summary.tsv'.format(species_tag, requested_alg),
    )
    summary_df.to_csv(out_path, sep='\t', index=False)
    return {
        'summary_path': out_path,
        'summary_df': summary_df,
    }


__all__ = [
    'BatchEffectResult',
    'annotate_metadata_with_batch_info',
    'build_batch_effect_summary_dataframe',
    'initialize_batch_info',
    'normalize_run_ids',
    'write_batch_effect_summary_tsv',
]
