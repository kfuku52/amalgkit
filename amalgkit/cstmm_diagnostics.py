"""Observed-only comparisons. These results never participate in normalization."""

import csv
import textwrap
from itertools import combinations

import numpy
import pandas

from amalgkit.normalization_tmm import calc_factor_tmm
from amalgkit.table_io import read_identifier_tsv


def write_observed_pair_diagnostics(counts, library_sizes, factors, path, global_reference=None):
    """Stream sample-pair audits without imputing or fitting a global scale.

    Ratios use the second sample relative to the first. A pair-specific TMM
    reference need not agree with CSTMM's common reference even without missing
    data. Differences are sensitivity diagnostics, not a significance test.
    """
    fields = ['reference_sample', 'sample', 'shared_observed', 'shared_positive',
              'retained_pairs', 'status', 'observed_factor_ratio', 'applied_factor_ratio',
              'log2_ratio_difference', 'purpose', 'global_reference_sample']
    values = counts.reindex(columns=factors.index).to_numpy(dtype=float)
    columns = list(factors.index)
    if global_reference is not None and global_reference not in columns:
        raise ValueError('The global TMM reference must be a retained sample.')
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for ref_idx, obs_idx in combinations(range(len(columns)), 2):
            # Compute directly against the actual global reference, including
            # when it occurs later in input order. Do not assume reciprocity at
            # floating-point trimming boundaries by inverting another estimate.
            if columns[obs_idx] == global_reference:
                ref_idx, obs_idx = obs_idx, ref_idx
            reference, sample = columns[ref_idx], columns[obs_idx]
            common = numpy.isfinite(values[:, ref_idx]) & numpy.isfinite(values[:, obs_idx])
            diagnostics = {}
            ratio = calc_factor_tmm(
                obs=values[common, obs_idx], ref=values[common, ref_idx],
                libsize_obs=float(library_sizes[sample]), libsize_ref=float(library_sizes[reference]),
                diagnostics=diagnostics,
            )
            estimable = diagnostics['status'] in {'estimated', 'near_identical_rates'}
            estimable = estimable and numpy.isfinite(ratio) and ratio > 0
            applied = float(factors[sample] / factors[reference])
            writer.writerow(dict(
                reference_sample=reference, sample=sample, shared_observed=int(common.sum()),
                shared_positive=diagnostics['positive_pairs'], retained_pairs=diagnostics['retained_pairs'],
                status=diagnostics['status'], observed_factor_ratio=ratio if estimable else '',
                applied_factor_ratio=applied,
                log2_ratio_difference=float(numpy.log2(ratio / applied)) if estimable else '',
                purpose='reference_only',
                global_reference_sample=global_reference if global_reference is not None else '',
            ))


def plot_observed_pair_diagnostics(table_path, output_path, global_reference):
    """Plot comparable factor ratios, not separately centered raw factors.

    Contract: a point is one sample pair; x=log2(f_sample/f_reference) from
    production, y=log2(observed-only pairwise factor). Both panels share axes
    and an identity line. The fixed-reference panel isolates the global anchor;
    the all-pairs panel also includes variation due to changing that anchor.
    Static PDF/PNG, one blue palette with a neutral identity line, no fitting
    or significance claim. Empty comparisons are explicitly counted.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    table = read_identifier_tsv(table_path, identifier_columns=('reference_sample', 'sample', 'global_reference_sample'))
    if not table.empty and not table['global_reference_sample'].eq(str(global_reference)).all():
        raise ValueError('The plot reference must match the global reference recorded in the pair table.')
    ratios = table[['applied_factor_ratio', 'observed_factor_ratio']].apply(pandas.to_numeric, errors='coerce').astype(float)
    valid = numpy.isfinite(ratios).all(axis=1) & ratios.gt(0).all(axis=1)
    valid &= table['status'].isin(['estimated', 'near_identical_rates'])
    fixed = table['reference_sample'].eq(str(global_reference))
    coordinates = numpy.log2(ratios.loc[valid].to_numpy(dtype=float))
    bound = max(0.5, float(numpy.abs(coordinates).max()) * 1.12) if coordinates.size else 0.5
    anchor = '\n'.join(textwrap.wrap(str(global_reference), width=90))
    anchor_lines = anchor.count('\n') + 1
    fig = Figure(figsize=(10.5, 6.0 + 0.18 * (anchor_lines - 1)), facecolor='white')
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)
    fig.subplots_adjust(left=0.10, right=0.97, bottom=0.23, top=0.78, wspace=0.32)
    fig.suptitle('TMM normalization factor comparison', fontsize=15, color='#252525', y=0.98)
    fig.text(0.5, 0.93, 'Global round-2 reference: ' + anchor, ha='center', va='top', fontsize=9)
    for ax, selected, title in zip(axes, [fixed, pandas.Series(True, index=table.index)],
                                    ['Same global reference', 'All sample pairs']):
        shown = valid & selected
        xy = numpy.log2(ratios.loc[shown].to_numpy(dtype=float))
        ax.plot([-bound, bound], [-bound, bound], '--', color='#555555', linewidth=1, label='Agreement (y = x)')
        if xy.size:
            ax.scatter(xy[:, 0], xy[:, 1], s=22, color='#3575A5', alpha=0.65,
                       edgecolors='none', rasterized=True)
        else:
            ax.text(0.5, 0.5, 'No estimable comparisons', transform=ax.transAxes, ha='center', fontsize=10)
        ax.set(xlim=(-bound, bound), ylim=(-bound, bound),
               xlabel='Global, with imputation\nlog2(f_sample / f_reference)',
               ylabel='Pairwise, observed only\nlog2(f_sample / f_reference)')
        ax.set_aspect('equal', adjustable='box')
        ax.set_title('{}\n{} plotted; {} not estimable'.format(
            title, int(shown.sum()), int((selected & ~valid).sum())), fontsize=11, pad=10)
        ax.grid(alpha=0.18, linewidth=0.5)
        ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(loc='upper left', frameon=False, fontsize=8)
    fig.text(0.10, 0.09, 'Both methods use the original library sizes. Points may overlap; pairs are not independent.', fontsize=9)
    fig.text(0.10, 0.055, 'Different gene sets / trimming can cause differences. Reference only: applied factors are unchanged.', fontsize=9)
    fig.savefig(output_path, dpi=160)
    return fig
