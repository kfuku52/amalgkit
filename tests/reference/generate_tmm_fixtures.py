"""Independent 60-digit Decimal/Fraction oracle; never imports AMALGKIT.

Regenerate with: python tests/reference/generate_tmm_fixtures.py
Use --all --output PATH for the extended 40-matrix numerical audit.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

import numpy as np


def average_ranks(values):
    positions = {}
    for position, index in enumerate(sorted(range(len(values)), key=values.__getitem__), 1):
        positions.setdefault(values[index], []).append(position)
    return [sum(positions[value]) / len(positions[value]) for value in values]


def pair_factor(obs, ref, obs_lib, ref_lib):
    pairs = [(Fraction(float(x)), Fraction(float(y))) for x, y in zip(obs, ref) if x > 0 and y > 0]
    if not pairs:
        return Decimal(1)
    m_ranks = average_ranks([x / y for x, y in pairs])
    a_ranks = average_ranks([x * y for x, y in pairs])
    n = len(pairs)
    m_low, a_low = int(n * .3) + 1, int(n * .05) + 1
    o_lib, r_lib = Decimal.from_float(float(obs_lib)), Decimal.from_float(float(ref_lib))
    ln2 = Decimal(2).ln()
    numerator, denominator = Decimal(0), Decimal(0)
    max_m = Decimal(0)
    for (x, y), m_rank, a_rank in zip(pairs, m_ranks, a_ranks):
        o = Decimal(x.numerator) / Decimal(x.denominator)
        r = Decimal(y.numerator) / Decimal(y.denominator)
        m = (o / o_lib / (r / r_lib)).ln() / ln2
        max_m = max(max_m, abs(m))
        if not (m_low <= m_rank <= n + 1 - m_low and a_low <= a_rank <= n + 1 - a_low):
            continue
        variance = (o_lib - o) / o_lib / o + (r_lib - r) / r_lib / r
        if variance == 0:
            continue
        numerator += m / variance
        denominator += 1 / variance
    if max_m < Decimal('1e-6') or denominator == 0:
        return Decimal(1)
    return (ln2 * numerator / denominator).exp()


def factors_for_reference(counts, libs, reference):
    raw = [pair_factor(counts[:, i], counts[:, reference], libs[i], libs[reference]) for i in range(len(libs))]
    center = (sum(value.ln() for value in raw) / len(raw)).exp()
    return [value / center for value in raw]


def fixture(index):
    rng = np.random.default_rng(8000 + index)
    genes = [12, 50, 250, 1000][index % 4]
    samples = [2, 3, 4, 8][(index // 4) % 4]
    values = rng.negative_binomial(3, .15, size=(genes, samples)).astype(float)
    if index % 3 == 0:
        values[rng.random(values.shape) < .65] = 0
    if index % 3 == 1:
        values *= rng.uniform(.1, 3., size=values.shape)
    if index % 3 == 2:
        values[0, :] *= 1000
    values[-1, :] = np.maximum(values[-1, :], 1)
    libs = values.sum(axis=0) * rng.uniform(1., 5., samples)
    nonzero = values[np.any(values > 0, axis=1)]
    quartiles = np.quantile(nonzero, .75, axis=0, method='linear') / libs
    if np.median(quartiles) < 1e-20:
        reference = int(np.argmax(np.sqrt(nonzero).sum(axis=0)))
    else:
        exact = [Fraction(float(q)) for q in quartiles]
        mean = sum(exact) / len(exact)
        reference = min(range(samples), key=lambda i: abs(exact[i] - mean))
    with localcontext() as context:
        context.prec = 60
        first = factors_for_reference(values, libs, reference)
        central = sorted(first)[(samples - 1) // 2 : samples // 2 + 1]
        median_reference = next(i for i, factor in enumerate(first) if factor in central)
        second = factors_for_reference(values, libs, median_reference)
        fixed = factors_for_reference(values, libs, 0)
    return {
        'case': index,
        'counts': values.tolist(),
        'library_sizes': libs.tolist(),
        'round1_reference': reference,
        'round2_reference': median_reference,
        'round1': [float(value) for value in first],
        'round2': [float(value) for value in second],
        'fixed_first': [float(value) for value in fixed],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--all', action='store_true')
    parser.add_argument('--edger', action='store_true', help='Also refresh fixed-reference edgeR fixtures using Rscript')
    parser.add_argument('--output', type=Path, default=Path(__file__).parents[1] / 'fixtures' / 'tmm.json')
    args = parser.parse_args()
    cases = range(40) if args.all else [0, 1, 3, 5, 9, 14, 15, 29]
    results = []
    for index in cases:
        results.append(fixture(index))
        print(f'Independent TMM oracle: case {index} complete', flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Keep the complete raw inputs in a compact, reproducible fixture.
    args.output.write_text(json.dumps({'precision': 60, 'cases': results}, separators=(',', ':')) + '\n', encoding='utf-8')
    if args.edger:
        # Preselected non-boundary examples, not a runtime filter on agreement.
        normal = [case for case in results if case['case'] in {0, 1, 3, 5, 9}]
        with tempfile.TemporaryDirectory(prefix='amalgkit-edger-reference-') as scratch:
            root = Path(scratch)
            for case in normal:
                with (root / f"counts-{case['case']}.tsv").open('w', newline='', encoding='utf-8') as handle:
                    csv.writer(handle, delimiter='\t').writerows(case['counts'])
                (root / f"libs-{case['case']}.txt").write_text(
                    '\n'.join(map(str, case['library_sizes'])), encoding='utf-8',
                )
            subprocess.run(  # noqa: S603 - checked-in reference script and temporary input paths
                ['Rscript', '--vanilla', str(Path(__file__).with_name('tmm_edger.R')), scratch,
                 *[str(case['case']) for case in normal]], check=True, timeout=60,
            )
            references = {
                'edgeR_version': (root / 'version.txt').read_text().strip(),
                'fixed_first': {
                    str(case['case']): np.loadtxt(root / f"factors-{case['case']}.txt").tolist()
                    for case in normal
                },
            }
            args.output.with_name('tmm-edger.json').write_text(json.dumps(references, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
