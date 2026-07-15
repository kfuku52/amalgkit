import os
import sys

import pandas


def read_getfastq_stats_row(output_dir, sra_id, warning_writer=None):
    stats_path = os.path.join(output_dir, 'getfastq_stats.tsv')
    if not os.path.isfile(stats_path):
        return None
    try:
        stats_df = pandas.read_csv(stats_path, sep='\t')
    except Exception as exc:
        if warning_writer is None:
            warning_writer = sys.stderr.write
        warning_writer('Failed to read getfastq stats file {}: {}\n'.format(stats_path, exc))
        return None
    if 'run' not in stats_df.columns:
        return None
    matched = stats_df.loc[stats_df['run'].fillna('').astype(str).str.strip() == str(sra_id)]
    if matched.empty:
        return None
    return matched.iloc[-1]
