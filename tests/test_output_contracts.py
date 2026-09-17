import json
import os

import pandas
import pytest

from amalgkit import output_contracts
from amalgkit.identifier_validation import TargetIdTracker
from amalgkit.table_io import read_identifier_tsv


def test_quant_validator_detects_cross_chunk_duplicate(tmp_path):
    run_id = "SRR001"
    target_ids = ["tx{}".format(index) for index in range(10_001)]
    target_ids[-1] = target_ids[0]
    pandas.DataFrame(
        {
            "target_id": target_ids,
            "length": [100] * len(target_ids),
            "eff_length": [90] * len(target_ids),
            "est_counts": [1] * len(target_ids),
            "tpm": [1] * len(target_ids),
        }
    ).to_csv(tmp_path / (run_id + "_abundance.tsv"), sep="\t", index=False)
    (tmp_path / (run_id + "_run_info.json")).write_text(
        json.dumps({"p_pseudoaligned": 50}),
        encoding="utf-8",
    )
    valid, error = output_contracts.validate_quant_output_files(run_id, str(tmp_path))

    assert not valid
    assert "duplicate target_id" in error


@pytest.mark.parametrize('last_id,valid', [('0001', False), ('001', True), ('NA', True)])
def test_duplicate_validation_remains_exact_after_spilling_to_disk(tmp_path, monkeypatch, last_id, valid):
    path = tmp_path / 'counts.tsv'
    path.write_text('target_id\tcount\n0001\t1\n2\t2\n3\t3\n' + last_id + '\t4\n', encoding='utf-8')
    trackers = []

    def small_tracker():
        tracker = TargetIdTracker(memory_limit=2)
        trackers.append(tracker)
        return tracker

    monkeypatch.setattr(output_contracts, 'TargetIdTracker', small_tracker)
    error = output_contracts.validate_nonempty_table(
        str(path), ['target_id', 'count'], 'counts', numeric_nonnegative_columns=['count'], chunk_size=1,
    )
    assert (error == '') is valid
    assert trackers[0]._scratch is not None
    assert not os.path.exists(trackers[0]._scratch.name)
    if not valid:
        assert 'duplicate target_id values: 0001' in error


def test_identifier_reader_preserves_na_ids_without_changing_numeric_na(tmp_path):
    path = tmp_path / 'counts.tsv'
    path.write_text('target_id\tR1\n0001\t1\nNA\tNA\n', encoding='utf-8')
    frame = read_identifier_tsv(path, index_col=0)
    assert frame.index.tolist() == ['0001', 'NA']
    assert pandas.isna(frame.loc['NA', 'R1'])


@pytest.mark.parametrize(
    'payload,error',
    [
        ('[]', 'must contain an object'),
        ('{', 'Failed to read'),
        ('{}', 'missing'),
        ('{"p_pseudoaligned": null}', 'invalid'),
        ('{"p_pseudoaligned": -1}', 'out-of-range'),
        ('{"p_pseudoaligned": NaN}', 'out-of-range'),
    ],
)
def test_run_info_reports_invalid_content(tmp_path, payload, error):
    path = tmp_path / 'run_info.json'
    path.write_text(payload)
    assert error in output_contracts.validate_quant_run_info_json(str(path))


def test_run_info_reports_unreadable_file(tmp_path, monkeypatch):
    def denied(*args, **kwargs):
        raise PermissionError('read denied')
    monkeypatch.setattr('builtins.open', denied)
    assert 'Failed to read quant run info JSON: read denied' == output_contracts.validate_quant_run_info_json(str(tmp_path / 'info.json'))
