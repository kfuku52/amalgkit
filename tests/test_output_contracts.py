import json
import os

import pandas
import pytest

from amalgkit import output_contracts
from amalgkit.identifier_validation import TargetIdTracker
from amalgkit.table_io import read_identifier_tsv


def test_quant_validator_streams_rows_and_detects_cross_chunk_duplicate(tmp_path, monkeypatch):
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
    real_read_csv = output_contracts.pandas.read_csv
    calls = []

    def recording_read_csv(*args, **kwargs):
        calls.append(kwargs.copy())
        return real_read_csv(*args, **kwargs)

    monkeypatch.setattr(output_contracts.pandas, "read_csv", recording_read_csv)

    valid, error = output_contracts.validate_quant_output_files(run_id, str(tmp_path))

    assert not valid
    assert "duplicate target_id" in error
    assert any(call.get("nrows") == 5 for call in calls)
    assert any(call.get("chunksize") == 10_000 for call in calls)
    assert all("nrows" in call or "chunksize" in call for call in calls)


def test_quant_run_info_rejects_nonstandard_nonfinite_number(tmp_path):
    run_info = tmp_path / "run_info.json"
    run_info.write_text('{"p_pseudoaligned": NaN}', encoding="utf-8")

    assert "out-of-range" in output_contracts.validate_quant_run_info_json(str(run_info))


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
