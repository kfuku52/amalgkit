import json

import pandas

from amalgkit import output_contracts


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
