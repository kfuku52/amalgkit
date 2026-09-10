"""Selection invariants and real FASTQ orchestration (no network or external binaries)."""

from collections import Counter
import gzip
import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from amalgkit import getfastq, getfastq_sampling as sampling
from amalgkit.metadata_utils import Metadata
from tests.support.gsa import manifest_row
from tests.test_gsa_workflow import native_args, install_native_inputs, write_metadata


def write_fastq(path, lengths, mate=None):
    content = b"".join(
        f"@r{i}{'/' + str(mate) if mate else ''}\n{'A' * length}\n+\n{'I' * length}\n".encode()
        for i, length in enumerate(lengths, 1)
    )
    path.write_bytes(gzip.compress(content) if str(path).endswith(".gz") else content)
    return content


def ids(path):
    with gzip.open(path, "rt") as handle:
        return [line.strip().split("/")[0] for i, line in enumerate(handle) if i % 4 == 0]


def test_prefix_nonreplacement_boundaries_and_golden_selection():
    first = sampling.selected_indices(101, 1, 20, 42, "run")
    assert sorted(first) == [7, 13, 23, 25, 32, 36, 39, 44, 45, 46, 47, 49, 50, 56, 60, 63, 71, 74, 83, 96]
    second = sampling.selected_indices(101, 21, 50, 42, "run")
    assert not first & second
    assert first | second == sampling.selected_indices(101, 1, 50, 42, "run")
    assert len(first) == 20 and min(first) >= 1 and max(first) <= 101
    assert first != sampling.selected_indices(101, 1, 20, 43, "run")
    assert first != sampling.selected_indices(101, 1, 20, 42, "other")
    assert list(sampling.selected_indices(1, 1, 1, 0, "run")) == [1]
    assert list(sampling.selected_indices(7, 1, 7, 0, "run")) == list(range(1, 8))
    for start, end in [(0, 1), (2, 1), (1, 102)]:
        with pytest.raises(ValueError):
            sampling.selected_indices(101, start, end, 0, "run")


def test_uniform_inclusion_and_block_composition():
    # Fixed seeds, theoretical hypergeometric mean/variance, no wall-time assertion.
    n, k, repeats = 40, 8, 4000
    frequencies = Counter()
    block_counts = []
    for seed in range(repeats):
        chosen = sampling.selected_indices(n, 1, k, seed, "synthetic")
        frequencies.update(chosen)
        block_counts.append(sum(index <= n // 2 for index in chosen))
    expected = repeats * k / n
    sigma = np.sqrt(repeats * k / n * (1 - k / n))
    assert all(abs(count - expected) < 6 * sigma for count in frequencies.values())
    variance = k * 0.5 * 0.5 * (n - k) / (n - 1)
    assert abs(np.mean(block_counts) - k / 2) < 6 * np.sqrt(variance / repeats)
    assert abs(np.var(block_counts) - variance) < 0.12 * variance


@pytest.mark.parametrize("paired", [False, True])
def test_sampling_filters_after_selection_and_preserves_sources(tmp_path, paired):
    paths = [tmp_path / "source1.fq.gz"]
    if paired:
        paths.append(tmp_path / "source2.fq.gz")
    lengths = [2, 10, 20, 3, 50, 10, 1, 10]
    before = [write_fastq(path, lengths, i + 1 if paired else None) for i, path in enumerate(paths)]
    outputs = [tmp_path / f"out{i}.fq.gz" for i in range(len(paths))]
    counts, manifest = sampling.sample_fastqs(
        paths, outputs, total=8, start=1, end=5, seed=3, run="x", min_length=5, run_dir=tmp_path
    )
    chosen = sampling.selected_indices(8, 1, 5, 3, "x")
    kept = [f"@r{i}" for i in sorted(chosen) if lengths[i - 1] >= 5]
    assert all(ids(output) == kept for output in outputs)
    assert counts["num_dumped"] == 5
    assert counts["num_written"] == len(kept)
    assert counts["bp_dumped"] == sum(lengths[i - 1] for i in chosen) * len(paths)
    assert manifest["candidate_bp"] == sum(lengths) * len(paths)
    assert [gzip.decompress(path.read_bytes()) for path in paths] == before
    assert all(not output.is_symlink() for output in outputs)


@pytest.mark.parametrize("fault", ["count", "id", "marker", "truncated", "total"])
def test_unselected_corruption_blocks_publication(tmp_path, fault):
    paths = [tmp_path / "r1.fq", tmp_path / "r2.fq"]
    for i, path in enumerate(paths, 1):
        write_fastq(path, [10] * 5, i)
    if fault == "count":
        write_fastq(paths[1], [10] * 4, 2)
    elif fault == "id":
        paths[1].write_bytes(paths[1].read_bytes().replace(b"@r5/", b"@wrong/"))
    elif fault == "marker":
        paths[1].write_bytes(paths[1].read_bytes().replace(b"/2", b"/1"))
    elif fault == "truncated":
        paths[1].write_bytes(paths[1].read_bytes()[:-5])
    outputs = [tmp_path / "out1.gz", tmp_path / "out2.gz"]
    for output in outputs:
        output.write_bytes(b"previous")
    with pytest.raises(ValueError):
        sampling.sample_fastqs(
            paths,
            outputs,
            total=6 if fault == "total" else 5,
            start=1,
            end=1,
            seed=0,
            run="x",
            min_length=1,
            run_dir=tmp_path,
        )
    assert all(output.read_bytes() == b"previous" for output in outputs)
    assert not (tmp_path / sampling.MANIFEST).exists()


def test_second_round_detects_source_change(tmp_path):
    source, output = tmp_path / "input.fq", tmp_path / "output.gz"
    write_fastq(source, [10] * 8)
    kwargs = dict(total=8, seed=0, run="x", min_length=1, run_dir=tmp_path)
    sampling.sample_fastqs([source], [output], start=1, end=3, **kwargs)
    original = output.read_bytes()
    source.write_bytes(source.read_bytes().replace(b"AAA", b"CCC"))
    with pytest.raises(ValueError, match="input changed"):
        sampling.sample_fastqs([source], [output], start=4, end=5, **kwargs)
    assert output.read_bytes() == original


@pytest.mark.parametrize("paired", [False, True])
def test_gsa_random_end_to_end_resume_and_manifest_tamper(tmp_path, monkeypatch, paired):
    row = manifest_row(paired, groups=2)
    write_metadata(tmp_path, [row])
    install_native_inputs(monkeypatch, row)
    args = native_args(tmp_path, "--sampling_method", "random", "--max_bp", "40", "--sampling_seed", "9")
    getfastq.getfastq_main(args)
    directory = tmp_path / "getfastq/CRR0001"
    output = directory / ("CRR0001_1.amalgkit.fastq.gz" if paired else "CRR0001.amalgkit.fastq.gz")
    original = gzip.decompress(output.read_bytes())
    stamp = output.stat().st_mtime_ns
    manifest_path = directory / sampling.MANIFEST
    manifest = json.loads(manifest_path.read_text())
    assert manifest["candidate_spots"] == 8
    assert manifest["rounds"][0]["num_dumped"] == (2 if paired else 4)
    stats = pd.read_csv(directory / "getfastq_stats.tsv", sep="\t").iloc[0]
    assert stats["spot_start_1st"] == stats["spot_end_1st"] == 0
    assert stats["getfastq_sampling_start_1st"] == 1
    getfastq.getfastq_main(args)
    assert output.stat().st_mtime_ns == stamp
    manifest_path.write_text("{}")
    getfastq.getfastq_main(args)
    assert gzip.decompress(output.read_bytes()) == original
    assert json.loads(manifest_path.read_text())["algorithm"] == sampling.ALGORITHM


def test_private_random_two_rounds_and_restart(tmp_path, monkeypatch):
    source = tmp_path / "private.fq"
    write_fastq(source, [10, 10, 2, 2, 10, 10, 2, 2] * 5)
    row = dict(
        run="private1",
        scientific_name="Arabidopsis thaliana",
        lib_layout="single",
        private_file="yes",
        read1_path=str(source),
        total_spots=40,
        total_bases=240,
        spot_length=6,
        exclusion="no",
    )
    write_metadata(tmp_path, [row])
    monkeypatch.setattr(getfastq, "check_getfastq_dependency", lambda args: None)
    args = native_args(
        tmp_path,
        "--sampling_method",
        "random",
        "--sampling_private",
        "yes",
        "--max_bp",
        "80",
        "--min_read_length",
        "5",
        "--sampling_seed",
        "5",
        "--tol",
        "0",
    )
    # Verify interruption after the first round recovers the same selection as a clean run.
    second = getfastq.sequence_extraction_2nd_round
    monkeypatch.setattr(
        getfastq, "sequence_extraction_2nd_round", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("interrupted"))
    )
    with pytest.raises(RuntimeError, match="interrupted"):
        getfastq.getfastq_main(args)
    monkeypatch.setattr(getfastq, "sequence_extraction_2nd_round", second)
    getfastq.getfastq_main(args)
    directory = tmp_path / "getfastq/private1"
    manifest = json.loads((directory / sampling.MANIFEST).read_text())
    assert len(manifest["rounds"]) == 2
    first, second_round = manifest["rounds"]
    assert first["end"] + 1 == second_round["start"]
    output = directory / "private1.amalgkit.fastq.gz"
    identifiers = ids(output)
    assert len(identifiers) == len(set(identifiers))
    result = gzip.decompress(output.read_bytes())
    stamp = output.stat().st_mtime_ns
    getfastq.getfastq_main(args)
    assert output.stat().st_mtime_ns == stamp
    args.redo = True
    getfastq.getfastq_main(args)
    assert gzip.decompress(output.read_bytes()) == result


def test_budget_excludes_unrequested_private_and_is_order_independent():
    rows = [
        dict(run="private", private_file="yes", total_bases=9999),
        dict(run="b", private_file="no", total_bases=1000),
        dict(run="a", private_file="no", total_bases=1000),
    ]
    args = SimpleNamespace(max_bp="200", sampling_method="random", sampling_private=False)
    for ordered in [rows, rows[::-1]]:
        params = getfastq.initialize_global_params(args, Metadata.from_DataFrame(pd.DataFrame(ordered)))
        assert params["num_bp_per_sra"] == 100
        assert params["total_sra_bp"] == 2000
        assert params["sampling_run_ids"] == ["a", "b"]


@pytest.mark.parametrize(
    "option,value", [("sampling_seed", -1), ("sampling_private", True), ("treat_identical_paired_as_single", True)]
)
def test_invalid_option_combinations(option, value):
    args = SimpleNamespace(
        sampling_method="contiguous" if option == "sampling_private" else "random", **{option: value}
    )
    with pytest.raises(ValueError):
        sampling.validate_options(args)


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("paired", [False, True])
def test_sra_raw_adapter_selects_before_min_length(tmp_path, monkeypatch, paired, compressed):
    from pathlib import Path

    row = dict(
        run="SRR1",
        lib_layout="paired" if paired else "single",
        total_spots=6,
        total_bases=60 * (2 if paired else 1),
        spot_length=10 * (2 if paired else 1),
    )
    metadata = Metadata.from_DataFrame(pd.DataFrame([row]))
    g = dict(num_bp_per_sra=30, max_bp=30)
    getfastq.initialize_columns(metadata, g)
    stat = dict(
        sra_id="SRR1",
        total_spot=6,
        layout=row["lib_layout"],
        spot_length=row["spot_length"],
        getfastq_sra_dir=str(tmp_path),
        metadata_idx=0,
    )

    def dump(raw_stat, args, table, start, end, return_file_state):
        assert (start, end) == (1, 6) and args.min_read_length == 0
        ext = ".fastq.gz" if compressed else ".fastq"
        for suffix in ["_1", "_2"] if paired else [""]:
            write_fastq(tmp_path / ("SRR1" + suffix + ext), [10, 2, 10, 2, 10, 10], int(suffix[-1]) if paired else None)
        getfastq.set_current_intermediate_extension(raw_stat, ext)
        return table, raw_stat, getfastq.RunFileState(work_dir=str(tmp_path))

    monkeypatch.setattr(getfastq, "run_fasterq_dump", dump)
    args = native_args(tmp_path, "--sampling_method", "random", "--min_read_length", "5")
    getfastq.extract_random_spots(args, stat, metadata, g, 1, 3)
    manifest = json.loads((tmp_path / sampling.MANIFEST).read_text())
    assert metadata.df.loc[0, "num_dumped"] == 3
    assert metadata.df.loc[0, "bp_dumped"] == manifest["rounds"][0]["bp_dumped"]
    assert not list(Path(tmp_path).glob("*.fastq"))
    if paired:
        assert ids(tmp_path / "SRR1_1.fastq.gz") == ids(tmp_path / "SRR1_2.fastq.gz")


def test_sampling_provenance_survives_merge(tmp_path):
    from amalgkit.merge import merge_fastp_stats_into_metadata

    directory = tmp_path / "getfastq/run1"
    directory.mkdir(parents=True)
    pd.DataFrame(
        [
            dict(
                run="run1",
                getfastq_sampling_method="random",
                getfastq_sampling_seed=0,
                getfastq_sampling_algorithm=sampling.ALGORITHM,
                getfastq_sampling_candidate_spots=100,
                getfastq_sampling_start_1st=1,
                getfastq_sampling_end_1st=10,
            )
        ]
    ).to_csv(directory / "getfastq_stats.tsv", sep="\t", index=False)
    metadata = Metadata.from_DataFrame(pd.DataFrame([dict(run="run1")]))
    result = merge_fastp_stats_into_metadata(metadata, str(tmp_path), max_workers=1)
    assert result.df.loc[0, "getfastq_sampling_algorithm"] == sampling.ALGORITHM
    assert result.df.loc[0, "getfastq_sampling_candidate_spots"] == 100


def test_missing_final_newline_is_safe_across_rounds(tmp_path):
    source = tmp_path / "source.fq"
    write_fastq(source, [10] * 4)
    source.write_bytes(source.read_bytes().rstrip(b"\n"))
    first, second = tmp_path / "first.gz", tmp_path / "second.gz"
    kwargs = dict(total=4, seed=0, run="x", min_length=1, run_dir=tmp_path)
    sampling.sample_fastqs([source], [first], start=1, end=2, **kwargs)
    sampling.sample_fastqs([source], [second], start=3, end=4, **kwargs)
    first.write_bytes(first.read_bytes() + second.read_bytes())
    assert len(ids(first)) == 4
    assert len(set(ids(first))) == 4


@pytest.mark.parametrize("paired", [False, True])
def test_empty_second_round_keeps_first_output(tmp_path, monkeypatch, paired):
    row = dict(
        run="private1",
        scientific_name="Arabidopsis thaliana",
        lib_layout="paired" if paired else "single",
        private_file="yes",
        total_spots=8,
        total_bases=80 * (2 if paired else 1),
        spot_length=10 * (2 if paired else 1),
        exclusion="no",
    )
    # Force the first two candidates to pass and all later candidates to fail.
    chosen = sampling.selected_indices(8, 1, 4, 0, "private1")
    passing = sorted(chosen)[:2]
    lengths = [10 if i in passing else 2 for i in range(1, 9)]
    for mate in range(1, 3 if paired else 2):
        source = tmp_path / f"source{mate}.fq"
        write_fastq(source, lengths, mate if paired else None)
        row[f"read{mate}_path"] = str(source)
    write_metadata(tmp_path, [row])
    monkeypatch.setattr(getfastq, "check_getfastq_dependency", lambda args: None)
    args = native_args(
        tmp_path,
        "--sampling_method",
        "random",
        "--sampling_private",
        "yes",
        "--max_bp",
        str(40 * (2 if paired else 1)),
        "--min_read_length",
        "5",
        "--tol",
        "0",
    )
    getfastq.getfastq_main(args)
    directory = tmp_path / "getfastq/private1"
    manifest = json.loads((directory / sampling.MANIFEST).read_text())
    assert len(manifest["rounds"]) == 2
    assert manifest["rounds"][1]["num_written"] == 0
    output = directory / ("private1_1.amalgkit.fastq.gz" if paired else "private1.amalgkit.fastq.gz")
    assert ids(output) == [f"@r{i}" for i in passing]
    stats = pd.read_csv(directory / "getfastq_stats.tsv", sep="\t").iloc[0]
    assert stats["bp_amalgkit"] == 20 * (2 if paired else 1)
    assert stats["bp_still_available"] == 0
    assert stats["bp_specified_for_extraction"] == sum(lengths) * (2 if paired else 1)
    getfastq.getfastq_main(args)


def test_parallel_and_metadata_reordering_keep_selection(tmp_path, monkeypatch):
    rows = []
    for run in ["private1", "private2"]:
        source = tmp_path / (run + ".fq")
        write_fastq(source, [10, 10, 2, 2, 10, 10, 2, 2] * 5)
        rows.append(
            dict(
                run=run,
                scientific_name="Arabidopsis thaliana",
                lib_layout="single",
                private_file="yes",
                read1_path=str(source),
                total_spots=40,
                total_bases=240,
                spot_length=6,
                exclusion="no",
                sampling_seed=1234,
            )
        )
    write_metadata(tmp_path, rows)
    monkeypatch.setattr(getfastq, "check_getfastq_dependency", lambda args: None)
    args = native_args(
        tmp_path,
        "--sampling_method",
        "random",
        "--sampling_private",
        "yes",
        "--max_bp",
        "160",
        "--min_read_length",
        "5",
        "--sampling_seed",
        "5",
        "--tol",
        "0",
    )
    getfastq.getfastq_main(args)
    outputs = [tmp_path / f"getfastq/{row['run']}/{row['run']}.amalgkit.fastq.gz" for row in rows]
    expected = [gzip.decompress(path.read_bytes()) for path in outputs]
    write_metadata(tmp_path, rows[::-1])
    args.threads, args.internal_jobs, args.redo = 2, 2, True
    getfastq.getfastq_main(args)
    assert [gzip.decompress(path.read_bytes()) for path in outputs] == expected
    # Existing run-selection seed is not overwritten by read-sampling initialization.
    metadata = Metadata.from_DataFrame(pd.DataFrame(rows))
    getfastq.initialize_columns(metadata, {"num_bp_per_sra": 80})
    assert metadata.df["sampling_seed"].tolist() == [1234, 1234]


def test_private_source_in_managed_directory_rejected_before_redo(tmp_path, monkeypatch):
    directory = tmp_path / "getfastq/private1"
    directory.mkdir(parents=True)
    source = directory / "private1.fastq.gz"
    content = write_fastq(source, [10] * 4)
    write_metadata(
        tmp_path,
        [
            dict(
                run="private1",
                scientific_name="Arabidopsis thaliana",
                lib_layout="single",
                private_file="yes",
                read1_path=str(source),
                total_spots=4,
                total_bases=40,
                spot_length=10,
                exclusion="no",
            )
        ],
    )
    monkeypatch.setattr(getfastq, "check_getfastq_dependency", lambda args: None)
    args = native_args(tmp_path, "--sampling_method", "random", "--sampling_private", "yes", "--redo", "yes")
    with pytest.raises(ValueError, match="outside its managed"):
        getfastq.getfastq_main(args)
    assert gzip.decompress(source.read_bytes()) == content


def test_private_content_change_invalidates_fingerprint(tmp_path):
    source = tmp_path / "source.fq"
    write_fastq(source, [10] * 4)
    table = Metadata.from_DataFrame(pd.DataFrame([dict(run="x", private_file="yes", read1_path=str(source))]))
    args = SimpleNamespace(sampling_method="random", sampling_private=True, sampling_seed=0)
    stat, params = dict(sra_id="x", layout="single", metadata_idx=0), dict(num_bp_per_sra=20)
    before = getfastq.build_getfastq_run_fingerprint(args, stat, params, table)
    source.write_bytes(source.read_bytes().replace(b"AAA", b"CCC"))
    assert getfastq.build_getfastq_run_fingerprint(args, stat, params, table) != before
