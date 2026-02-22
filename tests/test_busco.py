import gzip
import os
import pytest
import pandas

from types import SimpleNamespace

from amalgkit.busco import (
    normalize_busco_table,
    find_full_table,
    resolve_species_fasta,
    collect_species,
    select_tool,
    ensure_clean_dir,
    run_command,
    run_busco,
    busco_main,
)
from amalgkit.util import Metadata


def write_busco_table(path, rows):
    header = "# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n"
    with open(path, 'w') as f:
        f.write("# BUSCO version is: 6.0.0\n")
        f.write(header)
        for row in rows:
            f.write("\t".join(row) + "\n")


def test_normalize_busco_table(tmp_path):
    src = tmp_path / "full_table.tsv"
    rows = [
        ["BUSCO1", "Complete", "seq1", "100", "200", "url1", "desc1"],
        ["BUSCO2", "Missing", "-", "0", "0", "url2", "desc2"],
    ]
    write_busco_table(src, rows)
    dest = tmp_path / "normalized.tsv"
    normalize_busco_table(str(src), str(dest))
    content = dest.read_text()
    assert content.startswith("# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription")
    assert "BUSCO1" in content
    assert "BUSCO2" in content


def test_normalize_busco_table_ignores_non_header_busco_comment_rows(tmp_path):
    src = tmp_path / "full_table.tsv"
    src.write_text(
        "# BUSCO version is: 6.0.0\n"
        "# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n"
        "# BUSCO was run in mode: transcriptome\n"
        "BUSCO1\tComplete\tseq1\t100\t200\turl1\tdesc1\n"
    )
    dest = tmp_path / "normalized.tsv"

    normalize_busco_table(str(src), str(dest))

    content = dest.read_text()
    assert "BUSCO1" in content
    assert content.startswith("# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription")


def test_normalize_busco_table_without_header_preserves_first_row(tmp_path):
    src = tmp_path / "full_table.tsv"
    src.write_text(
        "BUSCO1\tComplete\tseq1\t100\t200\turl1\tdesc1\n"
        "BUSCO2\tMissing\t-\t0\t0\turl2\tdesc2\n"
    )
    dest = tmp_path / "normalized.tsv"

    normalize_busco_table(str(src), str(dest))

    content = dest.read_text()
    assert "BUSCO1" in content
    assert "BUSCO2" in content
    assert content.count("\n") == 3


def test_find_full_table_single(tmp_path):
    out_dir = tmp_path / "busco_out"
    out_dir.mkdir()
    table = out_dir / "full_table.tsv"
    table.write_text("Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
    found = find_full_table(str(out_dir))
    assert os.path.realpath(found) == os.path.realpath(str(table))


def test_find_full_table_multiple_raises(tmp_path):
    out_dir = tmp_path / "busco_out"
    out_dir.mkdir()
    (out_dir / "full_table.tsv").write_text("Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
    (out_dir / "full_table_copy.tsv").write_text("Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
    with pytest.raises(ValueError, match="Multiple BUSCO full_table files"):
        find_full_table(str(out_dir))


def test_find_full_table_nested_gz(tmp_path):
    out_dir = tmp_path / "busco_out"
    nested = out_dir / "run_busco" / "busco_output"
    nested.mkdir(parents=True)
    table = nested / "full_table.specific.busco.tsv.gz"
    with gzip.open(str(table), 'wt') as f:
        f.write("Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
    found = find_full_table(str(out_dir))
    assert os.path.realpath(found) == os.path.realpath(str(table))


def test_find_full_table_accepts_uppercase_gzip_extension(tmp_path):
    out_dir = tmp_path / "busco_out"
    out_dir.mkdir()
    table = out_dir / "full_table.TSV.GZ"
    with gzip.open(str(table), 'wt') as f:
        f.write("Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
    found = find_full_table(str(out_dir))
    assert os.path.realpath(found) == os.path.realpath(str(table))


def test_run_command_tolerates_non_utf8_output(monkeypatch):
    monkeypatch.setattr(
        'amalgkit.busco.subprocess.run',
        lambda *_args, **_kwargs: SimpleNamespace(returncode=0, stdout=b'\xff', stderr=b'\xfe'),
    )
    run_command(['dummy-cmd'])


def test_ensure_clean_dir_rejects_file_path(tmp_path):
    path_file = tmp_path / 'output_path'
    path_file.write_text('not a directory')

    with pytest.raises(NotADirectoryError, match='BUSCO output path exists but is not a directory'):
        ensure_clean_dir(str(path_file), redo=False)


def test_ensure_clean_dir_replaces_symlink_on_redo(tmp_path):
    real_output = tmp_path / 'real_output'
    real_output.mkdir()
    (real_output / 'keep.txt').write_text('keep')
    link_path = tmp_path / 'output_path'
    os.symlink(real_output, link_path)

    ensure_clean_dir(str(link_path), redo=True)

    assert os.path.isdir(str(link_path))
    assert not os.path.islink(str(link_path))
    assert os.path.exists(str(real_output / 'keep.txt'))


def test_ensure_clean_dir_replaces_broken_symlink_on_redo(tmp_path):
    missing_target = tmp_path / 'missing_target'
    link_path = tmp_path / 'output_path'
    os.symlink(missing_target, link_path)

    ensure_clean_dir(str(link_path), redo=True)

    assert os.path.isdir(str(link_path))
    assert not os.path.islink(str(link_path))


def test_ensure_clean_dir_rejects_broken_symlink_without_redo(tmp_path):
    missing_target = tmp_path / 'missing_target'
    link_path = tmp_path / 'output_path'
    os.symlink(missing_target, link_path)

    with pytest.raises(NotADirectoryError, match='BUSCO output path exists but is not a directory'):
        ensure_clean_dir(str(link_path), redo=False)


def test_normalize_busco_table_from_gzip(tmp_path):
    src = tmp_path / "full_table.tsv.gz"
    with gzip.open(str(src), 'wt') as f:
        f.write("# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
        f.write("BUSCO1\tComplete\tseq1\t100\t200\turl1\tdesc1\n")
    dest = tmp_path / "normalized.tsv"

    normalize_busco_table(str(src), str(dest))

    content = dest.read_text()
    assert content.startswith("# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription")
    assert "BUSCO1" in content


def test_normalize_busco_table_from_uppercase_gzip_extension(tmp_path):
    src = tmp_path / "full_table.TSV.GZ"
    with gzip.open(str(src), 'wt') as f:
        f.write("# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
        f.write("BUSCO1\tComplete\tseq1\t100\t200\turl1\tdesc1\n")
    dest = tmp_path / "normalized.tsv"

    normalize_busco_table(str(src), str(dest))

    content = dest.read_text()
    assert content.startswith("# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription")
    assert "BUSCO1" in content


def test_resolve_species_fasta(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    fasta_path = fasta_dir / "Homo_sapiens.fa"
    fasta_path.write_text(">seq\nATGC\n")
    result = resolve_species_fasta("Homo sapiens", str(fasta_dir))
    assert os.path.realpath(result) == os.path.realpath(str(fasta_path))


def test_resolve_species_fasta_accepts_uppercase_suffix(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    fasta_path = fasta_dir / "Homo_sapiens.FA"
    fasta_path.write_text(">seq\nATGC\n")
    result = resolve_species_fasta("Homo sapiens", str(fasta_dir))
    assert os.path.realpath(result) == os.path.realpath(str(fasta_path))


def test_resolve_species_fasta_multiple_raises(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    (fasta_dir / "Homo_sapiens.fa").write_text(">seq1\nATGC\n")
    (fasta_dir / "Homo_sapiens.v1.fa").write_text(">seq2\nATGC\n")
    with pytest.raises(ValueError, match="Found multiple reference fasta files"):
        resolve_species_fasta("Homo sapiens", str(fasta_dir))


def test_resolve_species_fasta_multiple_raises_with_unsorted_prefetched_names(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    (fasta_dir / "Homo_sapiens.fa").write_text(">seq1\nATGC\n")
    (fasta_dir / "Homo_sapiens.v1.fa").write_text(">seq2\nATGC\n")
    unsorted_names = [
        "Aardvark.fa",
        "Homo_sapiens.fa",
        "Mus_musculus.fa",
        "Homo_sapiens.v1.fa",
    ]
    with pytest.raises(ValueError, match="Found multiple reference fasta files"):
        resolve_species_fasta(
            "Homo sapiens",
            str(fasta_dir),
            fasta_filenames=unsorted_names,
        )


def test_resolve_species_fasta_ignores_similar_species_prefix(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    target = fasta_dir / "Homo_sapiens.fa"
    target.write_text(">seq\nATGC\n")
    (fasta_dir / "Homo_sapiens2.fa").write_text(">seq\nATGC\n")
    result = resolve_species_fasta("Homo sapiens", str(fasta_dir))
    assert os.path.realpath(result) == os.path.realpath(str(target))


def test_resolve_species_fasta_fallback_strips_dot_in_prefix(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    target = fasta_dir / "C_elegans.fa"
    target.write_text(">seq\nATGC\n")
    result = resolve_species_fasta("C. elegans", str(fasta_dir))
    assert os.path.realpath(result) == os.path.realpath(str(target))


def test_resolve_species_fasta_fallback_uses_genus_species_prefix(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    target = fasta_dir / "Canis_lupus.fa"
    target.write_text(">seq\nATGC\n")
    result = resolve_species_fasta("Canis lupus familiaris", str(fasta_dir))
    assert os.path.realpath(result) == os.path.realpath(str(target))


def test_resolve_species_fasta_fallback_handles_redundant_spaces(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    target = fasta_dir / "Canis_lupus.fa"
    target.write_text(">seq\nATGC\n")
    result = resolve_species_fasta("Canis   lupus familiaris", str(fasta_dir))
    assert os.path.realpath(result) == os.path.realpath(str(target))


def test_collect_species_reuses_prefetched_filenames(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    fasta_dir = out_dir / "fasta"
    fasta_dir.mkdir(parents=True)
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'scientific_name': ['Species A', 'Species B'],
        'run': ['R1', 'R2'],
        'exclusion': ['no', 'no'],
    }))
    args = SimpleNamespace(out_dir=str(out_dir), fasta_dir='inferred', fasta=None)

    prefetched = ['Species_A.fa', 'Species_B.fa']
    calls = []

    def fake_list_fasta_filenames(_fasta_dir):
        return prefetched

    def fake_resolve_species_fasta(sp, _fasta_dir, fasta_filenames=None):
        calls.append(fasta_filenames)
        return os.path.join(str(fasta_dir), sp.replace(' ', '_') + '.fa')

    monkeypatch.setattr('amalgkit.busco.list_fasta_filenames', fake_list_fasta_filenames)
    monkeypatch.setattr('amalgkit.busco.resolve_species_fasta', fake_resolve_species_fasta)

    species, fasta_map = collect_species(args, metadata)

    assert species == ['Species A', 'Species B']
    assert len(calls) == 2
    assert calls[0] is prefetched
    assert calls[1] is prefetched
    assert fasta_map['Species A'].endswith('Species_A.fa')
    assert fasta_map['Species B'].endswith('Species_B.fa')


def test_collect_species_ignores_blank_scientific_name_entries(tmp_path):
    out_dir = tmp_path / "out"
    fasta_dir = out_dir / "fasta"
    fasta_dir.mkdir(parents=True)
    (fasta_dir / "Species_A.fa").write_text(">a\nAAAA\n")
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'scientific_name': ['', 'Species A'],
        'run': ['R1', 'R2'],
        'exclusion': ['no', 'no'],
    }))
    args = SimpleNamespace(out_dir=str(out_dir), fasta_dir='inferred', fasta=None)

    species, fasta_map = collect_species(args, metadata)

    assert species == ['Species A']
    assert fasta_map == {'Species A': os.path.join(str(fasta_dir), 'Species_A.fa')}


def test_collect_species_raises_when_no_valid_scientific_name(tmp_path):
    out_dir = tmp_path / "out"
    fasta_dir = out_dir / "fasta"
    fasta_dir.mkdir(parents=True)
    (fasta_dir / "Species_A.fa").write_text(">a\nAAAA\n")
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'scientific_name': ['', float('nan')],
        'run': ['R1', 'R2'],
        'exclusion': ['no', 'no'],
    }))
    args = SimpleNamespace(out_dir=str(out_dir), fasta_dir='inferred', fasta=None)

    with pytest.raises(ValueError, match='No valid scientific_name entries were found in metadata'):
        collect_species(args, metadata)


def test_collect_species_rejects_missing_scientific_name_column(tmp_path):
    out_dir = tmp_path / "out"
    fasta_dir = out_dir / "fasta"
    fasta_dir.mkdir(parents=True)
    (fasta_dir / "Species_A.fa").write_text(">a\nAAAA\n")
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['R1'],
        'exclusion': ['no'],
    }))
    metadata.df = metadata.df.drop(columns=['scientific_name'])
    args = SimpleNamespace(out_dir=str(out_dir), fasta_dir='inferred', fasta=None)

    with pytest.raises(ValueError, match='Missing required metadata column\\(s\\) for busco: scientific_name'):
        collect_species(args, metadata)


def test_collect_species_with_explicit_fasta_strips_species_name(tmp_path):
    fasta_path = tmp_path / "input.fa"
    fasta_path.write_text(">a\nAAAA\n")
    args = SimpleNamespace(
        out_dir=str(tmp_path),
        fasta_dir='inferred',
        fasta=str(fasta_path),
        species=' Species A ',
    )

    species, fasta_map = collect_species(args, metadata=None)

    assert species == ['Species A']
    assert fasta_map == {'Species A': os.path.realpath(str(fasta_path))}


def test_collect_species_with_explicit_fasta_rejects_blank_species(tmp_path):
    fasta_path = tmp_path / "input.fa"
    fasta_path.write_text(">a\nAAAA\n")
    args = SimpleNamespace(
        out_dir=str(tmp_path),
        fasta_dir='inferred',
        fasta=str(fasta_path),
        species='   ',
    )

    with pytest.raises(ValueError, match='--species must not be empty'):
        collect_species(args, metadata=None)


def test_collect_species_with_explicit_fasta_rejects_missing_file(tmp_path):
    args = SimpleNamespace(
        out_dir=str(tmp_path),
        fasta_dir='inferred',
        fasta=str(tmp_path / "missing.fa"),
        species='Species A',
    )

    with pytest.raises(FileNotFoundError, match='FASTA file not found'):
        collect_species(args, metadata=None)


def test_collect_species_with_explicit_fasta_rejects_directory_path(tmp_path):
    fasta_dir = tmp_path / "fasta_dir"
    fasta_dir.mkdir()
    args = SimpleNamespace(
        out_dir=str(tmp_path),
        fasta_dir='inferred',
        fasta=str(fasta_dir),
        species='Species A',
    )

    with pytest.raises(IsADirectoryError, match='not a file'):
        collect_species(args, metadata=None)


def test_collect_species_rejects_fasta_dir_file_path(tmp_path):
    fasta_path = tmp_path / "fasta_path"
    fasta_path.write_text("not a directory")
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'scientific_name': ['Species A'],
        'run': ['R1'],
        'exclusion': ['no'],
    }))
    args = SimpleNamespace(
        out_dir=str(tmp_path),
        fasta_dir=str(fasta_path),
        fasta=None,
    )

    with pytest.raises(NotADirectoryError, match='not a directory'):
        collect_species(args, metadata)


def test_select_tool_auto_prefers_compleasm(monkeypatch):
    args = SimpleNamespace(tool='auto', compleasm_exe='compleasm', busco_exe='busco')

    def fake_which(name):
        if name == 'compleasm':
            return '/usr/bin/compleasm'
        return None

    monkeypatch.setattr('shutil.which', fake_which)
    assert select_tool(args) == 'compleasm'


def test_select_tool_auto_falls_back_to_busco(monkeypatch):
    args = SimpleNamespace(tool='auto', compleasm_exe='compleasm', busco_exe='busco')

    def fake_which(name):
        if name == 'busco':
            return '/usr/bin/busco'
        return None

    monkeypatch.setattr('shutil.which', fake_which)
    assert select_tool(args) == 'busco'


def test_run_busco_places_downloads_under_out_dir(tmp_path, monkeypatch):
    out_dir = tmp_path / 'out'
    busco_root = out_dir / 'busco'
    out_dir.mkdir()
    busco_root.mkdir()
    args = SimpleNamespace(
        busco_exe='busco',
        lineage='eukaryota_odb12',
        threads=4,
        redo=False,
        out_dir=str(out_dir),
    )
    captured = {}

    def fake_run_command(cmd):
        captured['cmd'] = cmd

    monkeypatch.setattr('amalgkit.busco.run_command', fake_run_command)
    output_dir = run_busco(
        fasta_path='/tmp/input.fa',
        sci_name='Species A',
        output_root=str(busco_root),
        args=args,
        extra_args=[],
    )

    assert output_dir == os.path.join(str(busco_root), 'Species_A')
    assert '--download_path' in captured['cmd']
    idx = captured['cmd'].index('--download_path')
    assert captured['cmd'][idx + 1] == os.path.join(str(out_dir), 'busco_downloads')
    assert (out_dir / 'busco_downloads').exists()


def test_run_busco_rejects_download_path_that_is_file(tmp_path):
    out_dir = tmp_path / 'out'
    busco_root = out_dir / 'busco'
    out_dir.mkdir()
    busco_root.mkdir()
    (out_dir / 'busco_downloads').write_text('not a directory')
    args = SimpleNamespace(
        busco_exe='busco',
        lineage='eukaryota_odb12',
        threads=4,
        redo=False,
        out_dir=str(out_dir),
    )

    with pytest.raises(NotADirectoryError, match='BUSCO download path exists but is not a directory'):
        run_busco(
            fasta_path='/tmp/input.fa',
            sci_name='Species A',
            output_root=str(busco_root),
            args=args,
            extra_args=[],
        )


def test_busco_main_rejects_nonpositive_species_jobs(tmp_path):
    args = SimpleNamespace(
        lineage='eukaryota_odb12',
        internal_jobs=0,
        tool='auto',
        tool_args=None,
        fasta=None,
        out_dir=str(tmp_path),
        metadata='inferred',
        fasta_dir='inferred',
        threads=1,
        redo=False,
        busco_exe='busco',
        compleasm_exe='compleasm',
    )
    with pytest.raises(ValueError, match='--internal_jobs must be > 0'):
        busco_main(args)


def test_busco_main_rejects_blank_lineage(tmp_path):
    args = SimpleNamespace(
        lineage='   ',
        internal_jobs=1,
        tool='auto',
        tool_args=None,
        fasta=None,
        out_dir=str(tmp_path),
        metadata='inferred',
        fasta_dir='inferred',
        threads=1,
        redo=False,
        busco_exe='busco',
        compleasm_exe='compleasm',
    )
    with pytest.raises(ValueError, match='--lineage is required'):
        busco_main(args)


def test_busco_main_rejects_out_dir_file_path(tmp_path):
    out_path = tmp_path / 'out_path'
    out_path.write_text('not a directory')
    args = SimpleNamespace(
        lineage='eukaryota_odb12',
        internal_jobs=1,
        tool='auto',
        tool_args=None,
        fasta=None,
        out_dir=str(out_path),
        metadata='inferred',
        fasta_dir='inferred',
        threads=1,
        redo=False,
        busco_exe='busco',
        compleasm_exe='compleasm',
    )
    with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
        busco_main(args)


def test_busco_main_rejects_busco_dir_file_path(tmp_path):
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    (out_dir / 'busco').write_text('not a directory')
    args = SimpleNamespace(
        lineage='eukaryota_odb12',
        internal_jobs=1,
        tool='auto',
        tool_args=None,
        fasta=None,
        out_dir=str(out_dir),
        metadata='inferred',
        fasta_dir='inferred',
        threads=1,
        redo=False,
        busco_exe='busco',
        compleasm_exe='compleasm',
    )
    with pytest.raises(NotADirectoryError, match='BUSCO path exists but is not a directory'):
        busco_main(args)


def test_busco_main_parallel_species_jobs(tmp_path, monkeypatch):
    out_dir = tmp_path / 'out'
    fasta_dir = out_dir / 'fasta'
    out_dir.mkdir()
    fasta_dir.mkdir(parents=True)
    (fasta_dir / 'Species_A.fa').write_text('>a\nAAAA\n')
    (fasta_dir / 'Species_B.fa').write_text('>b\nCCCC\n')

    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'scientific_name': ['Species A', 'Species B'],
        'run': ['R1', 'R2'],
        'exclusion': ['no', 'no'],
    }))
    args = SimpleNamespace(
        lineage='eukaryota_odb12',
        internal_jobs=2,
        tool='auto',
        tool_args=None,
        fasta=None,
        out_dir=str(out_dir),
        metadata='inferred',
        fasta_dir='inferred',
        threads=1,
        redo=False,
        busco_exe='busco',
        compleasm_exe='compleasm',
    )
    processed = []

    def fake_run_busco(fasta_path, sci_name, output_root, _args, _extra):
        processed.append(sci_name)
        out_species = os.path.join(output_root, sci_name.replace(' ', '_'))
        os.makedirs(out_species, exist_ok=True)
        table = os.path.join(out_species, 'full_table.tsv')
        with open(table, 'w') as f:
            f.write('Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n')
            f.write('BUSCO1\tComplete\tseq1\t100\t200\turl\tdesc\n')
        return out_species

    monkeypatch.setattr('amalgkit.busco.load_metadata', lambda _args: metadata)
    monkeypatch.setattr('amalgkit.busco.select_tool', lambda _args: 'busco')
    monkeypatch.setattr('amalgkit.busco.run_busco', fake_run_busco)

    busco_main(args)

    assert set(processed) == {'Species A', 'Species B'}
    assert (out_dir / 'busco' / 'Species_A_busco.tsv').exists()
    assert (out_dir / 'busco' / 'Species_B_busco.tsv').exists()


def test_busco_main_cpu_budget_caps_species_jobs_to_serial(tmp_path, monkeypatch):
    out_dir = tmp_path / 'out'
    fasta_dir = out_dir / 'fasta'
    out_dir.mkdir()
    fasta_dir.mkdir(parents=True)
    (fasta_dir / 'Species_A.fa').write_text('>a\nAAAA\n')
    (fasta_dir / 'Species_B.fa').write_text('>b\nCCCC\n')

    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'scientific_name': ['Species A', 'Species B'],
        'run': ['R1', 'R2'],
        'exclusion': ['no', 'no'],
    }))
    args = SimpleNamespace(
        lineage='eukaryota_odb12',
        internal_jobs=4,
        internal_cpu_budget=1,
        threads=2,
        tool='auto',
        tool_args=None,
        fasta=None,
        out_dir=str(out_dir),
        metadata='inferred',
        fasta_dir='inferred',
        redo=False,
        busco_exe='busco',
        compleasm_exe='compleasm',
    )
    processed = []

    def fake_run_busco(fasta_path, sci_name, output_root, _args, _extra):
        processed.append(sci_name)
        out_species = os.path.join(output_root, sci_name.replace(' ', '_'))
        os.makedirs(out_species, exist_ok=True)
        table = os.path.join(out_species, 'full_table.tsv')
        with open(table, 'w') as f:
            f.write('Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n')
            f.write('BUSCO1\tComplete\tseq1\t100\t200\turl\tdesc\n')
        return out_species

    def fail_if_called(*_args, **_kwargs):
        raise AssertionError('run_tasks_with_optional_threads should not be used when --internal_cpu_budget caps internal_jobs to 1.')

    monkeypatch.setattr('amalgkit.busco.load_metadata', lambda _args: metadata)
    monkeypatch.setattr('amalgkit.busco.select_tool', lambda _args: 'busco')
    monkeypatch.setattr('amalgkit.busco.run_busco', fake_run_busco)
    monkeypatch.setattr('amalgkit.busco.run_tasks_with_optional_threads', fail_if_called)

    busco_main(args)

    assert set(processed) == {'Species A', 'Species B'}
    assert args.threads == 1
