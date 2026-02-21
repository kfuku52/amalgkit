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
    table.write_text("Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
    found = find_full_table(str(out_dir))
    assert os.path.realpath(found) == os.path.realpath(str(table))


def test_resolve_species_fasta(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    fasta_path = fasta_dir / "Homo_sapiens.fa"
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


def test_busco_main_rejects_nonpositive_species_jobs(tmp_path):
    args = SimpleNamespace(
        lineage='eukaryota_odb12',
        species_jobs=0,
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
    with pytest.raises(ValueError, match='--species_jobs must be > 0'):
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
        species_jobs=2,
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
        species_jobs=4,
        cpu_budget=1,
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
        raise AssertionError('run_tasks_with_optional_threads should not be used when --cpu_budget caps species_jobs to 1.')

    monkeypatch.setattr('amalgkit.busco.load_metadata', lambda _args: metadata)
    monkeypatch.setattr('amalgkit.busco.select_tool', lambda _args: 'busco')
    monkeypatch.setattr('amalgkit.busco.run_busco', fake_run_busco)
    monkeypatch.setattr('amalgkit.busco.run_tasks_with_optional_threads', fail_if_called)

    busco_main(args)

    assert set(processed) == {'Species A', 'Species B'}
    assert args.threads == 1
