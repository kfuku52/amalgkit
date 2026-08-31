"""Execute documentation handoffs on tiny inputs; no SRA or taxonomy downloads."""

import json
import re
import runpy
import shlex
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas
import numpy
import pytest

from amalgkit.main import build_main_parser
from amalgkit.metadata_utils import Metadata, load_metadata


ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.integration


def _commands(page, command):
    shell_commands = runpy.run_path(str(ROOT / '.github/scripts/check_docs.py'))['shell_commands']
    return [shlex.split(source)[1:] for _, source in shell_commands((ROOT / '.wiki' / page).read_text())
            if shlex.split(source)[1] == command]


def _run(arguments):
    args = build_main_parser().parse_args(arguments)
    args.handler(args)


def _read(path):
    return pandas.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def test_yeast_tutorial_preserves_annotations_and_selects_reviewed_groups(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(_commands('Tutorial-1.md', 'dataset')[0])
    source = tmp_path / 'metadata' / 'metadata.tsv'
    original = (
        'run\tscientific_name\ttaxid\ttaxid_species\tsample_group\tgenotype\tnote\ttotal_spots\tbioproject\tbiosample\texclusion\n'
        '0001\tSaccharomyces cerevisiae\t4932\t4932\told\tWT\tNA\t2000000\tP1\tB1\tno\n'
        'NA\tSaccharomyces cerevisiae\t4932\t4932\theat_stress\t\tNA\t2000000\tP1\tB2\tno\n'
        'pombe1\tSchizosaccharomyces pombe\t4896\t4896\twild_type\t\tN/A\t2000000\tP2\tB3\tno\n'
    )
    source.write_text(original)
    tutorial = (ROOT / '.wiki/Tutorial-1.md').read_text()
    snippet = re.search(r"python - <<'PY'\n(.*?)\nPY", tutorial, re.S)[1]
    subprocess.run([sys.executable, '-c', snippet], cwd=tmp_path, check=True)
    assert source.read_text() == original
    reviewed = _read(tmp_path / 'metadata' / 'metadata_grouped.tsv').set_index('run')
    assert reviewed.loc['0001', 'sample_group'] == 'wt'
    assert reviewed.loc['NA', 'sample_group'] == 'heat_stress'
    assert reviewed.loc['0001', 'note'] == 'NA'
    assert reviewed.loc['pombe1', 'note'] == 'N/A'
    _run(_commands('Tutorial-1.md', 'select')[0])
    selected = _read(source).set_index('run')
    assert selected['sample_group'].to_dict() == {'0001': 'wt', 'NA': 'heat stress', 'pombe1': 'wild type'}
    assert selected['is_sampled'].eq('yes').all()
    assert selected['exclusion'].eq('no').all()
    assert selected['note'].to_dict() == {'0001': 'NA', 'NA': 'NA', 'pombe1': 'N/A'}


@pytest.mark.parametrize('merge_existing', [False, True])
def test_private_integration_docs_handoff_all_runs(tmp_path, monkeypatch, merge_existing):
    monkeypatch.chdir(tmp_path)
    fastq = tmp_path / 'private_fastq' / 'Saccharomyces_cerevisiae' / 'Private1.fastq'
    fastq.parent.mkdir(parents=True)
    fastq.write_text('@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nIIII\n')
    source = tmp_path / 'metadata' / 'metadata.tsv'
    original = 'run\tscientific_name\ttaxid\tsample_group\tis_sampled\texclusion\nPUBLIC\tSaccharomyces cerevisiae\t4932\twt\tyes\tno\n'
    if merge_existing:
        source.parent.mkdir()
        source.write_text(original)
    taxonomy = SimpleNamespace(
        get_name_translator=lambda names: {name: [4932] for name in names},
        get_lineage=lambda taxid: [taxid],
        get_rank=lambda taxids: dict.fromkeys(taxids, 'species'),
    )
    monkeypatch.setattr('amalgkit.integrate.get_ncbi_taxonomy', lambda args=None: taxonomy)
    parser = build_main_parser()
    args = parser.parse_args(_commands('amalgkit-integrate.md', 'integrate')[int(merge_existing)])
    # Exercise the actual Python FASTQ scanner on CI machines without SeqKit.
    args.seqkit_exe = str(tmp_path / 'uninstalled-seqkit')
    args.handler(args)
    expected = {'PUBLIC', 'Private1'} if merge_existing else {'Private1'}
    for command in ('getfastq', 'quant', 'merge'):
        consumer = parser.parse_args(_commands('amalgkit-integrate.md', command)[0])
        assert set(load_metadata(consumer).df['run']) == expected
    if merge_existing:
        assert source.read_text() == original
    else:
        assert not source.exists()


@pytest.mark.parametrize('instructions', ['readme', 'stdout'])
def test_generated_species_workspace_guide_connects_metadata_and_select(tmp_path, monkeypatch, capsys, instructions):
    monkeypatch.chdir(tmp_path)
    _run(['dataset', '--name', 'init', '--out_dir', './'])
    printed = capsys.readouterr().out
    (tmp_path / 'species.tsv').write_text('scientific_name\nArabidopsis thaliana\n')
    frame = Metadata.from_DataFrame(pandas.DataFrame([{
        'run': '0001', 'scientific_name': 'Arabidopsis thaliana', 'taxid': '3702', 'taxid_species': '3702',
        'sample_group': 'leaf', 'total_spots': 2000000, 'bioproject': 'P1', 'biosample': 'B1', 'exclusion': 'no',
    }]))

    def query_fixture(**kwargs):
        directory = Path(kwargs['out_dir'])
        metadata_path = directory / 'metadata' / 'metadata.tsv'
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        frame.df.to_csv(metadata_path, sep='\t', index=False)
        info = {'record_id_count': 1, 'used_cached_metadata': False}
        info_path = directory / 'query_info.json'
        info_path.write_text(json.dumps(info))
        return {'metadata': frame, 'query_label': kwargs['query_label'], 'query_info': info,
                'paths': {'query_info_path': str(info_path), 'metadata_path': str(metadata_path)}}

    monkeypatch.setattr('amalgkit.metadata._run_single_query', query_fixture)
    guide = (tmp_path / 'WORKSPACE_README.md').read_text()
    commands = (re.findall(r'\d+\. Run `(amalgkit [^`]+)`', guide) if instructions == 'readme'
                else re.findall(r'^  (amalgkit .+)$', printed, re.M))
    assert len(commands) == 2
    for command in commands:
        _run(shlex.split(command)[1:])
    manifest = _read(tmp_path / 'batch' / 'external_manifest.tsv')
    assert len(manifest) == 1
    selected_path = Path(manifest.loc[0, 'selected_metadata_path'])
    selected = _read(selected_path)
    assert selected.loc[0, 'run'] == '0001'
    assert selected.loc[0, 'is_sampled'] == 'yes'
    assert not (tmp_path / 'metadata' / 'metadata.tsv').exists()


def test_documented_long_read_chain_keeps_length_model_and_tmm_scale(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rows = []
    for species, offset in [('Species_A', 0), ('Species_B', 10)]:
        directory = tmp_path / 'merge' / species
        directory.mkdir(parents=True)
        genes = [f'{species}_g{i}' for i in range(6)]
        runs = [f'{species}_r{i}' for i in range(4)]
        counts = pandas.DataFrame({
            runs[0]: [100 + offset, 110, 120, 10, 20, 30],
            runs[1]: [110 + offset, 100, 130, 15, 25, 35],
            runs[2]: [10, 20, 30, 100 + offset, 110, 120],
            runs[3]: [15, 25, 35, 110 + offset, 100, 130],
        }, index=genes)
        counts.rename_axis('target_id').to_csv(directory / f'{species}_est_counts.tsv', sep='\t')
        (counts * 0 + 1).rename_axis('target_id').to_csv(directory / f'{species}_eff_length.tsv', sep='\t')
        pandas.DataFrame({'run': runs, 'backend': 'oarfish', 'length_model': 'none'}).to_csv(
            directory / f'{species}_quant_model.tsv', sep='\t', index=False)
        busco = tmp_path / 'busco'
        busco.mkdir(exist_ok=True)
        (busco / f'{species}_busco.tsv').write_text(''.join(
            f'OG{i}\tComplete\t{gene}\t50\t100\t-\t-\n' for i, gene in enumerate(genes)))
        rows += [{'run': run, 'scientific_name': species.replace('_', ' '), 'sample_group': group,
                  'bioproject': f'P{i % 2}', 'exclusion': 'no', 'mapping_rate': 100}
                 for i, (run, group) in enumerate(zip(runs, ['a', 'a', 'b', 'b']))]
    pandas.DataFrame(rows).to_csv(tmp_path / 'merge/metadata.tsv', sep='\t', index=False)
    for command in ('cstmm', 'wsfilter', 'csfilter', 'finalize'):
        _run(_commands('Metadata-and-normalization.md', command)[-1])
    for species in ('Species_A', 'Species_B'):
        original = pandas.read_csv(tmp_path / 'cstmm' / species / f'{species}_cstmm_counts.tsv', sep='\t', index_col=0)
        final = pandas.read_csv(tmp_path / 'finalize' / species / f'{species}_expression.tsv', sep='\t', index_col=0)
        assert len(final.columns) > 0
        numpy.testing.assert_allclose(final, numpy.log2(original.loc[final.index, final.columns] + 1))
        model = _read(tmp_path / 'cstmm' / species / f'{species}_quant_model.tsv')
        assert model['length_model'].eq('none').all()
