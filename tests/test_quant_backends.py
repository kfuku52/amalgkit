import os
import json
import pathlib
import subprocess

import pytest
import pandas
from types import SimpleNamespace

from amalgkit.quant import (
    build_kallisto_quant_command,
    build_oarfish_quant_command,
    check_layout_mismatch,
    resolve_input_fastq_files,
    find_species_index_files,
    find_species_fasta_files,
    check_kallisto_dependency,
    check_oarfish_dependency,
    resolve_quant_backend,
    resolve_oarfish_seq_tech,
)
from amalgkit.util import Metadata


def write_valid_quant_outputs(output_dir, sra_id='SRR001', target_id='tx1'):
    os.makedirs(output_dir, exist_ok=True)
    pandas.DataFrame([{
        'target_id': target_id,
        'length': 100,
        'eff_length': 90,
        'est_counts': 5,
        'tpm': 1.0,
    }]).to_csv(os.path.join(output_dir, sra_id + '_abundance.tsv'), sep='\t', index=False)
    with open(os.path.join(output_dir, sra_id + '_run_info.json'), 'w', encoding='utf-8') as handle:
        json.dump({'p_pseudoaligned': 50.0}, handle)


def write_safely_removed_fastq_marker(out_dir, sra_id='SRR001'):
    run_dir = pathlib.Path(out_dir) / 'getfastq' / sra_id
    run_dir.mkdir(parents=True, exist_ok=True)
    marker = run_dir / (sra_id + '.fastq.gz.safely_removed')
    marker.write_text(
        'This fastq file was safely removed after `amalgkit quant`.',
        encoding='utf-8',
    )
    return marker


def make_single_run_quant_metadata():
    return Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['SRR001'],
        'scientific_name': ['Species A'],
        'lib_layout': ['single'],
        'total_spots': [10],
        'spot_length': [100],
        'total_bases': [1000],
        'nominal_length': [200],
        'exclusion': ['no'],
    }))


# ---------------------------------------------------------------------------
# quant_output_exists (checks for abundance.tsv in output directory)
# ---------------------------------------------------------------------------


class TestCheckKallistoDependency:
    def test_raises_clear_error_when_kallisto_missing(self, monkeypatch):
        def fake_run(_cmd, stdout=None, stderr=None):
            raise FileNotFoundError('kallisto')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        with pytest.raises(FileNotFoundError, match='kallisto executable not found'):
            check_kallisto_dependency()

    def test_raises_when_kallisto_probe_returns_nonzero(self, monkeypatch):
        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 127, stdout=b'', stderr=b'')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        with pytest.raises(RuntimeError, match='kallisto dependency probe failed'):
            check_kallisto_dependency()


class TestCheckOarfishDependency:
    def test_raises_clear_error_when_oarfish_missing(self, monkeypatch):
        def fake_run(_cmd, stdout=None, stderr=None):
            raise FileNotFoundError('oarfish')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        with pytest.raises(FileNotFoundError, match='oarfish executable not found'):
            check_oarfish_dependency()

    def test_raises_when_oarfish_probe_returns_nonzero(self, monkeypatch):
        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 127, stdout=b'', stderr=b'')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        with pytest.raises(RuntimeError, match='oarfish dependency probe failed'):
            check_oarfish_dependency()


class TestQuantBackendAuto:
    def test_resolve_quant_backend_auto_uses_short_read_heuristic(self):
        args = SimpleNamespace(quant_backend='auto')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'lib_layout': ['single'],
            'spot_length': [150],
            'total_spots': [10],
            'total_bases': [1500],
            'exclusion': ['no'],
        }))

        observed = resolve_quant_backend(args, metadata, 'SRR001')

        assert observed == 'kallisto'

    def test_resolve_quant_backend_auto_detects_pacbio_platform(self):
        args = SimpleNamespace(quant_backend='auto')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'platform': ['PACBIO_SMRT'],
            'instrument': ['Sequel II'],
            'lib_layout': ['single'],
            'spot_length': [2500],
            'total_spots': [10],
            'total_bases': [25000],
            'exclusion': ['no'],
        }))

        observed = resolve_quant_backend(args, metadata, 'SRR001')

        assert observed == 'oarfish'

    def test_resolve_oarfish_seq_tech_auto_detects_ont_direct_rna(self):
        args = SimpleNamespace(quant_backend='auto', oarfish_seq_tech='auto')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'platform': ['OXFORD_NANOPORE'],
            'instrument': ['PromethION'],
            'sample_description': ['direct RNA sequencing library'],
            'lib_layout': ['single'],
            'spot_length': [1500],
            'total_spots': [10],
            'total_bases': [15000],
            'exclusion': ['no'],
        }))

        observed = resolve_oarfish_seq_tech(args, metadata, 'SRR001')

        assert observed == 'ont-drna'


class TestQuantOptionPassthrough:
    def test_build_kallisto_quant_command_appends_user_options_before_input(self):
        args = SimpleNamespace(threads=4, kallisto_options='--bias --seed 42')

        observed = build_kallisto_quant_command(
            args=args,
            in_files=['reads.fastq.gz'],
            lib_layout='single',
            output_dir='out_dir',
            index='ref.idx',
            nominal_length=250,
            fragment_sd=25,
        )

        assert observed == [
            'kallisto', 'quant', '--threads', '4', '--index', 'ref.idx', '-o', 'out_dir',
            '--single', '-l', '250', '-s', '25', '--bias', '--seed', '42', 'reads.fastq.gz',
        ]

    def test_build_oarfish_quant_command_appends_user_options_before_output(self):
        args = SimpleNamespace(threads=8, oarfish_options='--filter-group no-filters --best-n 25')

        observed = build_oarfish_quant_command(
            args=args,
            in_files=['reads.fastq.gz'],
            output_prefix='out_prefix',
            index='ref.mmi',
            seq_tech='ont-drna',
        )

        assert observed == [
            'oarfish', '-j', '8', '--reads', 'reads.fastq.gz', '--index', 'ref.mmi',
            '--seq-tech', 'ont-drna', '--filter-group', 'no-filters', '--best-n', '25',
            '-o', 'out_prefix',
        ]

    def test_build_kallisto_quant_command_rejects_invalid_option_string(self):
        args = SimpleNamespace(threads=1, kallisto_options='"unterminated')

        with pytest.raises(ValueError, match='Invalid --kallisto_options string'):
            build_kallisto_quant_command(
                args=args,
                in_files=['reads.fastq.gz'],
                lib_layout='single',
                output_dir='out_dir',
                index='ref.idx',
                nominal_length=200,
                fragment_sd=20,
            )

    def test_build_oarfish_quant_command_rejects_invalid_option_string(self):
        args = SimpleNamespace(threads=1, oarfish_options='"unterminated')

        with pytest.raises(ValueError, match='Invalid --oarfish_options string'):
            build_oarfish_quant_command(
                args=args,
                in_files=['reads.fastq.gz'],
                output_prefix='out_prefix',
                index='ref.mmi',
                seq_tech='pac-bio',
            )


# ---------------------------------------------------------------------------
# check_layout_mismatch (issue #80: corrects layout when files disagree)
# ---------------------------------------------------------------------------

class TestCheckLayoutMismatch:
    def test_paired_metadata_single_file(self, tmp_path):
        """Paired layout corrected to single when only one fastq found."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('data')
        result = check_layout_mismatch(sra_stat, str(tmp_path))
        assert result['layout'] == 'single'

    def test_paired_metadata_paired_files(self, tmp_path):
        """Paired layout unchanged when two fastq files found."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.amalgkit.fastq.gz').write_text('data')
        result = check_layout_mismatch(sra_stat, str(tmp_path))
        assert result['layout'] == 'paired'

    def test_single_metadata_no_change(self, tmp_path):
        """Single layout is not affected by mismatch check."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        result = check_layout_mismatch(sra_stat, str(tmp_path))
        assert result['layout'] == 'single'

    def test_uses_prefetched_files_without_rescan(self, tmp_path, monkeypatch):
        """When files are given, directory listing helper should not be called."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        files = {'SRR001.fastq.gz'}

        def fail_if_called(_output_dir):
            raise AssertionError('list_getfastq_run_files should not be called when files are provided.')

        monkeypatch.setattr('amalgkit.quant.list_getfastq_run_files', fail_if_called)
        result = check_layout_mismatch(sra_stat, str(tmp_path), files=files)
        assert result['layout'] == 'single'

    def test_ignores_similar_prefix_run_ids_when_counting_fastq_files(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        files = {'SRR0010.fastq.gz'}
        result = check_layout_mismatch(sra_stat, str(tmp_path), files=files)
        assert result['layout'] == 'paired'


class TestResolveInputFastqFiles:
    def test_paired_returns_deterministic_order(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        files = {'SRR001_2.fastq.gz', 'SRR001_1.fastq.gz'}
        result = resolve_input_fastq_files(
            sra_stat=sra_stat,
            output_dir_getfastq=str(tmp_path),
            ext='.fastq.gz',
            files=files,
        )
        assert result == [
            str(tmp_path / 'SRR001_1.fastq.gz'),
            str(tmp_path / 'SRR001_2.fastq.gz'),
        ]

    def test_single_returns_single_file(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        files = {'SRR001.fastq.gz'}
        result = resolve_input_fastq_files(
            sra_stat=sra_stat,
            output_dir_getfastq=str(tmp_path),
            ext='.fastq.gz',
            files=files,
        )
        assert result == [str(tmp_path / 'SRR001.fastq.gz')]

    def test_fallback_ignores_similar_prefix_run_ids(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        files = {'SRR0010.fastq.gz'}
        result = resolve_input_fastq_files(
            sra_stat=sra_stat,
            output_dir_getfastq=str(tmp_path),
            ext='.fastq.gz',
            files=files,
        )
        assert result == []


class TestIndexDiscoveryHelpers:
    def test_find_species_index_files_requires_exact_stem(self, tmp_path):
        (tmp_path / 'Homo_sapiens.idx').write_text('idx')
        (tmp_path / 'Homo_sapiens.v2.idx').write_text('idx')
        (tmp_path / 'Mus_musculus.idx').write_text('idx')
        result = find_species_index_files(str(tmp_path), 'Homo_sapiens')
        assert result == [
            str(tmp_path / 'Homo_sapiens.idx'),
        ]

    def test_find_species_fasta_files_filters_suffix(self, tmp_path):
        (tmp_path / 'Homo_sapiens.fa').write_text('>a\nAAAA\n')
        (tmp_path / 'Homo_sapiens.fasta.gz').write_text('gz')
        (tmp_path / 'Homo_sapiens.txt').write_text('no')
        (tmp_path / 'Mus_musculus.fa').write_text('>a\nAAAA\n')
        result = find_species_fasta_files(str(tmp_path), 'Homo_sapiens')
        assert result == [
            str(tmp_path / 'Homo_sapiens.fa'),
            str(tmp_path / 'Homo_sapiens.fasta.gz'),
        ]

    def test_find_species_fasta_files_accepts_uppercase_suffix(self, tmp_path):
        (tmp_path / 'Homo_sapiens.FA').write_text('>a\nAAAA\n')
        (tmp_path / 'Homo_sapiens.FASTA.GZ').write_text('gz')
        result = find_species_fasta_files(str(tmp_path), 'Homo_sapiens')
        assert result == [
            str(tmp_path / 'Homo_sapiens.FA'),
            str(tmp_path / 'Homo_sapiens.FASTA.GZ'),
        ]

    def test_find_species_fasta_files_normalizes_spaces_to_underscores(self, tmp_path):
        (tmp_path / 'Homo_sapiens.fa').write_text('>a\nAAAA\n')
        result = find_species_fasta_files(str(tmp_path), 'Homo sapiens')
        assert result == [str(tmp_path / 'Homo_sapiens.fa')]

    def test_find_species_index_files_rejects_similar_species_prefix(self, tmp_path):
        (tmp_path / 'Homo_sapiens.idx').write_text('idx')
        (tmp_path / 'Homo_sapiens2.idx').write_text('idx')
        result = find_species_index_files(str(tmp_path), 'Homo_sapiens')
        assert result == [str(tmp_path / 'Homo_sapiens.idx')]

    def test_find_species_index_files_ignores_directories(self, tmp_path):
        (tmp_path / 'Homo_sapiens.idx').mkdir()
        (tmp_path / 'Homo_sapiens.v2.idx').write_text('idx')
        result = find_species_index_files(str(tmp_path), 'Homo_sapiens')
        assert result == []

    def test_find_species_index_files_ignores_empty_index_files(self, tmp_path):
        (tmp_path / 'Homo_sapiens.idx').touch()
        result = find_species_index_files(str(tmp_path), 'Homo_sapiens')
        assert result == []

    def test_find_species_fasta_files_raises_for_file_path(self, tmp_path):
        fasta_path = tmp_path / 'fasta_path'
        fasta_path.write_text('not a directory')
        with pytest.raises(NotADirectoryError, match='Fasta path exists but is not a directory'):
            find_species_fasta_files(str(fasta_path), 'Homo_sapiens')

    def test_find_species_fasta_files_ignores_directories(self, tmp_path):
        (tmp_path / 'Homo_sapiens.fa').mkdir()
        (tmp_path / 'Homo_sapiens.fasta').write_text('>a\nAAAA\n')
        result = find_species_fasta_files(str(tmp_path), 'Homo_sapiens')
        assert result == [str(tmp_path / 'Homo_sapiens.fasta')]

    def test_find_species_fasta_files_requires_exact_species_stem(self, tmp_path):
        (tmp_path / 'Gorilla_gorilla.fa').write_text('>a\nAAAA\n')
        result = find_species_fasta_files(str(tmp_path), 'Gorilla_gorilla_gorilla')
        assert result == []

    def test_find_species_fasta_files_normalizes_redundant_underscores_without_shortening(self, tmp_path):
        (tmp_path / 'Canis_lupus_familiaris.fa').write_text('>a\nAAAA\n')
        result = find_species_fasta_files(str(tmp_path), 'Canis__lupus_familiaris')
        assert result == [str(tmp_path / 'Canis_lupus_familiaris.fa')]
