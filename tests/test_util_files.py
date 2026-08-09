import builtins
import json
import os
import warnings

import pytest
import pandas
import numpy

from amalgkit.util import (
    Metadata,
    detect_layout_from_file,
    is_there_unpaired_file,
    get_newest_intermediate_file_extension,
    get_mapping_rate,
    get_getfastq_run_dir,
    generate_multisp_busco_table,
)


# ---------------------------------------------------------------------------
# strtobool
# ---------------------------------------------------------------------------


class TestDetectLayoutFromFile:
    def test_paired_files_detected(self, tmp_path):
        """Paired-end files detected -> layout corrected to 'paired'."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.amalgkit.fastq.gz').write_text('data')
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'paired'

    def test_single_files_detected(self, tmp_path):
        """Single-end file detected -> layout corrected to 'single'."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('data')
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'single'

    def test_paired_plain_fastq_files_detected(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.fastq').write_text('data')
        (tmp_path / 'SRR001_2.fastq').write_text('data')
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'paired'

    def test_no_files_keeps_layout(self, tmp_path):
        """No fastq files -> layout unchanged."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'paired'

    def test_paired_safely_removed_files_detected(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.fastq.gz.safely_removed').write_text('')
        (tmp_path / 'SRR001_2.fastq.gz.safely_removed').write_text('')
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'paired'

    def test_single_safely_removed_file_detected(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.fastq.gz.safely_removed').write_text('')
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'single'


# ---------------------------------------------------------------------------
# is_there_unpaired_file
# ---------------------------------------------------------------------------

class TestIsThereUnpairedFile:
    def test_unpaired_file_present(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('data')
        assert is_there_unpaired_file(sra_stat, ['.amalgkit.fastq.gz']) is True

    def test_no_unpaired_file(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'getfastq_sra_dir': str(tmp_path),
        }
        assert is_there_unpaired_file(sra_stat, ['.amalgkit.fastq.gz']) is False


# ---------------------------------------------------------------------------
# get_newest_intermediate_file_extension
# ---------------------------------------------------------------------------

class TestGetNewestIntermediateFileExtension:
    def test_finds_amalgkit_extension(self, tmp_path):
        """Finds .amalgkit.fastq.gz as the newest extension."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.amalgkit.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_1.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.fastq.gz').write_text('data')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.amalgkit.fastq.gz'

    def test_finds_fastq_extension(self, tmp_path):
        """Falls back to .fastq.gz when no downstream extensions exist."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.fastq.gz').write_text('data')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.fastq.gz'

    def test_finds_plain_fastq_extension(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.fastq').write_text('data')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.fastq'

    def test_safely_removed(self, tmp_path):
        """Detects .safely_removed flag."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz.safely_removed').write_text('')
        (tmp_path / 'SRR001_2.amalgkit.fastq.gz.safely_removed').write_text('')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.safely_removed'

    def test_safely_removed_ignores_similar_run_prefix(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR0010_1.amalgkit.fastq.gz.safely_removed').write_text('')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == 'no_extension_found'

    def test_safely_removed_requires_both_mates_for_paired_layout(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz.safely_removed').write_text('')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == 'no_extension_found'

    def test_safely_removed_detected_despite_layout_mismatch_single_to_paired(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.fastq.gz.safely_removed').write_text('')
        (tmp_path / 'SRR001_2.fastq.gz.safely_removed').write_text('')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.safely_removed'

    def test_safely_removed_detected_despite_layout_mismatch_paired_to_single(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.fastq.gz.safely_removed').write_text('')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.safely_removed'

    def test_paired_extension_detection_requires_both_fastq_mates(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.fastq.gz').write_text('data')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == 'no_extension_found'

    def test_no_extension_found(self, tmp_path):
        """No matching files -> 'no_extension_found'."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == 'no_extension_found'

    def test_ignores_directory_named_as_fastq_file(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.fastq.gz').mkdir()
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == 'no_extension_found'


# ---------------------------------------------------------------------------
# get_mapping_rate (extracts p_pseudoaligned from quant run_info.json)
# ---------------------------------------------------------------------------

class TestGetMappingRate:
    def test_extracts_mapping_rate(self, tmp_path, sample_metadata):
        """Reads p_pseudoaligned from run_info.json into mapping_rate column."""
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        run_info = {'p_pseudoaligned': 85.5}
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps(run_info))
        m = get_mapping_rate(sample_metadata, str(quant_dir))
        assert m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0] == 85.5

    def test_reads_run_info_as_utf8_under_non_utf8_locale(self, tmp_path, sample_metadata, monkeypatch):
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        run_info_path = sra_dir / 'SRR001_run_info.json'
        run_info_path.write_text(
            json.dumps({'p_pseudoaligned': 85.5, 'note': '葉'}, ensure_ascii=False),
            encoding='utf-8',
        )
        original_open = builtins.open
        observed_encodings = []

        def non_utf8_locale_open(file, *args, **kwargs):
            if os.path.realpath(os.fspath(file)) == os.path.realpath(run_info_path):
                observed_encodings.append(kwargs.get('encoding'))
                if 'encoding' not in kwargs:
                    kwargs['encoding'] = 'ascii'
            return original_open(file, *args, **kwargs)

        monkeypatch.setattr(builtins, 'open', non_utf8_locale_open)

        metadata = get_mapping_rate(sample_metadata, str(quant_dir), max_workers=1)

        assert metadata.df.loc[metadata.df['run'] == 'SRR001', 'mapping_rate'].values[0] == 85.5
        assert observed_encodings == ['utf-8']

    def test_avoids_quant_root_listdir_scan(self, tmp_path, sample_metadata, monkeypatch):
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': 75.0}))

        def fail_if_listdir_called(_path):
            raise AssertionError('get_mapping_rate should not scan quant root with os.listdir.')

        monkeypatch.setattr('amalgkit.util.os.listdir', fail_if_listdir_called)

        m = get_mapping_rate(sample_metadata, str(quant_dir))

        assert m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0] == 75.0

    def test_skips_nan_run_ids_when_scanning_quant_subdirs(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', numpy.nan],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': 66.6}))

        m = get_mapping_rate(metadata, str(quant_dir))

        assert m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0] == 66.6

    def test_missing_quant_dir(self, sample_metadata, tmp_path):
        """Missing quant directory does not raise."""
        m = get_mapping_rate(sample_metadata, str(tmp_path / 'nonexistent'))
        assert isinstance(m, Metadata)

    def test_missing_run_info(self, tmp_path, sample_metadata):
        """Missing run_info.json for an SRA -> mapping_rate stays NaN."""
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        # No run_info.json
        m = get_mapping_rate(sample_metadata, str(quant_dir))
        assert numpy.isnan(m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0])

    def test_quant_path_is_file_does_not_crash(self, tmp_path, sample_metadata):
        """quant path that exists as a file should be treated as missing directory."""
        quant_file = tmp_path / 'quant'
        quant_file.write_text('not a directory')
        m = get_mapping_rate(sample_metadata, str(quant_file))
        assert isinstance(m, Metadata)

    def test_raises_when_run_column_missing(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        metadata.df = metadata.df.drop(columns=['run'])

        with pytest.raises(ValueError, match='Column \"run\" is required in metadata to compute mapping_rate'):
            get_mapping_rate(metadata, str(tmp_path / 'quant'))

    def test_duplicate_run_ids_are_all_updated(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR001'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': 42.0}))

        m = get_mapping_rate(metadata, str(quant_dir))

        assert m.df['mapping_rate'].tolist() == [42.0, 42.0]

    def test_metadata_run_whitespace_is_trimmed_for_matching(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [' SRR001 '],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': 33.3}))

        m = get_mapping_rate(metadata, str(quant_dir))

        assert m.df.loc[0, 'mapping_rate'] == 33.3

    def test_invalid_pseudoaligned_value_is_skipped(self, tmp_path, sample_metadata):
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': 'not-a-number'}))

        m = get_mapping_rate(sample_metadata, str(quant_dir))

        assert numpy.isnan(m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0])

    def test_numeric_string_pseudoaligned_value_is_parsed(self, tmp_path, sample_metadata):
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': '12.5'}))

        m = get_mapping_rate(sample_metadata, str(quant_dir))

        assert m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0] == 12.5

    def test_respects_max_workers_override(self, tmp_path, sample_metadata, monkeypatch):
        quant_dir = tmp_path / 'quant'
        for sra_id, value in [('SRR001', 12.5), ('SRR002', 33.3)]:
            sra_dir = quant_dir / sra_id
            sra_dir.mkdir(parents=True)
            (sra_dir / f'{sra_id}_run_info.json').write_text(json.dumps({'p_pseudoaligned': value}))
        observed = {}

        def fake_run_tasks(task_items, task_fn, max_workers):
            observed['max_workers'] = max_workers
            return {}, []

        monkeypatch.setattr('amalgkit.util.run_tasks_with_optional_threads', fake_run_tasks)
        get_mapping_rate(sample_metadata, str(quant_dir), max_workers=1)
        assert observed['max_workers'] == 1


# ---------------------------------------------------------------------------
# get_getfastq_run_dir (creates SRA-specific output directory)
# ---------------------------------------------------------------------------

class TestGetGetfastqRunDir:
    def test_creates_directory(self, tmp_path):
        """Creates getfastq/SRR001 directory and returns path."""
        class Args:
            out_dir = str(tmp_path)
        result = get_getfastq_run_dir(Args(), 'SRR001')
        assert os.path.isdir(result)
        assert result.endswith(os.path.join('getfastq', 'SRR001'))

    def test_existing_directory(self, tmp_path):
        """Returns existing directory without error."""
        class Args:
            out_dir = str(tmp_path)
        gf_dir = tmp_path / 'getfastq' / 'SRR001'
        gf_dir.mkdir(parents=True)
        result = get_getfastq_run_dir(Args(), 'SRR001')
        assert os.path.isdir(result)

    def test_rejects_out_dir_file_path(self, tmp_path):
        out_file = tmp_path / 'out_file'
        out_file.write_text('not a directory')

        class Args:
            out_dir = str(out_file)

        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            get_getfastq_run_dir(Args(), 'SRR001')

    def test_rejects_getfastq_root_file_path(self, tmp_path):
        (tmp_path / 'getfastq').write_text('not a directory')

        class Args:
            out_dir = str(tmp_path)

        with pytest.raises(NotADirectoryError, match='getfastq path exists but is not a directory'):
            get_getfastq_run_dir(Args(), 'SRR001')

    def test_rejects_getfastq_root_symlink(self, tmp_path):
        outside = tmp_path / 'outside'
        outside.mkdir()
        os.symlink(outside, tmp_path / 'getfastq')

        class Args:
            out_dir = str(tmp_path)

        with pytest.raises(NotADirectoryError, match='not a regular directory'):
            get_getfastq_run_dir(Args(), 'SRR001')

    def test_rejects_run_leaf_symlink(self, tmp_path):
        outside = tmp_path / 'outside'
        outside.mkdir()
        getfastq_root = tmp_path / 'getfastq'
        getfastq_root.mkdir()
        os.symlink(outside, getfastq_root / 'SRR001')

        class Args:
            out_dir = str(tmp_path)

        with pytest.raises(NotADirectoryError, match='Run path exists but is not a regular directory'):
            get_getfastq_run_dir(Args(), 'SRR001')

    @pytest.mark.parametrize('run_id', ['../escaped', '/tmp/escaped', 'nested/run', 'nested\\run'])
    def test_rejects_unsafe_run_component(self, tmp_path, run_id):
        class Args:
            out_dir = str(tmp_path)

        with pytest.raises(ValueError, match='Run ID'):
            get_getfastq_run_dir(Args(), run_id)


# ---------------------------------------------------------------------------
# generate_multisp_busco_table (merges BUSCO full_table.tsv files)
# ---------------------------------------------------------------------------

class TestGenerateMultispBuscoTable:
    def test_raises_when_busco_dir_missing(self, tmp_path):
        outfile = tmp_path / 'merged.tsv'
        missing_dir = tmp_path / 'busco_missing'

        with pytest.raises(FileNotFoundError, match='BUSCO directory not found'):
            generate_multisp_busco_table(str(missing_dir), str(outfile))

    def test_raises_when_busco_path_is_file(self, tmp_path):
        outfile = tmp_path / 'merged.tsv'
        busco_file = tmp_path / 'busco.tsv'
        busco_file.write_text('not a directory')

        with pytest.raises(NotADirectoryError, match='BUSCO path exists but is not a directory'):
            generate_multisp_busco_table(str(busco_file), str(outfile))

    def test_raises_when_no_busco_table_files_detected(self, tmp_path):
        outfile = tmp_path / 'merged.tsv'
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()

        with pytest.raises(FileNotFoundError, match='No BUSCO full table file'):
            generate_multisp_busco_table(str(busco_dir), str(outfile))

    def test_merges_two_species(self, tmp_path):
        """Merges BUSCO tables from two species into one output file."""
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        content_a = (
            '# comment line\n'
            'OG0001\tComplete\tgene1\t100\t200\thttp://odb\tgene desc\n'
            'OG0002\tComplete\tgene2\t90\t150\thttp://odb2\tgene desc2\n'
        )
        content_b = (
            '# comment line\n'
            'OG0001\tComplete\tgeneA\t95\t180\t-\t-\n'
            'OG0002\tMissing\t-\t0\t0\t-\t-\n'
        )
        (busco_dir / 'Species_A.tsv').write_text(content_a)
        (busco_dir / 'Species_B.tsv').write_text(content_b)
        outfile = tmp_path / 'merged.tsv'
        generate_multisp_busco_table(str(busco_dir), str(outfile))
        result = pandas.read_csv(str(outfile), sep='\t')
        assert 'Species_A' in result.columns
        assert 'Species_B' in result.columns
        assert result.shape[0] == 2

    def test_raises_when_all_busco_tables_fail_to_parse(self, tmp_path, monkeypatch):
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        (busco_dir / 'Species_A.tsv').write_text('x')
        outfile = tmp_path / 'merged.tsv'

        def fake_run_tasks(task_items, task_fn, max_workers):
            return {}, [(task_items[0], RuntimeError('bad table'))]

        monkeypatch.setattr('amalgkit.util.run_tasks_with_optional_threads', fake_run_tasks)

        with pytest.warns(UserWarning, match='Failed to parse BUSCO table'):
            with pytest.raises(ValueError, match='Failed to parse any BUSCO table'):
                generate_multisp_busco_table(str(busco_dir), str(outfile))

    def test_raises_on_duplicate_species_label_after_filename_parsing(self, tmp_path):
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        content = (
            '# comment line\n'
            'OG0001\tComplete\tgene1\t100\t200\thttp://odb\tgene desc\n'
        )
        (busco_dir / 'Homo_sapiens_strain1.tsv').write_text(content)
        (busco_dir / 'Homo_sapiens_strain2.tsv').write_text(content)
        outfile = tmp_path / 'merged.tsv'

        with pytest.raises(ValueError, match='Duplicate species label was detected across BUSCO tables'):
            generate_multisp_busco_table(str(busco_dir), str(outfile))

    def test_ignores_uncommented_busco_header_row(self, tmp_path):
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        content = (
            'Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n'
            'OG0001\tComplete\tgene1\t100\t200\thttp://odb\tgene desc\n'
        )
        (busco_dir / 'Species_A.tsv').write_text(content)
        (busco_dir / 'Species_B.tsv').write_text(content)
        outfile = tmp_path / 'merged.tsv'

        generate_multisp_busco_table(str(busco_dir), str(outfile))

        result = pandas.read_csv(str(outfile), sep='\t')
        assert result['busco_id'].tolist() == ['OG0001']

    def test_accepts_gzipped_tsv_inputs(self, tmp_path):
        import gzip
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        content = (
            '# comment line\n'
            'OG0001\tComplete\tgene1\t100\t200\thttp://odb\tgene desc\n'
        )
        with gzip.open(busco_dir / 'Species_A.tsv.gz', 'wt') as f:
            f.write(content)
        outfile = tmp_path / 'merged.tsv'

        generate_multisp_busco_table(str(busco_dir), str(outfile))

        result = pandas.read_csv(str(outfile), sep='\t')
        assert 'Species_A' in result.columns
        assert result.shape[0] == 1

    def test_ignores_directory_named_like_busco_table(self, tmp_path):
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        content = (
            '# comment line\n'
            'OG0001\tComplete\tgene1\t100\t200\thttp://odb\tgene desc\n'
        )
        (busco_dir / 'Species_A.tsv').write_text(content)
        (busco_dir / 'Species_B.tsv').mkdir()
        outfile = tmp_path / 'merged.tsv'

        with warnings.catch_warnings(record=True) as records:
            warnings.simplefilter('always')
            generate_multisp_busco_table(str(busco_dir), str(outfile))

        assert len(records) == 0
        result = pandas.read_csv(str(outfile), sep='\t')
        assert 'Species_A' in result.columns
        assert 'Species_B' not in result.columns
