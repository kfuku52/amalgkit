import pytest
import pandas
import gzip
import subprocess
import xml.etree.ElementTree as ET
from contextlib import contextmanager
from io import BytesIO

import os
import urllib.error
from types import SimpleNamespace

from defusedxml.common import EntitiesForbidden

from amalgkit.command_context import GetfastqRuntimeContext
from amalgkit.exceptions import AmalgkitExit
from amalgkit.getfastq import (
    getfastq_search_term,
    getfastq_getxml,
    getfastq_metadata,
    append_file_binary,
    run_mmseqs_rrna_filter,
    run_rrna_filter,
    run_mmseqs_contam_filter,
    run_mmseqs_easy_taxonomy_single_fastq,
    run_mmseqs_easy_search_single_fastq,
    resolve_mmseqs_dbtype,
    update_metadata_after_rrna_filter,
    build_mmseqs_query_input,
    append_mmseqs_sensitivity_option,
    append_mmseqs_positive_int_option,
    append_mmseqs_memory_limit_option,
    filter_fastq_by_core_set,
    iter_synchronized_mmseqs_query_chunks,
    normalize_paired_read_core,
    parse_mmseqs_search_matched_cores,
    parse_fasterq_dump_written_spots,
    parse_fasterq_dump_written_reads,
    infer_written_spots_from_written_reads,
    should_compress_fasterq_output_before_filters,
    normalize_fasterq_size_check,
    get_filter_execution_order,
    resolve_contam_filter_rank,
    resolve_run_target_rank_taxid,
    count_fastq_records_and_bases,
    summarize_layout_fastq_records_and_bases,
    ensure_contam_filter_metadata_rank_taxids,
    resolve_contam_filter_db_path,
    ensure_mmseqs_contam_taxonomy_db_exists,
    count_fastq_records,
    fetch_trace_run_xml_root,
)
from amalgkit.util import Metadata



class TestGetfastqSearchTerm:
    def test_id_only(self):
        assert getfastq_search_term('SRR123456') == 'SRR123456'

    def test_with_additional_term(self):
        result = getfastq_search_term('SRR123456', 'Homo sapiens[Organism]')
        assert result == 'SRR123456 AND Homo sapiens[Organism]'

    def test_none_additional(self):
        result = getfastq_search_term('PRJNA1', None)
        assert result == 'PRJNA1'

class TestNormalizeFasterqSizeCheck:
    @pytest.mark.parametrize('raw, expected', [
        ('on', 'on'),
        ('off', 'off'),
        ('only', 'only'),
        (' yes ', 'on'),
        ('NO', 'off'),
        ('unexpected', 'on'),
        (True, 'on'),
        (False, 'off'),
        (1, 'on'),
        (0, 'off'),
    ])
    def test_normalize(self, raw, expected):
        assert normalize_fasterq_size_check(raw) == expected


class TestShouldCompressFasterqOutputBeforeFilters:
    def test_compresses_when_no_filters(self):
        args = SimpleNamespace(fastp=False, rrna_filter=False)
        assert should_compress_fasterq_output_before_filters(args)

    def test_skips_compression_when_fastp_is_first_filter(self):
        args = SimpleNamespace(fastp=True, rrna_filter=False)
        assert not should_compress_fasterq_output_before_filters(args)

    def test_keeps_compression_when_rrna_is_first_filter(self):
        args = SimpleNamespace(fastp=True, rrna_filter=True, filter_order='rrna,fastp')
        assert not should_compress_fasterq_output_before_filters(args)

    def test_skips_compression_when_contam_filter_is_enabled(self):
        args = SimpleNamespace(fastp=False, rrna_filter=False, contam_filter=True)
        assert not should_compress_fasterq_output_before_filters(args)


class TestGetFilterExecutionOrder:
    def test_defaults_to_fastp_rrna_contam(self):
        args = SimpleNamespace(fastp=True, rrna_filter=True, contam_filter=True)
        assert get_filter_execution_order(args) == ['fastp', 'rrna', 'contam']

    def test_accepts_explicit_three_filter_order(self):
        args = SimpleNamespace(
            fastp=True,
            rrna_filter=True,
            contam_filter=True,
            filter_order='contam,rrna,fastp',
        )
        assert get_filter_execution_order(args) == ['contam', 'rrna', 'fastp']

    def test_accepts_enabled_subset_order(self):
        args = SimpleNamespace(
            fastp=False,
            rrna_filter=True,
            contam_filter=True,
            filter_order='contam,rrna',
        )
        assert get_filter_execution_order(args) == ['contam', 'rrna']

    def test_rejects_legacy_rrna_first_alias(self):
        args = SimpleNamespace(
            fastp=True,
            rrna_filter=True,
            contam_filter=True,
            filter_order='rrna_first',
        )
        with pytest.raises(ValueError, match='supports only'):
            get_filter_execution_order(args)

    def test_rejects_missing_enabled_filter(self):
        args = SimpleNamespace(
            fastp=True,
            rrna_filter=True,
            contam_filter=True,
            filter_order='rrna,contam',
        )
        with pytest.raises(ValueError, match='Missing enabled filter'):
            get_filter_execution_order(args)

    def test_rejects_duplicate_filter_name(self):
        args = SimpleNamespace(
            fastp=True,
            rrna_filter=True,
            contam_filter=False,
            filter_order='fastp,fastp,rrna',
        )
        with pytest.raises(ValueError, match='duplicate filter name'):
            get_filter_execution_order(args)


class TestResolveContamFilterRank:
    def test_defaults_to_superkingdom(self):
        args = SimpleNamespace()
        assert resolve_contam_filter_rank(args) == 'superkingdom'

    def test_normalizes_domain_alias_to_superkingdom(self):
        args = SimpleNamespace(contam_filter_rank='domain')
        assert resolve_contam_filter_rank(args) == 'superkingdom'


class TestResolveRunTargetRankTaxid:
    def test_accepts_taxid_domain_as_superkingdom_alias(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'taxid_domain': [2759],
        }))

        ind_sra, target_rank_taxid = resolve_run_target_rank_taxid(
            metadata=metadata,
            sra_stat={'sra_id': 'SRR001'},
            rank_name='superkingdom',
        )

        assert ind_sra == 0
        assert target_rank_taxid == 2759


class TestFasterqDumpWrittenCountParsing:
    def test_parse_written_spots(self):
        stderr_txt = '\n'.join([
            'spots read      : 10',
            'reads written   : 20',
            'spots written   : 10',
        ])
        assert parse_fasterq_dump_written_spots('', stderr_txt) == 10

    def test_parse_written_reads(self):
        stderr_txt = '\n'.join([
            'spots read      : 15,890,071',
            'reads written   : 31,780,142',
        ])
        assert parse_fasterq_dump_written_reads('', stderr_txt) == 31780142

    def test_infer_written_spots_single_layout(self):
        from amalgkit.getfastq import RunFileState
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        file_state = RunFileState(work_dir='/', files=set())
        assert infer_written_spots_from_written_reads(sra_stat, 123, file_state) == 123

    def test_infer_written_spots_paired_layout_without_singletons(self):
        from amalgkit.getfastq import RunFileState
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        file_state = RunFileState(work_dir='/', files={'SRR001_1.fastq', 'SRR001_2.fastq'})
        assert infer_written_spots_from_written_reads(sra_stat, 200, file_state) == 100

    def test_infer_written_spots_paired_layout_with_singletons_returns_none(self):
        from amalgkit.getfastq import RunFileState
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        file_state = RunFileState(work_dir='/', files={'SRR001.fastq', 'SRR001_1.fastq', 'SRR001_2.fastq'})
        assert infer_written_spots_from_written_reads(sra_stat, 201, file_state) is None


class TestCountFastqRecords:
    @staticmethod
    def _write_fastq(path, reads):
        with open(path, 'wt') as out:
            for i, seq in enumerate(reads):
                out.write('@r{}\n'.format(i))
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_counts_records_plain_and_gz(self, tmp_path):
        plain = tmp_path / 'x.fastq'
        self._write_fastq(str(plain), ['AAAA', 'CCCC', 'GGGG'])
        gz = tmp_path / 'x.fastq.gz'
        with open(str(plain), 'rb') as fin, gzip.open(str(gz), 'wb') as fout:
            fout.write(fin.read())

        assert count_fastq_records(str(plain)) == 3
        assert count_fastq_records(str(gz)) == 3

    def test_counts_records_and_bases_plain_and_gz(self, tmp_path):
        plain = tmp_path / 'x.fastq'
        self._write_fastq(str(plain), ['AAAA', 'CC', 'GGGGG'])
        gz = tmp_path / 'x.fastq.gz'
        with open(str(plain), 'rb') as fin, gzip.open(str(gz), 'wb') as fout:
            fout.write(fin.read())

        num_plain, bp_plain = count_fastq_records_and_bases(str(plain))
        num_gz, bp_gz = count_fastq_records_and_bases(str(gz))
        assert num_plain == 3
        assert bp_plain == (4 + 2 + 5)
        assert num_gz == 3
        assert bp_gz == (4 + 2 + 5)

    def test_counts_valid_record_without_terminal_newline(self, tmp_path):
        payload = b'@r0\nAAAA\n+\nIIII'
        plain = tmp_path / 'no-final-newline.fastq'
        plain.write_bytes(payload)
        gz = tmp_path / 'no-final-newline.fastq.gz'
        with gzip.open(gz, 'wb') as handle:
            handle.write(payload)

        assert count_fastq_records(str(plain)) == 1
        assert count_fastq_records(str(gz)) == 1
        assert count_fastq_records_and_bases(str(plain)) == (1, 4)
        assert count_fastq_records_and_bases(str(gz)) == (1, 4)

    def test_warns_on_truncated_fastq(self, tmp_path, capsys):
        path = tmp_path / 'bad.fastq'
        with open(path, 'wt') as out:
            out.write('@r0\nAAAA\n+\nIIII\n')
            out.write('@r1\nCCCC\n+\n')  # truncated

        assert count_fastq_records(str(path)) == 1
        assert 'not divisible by 4' in capsys.readouterr().err


class TestAppendFileBinary:
    def test_appends_source_to_destination(self, tmp_path):
        src = tmp_path / 'src.bin'
        dst = tmp_path / 'dst.bin'
        src.write_bytes(b'BBBB')
        dst.write_bytes(b'AAAA')

        append_file_binary(str(src), str(dst))

        assert dst.read_bytes() == b'AAAABBBB'

    def test_falls_back_to_copyfileobj_when_sendfile_raises(self, tmp_path, monkeypatch):
        src = tmp_path / 'src.bin'
        dst = tmp_path / 'dst.bin'
        src.write_bytes(b'12345')
        dst.write_bytes(b'abc')

        if not hasattr(os, 'sendfile'):
            append_file_binary(str(src), str(dst))
            assert dst.read_bytes() == b'abc12345'
            return

        def fake_sendfile(*_args, **_kwargs):
            raise OSError('x')

        monkeypatch.setattr('amalgkit.getfastq.os.sendfile', fake_sendfile)
        append_file_binary(str(src), str(dst))
        assert dst.read_bytes() == b'abc12345'


class TestFastqStatsWithSeqkit:
    @staticmethod
    def _write_fastq(path, reads):
        with open(path, 'wt') as out:
            for i, seq in enumerate(reads):
                out.write('@r{}\n'.format(i))
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_summarize_layout_uses_seqkit_stats_when_args_provided(self, tmp_path, monkeypatch):
        pair1 = tmp_path / 'x_1.fastq'
        pair2 = tmp_path / 'x_2.fastq'
        self._write_fastq(str(pair1), ['AAAA', 'CC'])
        self._write_fastq(str(pair2), ['TT', 'GG'])

        class Args:
            threads = 2
            dump_print = False
            seqkit_exe = 'seqkit'

        def fake_seqkit_run(cmd, stdout=None, stderr=None):
            assert cmd[0] == 'seqkit'
            assert cmd[1] == 'stats'
            rows = [
                'file\tformat\ttype\tnum_seqs\tsum_len\tavg_len',
                '{}\tFASTQ\tDNA\t2\t6\t3'.format(os.path.abspath(str(pair1))),
                '{}\tFASTQ\tDNA\t2\t4\t2'.format(os.path.abspath(str(pair2))),
            ]
            return subprocess.CompletedProcess(cmd, 0, stdout=('\n'.join(rows) + '\n').encode('utf8'), stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_seqkit_run)

        num_spots, bp_total = summarize_layout_fastq_records_and_bases(
            layout='paired',
            pair1_path=str(pair1),
            pair2_path=str(pair2),
            args=Args(),
        )
        assert num_spots == 2
        assert bp_total == 10

    def test_summarize_layout_falls_back_to_python_parser_when_seqkit_fails(self, tmp_path, monkeypatch):
        pair1 = tmp_path / 'x_1.fastq'
        pair2 = tmp_path / 'x_2.fastq'
        self._write_fastq(str(pair1), ['AAAA', 'CC'])
        self._write_fastq(str(pair2), ['TT', 'GGG'])

        class Args:
            threads = 2
            dump_print = False
            seqkit_exe = 'seqkit'

        monkeypatch.setattr(
            'amalgkit.getfastq.subprocess.run',
            lambda cmd, stdout=None, stderr=None: subprocess.CompletedProcess(cmd, 1, stdout=b'', stderr=b'failed'),
        )

        num_spots, bp_total = summarize_layout_fastq_records_and_bases(
            layout='paired',
            pair1_path=str(pair1),
            pair2_path=str(pair2),
            args=Args(),
        )
        assert num_spots == 2
        assert bp_total == (4 + 2 + 2 + 3)


class TestRrnaMetadataUpdate:
    def test_update_metadata_uses_known_input_counts_without_rescanning_input(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'num_rrna_in': [0],
            'num_rrna_out': [0],
            'bp_rrna_in': [0],
            'bp_rrna_out': [0],
        }))
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'metadata_idx': 0,
        }
        observed = {'calls': []}

        def fake_summarize(layout, single_path=None, pair1_path=None, pair2_path=None, args=None):
            observed['calls'].append((layout, single_path, pair1_path, pair2_path))
            return 7, 70

        monkeypatch.setattr('amalgkit.getfastq.summarize_layout_fastq_records_and_bases', fake_summarize)
        metadata = update_metadata_after_rrna_filter(
            metadata=metadata,
            sra_stat=sra_stat,
            input_paths={'': '/tmp/in.fastq.gz'},
            output_paths={'': '/tmp/out.fastq.gz'},
            args=None,
            known_input_counts={'num_spots': 11, 'bp_total': 110},
        )
        assert metadata.df.loc[0, 'num_rrna_in'] == 11
        assert metadata.df.loc[0, 'bp_rrna_in'] == 110
        assert metadata.df.loc[0, 'num_rrna_out'] == 7
        assert metadata.df.loc[0, 'bp_rrna_out'] == 70
        assert len(observed['calls']) == 1
        assert observed['calls'][0][1] == '/tmp/out.fastq.gz'

    def test_update_metadata_uses_known_input_and_output_counts_without_fastq_scan(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'num_rrna_in': [0],
            'num_rrna_out': [0],
            'bp_rrna_in': [0],
            'bp_rrna_out': [0],
        }))
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'metadata_idx': 0,
        }

        def fail_summarize(*_args, **_kwargs):
            raise AssertionError('FASTQ scan should not run when known_input_counts and known_output_counts are provided.')

        monkeypatch.setattr('amalgkit.getfastq.summarize_layout_fastq_records_and_bases', fail_summarize)
        metadata = update_metadata_after_rrna_filter(
            metadata=metadata,
            sra_stat=sra_stat,
            input_paths={'': '/tmp/in.fastq.gz'},
            output_paths={'': '/tmp/out.fastq.gz'},
            args=None,
            known_input_counts={'num_spots': 11, 'bp_total': 110},
            known_output_counts={'num_spots': 7, 'bp_total': 70},
        )
        assert metadata.df.loc[0, 'num_rrna_in'] == 11
        assert metadata.df.loc[0, 'bp_rrna_in'] == 110
        assert metadata.df.loc[0, 'num_rrna_out'] == 7
        assert metadata.df.loc[0, 'bp_rrna_out'] == 70


class TestRrnaReadIdHandling:
    def test_preserves_sra_accession_spot_suffixes(self):
        assert normalize_paired_read_core('SRR000001.1') == 'SRR000001.1'
        assert normalize_paired_read_core('SRR000001.2') == 'SRR000001.2'
        assert normalize_paired_read_core('SRR000001.1/1') == 'SRR000001.1'
        assert normalize_paired_read_core('SRR000001.1/2') == 'SRR000001.1'

    def test_one_spot_hit_does_not_remove_the_next_spot(self, tmp_path):
        result_tsv = tmp_path / 'result.tsv'
        result_tsv.write_text('SRR000001.1\ttarget\n')
        input_fastq = tmp_path / 'input.fastq.gz'
        output_fastq = tmp_path / 'output.fastq.gz'
        with gzip.open(input_fastq, 'wt') as handle:
            handle.write('@SRR000001.1\nAAAA\n+\nIIII\n')
            handle.write('@SRR000001.2\nCCCC\n+\nIIII\n')

        remove_cores = parse_mmseqs_search_matched_cores(str(result_tsv))
        counts = filter_fastq_by_core_set(
            input_path=str(input_fastq),
            output_path=str(output_fastq),
            remove_cores=remove_cores,
        )

        assert remove_cores == {'SRR000001.1'}
        assert counts[:2] == (2, 1)
        with gzip.open(output_fastq, 'rt') as handle:
            assert handle.readline().strip() == '@SRR000001.2'

    def test_rejects_paired_mate_order_mismatch(self, tmp_path):
        in1 = tmp_path / 'r1.fastq.gz'
        in2 = tmp_path / 'r2.fastq.gz'
        query_root = tmp_path / 'queries'
        query_root.mkdir()
        with gzip.open(in1, 'wt') as handle:
            handle.write('@A/1\nAAAA\n+\nIIII\n@B/1\nCCCC\n+\nIIII\n')
        with gzip.open(in2, 'wt') as handle:
            handle.write('@A/2\nTTTT\n+\nIIII\n@C/2\nGGGG\n+\nIIII\n')

        with pytest.raises(ValueError, match='mate IDs are out of sync'):
            list(iter_synchronized_mmseqs_query_chunks(
                input_path_by_suffix={'_1': str(in1), '_2': str(in2)},
                query_root=str(query_root),
                query_tag='SRR001',
                chunk_spots=10,
            ))

    def test_missing_mmseqs_result_is_not_treated_as_no_hits(self, tmp_path):
        with pytest.raises(FileNotFoundError, match='was not generated'):
            parse_mmseqs_search_matched_cores(str(tmp_path / 'missing.tsv'))


class TestRunMmseqsRrnaFilter:
    @staticmethod
    def _write_fastq_gz(path, read_id_seq_pairs):
        with gzip.open(path, 'wt') as out:
            for read_id, seq in read_id_seq_pairs:
                out.write('@{}\n'.format(read_id))
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_paired_strict_pair_removal(self, tmp_path, monkeypatch, capsys):
        sra_id = 'SRR001'
        in1 = tmp_path / '{}_1.fastq.gz'.format(sra_id)
        in2 = tmp_path / '{}_2.fastq.gz'.format(sra_id)
        self._write_fastq_gz(str(in1), [('rA/1', 'AAAA'), ('rB/1', 'CCCC'), ('rC/1', 'GGGG')])
        self._write_fastq_gz(str(in2), [('rA/2', 'TTTT'), ('rB/2', 'GGGG'), ('rC/2', 'AAAA')])

        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'num_rrna_in': [0],
            'num_rrna_out': [0],
            'bp_rrna_in': [0],
            'bp_rrna_out': [0],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'metadata_idx': 0,
            'current_ext': '.fastq.gz',
            'getfastq_sra_dir': str(tmp_path),
        }

        args = SimpleNamespace(
            rrna_filter=True,
            mmseqs_exe='mmseqs',
            remove_tmp=False,
            dump_print=False,
            threads=1,
            out_dir=str(tmp_path),
            download_dir='inferred',
            download_lock_dir='inferred',
            rrna_filter_jobs=1,
            rrna_filter_chunk_spots=2,
            rrna_filter_memory_limit='32G',
        )

        monkeypatch.setattr('amalgkit.getfastq.ensure_mmseqs_rrna_reference_db_exists', lambda _args: '/tmp/mock_rrna_db')

        observed = {'input_paths': [], 'semaphore_timeout_seconds': 'unset'}

        @contextmanager
        def capture_rrna_semaphore(**kwargs):
            observed['semaphore_timeout_seconds'] = kwargs.get('timeout_seconds')
            yield None

        monkeypatch.setattr(
            'amalgkit.getfastq.maybe_acquire_download_semaphore',
            capture_rrna_semaphore,
        )

        def fake_run_mmseqs_search(args, input_path, target_db, result_tsv, tmp_dir):
            _ = (args, target_db, tmp_dir)
            observed['input_paths'].append(input_path)
            os.makedirs(os.path.dirname(result_tsv), exist_ok=True)
            lines = []
            with gzip.open(input_path, 'rt') as fin:
                while True:
                    header = fin.readline()
                    if header == '':
                        break
                    seq = fin.readline()
                    plus = fin.readline()
                    qual = fin.readline()
                    _ = (seq, plus, qual)
                    read_id = header.strip().lstrip('@')
                    if read_id.startswith('rB/'):
                        lines.append(read_id)
            with open(result_tsv, 'wt') as fout:
                fout.write('\n'.join(lines) + '\n')
            return result_tsv

        monkeypatch.setattr('amalgkit.getfastq.run_mmseqs_easy_search_single_fastq', fake_run_mmseqs_search)

        metadata, run_file_state = run_mmseqs_rrna_filter(
            sra_stat=sra_stat,
            args=args,
            output_dir=str(tmp_path),
            metadata=metadata,
            files={'{}_1.fastq.gz'.format(sra_id), '{}_2.fastq.gz'.format(sra_id)},
            return_file_state=True,
        )

        out1 = tmp_path / '{}_1.rrna-filtered.fastq.gz'.format(sra_id)
        out2 = tmp_path / '{}_2.rrna-filtered.fastq.gz'.format(sra_id)
        assert out1.exists()
        assert out2.exists()
        assert count_fastq_records(str(out1)) == 2
        assert count_fastq_records(str(out2)) == 2
        assert metadata.df.loc[0, 'num_rrna_in'] == 3
        assert metadata.df.loc[0, 'num_rrna_out'] == 2
        assert metadata.df.loc[0, 'bp_rrna_in'] == 24
        assert metadata.df.loc[0, 'bp_rrna_out'] == 16
        assert metadata.df.loc[0, 'sec_rrna_search'] >= 0.0
        assert metadata.df.loc[0, 'sec_rrna_rewrite'] >= 0.0
        assert metadata.df.loc[0, 'sec_rrna_filter'] >= 0.0
        assert len(observed['input_paths']) == 2
        assert all(path.endswith('.fastq.gz') for path in observed['input_paths'])
        work_roots = {
            os.path.relpath(path, str(tmp_path)).split(os.sep, 1)[0]
            for path in observed['input_paths']
        }
        assert len(work_roots) == 1
        assert next(iter(work_roots)).startswith('mmseqs_rrna_work.')
        assert not (tmp_path / 'mmseqs_rrna_work').exists()
        assert observed['semaphore_timeout_seconds'] is None
        assert run_file_state.has('{}_1.rrna-filtered.fastq.gz'.format(sra_id))
        assert run_file_state.has('{}_2.rrna-filtered.fastq.gz'.format(sra_id))
        out = capsys.readouterr().out
        assert 'Time elapsed for MMseqs rRNA search ({}):'.format(sra_id) in out
        assert 'Time elapsed for rRNA FASTQ rewrite ({}):'.format(sra_id) in out
        assert 'Time elapsed for MMseqs rRNA filter ({}):'.format(sra_id) in out


class TestRrnaFilterDispatch:
    def test_dispatches_to_mmseqs(self, monkeypatch):
        args = SimpleNamespace(rrna_filter=True)
        expected = ('metadata_obj', 'file_state_obj')

        monkeypatch.setattr(
            'amalgkit.getfastq.run_mmseqs_rrna_filter',
            lambda **kwargs: expected,
        )

        observed = run_rrna_filter(
            sra_stat={'sra_id': 'SRR001', 'layout': 'single'},
            args=args,
            output_dir='/tmp',
            metadata='metadata_obj',
            return_file_state=True,
        )
        assert observed == expected


class TestRunMmseqsContamFilter:
    @staticmethod
    def _write_fastq_gz(path, read_id_seq_pairs):
        with gzip.open(path, 'wt') as out:
            for read_id, seq in read_id_seq_pairs:
                out.write('@{}\n'.format(read_id))
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_paired_strict_pair_removal_and_unclassified_keep(self, tmp_path, monkeypatch, capsys):
        sra_id = 'SRR001'
        in1 = tmp_path / '{}_1.fastq.gz'.format(sra_id)
        in2 = tmp_path / '{}_2.fastq.gz'.format(sra_id)
        self._write_fastq_gz(str(in1), [('rA/1', 'AAAA'), ('rB/1', 'CCCC'), ('rC/1', 'GGGG')])
        self._write_fastq_gz(str(in2), [('rA/2', 'TTTT'), ('rB/2', 'GGGG'), ('rC/2', 'AAAA')])

        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'taxid_phylum': [500],
            'num_contam_in': [0],
            'num_contam_out': [0],
            'bp_contam_in': [0],
            'bp_contam_out': [0],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'metadata_idx': 0,
            'current_ext': '.fastq.gz',
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            contam_filter = True
            contam_filter_rank = 'phylum'
            contam_filter_db = 'inferred'
            contam_filter_db_name = 'UniRef90'
            mmseqs_exe = 'mmseqs'
            remove_tmp = False
            dump_print = False
            threads = 1
            out_dir = str(tmp_path)
            download_dir = 'inferred'

        monkeypatch.setattr('amalgkit.getfastq.ensure_mmseqs_contam_taxonomy_db_exists', lambda _args: '/tmp/mockdb')

        observed = {'input_paths': []}

        def fake_run_mmseqs(args, input_path, target_db, result_prefix, tmp_dir, runtime_context=None):
            _ = (args, target_db, tmp_dir, runtime_context)
            observed['input_paths'].append(input_path)
            os.makedirs(os.path.dirname(result_prefix), exist_ok=True)
            lca_path = result_prefix + '_lca.tsv'
            lines = []
            with gzip.open(input_path, 'rt') as fin:
                while True:
                    header = fin.readline()
                    if header == '':
                        break
                    seq = fin.readline()
                    plus = fin.readline()
                    qual = fin.readline()
                    _ = (seq, plus, qual)
                    read_id = header.strip().lstrip('@')
                    if read_id == 'rA/1' or read_id == 'rA/2':
                        lines.append('{}\t11\tphylum\tmatch'.format(read_id))
                    elif read_id == 'rB/1':
                        lines.append('{}\t22\tphylum\tmismatch'.format(read_id))
                    elif read_id == 'rB/2':
                        lines.append('{}\t11\tphylum\tmatch'.format(read_id))
                    elif read_id == 'rC/1' or read_id == 'rC/2':
                        lines.append('{}\t0\tno rank\tunclassified'.format(read_id))
            with open(lca_path, 'wt') as fout:
                fout.write('\n'.join(lines) + '\n')
            return lca_path

        # Map assigned taxid to phylum-level taxid.
        monkeypatch.setattr(
            'amalgkit.getfastq.resolve_taxid_at_rank',
            lambda taxid, rank_name, ncbi, rank_cache: 500 if int(taxid) == 11 else (600 if int(taxid) == 22 else None),
        )
        monkeypatch.setattr('amalgkit.getfastq.run_mmseqs_easy_taxonomy_single_fastq', fake_run_mmseqs)
        monkeypatch.setattr('amalgkit.getfastq.get_ncbi_taxonomy', lambda args=None: object())
        perf_counter_values = iter([20.0, 21.0, 23.5, 24.5])
        monkeypatch.setattr('amalgkit.getfastq.time.perf_counter', lambda: next(perf_counter_values))

        metadata, run_file_state = run_mmseqs_contam_filter(
            sra_stat=sra_stat,
            args=Args(),
            output_dir=str(tmp_path),
            metadata=metadata,
            files={'{}_1.fastq.gz'.format(sra_id), '{}_2.fastq.gz'.format(sra_id)},
            return_file_state=True,
        )

        out1 = tmp_path / '{}_1.contam-filtered.fastq.gz'.format(sra_id)
        out2 = tmp_path / '{}_2.contam-filtered.fastq.gz'.format(sra_id)
        assert out1.exists()
        assert out2.exists()
        assert count_fastq_records(str(out1)) == 2
        assert count_fastq_records(str(out2)) == 2
        assert metadata.df.loc[0, 'num_contam_in'] == 3
        assert metadata.df.loc[0, 'num_contam_out'] == 2
        assert metadata.df.loc[0, 'bp_contam_in'] == 24
        assert metadata.df.loc[0, 'bp_contam_out'] == 16
        assert metadata.df.loc[0, 'sec_ete_taxonomy'] == pytest.approx(2.5)
        assert metadata.df.loc[0, 'sec_contam_filter'] == pytest.approx(4.5)
        assert len(observed['input_paths']) == 1
        assert observed['input_paths'][0].endswith('_combined.fastq.gz')
        assert run_file_state.has('{}_1.contam-filtered.fastq.gz'.format(sra_id))
        assert run_file_state.has('{}_2.contam-filtered.fastq.gz'.format(sra_id))
        out = capsys.readouterr().out
        assert 'Time elapsed for NCBI taxonomy initialization ({}):'.format(sra_id) in out
        assert 'Time elapsed for contaminant filter ({}):'.format(sra_id) in out


class TestBuildMmseqsQueryInput:
    def test_concatenates_paired_gzip_queries_once(self, tmp_path):
        query_root = tmp_path / 'mmseqs'
        query_root.mkdir()
        in1 = tmp_path / 'SRR001_1.fastq.gz'
        in2 = tmp_path / 'SRR001_2.fastq.gz'
        with gzip.open(in1, 'wt') as out:
            out.write('@rA/1\nAAAA\n+\nIIII\n')
        with gzip.open(in2, 'wt') as out:
            out.write('@rA/2\nTTTT\n+\nIIII\n')

        combined_path = build_mmseqs_query_input(
            query_root=str(query_root),
            query_tag='SRR001',
            input_path_by_suffix={'_1': str(in1), '_2': str(in2)},
        )

        assert combined_path.endswith('_combined.fastq.gz')
        with gzip.open(combined_path, 'rt') as fin:
            assert fin.read() == '@rA/1\nAAAA\n+\nIIII\n@rA/2\nTTTT\n+\nIIII\n'


class TestEnsureContamFilterMetadataRankTaxids:
    def test_accepts_existing_taxid_domain_for_superkingdom(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['NuKS2-prey'],
            'scientific_name': ['Arabidopsis thaliana'],
            'taxid': [3702],
            'taxid_domain': [2759],
            'taxid_kingdom': [33090],
            'taxid_phylum': [35493],
            'taxid_class': [3398],
            'taxid_order': [3699],
            'taxid_family': [3700],
            'taxid_genus': [3701],
            'taxid_species': [3702],
        }))
        metadata.df['taxid'] = metadata.df['taxid'].astype('Int64')
        for col in [col for col in metadata.df.columns if col.startswith('taxid_') and col != 'taxid']:
            metadata.df[col] = pandas.to_numeric(metadata.df[col], errors='coerce').astype('Int64')

        captured = {'called': 0}

        def fake_add_standard_rank_taxids(self, args=None):
            _ = args
            captured['called'] += 1

        monkeypatch.setattr(Metadata, 'add_standard_rank_taxids', fake_add_standard_rank_taxids)

        class Args:
            contam_filter = True
            contam_filter_rank = 'superkingdom'

        out = ensure_contam_filter_metadata_rank_taxids(metadata, Args())

        assert captured['called'] == 0
        assert out.df.loc[0, 'taxid_domain'] == 2759

    def test_rebuilds_required_rank_when_column_exists_but_value_is_missing(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['NuKS2-prey'],
            'scientific_name': ['Arabidopsis thaliana'],
            'taxid': [3702],
            'taxid_domain': [2759],
            'taxid_kingdom': [33090],
            'taxid_phylum': [pandas.NA],
            'taxid_class': [pandas.NA],
            'taxid_order': [pandas.NA],
            'taxid_family': [pandas.NA],
            'taxid_genus': [pandas.NA],
            'taxid_species': [pandas.NA],
        }))
        metadata.df['taxid'] = metadata.df['taxid'].astype('Int64')
        for col in [col for col in metadata.df.columns if col.startswith('taxid_') and col != 'taxid']:
            metadata.df[col] = pandas.to_numeric(metadata.df[col], errors='coerce').astype('Int64')

        captured = {'called': 0}

        def fake_add_standard_rank_taxids(self, args=None):
            _ = args
            captured['called'] += 1
            self.df['taxid_domain'] = pandas.Series([2759], dtype='Int64')
            self.df['taxid_kingdom'] = pandas.Series([33090], dtype='Int64')
            self.df['taxid_phylum'] = pandas.Series([35493], dtype='Int64')
            self.df['taxid_class'] = pandas.Series([3398], dtype='Int64')
            self.df['taxid_order'] = pandas.Series([3699], dtype='Int64')
            self.df['taxid_family'] = pandas.Series([3700], dtype='Int64')
            self.df['taxid_genus'] = pandas.Series([3701], dtype='Int64')
            self.df['taxid_species'] = pandas.Series([3702], dtype='Int64')

        monkeypatch.setattr(Metadata, 'add_standard_rank_taxids', fake_add_standard_rank_taxids)

        class Args:
            contam_filter = True
            contam_filter_rank = 'phylum'

        out = ensure_contam_filter_metadata_rank_taxids(metadata, Args())

        assert captured['called'] == 1
        assert out.df.loc[0, 'taxid_phylum'] == 35493

    def test_reports_unresolved_required_rank_after_rebuild(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR-unresolved'],
            'scientific_name': ['Unknown species'],
            'taxid': [999999999],
            'taxid_phylum': [pandas.NA],
        }))

        def fake_add_standard_rank_taxids(self, args=None):
            _ = args
            self.df['taxid_phylum'] = pandas.Series([pandas.NA], dtype='Int64')

        monkeypatch.setattr(Metadata, 'add_standard_rank_taxids', fake_add_standard_rank_taxids)

        class Args:
            contam_filter = True
            contam_filter_rank = 'phylum'

        with pytest.raises(
            ValueError,
            match=r'Missing taxid_phylum in metadata for 1 run\(s\): SRR-unresolved',
        ):
            ensure_contam_filter_metadata_rank_taxids(metadata, Args())


class TestRunMmseqsEasyTaxonomy:
    def test_adds_search_type_for_nucleotide_db(self, monkeypatch):
        class Args:
            mmseqs_exe = 'mmseqs'
            threads = 2
            dump_print = False
            contam_filter_sensitivity = 'auto'

        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.resolve_mmseqs_easy_taxonomy_search_type', lambda **_kwargs: '3')
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        out_path = run_mmseqs_easy_taxonomy_single_fastq(
            args=Args(),
            input_path='/tmp/in.fastq.gz',
            target_db='/tmp/db',
            result_prefix='/tmp/out/result',
            tmp_dir='/tmp/out/tmp',
        )

        assert out_path == '/tmp/out/result_lca.tsv'
        assert '--search-type' in observed['cmd']
        assert observed['cmd'][observed['cmd'].index('--search-type') + 1] == '3'
        assert '-s' not in observed['cmd']
        assert '--report-mode' in observed['cmd']
        assert observed['cmd'][observed['cmd'].index('--report-mode') + 1] == '2'

    def test_omits_search_type_for_aminoacid_db(self, monkeypatch):
        class Args:
            mmseqs_exe = 'mmseqs'
            threads = 2
            dump_print = False
            contam_filter_sensitivity = 2.5
            contam_filter_max_seqs = 'auto'

        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.resolve_mmseqs_easy_taxonomy_search_type', lambda **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_mmseqs_easy_taxonomy_single_fastq(
            args=Args(),
            input_path='/tmp/in.fastq.gz',
            target_db='/tmp/db',
            result_prefix='/tmp/out/result',
            tmp_dir='/tmp/out/tmp',
        )
        assert '--search-type' not in observed['cmd']
        assert '-s' in observed['cmd']
        assert observed['cmd'][observed['cmd'].index('-s') + 1] == '2.5'
        assert '--report-mode' in observed['cmd']
        assert observed['cmd'][observed['cmd'].index('--report-mode') + 1] == '2'

    def test_appends_contam_max_seqs_when_configured(self, monkeypatch):
        class Args:
            mmseqs_exe = 'mmseqs'
            threads = 2
            dump_print = False
            contam_filter_sensitivity = 'auto'
            contam_filter_max_seqs = 20

        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.resolve_mmseqs_easy_taxonomy_search_type', lambda **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_mmseqs_easy_taxonomy_single_fastq(
            args=Args(),
            input_path='/tmp/in.fastq.gz',
            target_db='/tmp/db',
            result_prefix='/tmp/out/result',
            tmp_dir='/tmp/out/tmp',
        )

        assert '--max-seqs' in observed['cmd']
        assert observed['cmd'][observed['cmd'].index('--max-seqs') + 1] == '20'


class TestRunMmseqsEasySearch:
    def test_appends_rrna_sensitivity_when_configured(self, monkeypatch):
        args = SimpleNamespace(
            mmseqs_exe='mmseqs',
            threads=2,
            dump_print=False,
            rrna_filter_sensitivity=1.5,
            rrna_filter_max_seqs='auto',
        )
        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            with open(cmd[4], 'wt') as fout:
                fout.write('')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr(
            'amalgkit.getfastq.ensure_mmseqs_rrna_search_index_exists',
            lambda **kwargs: kwargs['db_path'],
        )
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_mmseqs_easy_search_single_fastq(
            args=args,
            input_path='/tmp/in.fastq.gz',
            target_db='/tmp/db',
            result_tsv='/tmp/out.tsv',
            tmp_dir='/tmp/tmp',
        )

        assert observed['cmd'][0:2] == ['mmseqs', 'easy-search']
        assert '-s' in observed['cmd']
        assert observed['cmd'][observed['cmd'].index('-s') + 1] == '1.5'
        assert '--max-accept' in observed['cmd']
        assert observed['cmd'][observed['cmd'].index('--max-accept') + 1] == '1'
        assert observed['cmd'][observed['cmd'].index('--split-memory-limit') + 1] == '32G'
        assert observed['cmd'][observed['cmd'].index('--db-load-mode') + 1] == '2'

    def test_omits_rrna_sensitivity_when_auto(self, monkeypatch):
        args = SimpleNamespace(
            mmseqs_exe='mmseqs',
            threads=2,
            dump_print=False,
            rrna_filter_sensitivity='auto',
            rrna_filter_max_seqs='auto',
        )
        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            with open(cmd[4], 'wt') as fout:
                fout.write('')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr(
            'amalgkit.getfastq.ensure_mmseqs_rrna_search_index_exists',
            lambda **kwargs: kwargs['db_path'],
        )
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_mmseqs_easy_search_single_fastq(
            args=args,
            input_path='/tmp/in.fastq.gz',
            target_db='/tmp/db',
            result_tsv='/tmp/out.tsv',
            tmp_dir='/tmp/tmp',
        )

        assert '-s' not in observed['cmd']
        assert '--max-accept' in observed['cmd']
        assert observed['cmd'][observed['cmd'].index('--max-accept') + 1] == '1'

    def test_appends_rrna_max_seqs_when_configured(self, monkeypatch):
        args = SimpleNamespace(
            mmseqs_exe='mmseqs',
            threads=2,
            dump_print=False,
            rrna_filter_sensitivity='auto',
            rrna_filter_max_seqs=20,
        )
        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            with open(cmd[4], 'wt') as fout:
                fout.write('')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr(
            'amalgkit.getfastq.ensure_mmseqs_rrna_search_index_exists',
            lambda **kwargs: kwargs['db_path'],
        )
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_mmseqs_easy_search_single_fastq(
            args=args,
            input_path='/tmp/in.fastq.gz',
            target_db='/tmp/db',
            result_tsv='/tmp/out.tsv',
            tmp_dir='/tmp/tmp',
        )

        assert '--max-seqs' in observed['cmd']
        assert observed['cmd'][observed['cmd'].index('--max-seqs') + 1] == '20'

    def test_rejects_fatal_child_failure_even_with_zero_exit_code(self, monkeypatch):
        args = SimpleNamespace(
            mmseqs_exe='mmseqs',
            threads=1,
            dump_print=False,
        )

        monkeypatch.setattr(
            'amalgkit.getfastq.ensure_mmseqs_rrna_search_index_exists',
            lambda **kwargs: kwargs['db_path'],
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.subprocess.run',
            lambda cmd, stdout=None, stderr=None: subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'',
                stderr=b'Error: search died\n',
            ),
        )

        with pytest.raises(RuntimeError, match='fatal child-process failure'):
            run_mmseqs_easy_search_single_fastq(
                args=args,
                input_path='/tmp/in.fastq.gz',
                target_db='/tmp/db',
                result_tsv='/tmp/out.tsv',
                tmp_dir='/tmp/tmp',
            )

    def test_appends_rrna_split_memory_limit(self, monkeypatch):
        args = SimpleNamespace(
            mmseqs_exe='mmseqs',
            threads=2,
            dump_print=False,
            rrna_filter_sensitivity='auto',
            rrna_filter_max_seqs='auto',
            rrna_filter_memory_limit='32G',
        )
        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            with open(cmd[4], 'wt') as fout:
                fout.write('')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr(
            'amalgkit.getfastq.ensure_mmseqs_rrna_search_index_exists',
            lambda **kwargs: kwargs['db_path'],
        )
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_mmseqs_easy_search_single_fastq(
            args=args,
            input_path='/tmp/in.fastq.gz',
            target_db='/tmp/db',
            result_tsv='/tmp/out.tsv',
            tmp_dir='/tmp/tmp',
        )

        assert observed['cmd'][observed['cmd'].index('--split-memory-limit') + 1] == '32G'


class TestAppendMmseqsSensitivityOption:
    def test_leaves_command_unchanged_for_auto(self):
        command = ['mmseqs', 'easy-search']
        out = append_mmseqs_sensitivity_option(command[:], 'auto')
        assert out == ['mmseqs', 'easy-search']

    def test_appends_float_for_numeric_value(self):
        command = ['mmseqs', 'easy-search']
        out = append_mmseqs_sensitivity_option(command[:], 4.0)
        assert out == ['mmseqs', 'easy-search', '-s', '4']


class TestAppendMmseqsPositiveIntOption:
    def test_leaves_command_unchanged_for_auto(self):
        command = ['mmseqs', 'easy-search']
        out = append_mmseqs_positive_int_option(command[:], '--max-seqs', 'auto')
        assert out == ['mmseqs', 'easy-search']

    def test_appends_integer_for_numeric_value(self):
        command = ['mmseqs', 'easy-search']
        out = append_mmseqs_positive_int_option(command[:], '--max-seqs', 20)
        assert out == ['mmseqs', 'easy-search', '--max-seqs', '20']

    @pytest.mark.parametrize('invalid_value', ['1.5G', '1g', '1P', '0G', '-1G'])
    def test_rejects_memory_values_not_accepted_by_mmseqs(self, invalid_value):
        with pytest.raises(ValueError, match='split-memory-limit'):
            append_mmseqs_memory_limit_option(
                ['mmseqs', 'easy-search'],
                invalid_value,
            )

    @pytest.mark.parametrize('valid_value', ['800B', '5K', '10M', '32G', '1T'])
    def test_accepts_mmseqs_memory_grammar(self, valid_value):
        command = append_mmseqs_memory_limit_option(
            ['mmseqs', 'easy-search'],
            valid_value,
        )
        assert command[-2:] == ['--split-memory-limit', valid_value]

    def test_resolve_mmseqs_dbtype_uses_runtime_context_without_mutating_args(self, monkeypatch):
        args = SimpleNamespace(mmseqs_exe='mmseqs')
        runtime_context = GetfastqRuntimeContext()
        calls = {'n': 0}

        def fake_run(cmd, stdout=None, stderr=None):
            calls['n'] += 1
            return subprocess.CompletedProcess(cmd, 0, stdout=b'nucleotide\n', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        first = resolve_mmseqs_dbtype(args, '/tmp/db', runtime_context=runtime_context)
        second = resolve_mmseqs_dbtype(args, '/tmp/db', runtime_context=runtime_context)

        assert first == 'nucleotide'
        assert second == 'nucleotide'
        assert calls['n'] == 1
        assert not hasattr(args, '_mmseqs_dbtype_cache')


class TestContamFilterDbPathResolution:
    def test_inferred_db_path_uses_out_dir_downloads_when_download_dir_is_inferred(self, tmp_path):
        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir='inferred',
            contam_filter_db='inferred',
            contam_filter_db_name='UniRef90',
        )
        observed = resolve_contam_filter_db_path(args)
        expected = os.path.join(os.path.realpath(str(tmp_path / 'out')), 'downloads', 'mmseqs2', 'UniRef90_DB')
        assert observed == expected

    def test_inferred_db_path_uses_custom_download_dir(self, tmp_path):
        custom_download_dir = tmp_path / 'shared_downloads'
        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(custom_download_dir),
            contam_filter_db='inferred',
            contam_filter_db_name='UniRef90',
        )
        observed = resolve_contam_filter_db_path(args)
        expected = os.path.join(os.path.realpath(str(custom_download_dir)), 'mmseqs2', 'UniRef90_DB')
        assert observed == expected

    def test_mmseqs_db_download_uses_lock_file(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        args = SimpleNamespace(
            out_dir=str(out_dir),
            download_dir='inferred',
            contam_filter_db='inferred',
            contam_filter_db_name='UniRef90',
            mmseqs_exe='mmseqs',
            threads=1,
            dump_print=False,
        )
        captured = {'lock_path': None, 'cmd': None}

        class DummyLock:
            def __init__(self, lock_path, lock_label='Lock', poll_seconds=5, timeout_seconds=3600):
                _ = (lock_label, poll_seconds, timeout_seconds)
                captured['lock_path'] = lock_path

            def __enter__(self):
                return None

            def __exit__(self, exc_type, exc, tb):
                return False

        def fake_run(cmd, stdout=None, stderr=None):
            _ = (stdout, stderr)
            captured['cmd'] = cmd
            db_path = cmd[3]
            with open(db_path, 'wt') as fout:
                fout.write('mockdb')
            with open(db_path + '.dbtype', 'wt') as fout:
                fout.write('mock')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.acquire_exclusive_lock', DummyLock)
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        db_path = ensure_mmseqs_contam_taxonomy_db_exists(args)
        assert captured['cmd'] is not None
        assert captured['lock_path'] == os.path.join(
            os.path.realpath(str(out_dir)),
            'downloads',
            'locks',
            'mmseqs_db_UniRef90_DB.lock',
        )
        assert captured['cmd'][4] == os.path.dirname(db_path)
        assert os.path.isfile(db_path + '.ready')

    def test_mmseqs_db_existing_without_ready_marker_is_rebuilt(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        args = SimpleNamespace(
            out_dir=str(out_dir),
            download_dir='inferred',
            contam_filter_db='inferred',
            contam_filter_db_name='UniRef90',
            mmseqs_exe='mmseqs',
            threads=1,
            dump_print=False,
        )
        db_path = resolve_contam_filter_db_path(args)
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        with open(db_path, 'wt') as fout:
            fout.write('mockdb')
        with open(db_path + '.dbtype', 'wt') as fout:
            fout.write('mock')

        calls = []

        def fake_run(cmd, stdout=None, stderr=None):
            _ = (stdout, stderr)
            calls.append(cmd)
            with open(cmd[3], 'wt') as fout:
                fout.write('rebuilt')
            with open(cmd[3] + '.dbtype', 'wt') as fout:
                fout.write('rebuilt-type')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        observed = ensure_mmseqs_contam_taxonomy_db_exists(args)
        assert observed == db_path
        assert len(calls) == 1
        assert open(db_path).read() == 'rebuilt'
        assert os.path.isfile(db_path + '.ready')

class TestGetfastqXmlRetrieval:
    def test_trace_run_xml_rejects_entities(self, monkeypatch):
        malicious_xml = b'<!DOCTYPE root [<!ENTITY xxe SYSTEM "file:///etc/passwd">]><root>&xxe;</root>'
        monkeypatch.setattr(
            'amalgkit.getfastq.urllib.request.urlopen',
            lambda *_args, **_kwargs: BytesIO(malicious_xml),
        )

        with pytest.raises(EntitiesForbidden):
            fetch_trace_run_xml_root('SRR000000')

    class _DummyTree:
        def __init__(self, root):
            self._root = root

        def getroot(self):
            return self._root

    def test_returns_empty_root_when_no_records(self, monkeypatch):
        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': []})
        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', lambda **kwargs: (_ for _ in ()).throw(AssertionError('efetch should not be called')))

        root = getfastq_getxml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE_SET'

    def test_retries_esearch_once_on_urlerror(self, monkeypatch):
        esearch_calls = {'n': 0}

        def flaky_esearch(**_kwargs):
            esearch_calls['n'] += 1
            if esearch_calls['n'] == 1:
                raise urllib.error.URLError('temporary network failure')
            return object()

        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', flaky_esearch)
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': []})
        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', lambda **kwargs: object())

        root = getfastq_getxml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE_SET'
        assert esearch_calls['n'] == 2

    def test_retries_efetch_once_on_urlerror(self, monkeypatch):
        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': ['ID1']})
        efetch_calls = {'n': 0}

        def flaky_efetch(**_kwargs):
            efetch_calls['n'] += 1
            if efetch_calls['n'] == 1:
                raise urllib.error.URLError('temporary network failure')
            return object()

        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', flaky_efetch)
        monkeypatch.setattr(
            'amalgkit.sra.parse_untrusted_xml',
            lambda handle: self._DummyTree(ET.Element('EXPERIMENT_PACKAGE'))
        )
        monkeypatch.setattr('amalgkit.getfastq.time.sleep', lambda *_args, **_kwargs: None)

        root = getfastq_getxml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE'
        assert efetch_calls['n'] == 2

    def test_batches_without_extra_request_on_exact_multiple(self, monkeypatch):
        id_list = ['ID{}'.format(i) for i in range(2000)]
        efetch_calls = []

        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': id_list})

        def fake_efetch(**kwargs):
            efetch_calls.append(list(kwargs['id']))
            return object()

        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', fake_efetch)
        monkeypatch.setattr(
            'amalgkit.sra.parse_untrusted_xml',
            lambda handle: self._DummyTree(ET.Element('EXPERIMENT_PACKAGE'))
        )

        root = getfastq_getxml(search_term='SRR000000', retmax=1000)

        assert root is not None
        assert len(efetch_calls) == 2
        assert [len(c) for c in efetch_calls] == [1000, 1000]

    def test_retries_when_xml_chunk_parse_fails_once(self, monkeypatch):
        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': ['ID1']})
        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', lambda **kwargs: object())
        parse_calls = {'n': 0}

        def flaky_parse(_handle):
            parse_calls['n'] += 1
            if parse_calls['n'] == 1:
                raise ET.ParseError('truncated xml')
            return self._DummyTree(ET.Element('EXPERIMENT_PACKAGE'))

        monkeypatch.setattr('amalgkit.sra.parse_untrusted_xml', flaky_parse)

        root = getfastq_getxml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE'
        assert parse_calls['n'] == 2

    def test_raises_when_xml_chunk_parsing_never_succeeds(self, monkeypatch):
        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': ['ID1']})
        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', lambda **kwargs: object())
        monkeypatch.setattr(
            'amalgkit.sra.parse_untrusted_xml',
            lambda _handle: (_ for _ in ()).throw(ET.ParseError('broken xml')),
        )

        with pytest.raises(RuntimeError, match='Failed to parse Entrez XML chunk'):
            getfastq_getxml(search_term='SRR000000', retmax=1000)

    def test_merges_package_set_chunks_without_nested_container(self, monkeypatch):
        id_list = ['ID{}'.format(i) for i in range(2000)]
        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': id_list})
        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', lambda **kwargs: object())

        def fake_parse(_handle):
            root = ET.Element('EXPERIMENT_PACKAGE_SET')
            ET.SubElement(root, 'EXPERIMENT_PACKAGE')
            return self._DummyTree(root)

        monkeypatch.setattr('amalgkit.sra.parse_untrusted_xml', fake_parse)

        root = getfastq_getxml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE_SET'
        assert len(root.findall('./EXPERIMENT_PACKAGE')) == 2
        assert len(root.findall('./EXPERIMENT_PACKAGE_SET')) == 0

    def test_wraps_non_set_chunks_to_preserve_multiple_records(self, monkeypatch):
        id_list = ['ID{}'.format(i) for i in range(2000)]
        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': id_list})
        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', lambda **kwargs: object())

        def fake_parse(_handle):
            root = ET.Element('EXPERIMENT_PACKAGE')
            ET.SubElement(root, 'RUN_SET')
            return self._DummyTree(root)

        monkeypatch.setattr('amalgkit.sra.parse_untrusted_xml', fake_parse)

        root = getfastq_getxml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE_SET'
        assert len(root.findall('./EXPERIMENT_PACKAGE')) == 2

    def test_raises_when_error_tag_present(self, monkeypatch):
        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': ['ID1']})
        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', lambda **kwargs: object())
        err_root = ET.Element('EXPERIMENT_PACKAGE')
        err = ET.SubElement(err_root, 'Error')
        err.text = 'SRA error'
        monkeypatch.setattr(
            'amalgkit.sra.parse_untrusted_xml',
            lambda handle: self._DummyTree(err_root)
        )

        with pytest.raises(RuntimeError, match='Error found in Entrez XML response'):
            getfastq_getxml(search_term='SRR000000', retmax=1000)


class TestGetfastqMetadataIdFiltering:
    def test_rejects_simultaneous_id_and_id_list(self, tmp_path):
        id_list_path = tmp_path / 'ids.txt'
        id_list_path.write_text('SRR001\n')

        class Args:
            id = 'SRR001'
            id_list = str(id_list_path)
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)

        with pytest.raises(ValueError, match='mutually exclusive'):
            getfastq_metadata(Args())

    def test_rejects_missing_id_list_path(self, tmp_path):
        class Args:
            id = None
            id_list = str(tmp_path / 'missing_ids.txt')
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)

        with pytest.raises(FileNotFoundError, match='SRA ID list file not found'):
            getfastq_metadata(Args())

    def test_rejects_directory_id_list_path(self, tmp_path):
        id_list_dir = tmp_path / 'ids_dir'
        id_list_dir.mkdir()

        class Args:
            id = None
            id_list = str(id_list_dir)
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)

        with pytest.raises(IsADirectoryError, match='SRA ID list path exists but is not a file'):
            getfastq_metadata(Args())

    def test_id_mode_layout_filter_normalizes_metadata_layout_values(self, tmp_path, monkeypatch):
        def fake_getxml(search_term, retmax=1000):
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root=None):
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': ['SRR001', 'SRR002'],
                'lib_layout': ['  PAIRED ', 'single'],
                'scientific_name': ['Sp1', 'Sp2'],
                'total_bases': ['100', '100'],
                'spot_length': ['50', '50'],
                'exclusion': ['no', 'no'],
            }))

        class Args:
            id = 'SRR001'
            id_list = None
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'paired'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)

        metadata = getfastq_metadata(Args())

        assert metadata.df['run'].tolist() == ['SRR001']

    def test_id_mode_invalid_explicit_layout_raises(self, tmp_path, monkeypatch):
        def fake_getxml(search_term, retmax=1000):
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root=None):
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': ['SRR001'],
                'lib_layout': ['single'],
                'scientific_name': ['Sp1'],
                'total_bases': ['100'],
                'spot_length': ['50'],
                'exclusion': ['no'],
            }))

        class Args:
            id = 'SRR001'
            id_list = None
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'paired_end'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)

        with pytest.raises(ValueError, match='--layout must be one of'):
            getfastq_metadata(Args())

    def test_id_mode_sciname_filter_strips_whitespace(self, tmp_path, monkeypatch):
        def fake_getxml(search_term, retmax=1000):
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root=None):
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': ['SRR001', 'SRR002'],
                'lib_layout': ['single', 'single'],
                'scientific_name': [' Target species ', 'Other species'],
                'total_bases': ['100', '100'],
                'spot_length': ['50', '50'],
                'exclusion': ['no', 'no'],
            }))

        class Args:
            id = 'SRR001'
            id_list = None
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = 'Target species'
            metadata = 'unused'
            out_dir = str(tmp_path)

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)

        metadata = getfastq_metadata(Args())

        assert metadata.df['run'].tolist() == ['SRR001']


class TestGetfastqMetadataIdListParsing:
    def test_ignores_blank_and_comment_lines(self, tmp_path, monkeypatch):
        id_list_path = tmp_path / 'ids.txt'
        id_list_path.write_text('\n'.join([
            '# comment',
            'SRR100',
            '',
            '   # indented comment',
            'SRR200  ',
            '   ',
        ]))

        called_terms = []

        def fake_getxml(search_term, retmax=1000):
            called_terms.append(search_term)
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root):
            sid = xml_root.get('sid')
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': [sid],
                'lib_layout': ['single'],
                'scientific_name': ['Testus species'],
                'total_bases': ['100'],
                'spot_length': ['50'],
                'exclusion': ['no'],
            }))

        class Args:
            id = None
            id_list = str(id_list_path)
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)
            internal_jobs = 1

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)

        metadata = getfastq_metadata(Args())

        assert called_terms == ['SRR100', 'SRR200']
        assert set(metadata.df['run'].tolist()) == {'SRR100', 'SRR200'}

    def test_parallel_id_list_fetch_uses_capped_workers(self, tmp_path, monkeypatch):
        id_list_path = tmp_path / 'ids.txt'
        id_list_path.write_text('\n'.join(['SRR100', 'SRR200', 'SRR300']))

        called_terms = []
        captured = {}

        def fake_getxml(search_term, retmax=1000):
            called_terms.append(search_term)
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root):
            sid = xml_root.get('sid')
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': [sid],
                'lib_layout': ['single'],
                'scientific_name': ['Testus species'],
                'total_bases': ['100'],
                'spot_length': ['50'],
                'exclusion': ['no'],
            }))

        def fake_run_tasks(task_items, task_fn, max_workers=1):
            captured['max_workers'] = max_workers
            results = {}
            failures = []
            for item in task_items:
                results[item] = task_fn(item)
            return results, failures

        class Args:
            id = None
            id_list = str(id_list_path)
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)
            internal_jobs=8

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)
        monkeypatch.setattr('amalgkit.getfastq.run_tasks_with_optional_threads', fake_run_tasks)

        metadata = getfastq_metadata(Args())

        assert captured['max_workers'] == 3
        assert called_terms == ['SRR100', 'SRR200', 'SRR300']
        assert metadata.df['run'].tolist() == ['SRR100', 'SRR200', 'SRR300']

    def test_id_list_fetch_respects_total_core_budget(self, tmp_path, monkeypatch):
        id_list_path = tmp_path / 'ids.txt'
        id_list_path.write_text('\n'.join(['SRR100', 'SRR200', 'SRR300']))

        captured = {}

        def fake_getxml(search_term, retmax=1000):
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root):
            sid = xml_root.get('sid')
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': [sid],
                'lib_layout': ['single'],
                'scientific_name': ['Testus species'],
                'total_bases': ['100'],
                'spot_length': ['50'],
                'exclusion': ['no'],
            }))

        def fake_run_tasks(task_items, task_fn, max_workers=1):
            captured['max_workers'] = max_workers
            results = {}
            for item in task_items:
                results[item] = task_fn(item)
            return results, []

        class Args:
            id = None
            id_list = str(id_list_path)
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)
            internal_jobs = 8
            threads = 1
            internal_cpu_budget = 4

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)
        monkeypatch.setattr('amalgkit.getfastq.run_tasks_with_optional_threads', fake_run_tasks)

        _ = getfastq_metadata(Args())

        assert captured['max_workers'] == 1

    def test_id_list_keeps_duplicate_ids_in_output_order(self, tmp_path, monkeypatch):
        id_list_path = tmp_path / 'ids.txt'
        id_list_path.write_text('\n'.join(['SRR100', 'SRR100', 'SRR200']))

        def fake_getxml(search_term, retmax=1000):
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root):
            sid = xml_root.get('sid')
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': [sid],
                'lib_layout': ['single'],
                'scientific_name': ['Testus species'],
                'total_bases': ['100'],
                'spot_length': ['50'],
                'exclusion': ['no'],
            }))

        class Args:
            id = None
            id_list = str(id_list_path)
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)
            internal_jobs=4

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)

        metadata = getfastq_metadata(Args())

        assert metadata.df['run'].tolist() == ['SRR100', 'SRR100', 'SRR200']

    def test_id_list_deduplicates_entrez_fetches(self, tmp_path, monkeypatch):
        id_list_path = tmp_path / 'ids.txt'
        id_list_path.write_text('\n'.join(['SRR100', 'SRR100', 'SRR200']))
        called_terms = []

        def fake_getxml(search_term, retmax=1000):
            called_terms.append(search_term)
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root):
            sid = xml_root.get('sid')
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': [sid],
                'lib_layout': ['single'],
                'scientific_name': ['Testus species'],
                'total_bases': ['100'],
                'spot_length': ['50'],
                'exclusion': ['no'],
            }))

        class Args:
            id = None
            id_list = str(id_list_path)
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)
            internal_jobs=4

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)

        metadata = getfastq_metadata(Args())

        assert set(called_terms) == {'SRR100', 'SRR200'}
        assert len(called_terms) == 2
        assert metadata.df['run'].tolist() == ['SRR100', 'SRR100', 'SRR200']

    def test_id_list_applies_layout_and_sciname_filters(self, tmp_path, monkeypatch):
        id_list_path = tmp_path / 'ids.txt'
        id_list_path.write_text('\n'.join(['SRR100', 'SRR200']))

        def fake_getxml(search_term, retmax=1000):
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root=None):
            sid = xml_root.get('sid')
            if sid == 'SRR100':
                return Metadata.from_DataFrame(pandas.DataFrame({
                    'run': [sid],
                    'lib_layout': ['  PAIRED '],
                    'scientific_name': ['Target species'],
                    'total_bases': ['100'],
                    'spot_length': ['50'],
                    'exclusion': ['no'],
                }))
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': [sid],
                'lib_layout': ['single'],
                'scientific_name': ['Other species'],
                'total_bases': ['100'],
                'spot_length': ['50'],
                'exclusion': ['no'],
            }))

        class Args:
            id = None
            id_list = str(id_list_path)
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'paired'
            sci_name = 'Target species'
            metadata = 'unused'
            out_dir = str(tmp_path)
            internal_jobs=2

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)

        metadata = getfastq_metadata(Args())

        assert metadata.df['run'].tolist() == ['SRR100']

    def test_id_list_batch_selects_single_row_without_is_sampled(self, tmp_path, monkeypatch):
        id_list_path = tmp_path / 'ids.txt'
        id_list_path.write_text('\n'.join(['SRR100', 'SRR200']))

        def fake_getxml(search_term, retmax=1000):
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root=None):
            sid = xml_root.get('sid')
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': [sid],
                'lib_layout': ['single'],
                'scientific_name': ['Target species'],
                'total_bases': ['100'],
                'spot_length': ['50'],
                'exclusion': ['no'],
            }))

        class Args:
            id = None
            id_list = str(id_list_path)
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)
            internal_jobs = 2
            batch = 2

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)

        metadata = getfastq_metadata(Args())

        assert metadata.df['run'].tolist() == ['SRR200']

    def test_id_list_batch_too_large_exits_cleanly(self, tmp_path, monkeypatch):
        id_list_path = tmp_path / 'ids.txt'
        id_list_path.write_text('SRR100\n')

        def fake_getxml(search_term, retmax=1000):
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root=None):
            sid = xml_root.get('sid')
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': [sid],
                'lib_layout': ['single'],
                'scientific_name': ['Target species'],
                'total_bases': ['100'],
                'spot_length': ['50'],
                'exclusion': ['no'],
            }))

        class Args:
            id = None
            id_list = str(id_list_path)
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)
            internal_jobs = 1
            batch = 2

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)

        with pytest.raises(AmalgkitExit) as exit_info:
            getfastq_metadata(Args())

        assert exit_info.value.exit_code == 0
