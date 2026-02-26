import pytest
import pandas
import numpy
import time
import gzip
import subprocess
import xml.etree.ElementTree as ET

import os
import urllib.error
from types import SimpleNamespace

from amalgkit.getfastq import (
    getfastq_search_term,
    getfastq_getxml,
    getfastq_metadata,
    get_range,
    get_layout,
    detect_concat_input_files,
    concat_fastq,
    append_file_binary,
    remove_experiment_without_run,
    check_metadata_validity,
    initialize_global_params,
    rename_fastq,
    rename_reads,
    remove_old_intermediate_files,
    remove_intermediate_files,
    is_getfastq_output_present,
    initialize_columns,
    calc_2nd_ranges,
    has_remaining_spots_after_first_round,
    is_2nd_round_needed,
    get_identical_paired_ratio,
    maybe_treat_paired_as_single,
    parse_fastp_metrics,
    parse_fastp_summary_counts,
    update_fastp_metrics,
    write_fastp_stats,
    write_getfastq_stats,
    run_fastp,
    run_sortmerna_rrna_filter,
    run_mmseqs_contam_filter,
    run_mmseqs_easy_taxonomy_single_fastq,
    parse_sortmerna_output_counts,
    update_metadata_after_rrna_filter,
    parse_fasterq_dump_written_spots,
    parse_fasterq_dump_written_reads,
    infer_written_spots_from_written_reads,
    sequence_extraction,
    sequence_extraction_2nd_round,
    sequence_extraction_private,
    estimate_num_written_spots_from_fastq,
    compress_fasterq_output_files,
    should_compress_fasterq_output_before_filters,
    normalize_fasterq_size_check,
    count_fastq_records_and_bases,
    summarize_layout_fastq_records_and_bases,
    download_sra,
    run_fasterq_dump,
    getfastq_main,
    check_getfastq_dependency,
    resolve_sortmerna_refs,
    count_fastq_records,
    collect_valid_run_ids,
    remove_sra_files,
    remove_sra_path,
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
        args = SimpleNamespace(fastp=True, rrna_filter=True, filter_order='rrna_first')
        assert not should_compress_fasterq_output_before_filters(args)

    def test_skips_compression_when_contam_filter_is_enabled(self):
        args = SimpleNamespace(fastp=False, rrna_filter=False, contam_filter=True)
        assert not should_compress_fasterq_output_before_filters(args)


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


class TestSortmernaLogCounts:
    def test_parse_sortmerna_output_counts_single(self):
        stdout_txt = '\n'.join([
            'all_reads_count= 1000 all_reads_len= 100000',
            '[run:11] && Reads added: 1000 Num aligned reads (passing E-value): 240',
        ])
        out = parse_sortmerna_output_counts(
            stdout_txt=stdout_txt,
            stderr_txt='',
            layout='single',
            num_input_files=1,
            known_input_counts=None,
        )
        assert out == {'num_spots': 760, 'bp_total': 76000}

    def test_parse_sortmerna_output_counts_paired_uses_known_input_for_divisor_and_bp(self):
        stdout_txt = '\n'.join([
            'all_reads_count= 10000 all_reads_len= 1000000',
            '[run:22] && Reads added: 10000 Num aligned reads (passing E-value): 5944',
        ])
        out = parse_sortmerna_output_counts(
            stdout_txt=stdout_txt,
            stderr_txt='',
            layout='paired',
            num_input_files=2,
            known_input_counts={'num_spots': 5000, 'bp_total': 1000000},
        )
        assert out == {'num_spots': 2028, 'bp_total': 405600}


class TestRunSortmernaWithLogCounts:
    @staticmethod
    def _write_fastq_gz(path, reads):
        with gzip.open(path, 'wt') as out:
            for i, seq in enumerate(reads):
                out.write('@r{}\n'.format(i))
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_uses_sortmerna_log_counts_without_fastq_rescan(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        in_path = tmp_path / '{}.fastp.fastq.gz'.format(sra_id)
        self._write_fastq_gz(str(in_path), ['AAAA'])
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'num_rrna_in': [0],
            'num_rrna_out': [0],
            'bp_rrna_in': [0],
            'bp_rrna_out': [0],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'metadata_idx': 0,
            'current_ext': '.fastp.fastq.gz',
            'getfastq_sra_dir': str(tmp_path),
        }
        args = SimpleNamespace(
            rrna_filter=True,
            threads=1,
            sortmerna_exe='sortmerna',
            sortmerna_option='',
            remove_tmp=False,
            dump_print=False,
            seqkit_exe='seqkit',
            out_dir=str(tmp_path),
            download_dir='inferred',
        )

        monkeypatch.setattr('amalgkit.getfastq.resolve_sortmerna_refs', lambda _args: ['/tmp/ref.fa'])

        def fake_execute_sortmerna(sortmerna_command, args):
            other_prefix = sortmerna_command[sortmerna_command.index('--other') + 1]
            out_path = other_prefix + '.fastq.gz'
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            self._write_fastq_gz(out_path, ['AAAA'])
            stdout_txt = '\n'.join([
                'all_reads_count= 10 all_reads_len= 1000',
                '[run:11] && Reads added: 10 Num aligned reads (passing E-value): 3',
            ]) + '\n'
            return subprocess.CompletedProcess(sortmerna_command, 0, stdout=stdout_txt.encode('utf8'), stderr=b'')

        def fail_summarize(*_args, **_kwargs):
            raise AssertionError('FASTQ scan should not be called when SortMeRNA log counts are available.')

        monkeypatch.setattr('amalgkit.getfastq.execute_sortmerna_command', fake_execute_sortmerna)
        monkeypatch.setattr('amalgkit.getfastq.summarize_layout_fastq_records_and_bases', fail_summarize)

        metadata, run_file_state = run_sortmerna_rrna_filter(
            sra_stat=sra_stat,
            args=args,
            output_dir=str(tmp_path),
            metadata=metadata,
            files={'{}.fastp.fastq.gz'.format(sra_id)},
            known_input_counts={'num_spots': 10, 'bp_total': 1000},
            return_file_state=True,
        )

        assert metadata.df.loc[0, 'num_rrna_in'] == 10
        assert metadata.df.loc[0, 'bp_rrna_in'] == 1000
        assert metadata.df.loc[0, 'num_rrna_out'] == 7
        assert metadata.df.loc[0, 'bp_rrna_out'] == 700
        assert run_file_state.has('{}.rrna.fastq.gz'.format(sra_id))


class TestRunMmseqsContamFilter:
    @staticmethod
    def _write_fastq_gz(path, read_id_seq_pairs):
        with gzip.open(path, 'wt') as out:
            for read_id, seq in read_id_seq_pairs:
                out.write('@{}\n'.format(read_id))
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_paired_strict_pair_removal_and_unclassified_keep(self, tmp_path, monkeypatch):
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

        def fake_run_mmseqs(args, input_path, target_db, result_prefix, tmp_dir):
            _ = (args, target_db, tmp_dir)
            os.makedirs(os.path.dirname(result_prefix), exist_ok=True)
            lca_path = result_prefix + '_lca.tsv'
            if input_path.endswith('_1.fastq.gz'):
                # rA -> match, rB -> mismatch, rC -> unclassified(keep)
                lines = [
                    'rA/1\t11\tphylum\tmatch',
                    'rB/1\t22\tphylum\tmismatch',
                    'rC/1\t0\tno rank\tunclassified',
                ]
            else:
                lines = [
                    'rA/2\t11\tphylum\tmatch',
                    'rB/2\t11\tphylum\tmatch',
                    'rC/2\t0\tno rank\tunclassified',
                ]
            with open(lca_path, 'wt') as fout:
                fout.write('\n'.join(lines) + '\n')
            return lca_path

        # Map assigned taxid to phylum-level taxid.
        monkeypatch.setattr(
            'amalgkit.getfastq.resolve_taxid_at_rank',
            lambda taxid, rank_name, ncbi, rank_cache: 500 if int(taxid) == 11 else (600 if int(taxid) == 22 else None),
        )
        monkeypatch.setattr('amalgkit.getfastq.run_mmseqs_easy_taxonomy_single_fastq', fake_run_mmseqs)
        monkeypatch.setattr('amalgkit.getfastq.ete4.NCBITaxa', lambda: object())

        metadata, run_file_state = run_mmseqs_contam_filter(
            sra_stat=sra_stat,
            args=Args(),
            output_dir=str(tmp_path),
            metadata=metadata,
            files={'{}_1.fastq.gz'.format(sra_id), '{}_2.fastq.gz'.format(sra_id)},
            return_file_state=True,
        )

        out1 = tmp_path / '{}_1.contam.fastq.gz'.format(sra_id)
        out2 = tmp_path / '{}_2.contam.fastq.gz'.format(sra_id)
        assert out1.exists()
        assert out2.exists()
        assert count_fastq_records(str(out1)) == 2
        assert count_fastq_records(str(out2)) == 2
        assert metadata.df.loc[0, 'num_contam_in'] == 3
        assert metadata.df.loc[0, 'num_contam_out'] == 2
        assert metadata.df.loc[0, 'bp_contam_in'] == 24
        assert metadata.df.loc[0, 'bp_contam_out'] == 16
        assert run_file_state.has('{}_1.contam.fastq.gz'.format(sra_id))
        assert run_file_state.has('{}_2.contam.fastq.gz'.format(sra_id))


class TestRunMmseqsEasyTaxonomy:
    def test_adds_search_type_for_nucleotide_db(self, monkeypatch):
        class Args:
            mmseqs_exe = 'mmseqs'
            threads = 2
            dump_print = False

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

    def test_omits_search_type_for_aminoacid_db(self, monkeypatch):
        class Args:
            mmseqs_exe = 'mmseqs'
            threads = 2
            dump_print = False

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


class TestGetfastqXmlRetrieval:
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
            'amalgkit.getfastq.ET.parse',
            lambda handle: self._DummyTree(ET.Element('EXPERIMENT_PACKAGE'))
        )

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
            'amalgkit.getfastq.ET.parse',
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

        monkeypatch.setattr('amalgkit.getfastq.ET.parse', flaky_parse)

        root = getfastq_getxml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE'
        assert parse_calls['n'] == 2

    def test_raises_when_xml_chunk_parsing_never_succeeds(self, monkeypatch):
        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': ['ID1']})
        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', lambda **kwargs: object())
        monkeypatch.setattr(
            'amalgkit.getfastq.ET.parse',
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

        monkeypatch.setattr('amalgkit.getfastq.ET.parse', fake_parse)

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

        monkeypatch.setattr('amalgkit.getfastq.ET.parse', fake_parse)

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
            'amalgkit.getfastq.ET.parse',
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

        with pytest.raises(SystemExit) as exit_info:
            getfastq_metadata(Args())

        assert exit_info.value.code == 0


class TestGetRange:
    def test_total_within_max(self):
        sra_stat = {'total_spot': 1000, 'num_read_per_sra': 500}
        start, end = get_range(sra_stat, offset=0, total_sra_bp=100, max_bp=200)
        assert start == 1
        assert end == 1000

    def test_total_exceeds_max_with_offset(self):
        sra_stat = {'total_spot': 10000, 'num_read_per_sra': 5000}
        start, end = get_range(sra_stat, offset=100, total_sra_bp=2000, max_bp=1000)
        assert start == 100
        assert end == 5100

    def test_total_exceeds_max_offset_too_large(self):
        sra_stat = {'total_spot': 6000, 'num_read_per_sra': 5000}
        start, end = get_range(sra_stat, offset=2000, total_sra_bp=2000, max_bp=1000)
        # total_spot > num_read_per_sra but total_spot <= num_read_per_sra + offset
        assert start == 1000  # total_spot - num_read_per_sra
        assert end == 6000

    def test_total_spot_less_than_num_reads(self):
        sra_stat = {'total_spot': 3000, 'num_read_per_sra': 5000}
        start, end = get_range(sra_stat, offset=0, total_sra_bp=2000, max_bp=1000)
        assert start == 1
        assert end == 3000


# ---------------------------------------------------------------------------
# get_layout (wiki: auto-detects paired/single from metadata)
# ---------------------------------------------------------------------------

class TestGetLayout:
    def test_auto_prefers_paired(self):
        """Wiki: auto layout prefers paired when multiple layouts exist."""
        class Args:
            layout = 'auto'
        data = {'lib_layout': ['paired', 'single', 'paired']}
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2', 'R3'],
            'lib_layout': ['paired', 'single', 'paired'],
            'exclusion': ['no', 'no', 'no'],
        }))
        result = get_layout(Args(), m)
        assert result == 'paired'

    def test_auto_single_only(self):
        """When all samples are single-end, auto returns single."""
        class Args:
            layout = 'auto'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2'],
            'lib_layout': ['single', 'single'],
            'exclusion': ['no', 'no'],
        }))
        result = get_layout(Args(), m)
        assert result == 'single'

    def test_explicit_override(self):
        """Explicit layout overrides metadata."""
        class Args:
            layout = 'single'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'lib_layout': ['paired'],
            'exclusion': ['no'],
        }))
        result = get_layout(Args(), m)
        assert result == 'single'

    def test_auto_normalizes_case_and_whitespace(self):
        class Args:
            layout = 'auto'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2'],
            'lib_layout': ['  PAIRED ', 'single'],
            'exclusion': ['no', 'no'],
        }))
        result = get_layout(Args(), m)
        assert result == 'paired'

    def test_auto_raises_when_no_valid_layout(self):
        class Args:
            layout = 'auto'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2'],
            'lib_layout': [None, ''],
            'exclusion': ['no', 'no'],
        }))
        with pytest.raises(ValueError, match='No valid lib_layout'):
            get_layout(Args(), m)

    def test_auto_raises_when_lib_layout_column_missing(self):
        class Args:
            layout = 'auto'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'exclusion': ['no'],
        }))
        m.df = m.df.drop(columns=['lib_layout'])
        with pytest.raises(ValueError, match='Column \"lib_layout\" is required'):
            get_layout(Args(), m)

    def test_explicit_invalid_layout_raises(self):
        class Args:
            layout = 'paired_end'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'lib_layout': ['paired'],
            'exclusion': ['no'],
        }))
        with pytest.raises(ValueError, match='--layout must be one of'):
            get_layout(Args(), m)


# ---------------------------------------------------------------------------
# concat_fastq
# ---------------------------------------------------------------------------

class TestConcatFastq:
    def _metadata_single(self, runs):
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': runs,
            'lib_layout': ['single'] * len(runs),
            'total_spots': [1] * len(runs),
            'spot_length': [4] * len(runs),
            'total_bases': [4] * len(runs),
            'scientific_name': ['Sp'] * len(runs),
            'exclusion': ['no'] * len(runs),
        }))

    def _metadata_paired(self, runs):
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': runs,
            'lib_layout': ['paired'] * len(runs),
            'total_spots': [1] * len(runs),
            'spot_length': [4] * len(runs),
            'total_bases': [4] * len(runs),
            'scientific_name': ['Sp'] * len(runs),
            'exclusion': ['no'] * len(runs),
        }))

    def test_single_file_uses_single_directory_scan(self, tmp_path, monkeypatch):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('ACGT\n')
        metadata = self._metadata_single(['SRR001'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'NEWID_',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}
        calls = {'num': 0}

        def fake_list_run_dir_files(work_dir):
            calls['num'] += 1
            return set(os.listdir(work_dir))

        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fake_list_run_dir_files)

        concat_fastq(args, metadata, str(tmp_path), g)

        assert calls['num'] == 1
        assert (tmp_path / 'NEWID_SRR001.amalgkit.fastq.gz').exists()

    def test_single_file_with_run_id_does_not_duplicate_prefix(self, tmp_path):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('ACGT\n')
        metadata = self._metadata_single(['SRR001'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'SRR001',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        concat_fastq(args, metadata, str(tmp_path), g)

        assert (tmp_path / 'SRR001.amalgkit.fastq.gz').exists()
        assert not (tmp_path / 'SRR001SRR001.amalgkit.fastq.gz').exists()

    def test_remove_tmp_reuses_prefetched_file_set_for_extension_lookup(self, tmp_path, monkeypatch):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        (tmp_path / 'SRR002.amalgkit.fastq.gz').write_text('CCCC\n')
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': True,
        })()
        g = {'num_bp_per_sra': 4}
        seen = []

        def fake_get_newest_intermediate_file_extension(sra_stat, work_dir, files=None):
            assert files is not None
            seen.append((sra_stat['sra_id'], files))
            return '.amalgkit.fastq.gz'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            fake_get_newest_intermediate_file_extension,
        )
        monkeypatch.setattr('amalgkit.getfastq.remove_intermediate_files', lambda sra_stat, ext, work_dir: None)

        concat_fastq(args, metadata, str(tmp_path), g)

        assert (tmp_path / 'MERGED.amalgkit.fastq.gz').exists()
        assert [run_id for run_id, _ in seen] == ['SRR001', 'SRR002']
        assert seen[0][1] is seen[1][1]

    def test_concat_uses_system_cat_when_available(self, tmp_path, monkeypatch):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_bytes(b'AAAA\n')
        (tmp_path / 'SRR002.amalgkit.fastq.gz').write_bytes(b'CCCC\n')
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}
        captured = {}

        monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda name: '/bin/cat' if name == 'cat' else None)

        def fake_run(cmd, stdout=None, stderr=None):
            captured['cmd'] = cmd
            stdout.write(b'AAAA\nCCCC\n')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('append_file_binary should not be used when system cat succeeds.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.append_file_binary', fail_if_called)

        concat_fastq(args, metadata, str(tmp_path), g)

        assert captured['cmd'][0] == '/bin/cat'
        assert (tmp_path / 'MERGED.amalgkit.fastq.gz').exists()
        assert (tmp_path / 'MERGED.amalgkit.fastq.gz').read_bytes() == b'AAAA\nCCCC\n'

    def test_skips_missing_run_ids_while_concatenating(self, tmp_path):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        (tmp_path / 'SRR002.amalgkit.fastq.gz').write_text('CCCC\n')
        metadata = self._metadata_single(['SRR001', numpy.nan, 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        concat_fastq(args, metadata, str(tmp_path), g)

        merged = tmp_path / 'MERGED.amalgkit.fastq.gz'
        assert merged.exists()
        assert merged.read_text() == 'AAAA\nCCCC\n'

    def test_raises_when_expected_input_fastq_is_missing(self, tmp_path, monkeypatch):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        # Simulate a race where a second input disappears after an initial directory snapshot.
        monkeypatch.setattr(
            'amalgkit.getfastq.list_run_dir_files',
            lambda _work_dir: {'SRR001.amalgkit.fastq.gz', 'SRR002.amalgkit.fastq.gz'},
        )

        with pytest.raises(FileNotFoundError, match='Dumped fastq not found'):
            concat_fastq(args, metadata, str(tmp_path), g)

    def test_single_shortcut_does_not_hide_missing_runs(self, tmp_path):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        with pytest.raises(FileNotFoundError, match='Dumped fastq not found'):
            concat_fastq(args, metadata, str(tmp_path), g)

    def test_paired_shortcut_does_not_hide_missing_mates_across_runs(self, tmp_path):
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz').write_text('AAAA\n')
        (tmp_path / 'SRR002_2.amalgkit.fastq.gz').write_text('CCCC\n')
        metadata = self._metadata_paired(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        with pytest.raises(FileNotFoundError, match='Dumped fastq not found'):
            concat_fastq(args, metadata, str(tmp_path), g)

    def test_ignores_duplicate_run_ids_in_metadata(self, tmp_path):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        metadata = self._metadata_single(['SRR001', 'SRR001'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED_',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        concat_fastq(args, metadata, str(tmp_path), g)

        merged = tmp_path / 'MERGED_SRR001.amalgkit.fastq.gz'
        assert merged.exists()
        assert merged.read_text() == 'AAAA\n'

    def test_raises_when_concat_output_path_exists_as_directory(self, tmp_path):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        (tmp_path / 'SRR002.amalgkit.fastq.gz').write_text('CCCC\n')
        (tmp_path / 'MERGED.amalgkit.fastq.gz').mkdir()
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        with pytest.raises(IsADirectoryError, match='Concatenation output path exists but is not a file'):
            concat_fastq(args, metadata, str(tmp_path), g)

    def test_raises_when_concat_input_path_is_directory(self, tmp_path, monkeypatch):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        (tmp_path / 'SRR002.amalgkit.fastq.gz').mkdir()
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        monkeypatch.setattr(
            'amalgkit.getfastq.list_run_dir_files',
            lambda _work_dir: {'SRR001.amalgkit.fastq.gz', 'SRR002.amalgkit.fastq.gz'},
        )

        with pytest.raises(IsADirectoryError, match='Concatenation input path exists but is not a file'):
            concat_fastq(args, metadata, str(tmp_path), g)

    def test_paired_concat_uses_parallel_workers_when_threads_gt_one(self, tmp_path, monkeypatch):
        metadata = self._metadata_paired(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
            'threads': 2,
        })()
        g = {'num_bp_per_sra': 4}
        observed = {'max_workers': None, 'subexts': []}

        def fake_concat_fastq_files_for_subext(run_ids, subext, inext, output_dir, outfile_path):
            assert run_ids == ['SRR001', 'SRR002']
            assert inext == '.amalgkit.fastq.gz'
            observed['subexts'].append(subext)
            with open(outfile_path, 'wt') as out:
                out.write('dummy-{}\n'.format(subext))

        def fake_run_tasks(task_items, task_fn, max_workers=1):
            observed['max_workers'] = max_workers
            results = {}
            for item in task_items:
                results[item] = task_fn(item)
            return results, []

        monkeypatch.setattr('amalgkit.getfastq.concat_fastq_files_for_subext', fake_concat_fastq_files_for_subext)
        monkeypatch.setattr('amalgkit.getfastq.run_tasks_with_optional_threads', fake_run_tasks)
        monkeypatch.setattr(
            'amalgkit.getfastq.list_run_dir_files',
            lambda _work_dir: {
                'SRR001_1.amalgkit.fastq.gz',
                'SRR001_2.amalgkit.fastq.gz',
                'SRR002_1.amalgkit.fastq.gz',
                'SRR002_2.amalgkit.fastq.gz',
            },
        )

        concat_fastq(args, metadata, str(tmp_path), g)

        assert observed['max_workers'] == 2
        assert set(observed['subexts']) == {'_1', '_2'}
        assert (tmp_path / 'MERGED_1.amalgkit.fastq.gz').exists()
        assert (tmp_path / 'MERGED_2.amalgkit.fastq.gz').exists()


class TestDetectConcatInputFiles:
    def test_does_not_match_prefix_only_run_ids(self):
        output_files = {
            'SRR1.amalgkit.fastq.gz',
            'SRR10.amalgkit.fastq.gz',
            'SRR100_1.amalgkit.fastq.gz',
            'SRR100_2.amalgkit.fastq.gz',
        }
        detected = detect_concat_input_files(
            output_files=output_files,
            run_ids=['SRR1', 'SRR100'],
            inext='.amalgkit.fastq.gz',
        )
        assert detected == [
            'SRR1.amalgkit.fastq.gz',
            'SRR100_1.amalgkit.fastq.gz',
            'SRR100_2.amalgkit.fastq.gz',
        ]


class TestCollectValidRunIds:
    def test_filters_blank_nan_and_duplicates(self):
        values = ['SRR001', '', '   ', numpy.nan, None, 'SRR001', 'SRR002  ']
        assert collect_valid_run_ids(values, unique=True) == ['SRR001', 'SRR002']


# ---------------------------------------------------------------------------
# remove_experiment_without_run
# ---------------------------------------------------------------------------

class TestRemoveExperimentWithoutRun:
    def test_removes_empty_run(self):
        """Experiments without run IDs should be filtered out."""
        data = {
            'run': ['SRR001', '', 'SRR003'],
            'scientific_name': ['Sp1', 'Sp1', 'Sp1'],
            'exclusion': ['no', 'no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = remove_experiment_without_run(m)
        assert m.df.shape[0] == 2
        assert '' not in m.df['run'].values

    def test_no_removal_needed(self):
        """All experiments have runs, nothing removed."""
        data = {
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = remove_experiment_without_run(m)
        assert m.df.shape[0] == 2

    def test_removes_nan_and_whitespace_runs(self):
        data = {
            'run': ['SRR001', numpy.nan, None, '   ', 'SRR002  '],
            'scientific_name': ['Sp1', 'Sp1', 'Sp1', 'Sp1', 'Sp1'],
            'exclusion': ['no', 'no', 'no', 'no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = remove_experiment_without_run(m)
        assert m.df['run'].tolist() == ['SRR001', 'SRR002']

    def test_handles_all_nan_run_column_without_typeerror(self):
        data = {
            'run': [numpy.nan, numpy.nan],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = remove_experiment_without_run(m)
        assert m.df.shape[0] == 0


# ---------------------------------------------------------------------------
# check_metadata_validity (issues #96, #110: empty total_bases/total_spots)
# ---------------------------------------------------------------------------

class TestCheckMetadataValidity:
    def test_rejects_empty_metadata_table(self):
        data = {
            'run': [],
            'scientific_name': [],
            'exclusion': [],
            'total_bases': [],
            'total_spots': [],
            'spot_length': [],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='No SRA entry found'):
            check_metadata_validity(m)

    def test_rejects_missing_run_ids(self):
        data = {
            'run': ['SRR001', '   '],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
            'total_bases': [1000000, 1000000],
            'total_spots': [5000, 5000],
            'spot_length': [200, 200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='Missing Run ID\\(s\\) were detected'):
            check_metadata_validity(m)

    def test_fills_missing_total_bases(self):
        """Issue #96: Empty total_bases should be filled with placeholder 999999999999."""
        data = {
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
            'total_bases': [1000000, 0],
            'total_spots': [5000, 5000],
            'spot_length': [200, 200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        assert m.df.loc[m.df['run'] == 'SRR002', 'total_bases'].values[0] == 999999999999

    def test_fills_missing_total_spots(self):
        """Issue #110: Empty total_spots should be filled with placeholder."""
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000000],
            'total_spots': [0],
            'spot_length': [200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        # total_spots was 0, should be filled with a value
        assert m.df.loc[0, 'total_spots'] > 0

    def test_valid_metadata_unchanged(self):
        """Valid metadata should pass through unchanged."""
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [2000000000],
            'total_spots': [10000000],
            'spot_length': [200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        assert m.df.loc[0, 'total_bases'] == 2000000000
        assert m.df.loc[0, 'total_spots'] == 10000000

    def test_uses_placeholder_when_spot_length_is_zero(self):
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000000],
            'total_spots': [0],
            'spot_length': [0],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        assert int(m.df.loc[0, 'total_spots']) == 999999999999

    def test_fills_non_numeric_total_bases(self):
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': ['abc'],
            'total_spots': [100],
            'spot_length': [100],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        assert int(m.df.loc[0, 'total_bases']) == 999999999999

    def test_estimates_total_spots_from_non_numeric_value(self):
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000],
            'total_spots': ['abc'],
            'spot_length': [100],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        assert int(m.df.loc[0, 'total_spots']) == 10

    def test_rejects_duplicate_run_ids(self):
        data = {
            'run': ['SRR001', ' SRR001 '],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
            'total_bases': [1000000, 1000000],
            'total_spots': [5000, 5000],
            'spot_length': [200, 200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='Duplicate Run ID\\(s\\) were detected'):
            check_metadata_validity(m)

    def test_rejects_missing_required_columns(self):
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000000],
            'total_spots': [5000],
            'spot_length': [200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.df = m.df.drop(columns=['total_spots'])

        with pytest.raises(ValueError, match='Missing required metadata column\\(s\\) for getfastq: total_spots'):
            check_metadata_validity(m)


# ---------------------------------------------------------------------------
# initialize_global_params (wiki: calculates per-SRA bp targets)
# ---------------------------------------------------------------------------

class TestInitializeGlobalParams:
    def test_basic_calculation(self):
        """Calculates max_bp, num_sra, num_bp_per_sra, total_sra_bp."""
        class Args:
            max_bp = '1,000,000'
        data = {
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
            'total_bases': [500000, 500000],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        g = initialize_global_params(Args(), m)
        assert g['max_bp'] == 1000000
        assert g['num_sra'] == 2
        assert g['num_bp_per_sra'] == 500000
        assert g['total_sra_bp'] == 1000000

    def test_comma_removal(self):
        """Commas in max_bp string should be removed."""
        class Args:
            max_bp = '999,999,999,999,999'
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        g = initialize_global_params(Args(), m)
        assert g['max_bp'] == 999999999999999

    def test_rejects_non_positive_max_bp(self):
        class Args:
            max_bp = '0'

        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='--max_bp must be > 0'):
            initialize_global_params(Args(), m)

    def test_rejects_max_bp_too_small_for_num_sra(self):
        class Args:
            max_bp = '1'

        data = {
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
            'total_bases': [500, 500],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='--max_bp \\(1\\) is too small for 2 SRA runs'):
            initialize_global_params(Args(), m)


# ---------------------------------------------------------------------------
# rename_fastq (renames fastq files by extension)
# ---------------------------------------------------------------------------

class TestRenameFastq:
    def test_rename_single(self, tmp_path):
        """Renames single-end fastq file."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        (tmp_path / 'SRR001.fastq.gz').write_text('data')
        rename_fastq(sra_stat, str(tmp_path), '.fastq.gz', '.amalgkit.fastq.gz')
        assert os.path.exists(str(tmp_path / 'SRR001.amalgkit.fastq.gz'))
        assert not os.path.exists(str(tmp_path / 'SRR001.fastq.gz'))
        assert sra_stat['current_ext'] == '.amalgkit.fastq.gz'

    def test_rename_paired(self, tmp_path):
        """Renames paired-end fastq files."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        (tmp_path / 'SRR001_1.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.fastq.gz').write_text('data')
        rename_fastq(sra_stat, str(tmp_path), '.fastq.gz', '.amalgkit.fastq.gz')
        assert os.path.exists(str(tmp_path / 'SRR001_1.amalgkit.fastq.gz'))
        assert os.path.exists(str(tmp_path / 'SRR001_2.amalgkit.fastq.gz'))
        assert sra_stat['current_ext'] == '.amalgkit.fastq.gz'

    def test_rejects_directory_as_source_path(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        (tmp_path / 'SRR001.fastq.gz').mkdir()

        with pytest.raises(IsADirectoryError, match='Intermediate path exists but is not a file'):
            rename_fastq(sra_stat, str(tmp_path), '.fastq.gz', '.amalgkit.fastq.gz')


# ---------------------------------------------------------------------------
# remove_old_intermediate_files (removes old files but keeps .sra)
# ---------------------------------------------------------------------------

class TestRemoveOldIntermediateFiles:
    def test_removes_intermediate_keeps_sra(self, tmp_path):
        """Removes intermediate files but keeps .sra files."""
        (tmp_path / 'SRR001.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_1.fastp.fastq.gz').write_text('data')
        (tmp_path / 'SRR001.sra').write_text('data')
        remove_old_intermediate_files('SRR001', str(tmp_path))
        assert not os.path.exists(str(tmp_path / 'SRR001.fastq.gz'))
        assert not os.path.exists(str(tmp_path / 'SRR001_1.fastp.fastq.gz'))
        assert os.path.exists(str(tmp_path / 'SRR001.sra'))

    def test_no_files_to_remove(self, tmp_path):
        """No matching intermediate files -> no error."""
        (tmp_path / 'other_file.txt').write_text('data')
        remove_old_intermediate_files('SRR001', str(tmp_path))

    def test_ignores_similar_prefix_files(self, tmp_path):
        (tmp_path / 'SRR001.fastq.gz').write_text('data')
        (tmp_path / 'SRR0010.fastq.gz').write_text('data')
        remove_old_intermediate_files('SRR001', str(tmp_path))
        assert not os.path.exists(str(tmp_path / 'SRR001.fastq.gz'))
        assert os.path.exists(str(tmp_path / 'SRR0010.fastq.gz'))


# ---------------------------------------------------------------------------
# remove_intermediate_files (removes single/paired intermediate files)
# ---------------------------------------------------------------------------

class TestRemoveIntermediateFiles:
    def test_removes_single(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        (tmp_path / 'SRR001.fastq.gz').write_text('data')
        remove_intermediate_files(sra_stat, '.fastq.gz', str(tmp_path))
        assert not os.path.exists(str(tmp_path / 'SRR001.fastq.gz'))

    def test_removes_paired(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        (tmp_path / 'SRR001_1.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.fastq.gz').write_text('data')
        remove_intermediate_files(sra_stat, '.fastq.gz', str(tmp_path))
        assert not os.path.exists(str(tmp_path / 'SRR001_1.fastq.gz'))
        assert not os.path.exists(str(tmp_path / 'SRR001_2.fastq.gz'))

    def test_missing_files_no_error(self, tmp_path):
        """Missing files should not raise."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        remove_intermediate_files(sra_stat, '.fastq.gz', str(tmp_path))

    def test_raises_when_intermediate_path_is_directory(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        (tmp_path / 'SRR001.fastq.gz').mkdir()

        with pytest.raises(IsADirectoryError, match='Intermediate path exists but is not a file'):
            remove_intermediate_files(sra_stat, '.fastq.gz', str(tmp_path))


# ---------------------------------------------------------------------------
# is_getfastq_output_present (checks for output files)
# ---------------------------------------------------------------------------

class TestIsGetfastqOutputPresent:
    def test_output_present_single(self, tmp_path):
        """Single-end output detected."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('data')
        observed = is_getfastq_output_present(sra_stat)
        assert observed
        assert isinstance(observed, bool)

    def test_output_present_paired(self, tmp_path):
        """Paired-end output detected."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.amalgkit.fastq.gz').write_text('data')
        assert is_getfastq_output_present(sra_stat)

    def test_safely_removed_counts(self, tmp_path):
        """Safely removed output counts as present."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.amalgkit.fastq.gz.safely_removed').write_text('')
        assert is_getfastq_output_present(sra_stat)

    def test_output_missing(self, tmp_path):
        """No output files -> False."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        assert not is_getfastq_output_present(sra_stat)

    def test_uses_prefetched_file_set(self, tmp_path, monkeypatch):
        """When files set is provided, no directory re-scan is needed."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        files = {'SRR001.amalgkit.fastq.gz'}

        def fail_if_called(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when files are provided.')

        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_if_called)
        assert is_getfastq_output_present(sra_stat, files=files)


# ---------------------------------------------------------------------------
# initialize_columns (initializes tracking columns for getfastq metrics)
# ---------------------------------------------------------------------------

class TestInitializeColumns:
    def test_initializes_columns(self):
        """Adds all tracking columns to metadata DataFrame."""
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_spots': [10000000],
            'total_bases': [2000000000],
            'size': [500000000],
            'nominal_length': [200],
            'nominal_sdev': [0],
            'spot_length': [200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        g = {'num_bp_per_sra': 1000000}
        m = initialize_columns(m, g)
        assert 'bp_amalgkit' in m.df.columns
        assert 'bp_dumped' in m.df.columns
        assert 'rate_obtained' in m.df.columns
        assert 'layout_amalgkit' in m.df.columns
        assert 'fastp_duplication_rate' in m.df.columns
        assert 'fastp_insert_size_peak' in m.df.columns
        assert 'num_rrna_in' in m.df.columns
        assert 'num_rrna_out' in m.df.columns
        assert 'bp_rrna_in' in m.df.columns
        assert 'bp_rrna_out' in m.df.columns
        assert 'spot_length_amalgkit' in m.df.columns
        assert m.df.loc[0, 'bp_until_target_size'] == 1000000
        assert m.df.loc[0, 'bp_dumped'] == 0
        assert m.df.loc[0, 'spot_length_amalgkit'] == 0
        assert m.df.loc[0, 'layout_amalgkit'] == ''
        assert numpy.isnan(m.df.loc[0, 'fastp_duplication_rate'])
        assert numpy.isnan(m.df.loc[0, 'fastp_insert_size_peak'])


# ---------------------------------------------------------------------------
# is_2nd_round_needed (triggers compensatory extraction based on tolerated loss)
# ---------------------------------------------------------------------------

class TestIs2ndRoundNeeded:
    def test_zero_obtained_requires_second_round(self):
        assert is_2nd_round_needed(rate_obtained_1st=0.0, tol=1.0)

    def test_half_obtained_requires_second_round_at_default_tol(self):
        assert is_2nd_round_needed(rate_obtained_1st=0.5, tol=1.0)

    def test_within_default_tolerance_skips_second_round(self):
        assert not is_2nd_round_needed(rate_obtained_1st=0.99, tol=1.0)

    def test_custom_tolerance_boundary(self):
        assert not is_2nd_round_needed(rate_obtained_1st=0.95, tol=5.0)
        assert is_2nd_round_needed(rate_obtained_1st=0.9499, tol=5.0)


class TestCalc2ndRanges:
    def test_non_zero_index_uses_label_based_access(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'bp_until_target_size': [300.0, 600.0, 900.0],
            'rate_obtained': [0.5, numpy.nan, 1.0],
            'spot_length_amalgkit': [100.0, 200.0, 300.0],
            'total_spots': [1000, 1000, 1000],
            'spot_end_1st': [100, 200, 300],
        }))
        metadata.df.index = [10, 20, 30]

        out = calc_2nd_ranges(metadata)

        assert list(out.df.index) == [10, 20, 30]
        assert out.df.loc[10, 'spot_start_2nd'] == 101
        assert out.df.loc[20, 'spot_start_2nd'] == 201
        assert out.df.loc[30, 'spot_start_2nd'] == 301
        assert out.df.loc[10, 'spot_end_2nd'] == 108
        assert out.df.loc[20, 'spot_end_2nd'] == 205
        assert out.df.loc[30, 'spot_end_2nd'] == 305

    def test_redistributes_missing_reads_when_other_sra_has_capacity(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'bp_until_target_size': [50.0, 0.0],
            'rate_obtained': [1.0, 1.0],
            'spot_length_amalgkit': [1.0, 1.0],
            'total_spots': [10, 100],
            'spot_end_1st': [0, 0],
        }))

        out = calc_2nd_ranges(metadata)

        assert out.df.loc[0, 'spot_start_2nd'] == 1
        assert out.df.loc[1, 'spot_start_2nd'] == 1
        assert out.df.loc[0, 'spot_end_2nd'] == 10
        assert out.df.loc[1, 'spot_end_2nd'] == 44

    def test_zero_rate_obtained_uses_base_target_without_inf(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'bp_until_target_size': [100.0],
            'rate_obtained': [0.0],
            'spot_length_amalgkit': [10.0],
            'total_spots': [1000],
            'spot_end_1st': [0],
        }))

        out = calc_2nd_ranges(metadata)

        assert out.df.loc[0, 'spot_start_2nd'] == 1
        assert out.df.loc[0, 'spot_end_2nd'] == 12

    def test_never_emits_end_before_start_when_total_spots_is_smaller(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'bp_until_target_size': [100.0],
            'rate_obtained': [1.0],
            'spot_length_amalgkit': [10.0],
            'total_spots': [5],
            'spot_end_1st': [10],
        }))

        out = calc_2nd_ranges(metadata)

        assert out.df.loc[0, 'spot_start_2nd'] == 11
        assert out.df.loc[0, 'spot_end_2nd'] == 11

    def test_uses_spot_length_when_spot_length_amalgkit_is_missing(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'bp_until_target_size': [100.0],
            'rate_obtained': [1.0],
            'spot_length': [10.0],
            'total_spots': [1000],
            'spot_end_1st': [0],
        }))

        out = calc_2nd_ranges(metadata)

        assert out.df.loc[0, 'spot_start_2nd'] == 1
        assert out.df.loc[0, 'spot_end_2nd'] == 12


class TestHasRemainingSpotsAfterFirstRound:
    def test_returns_true_when_any_spot_remains(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'total_spots': [10, 20],
            'spot_end_1st': [10, 19],
        }))
        assert has_remaining_spots_after_first_round(metadata)

    def test_returns_false_when_all_spots_are_exhausted(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'total_spots': [10, 20],
            'spot_end_1st': [10, 20],
        }))
        assert not has_remaining_spots_after_first_round(metadata)

    def test_returns_true_when_required_columns_are_missing(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'total_spots': [10],
        }))
        assert has_remaining_spots_after_first_round(metadata)


class TestFastpMetrics:
    def test_parse_fastp_metrics_extracts_duplication_and_insert_size(self):
        stderr_txt = '\n'.join([
            'Filtering result:',
            'reads passed filter: 274442354',
            'Duplication rate: 22.9565%',
            'Insert size peak (evaluated by paired-end reads): 138',
            'JSON report: /dev/null',
        ])
        duplication_rate, insert_size_peak = parse_fastp_metrics(stderr_txt)
        assert duplication_rate == pytest.approx(22.9565)
        assert insert_size_peak == pytest.approx(138.0)

    def test_parse_fastp_metrics_returns_nan_when_not_present(self):
        duplication_rate, insert_size_peak = parse_fastp_metrics('no matching lines')
        assert numpy.isnan(duplication_rate)
        assert numpy.isnan(insert_size_peak)

    def test_parse_fastp_summary_counts_raises_on_truncated_section(self):
        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            # missing total bases line
        ])
        with pytest.raises(RuntimeError, match='Unexpected fastp stderr format'):
            parse_fastp_summary_counts(stderr_txt)

    def test_update_fastp_metrics_weighted_average(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'num_fastp_in': [100],
            'fastp_duplication_rate': [10.0],
            'fastp_insert_size_peak': [200.0],
        }))
        update_fastp_metrics(
            metadata=metadata,
            ind_sra=0,
            current_num_in=300,
            duplication_rate=20.0,
            insert_size_peak=100.0,
        )
        assert metadata.df.loc[0, 'fastp_duplication_rate'] == pytest.approx(17.5)
        assert metadata.df.loc[0, 'fastp_insert_size_peak'] == pytest.approx(125.0)

    def test_update_fastp_metrics_sets_initial_values(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'num_fastp_in': [0],
            'fastp_duplication_rate': [numpy.nan],
            'fastp_insert_size_peak': [numpy.nan],
        }))
        update_fastp_metrics(
            metadata=metadata,
            ind_sra=0,
            current_num_in=50,
            duplication_rate=33.0,
            insert_size_peak=222.0,
        )
        assert metadata.df.loc[0, 'fastp_duplication_rate'] == pytest.approx(33.0)
        assert metadata.df.loc[0, 'fastp_insert_size_peak'] == pytest.approx(222.0)

    def test_write_fastp_stats_writes_tsv(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'fastp_duplication_rate': [22.9565],
            'fastp_insert_size_peak': [138.0],
            'num_fastp_in': [1000],
            'num_fastp_out': [900],
            'bp_fastp_in': [100000],
            'bp_fastp_out': [90000],
        }))
        sra_stat = {'sra_id': 'SRR001'}
        write_fastp_stats(sra_stat=sra_stat, metadata=metadata, output_dir=str(tmp_path))
        out_path = tmp_path / 'fastp_stats.tsv'
        assert out_path.exists()
        out = pandas.read_csv(out_path, sep='\t')
        assert out.loc[0, 'run'] == 'SRR001'
        assert out.loc[0, 'fastp_duplication_rate'] == pytest.approx(22.9565)
        assert out.loc[0, 'fastp_insert_size_peak'] == pytest.approx(138.0)
        assert out.loc[0, 'num_fastp_in'] == 1000
        assert out.loc[0, 'num_fastp_out'] == 900
        assert out.loc[0, 'bp_fastp_in'] == 100000
        assert out.loc[0, 'bp_fastp_out'] == 90000

    def test_write_getfastq_stats_writes_tsv(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'num_dumped': [10],
            'num_rejected': [3],
            'num_written': [7],
            'num_fastp_in': [7],
            'num_fastp_out': [6],
            'num_rrna_in': [6],
            'num_rrna_out': [5],
            'bp_dumped': [1000],
            'bp_rejected': [300],
            'bp_written': [700],
            'bp_fastp_in': [700],
            'bp_fastp_out': [600],
            'bp_rrna_in': [600],
            'bp_rrna_out': [500],
            'bp_discarded': [500],
            'fastp_duplication_rate': [12.0],
            'fastp_insert_size_peak': [250.0],
        }))
        sra_stat = {'sra_id': 'SRR001'}
        write_getfastq_stats(sra_stat=sra_stat, metadata=metadata, output_dir=str(tmp_path))
        out_path = tmp_path / 'getfastq_stats.tsv'
        assert out_path.exists()
        out = pandas.read_csv(out_path, sep='\t')
        assert out.loc[0, 'run'] == 'SRR001'
        assert out.loc[0, 'num_dumped'] == 10
        assert out.loc[0, 'num_rejected'] == 3
        assert out.loc[0, 'bp_rejected'] == 300
        assert out.loc[0, 'bp_discarded'] == 500
        assert out.loc[0, 'fastp_duplication_rate'] == pytest.approx(12.0)
        assert out.loc[0, 'fastp_insert_size_peak'] == pytest.approx(250.0)


class TestRunFastp:
    @staticmethod
    def _build_metadata():
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'num_fastp_in': [0],
            'num_fastp_out': [0],
            'bp_fastp_in': [0],
            'bp_fastp_out': [0],
            'fastp_duplication_rate': [numpy.nan],
            'fastp_insert_size_peak': [numpy.nan],
        }))

    @staticmethod
    def _materialize_fastp_outputs(cmd):
        for flag in ['--out1', '--out2']:
            if flag not in cmd:
                continue
            out_path = cmd[cmd.index(flag) + 1]
            with open(out_path, 'wb') as fout:
                fout.write(b'@r0\nAAAA\n+\nIIII\n')

    def test_uses_configured_fastp_exe_and_shlex_parsing(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        captured = {}

        class Args:
            threads = 4
            min_read_length = 25
            fastp_option = '--adapter_sequence "A B"'
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8')

        def fake_run(cmd, stdout=None, stderr=None):
            captured['cmd'] = cmd
            self._materialize_fastp_outputs(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

        assert captured['cmd'][0] == 'fastp-custom'
        assert '--adapter_sequence' in captured['cmd']
        assert 'A B' in captured['cmd']
        assert metadata.df.loc[0, 'num_fastp_in'] == 10
        assert metadata.df.loc[0, 'num_fastp_out'] == 8
        assert metadata.df.loc[0, 'bp_fastp_in'] == 100
        assert metadata.df.loc[0, 'bp_fastp_out'] == 80

    def test_raises_when_fastp_exits_nonzero(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 1, stdout=b'', stderr=b'fastp failed')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(RuntimeError, match='fastp did not finish safely'):
            run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

    def test_raises_when_fastp_stderr_is_truncated(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        def fake_run(cmd, stdout=None, stderr=None):
            self._materialize_fastp_outputs(cmd)
            bad_stderr = '\n'.join([
                ' before filtering:',
                'total reads: 10',
            ]).encode('utf8')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=bad_stderr)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(RuntimeError, match='Unexpected fastp stderr format'):
            run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

    def test_uses_cached_extension_without_redetection(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single', 'current_ext': '.fastq.gz'}
        captured = {}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('get_newest_intermediate_file_extension should not be called with cached current_ext.')

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            fail_if_called,
        )

        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8')

        def fake_run(cmd, stdout=None, stderr=None):
            captured['cmd'] = cmd
            self._materialize_fastp_outputs(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

        assert '--in1' in captured['cmd']
        in1_arg = captured['cmd'][captured['cmd'].index('--in1') + 1]
        assert in1_arg.endswith('SRR001.fastq.gz')
        assert sra_stat['current_ext'] == '.fastp.fastq.gz'

    def test_uses_prefetched_files_for_extension_detection_without_rescan(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        captured = {'seen_files': None}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        def fake_get_newest_intermediate_file_extension(sra_stat, work_dir, files=None):
            captured['seen_files'] = set(files)
            return '.fastq.gz'

        def fail_if_listed(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when files are provided.')

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            fake_get_newest_intermediate_file_extension,
        )
        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_if_listed)

        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8')

        def fake_run(cmd, stdout=None, stderr=None):
            self._materialize_fastp_outputs(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        metadata, files_out = run_fastp(
            sra_stat=sra_stat,
            args=Args(),
            output_dir=str(tmp_path),
            metadata=metadata,
            files={'SRR001.fastq.gz'},
            return_files=True,
        )

        assert captured['seen_files'] == {'SRR001.fastq.gz'}
        assert 'SRR001.fastp.fastq.gz' in files_out
        assert metadata.df.loc[0, 'num_fastp_out'] == 8

    def test_tolerates_non_utf8_fastp_stderr(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = True
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        stderr_bytes = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8') + b'\xff'

        def fake_run(cmd, stdout=None, stderr=None):
            self._materialize_fastp_outputs(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'\xff', stderr=stderr_bytes)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

        assert metadata.df.loc[0, 'num_fastp_in'] == 10
        assert metadata.df.loc[0, 'num_fastp_out'] == 8

    def test_raises_when_fastp_output_file_is_missing(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8')

        def fake_run(cmd, stdout=None, stderr=None):
            # Intentionally do not materialize --out1.
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(FileNotFoundError, match='fastp output file was not generated'):
            run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

    def test_raises_when_fastp_output_path_is_directory(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8')

        def fake_run(cmd, stdout=None, stderr=None):
            out1 = cmd[cmd.index('--out1') + 1]
            os.makedirs(out1, exist_ok=True)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(IsADirectoryError, match='fastp output path exists but is not a file'):
            run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)


class TestIdenticalPairedReads:
    @staticmethod
    def _write_fastq_gz(path, seqs):
        with gzip.open(path, 'wt') as out:
            for i, seq in enumerate(seqs, start=1):
                out.write('@read{}\n{}\n+\n{}\n'.format(i, seq, 'I' * len(seq)))

    @staticmethod
    def _write_fastq_plain(path, seqs):
        with open(path, 'wt') as out:
            for i, seq in enumerate(seqs, start=1):
                out.write('@read{}\n{}\n+\n{}\n'.format(i, seq, 'I' * len(seq)))

    def test_get_identical_paired_ratio(self, tmp_path):
        read1 = tmp_path / 'read1.fastq.gz'
        read2 = tmp_path / 'read2.fastq.gz'
        self._write_fastq_gz(read1, ['AAAA', 'CCCC', 'GGGG'])
        self._write_fastq_gz(read2, ['AAAA', 'TTTT', 'GGGG'])
        ratio, num_checked, read_length = get_identical_paired_ratio(str(read1), str(read2), num_checked_reads=3)
        assert ratio == pytest.approx(2 / 3)
        assert num_checked == 3
        assert read_length == 4

    def test_get_identical_paired_ratio_plain_fastq(self, tmp_path):
        read1 = tmp_path / 'read1.fastq'
        read2 = tmp_path / 'read2.fastq'
        self._write_fastq_plain(read1, ['AAAA', 'CCCC', 'GGGG'])
        self._write_fastq_plain(read2, ['AAAA', 'TTTT', 'GGGG'])
        ratio, num_checked, read_length = get_identical_paired_ratio(str(read1), str(read2), num_checked_reads=3)
        assert ratio == pytest.approx(2 / 3)
        assert num_checked == 3
        assert read_length == 4

    def test_maybe_treat_paired_as_single_converts_files(self, tmp_path):
        sra_id = 'SRR001'
        read1 = tmp_path / '{}_1.fastq.gz'.format(sra_id)
        read2 = tmp_path / '{}_2.fastq.gz'.format(sra_id)
        self._write_fastq_gz(read1, ['AAAA', 'CCCC', 'GGGG'])
        self._write_fastq_gz(read2, ['AAAA', 'CCCC', 'GGGG'])
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['paired'],
            'layout_amalgkit': ['paired'],
            'spot_length': [8],
            'total_spots': [3],
            'total_bases': [24],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 8,
        }
        metadata, sra_stat = maybe_treat_paired_as_single(
            sra_stat=sra_stat,
            metadata=metadata,
            work_dir=str(tmp_path),
            threshold=0.99,
            num_checked_reads=3,
        )
        assert sra_stat['layout'] == 'single'
        assert os.path.exists(str(tmp_path / '{}.fastq.gz'.format(sra_id)))
        assert not os.path.exists(str(read1))
        assert not os.path.exists(str(read2))
        assert metadata.df.loc[0, 'layout_amalgkit'] == 'single'
        assert metadata.df.loc[0, 'spot_length'] == 4

    def test_maybe_treat_paired_as_single_keeps_paired_when_ratio_low(self, tmp_path):
        sra_id = 'SRR001'
        read1 = tmp_path / '{}_1.fastq.gz'.format(sra_id)
        read2 = tmp_path / '{}_2.fastq.gz'.format(sra_id)
        self._write_fastq_gz(read1, ['AAAA', 'CCCC', 'GGGG'])
        self._write_fastq_gz(read2, ['TTTT', 'CCCC', 'AAAA'])
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['paired'],
            'layout_amalgkit': ['paired'],
            'spot_length': [8],
            'total_spots': [3],
            'total_bases': [24],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 8,
        }
        metadata, sra_stat = maybe_treat_paired_as_single(
            sra_stat=sra_stat,
            metadata=metadata,
            work_dir=str(tmp_path),
            threshold=0.99,
            num_checked_reads=3,
        )
        assert sra_stat['layout'] == 'paired'
        assert os.path.exists(str(read1))
        assert os.path.exists(str(read2))
        assert metadata.df.loc[0, 'layout_amalgkit'] == 'paired'

    def test_maybe_treat_paired_as_single_uses_prefetched_files_without_rescan(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        read1 = tmp_path / '{}_1.fastq.gz'.format(sra_id)
        read2 = tmp_path / '{}_2.fastq.gz'.format(sra_id)
        self._write_fastq_gz(read1, ['AAAA', 'CCCC'])
        self._write_fastq_gz(read2, ['AAAA', 'CCCC'])
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['paired'],
            'layout_amalgkit': ['paired'],
            'spot_length': [8],
            'total_spots': [2],
            'total_bases': [16],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 8,
        }

        def fail_if_called(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when files are provided.')

        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_if_called)
        metadata, sra_stat = maybe_treat_paired_as_single(
            sra_stat=sra_stat,
            metadata=metadata,
            work_dir=str(tmp_path),
            threshold=0.99,
            num_checked_reads=2,
            files={'{}_1.fastq.gz'.format(sra_id), '{}_2.fastq.gz'.format(sra_id)},
        )
        assert sra_stat['layout'] == 'single'


class TestWrittenSpotEstimation:
    @staticmethod
    def _write_fastq(path, seqs):
        with open(path, 'wt') as out:
            for i, seq in enumerate(seqs):
                out.write('@r{}\n{}\n+\n{}\n'.format(i, seq, 'I' * len(seq)))

    def test_uses_prefetched_files_without_rescan(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        fastq_path = tmp_path / '{}.fastq'.format(sra_id)
        self._write_fastq(fastq_path, ['AAAA', 'CCCC', 'GGGG'])
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
        }

        def fail_if_called(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when files are provided.')

        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_if_called)
        observed = estimate_num_written_spots_from_fastq(sra_stat, files={'{}.fastq'.format(sra_id)})
        assert observed == 3


class TestSraRecovery:
    @staticmethod
    def _metadata_for_extraction(sra_id):
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'num_dumped': [0],
            'num_rejected': [0],
            'num_written': [0],
            'bp_dumped': [0],
            'bp_rejected': [0],
            'bp_written': [0],
            'layout_amalgkit': ['paired'],
        }))

    @staticmethod
    def _args_for_fasterq_dump():
        class Args:
            threads = 2
            min_read_length = 25
            dump_print = False
            fasterq_dump_exe = 'fasterq-dump'
            fasterq_size_check = True
            fasterq_disk_limit = None
            fasterq_disk_limit_tmp = None
        return Args()

    @pytest.fixture(autouse=True)
    def _disable_fasterq_output_validation(self, monkeypatch):
        monkeypatch.setattr('amalgkit.getfastq.ensure_fasterq_output_files_exist', lambda **_kwargs: None)

    def test_remove_sra_files_deletes_matching_sra_files(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001.sra').write_text('a')
        (sra_dir / 'SRR001.sra.vdbcache').write_text('b')
        (sra_dir / 'other.txt').write_text('keep')

        remove_sra_files(metadata, str(tmp_path))

        assert not (sra_dir / 'SRR001.sra').exists()
        assert not (sra_dir / 'SRR001.sra.vdbcache').exists()
        assert (sra_dir / 'other.txt').exists()

    def test_remove_sra_files_avoids_root_listdir_scan(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001.sra').write_text('a')

        def fail_if_listdir_called(_path):
            raise AssertionError('remove_sra_files should not call os.listdir on getfastq root.')

        monkeypatch.setattr('amalgkit.getfastq.os.listdir', fail_if_listdir_called)
        remove_sra_files(metadata, str(tmp_path))

        assert not (sra_dir / 'SRR001.sra').exists()

    def test_remove_sra_files_ignores_non_directory_entries(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        getfastq_root = tmp_path / 'getfastq'
        getfastq_root.mkdir(parents=True)
        (getfastq_root / 'SRR001').write_text('not a directory')

        remove_sra_files(metadata, str(tmp_path))

        assert (getfastq_root / 'SRR001').exists()

    def test_remove_sra_files_handles_missing_getfastq_root(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))

        remove_sra_files(metadata, str(tmp_path))

    def test_remove_sra_files_skips_missing_run_ids(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', numpy.nan, ''],
            'scientific_name': ['Sp1', 'Sp1', 'Sp1'],
            'exclusion': ['no', 'no', 'no'],
        }))
        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001.sra').write_text('a')

        remove_sra_files(metadata, str(tmp_path))

        assert not (sra_dir / 'SRR001.sra').exists()

    def test_remove_sra_files_keeps_non_artifact_suffix(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001.sra').write_text('a')
        (sra_dir / 'SRR001.sra.vdbcache').write_text('b')
        (sra_dir / 'SRR001.sra2').write_text('keep')

        remove_sra_files(metadata, str(tmp_path))

        assert not (sra_dir / 'SRR001.sra').exists()
        assert not (sra_dir / 'SRR001.sra.vdbcache').exists()
        assert (sra_dir / 'SRR001.sra2').exists()

    def test_remove_sra_path_file_and_directory(self, tmp_path):
        file_path = tmp_path / 'SRR001.sra'
        file_path.write_text('dummy')
        remove_sra_path(str(file_path))
        assert not file_path.exists()

        dir_path = tmp_path / 'SRR002.sra'
        (dir_path / 'tbl').mkdir(parents=True)
        (dir_path / 'tbl' / 'x').write_text('dummy')
        remove_sra_path(str(dir_path))
        assert not dir_path.exists()

    def test_run_fasterq_dump_retries_once_after_redownload(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        sra_path = tmp_path / '{}.sra'.format(sra_id)
        sra_path.mkdir()
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        run_calls = {'count': 0}

        def fake_run(cmd, stdout=None, stderr=None):
            run_calls['count'] += 1
            if run_calls['count'] == 1:
                return subprocess.CompletedProcess(cmd, 1, stdout=b'', stderr=b'fasterq-dump failed')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        redownload_calls = []

        def fake_download_sra(metadata, sra_stat, args, work_dir, overwrite=False):
            redownload_calls.append(overwrite)
            with open(os.path.join(work_dir, sra_stat['sra_id'] + '.sra'), 'w') as fh:
                fh.write('fresh')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.download_sra', fake_download_sra)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 4)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, sra_stat_out = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert run_calls['count'] == 2
        assert redownload_calls == [True]
        assert sra_path.is_file()
        assert metadata.df.loc[0, 'num_written'] == 4
        assert metadata.df.loc[0, 'num_dumped'] == 10
        assert metadata.df.loc[0, 'num_rejected'] == 6
        assert metadata.df.loc[0, 'bp_written'] == 400
        assert metadata.df.loc[0, 'bp_dumped'] == 1000
        assert metadata.df.loc[0, 'bp_rejected'] == 600
        assert sra_stat_out['layout'] == 'paired'

    def test_run_fasterq_dump_exits_when_retry_fails(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        sra_path = tmp_path / '{}.sra'.format(sra_id)
        sra_path.write_text('broken')
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        run_calls = {'count': 0}

        def fake_run(cmd, stdout=None, stderr=None):
            run_calls['count'] += 1
            return subprocess.CompletedProcess(cmd, 1, stdout=b'', stderr=b'fasterq-dump failed')

        redownload_calls = []

        def fake_download_sra(metadata, sra_stat, args, work_dir, overwrite=False):
            redownload_calls.append(overwrite)
            with open(os.path.join(work_dir, sra_stat['sra_id'] + '.sra'), 'w') as fh:
                fh.write('fresh')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.download_sra', fake_download_sra)

        with pytest.raises(SystemExit):
            run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert run_calls['count'] == 2
        assert redownload_calls == [True]

    def test_run_fasterq_dump_no_redownload_when_first_attempt_succeeds(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        run_calls = {'count': 0}

        def fake_run(cmd, stdout=None, stderr=None):
            run_calls['count'] += 1
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        def fail_download(*_args, **_kwargs):
            raise AssertionError('download_sra should not be called when first extraction succeeds.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.download_sra', fail_download)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 5)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert run_calls['count'] == 1
        assert metadata.df.loc[0, 'num_written'] == 5
        assert metadata.df.loc[0, 'num_dumped'] == 10
        assert metadata.df.loc[0, 'num_rejected'] == 5
        assert metadata.df.loc[0, 'bp_written'] == 500
        assert metadata.df.loc[0, 'bp_dumped'] == 1000

    def test_run_fasterq_dump_uses_reported_spots_without_recount_for_partial_range(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'spots written      : 7\n',
                stderr=b'',
            )

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when spots written is reported.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=2, end=8)

        assert metadata.df.loc[0, 'num_written'] == 7
        assert metadata.df.loc[0, 'num_dumped'] == 7
        assert metadata.df.loc[0, 'num_rejected'] == 0
        assert metadata.df.loc[0, 'bp_written'] == 700
        assert metadata.df.loc[0, 'bp_dumped'] == 700
        assert metadata.df.loc[0, 'bp_rejected'] == 0

    def test_run_fasterq_dump_ignores_written_total_line(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 0, stdout=b'written markers should be ignored', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 6667)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert metadata.df.loc[0, 'num_written'] == 6667
        assert metadata.df.loc[0, 'num_dumped'] == 6667
        assert metadata.df.loc[0, 'num_rejected'] == 0

    def test_run_fasterq_dump_infers_spots_from_reported_written_reads_for_paired_layout(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'',
                stderr=b'reads written   : 20\n',
            )

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when spots are inferred from reads written.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert metadata.df.loc[0, 'num_written'] == 10
        assert metadata.df.loc[0, 'num_dumped'] == 10
        assert metadata.df.loc[0, 'num_rejected'] == 0

    def test_run_fasterq_dump_falls_back_to_fastq_count_when_singletons_exist(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        (tmp_path / '{}.fastq'.format(sra_id)).write_text('@r0\nA\n+\nI\n')
        observed = {'estimate_called': False}

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'',
                stderr=b'reads written   : 21\n',
            )

        def fake_estimate(*_args, **_kwargs):
            observed['estimate_called'] = True
            return 3

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fake_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert observed['estimate_called']
        assert metadata.df.loc[0, 'num_written'] == 3

    def test_run_fasterq_dump_skips_pre_fastp_compression(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        args.fastp = True
        args.rrna_filter = False

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'',
                stderr=b'reads written   : 20\n',
            )

        def fail_compress(*_args, **_kwargs):
            raise AssertionError('compress_fasterq_output_files should be skipped when fastp is the next filter.')

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when reads written is usable.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', fail_compress)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, sra_stat_out = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert sra_stat_out['current_ext'] == '.fastq'
        assert metadata.df.loc[0, 'num_written'] == 10

    def test_run_fasterq_dump_skips_pre_rrna_compression(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        args.fastp = False
        args.rrna_filter = True

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'',
                stderr=b'reads written   : 20\n',
            )

        def fail_compress(*_args, **_kwargs):
            raise AssertionError('compress_fasterq_output_files should be skipped when a downstream filter exists.')

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when reads written is usable.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', fail_compress)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, sra_stat_out = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert sra_stat_out['current_ext'] == '.fastq'
        assert metadata.df.loc[0, 'num_written'] == 10

    def test_run_fasterq_dump_skips_trim_for_full_range(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
            'total_spot': 10,
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        def fail_trim(*_args, **_kwargs):
            raise AssertionError('trim_fasterq_output_files should be skipped for full-range extraction.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', fail_trim)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 10)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert metadata.df.loc[0, 'num_written'] == 10

    def test_run_fasterq_dump_skips_trim_for_partial_range(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
            'total_spot': 20,
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'spots written      : 7\n',
                stderr=b'',
            )

        def fail_trim(*_args, **_kwargs):
            raise AssertionError('trim_fasterq_output_files should not be called when fasterq-dump is spot-range limited.')

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when spots written is reported.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', fail_trim)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=2, end=8)

        assert metadata.df.loc[0, 'num_written'] == 7

    def test_run_fasterq_dump_trims_when_spot_range_flags_are_not_supported(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
            'total_spot': 20,
        }
        args = self._args_for_fasterq_dump()
        args.fastp = True
        args._fasterq_supports_spot_range = False
        observed = {'cmd': None, 'trim': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            return subprocess.CompletedProcess(cmd, 0, stdout=b'spots written : 100\n', stderr=b'')

        def fake_trim(*_args, **kwargs):
            observed['trim'] = (kwargs.get('start'), kwargs.get('end'), kwargs.get('compress_to_gz'))
            return {'': 0, '_1': 3, '_2': 3}, kwargs.get('file_state')

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when trim counts are available.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', fake_trim)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=2, end=8)

        cmd = observed['cmd']
        assert '-N' not in cmd
        assert '-X' not in cmd
        assert observed['trim'] == (2, 8, False)
        assert metadata.df.loc[0, 'num_written'] == 3
        assert metadata.df.loc[0, 'num_dumped'] == 7
        assert metadata.df.loc[0, 'num_rejected'] == 4

    def test_run_fasterq_dump_full_range_uses_reported_spots_without_recount(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
            'total_spot': 10,
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'spots written      : 7\n',
                stderr=b'',
            )

        def fail_trim(*_args, **_kwargs):
            raise AssertionError('trim_fasterq_output_files should be skipped for full-range extraction.')

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when spots written is reported.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', fail_trim)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert metadata.df.loc[0, 'num_written'] == 7
        assert metadata.df.loc[0, 'num_dumped'] == 10
        assert metadata.df.loc[0, 'num_rejected'] == 3

    def test_run_fasterq_dump_passes_size_check_and_disk_limits(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()
        args.fasterq_size_check = False
        args.fasterq_disk_limit = '10G'
        args.fasterq_disk_limit_tmp = '20G'
        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 1)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        run_fasterq_dump(sra_stat, args, metadata, start=1, end=1)

        cmd = observed['cmd']
        assert '--size-check' in cmd
        assert cmd[cmd.index('--size-check') + 1] == 'off'
        assert '--disk-limit' in cmd
        assert cmd[cmd.index('--disk-limit') + 1] == '10G'
        assert '--disk-limit-tmp' in cmd
        assert cmd[cmd.index('--disk-limit-tmp') + 1] == '20G'
        assert '-N' in cmd
        assert cmd[cmd.index('-N') + 1] == '1'
        assert '-X' in cmd
        assert cmd[cmd.index('-X') + 1] == '1'

    def test_run_fasterq_dump_passes_requested_spot_range(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()
        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 8)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        run_fasterq_dump(sra_stat, args, metadata, start=5, end=12)

        cmd = observed['cmd']
        assert '-N' in cmd
        assert cmd[cmd.index('-N') + 1] == '5'
        assert '-X' in cmd
        assert cmd[cmd.index('-X') + 1] == '12'

    def test_run_fasterq_dump_skips_spot_range_flags_for_full_range(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
            'total_spot': 10,
        }
        args = self._args_for_fasterq_dump()
        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'spots written : 10')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 10)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        cmd = observed['cmd']
        assert '-N' not in cmd
        assert '-X' not in cmd

    def test_sequence_extraction_reuses_run_files_without_rescan_when_fastp_disabled(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = False
            read_name = 'trinity'

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            if return_files:
                return metadata, sra_stat, {'{}.fastq.gz'.format(sra_id)}
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if file_state is not None:
                files = file_state.to_set()
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=files)
            if return_files:
                return metadata, sra_stat, set(files)
            return metadata, sra_stat

        captured = {'files': None}

        def fake_rename_reads(sra_stat, args, output_dir, files=None, file_state=None, return_file_state=False):
            if file_state is not None:
                files = file_state.to_set()
            captured['files'] = set(files)
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return RunFileState(work_dir=output_dir, files=files)
            return set(files)

        def fail_list_run_dir_files(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when run file snapshot is already available.')

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.rename_reads', fake_rename_reads)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_list_run_dir_files)

        metadata = sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert captured['files'] == {'{}.fastq.gz'.format(sra_id)}
        assert metadata.df.loc[0, 'num_written'] == 3
        assert metadata.df.loc[0, 'num_dumped'] == 3

    def test_sequence_extraction_fastp_trinity_reuses_run_files_without_rescan(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = True
            read_name = 'trinity'

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            if return_files:
                return metadata, sra_stat, {'{}.fastq.gz'.format(sra_id)}
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if file_state is not None:
                files = file_state.to_set()
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=files)
            if return_files:
                return metadata, sra_stat, set(files)
            return metadata, sra_stat

        def fake_run_fastp(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if file_state is not None:
                files = file_state.to_set()
            assert set(files) == {'{}.fastq.gz'.format(sra_id)}
            idx = 0
            metadata.df.at[idx, 'bp_fastp_out'] += 250
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.fastp.fastq.gz'.format(sra_id)})
            if return_files:
                return metadata, {'{}.fastp.fastq.gz'.format(sra_id)}
            return metadata

        captured = {'files': None}

        def fake_rename_reads(sra_stat, args, output_dir, files=None, file_state=None, return_file_state=False):
            if file_state is not None:
                files = file_state.to_set()
            captured['files'] = set(files)
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return RunFileState(work_dir=output_dir, files=files)
            return set(files)

        def fail_list_run_dir_files(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when run file snapshot is already available.')

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.rename_reads', fake_rename_reads)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_list_run_dir_files)

        metadata = sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert captured['files'] == {'{}.fastp.fastq.gz'.format(sra_id)}
        assert metadata.df.loc[0, 'bp_fastp_out'] == 250
        assert metadata.df.loc[0, 'bp_amalgkit'] == 250

    def test_sequence_extraction_uses_rrna_filtered_bp_for_target_size(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = True
            rrna_filter = True
            read_name = 'default'

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=file_state.to_set())
            return metadata, sra_stat

        def fake_run_fastp(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            idx = 0
            metadata.df.at[idx, 'bp_fastp_out'] += 250
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.fastp.fastq.gz'.format(sra_id)})
            return metadata

        def fake_run_sortmerna_rrna_filter(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            known_input_counts=None,
            return_files=False,
            return_file_state=False,
        ):
            idx = 0
            metadata.df.at[idx, 'bp_rrna_in'] += 250
            metadata.df.at[idx, 'bp_rrna_out'] += 180
            sra_stat['current_ext'] = '.rrna.fastq.gz'
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.rrna.fastq.gz'.format(sra_id)})
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.run_sortmerna_rrna_filter', fake_run_sortmerna_rrna_filter)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)

        metadata = sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert metadata.df.loc[0, 'bp_fastp_out'] == 250
        assert metadata.df.loc[0, 'bp_rrna_out'] == 180
        assert metadata.df.loc[0, 'bp_amalgkit'] == 180

    def test_sequence_extraction_passes_fastp_output_counts_to_rrna(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = True
            rrna_filter = True
            read_name = 'default'

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=file_state.to_set())
            return metadata, sra_stat

        def fake_run_fastp(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            idx = 0
            metadata.df.at[idx, 'num_fastp_out'] += 2
            metadata.df.at[idx, 'bp_fastp_out'] += 250
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.fastp.fastq.gz'.format(sra_id)})
            return metadata

        observed = {'known_input_counts': None}

        def fake_run_sortmerna_rrna_filter(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            known_input_counts=None,
            return_files=False,
            return_file_state=False,
        ):
            observed['known_input_counts'] = known_input_counts
            idx = 0
            metadata.df.at[idx, 'num_rrna_in'] += 2
            metadata.df.at[idx, 'num_rrna_out'] += 1
            metadata.df.at[idx, 'bp_rrna_in'] += 250
            metadata.df.at[idx, 'bp_rrna_out'] += 180
            sra_stat['current_ext'] = '.rrna.fastq.gz'
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.rrna.fastq.gz'.format(sra_id)})
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.run_sortmerna_rrna_filter', fake_run_sortmerna_rrna_filter)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)

        sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert observed['known_input_counts'] == {'num_spots': 2, 'bp_total': 250}

    def test_sequence_extraction_respects_rrna_first_order_for_bp_target(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = True
            rrna_filter = True
            filter_order = 'rrna_first'
            read_name = 'default'

        observed_order = []

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=file_state.to_set())
            return metadata, sra_stat

        def fake_run_sortmerna_rrna_filter(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            known_input_counts=None,
            return_files=False,
            return_file_state=False,
        ):
            observed_order.append('rrna')
            idx = 0
            metadata.df.at[idx, 'bp_rrna_in'] += 300
            metadata.df.at[idx, 'bp_rrna_out'] += 210
            sra_stat['current_ext'] = '.rrna.fastq.gz'
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.rrna.fastq.gz'.format(sra_id)})
            return metadata

        def fake_run_fastp(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            observed_order.append('fastp')
            idx = 0
            metadata.df.at[idx, 'bp_fastp_in'] += 210
            metadata.df.at[idx, 'bp_fastp_out'] += 180
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.fastp.fastq.gz'.format(sra_id)})
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.run_sortmerna_rrna_filter', fake_run_sortmerna_rrna_filter)
        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)

        metadata = sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert observed_order == ['rrna', 'fastp']
        assert metadata.df.loc[0, 'bp_rrna_out'] == 210
        assert metadata.df.loc[0, 'bp_fastp_out'] == 180
        assert metadata.df.loc[0, 'bp_amalgkit'] == 180


class TestRunFasterqDumpOutputValidation:
    @staticmethod
    def _metadata_for_extraction(sra_id):
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'num_dumped': [0],
            'num_rejected': [0],
            'num_written': [0],
            'bp_dumped': [0],
            'bp_rejected': [0],
            'bp_written': [0],
            'layout_amalgkit': ['single'],
        }))

    @staticmethod
    def _args_for_fasterq_dump():
        class Args:
            threads = 1
            min_read_length = 25
            dump_print = False
            fasterq_dump_exe = 'fasterq-dump'
            fasterq_size_check = True
            fasterq_disk_limit = None
            fasterq_disk_limit_tmp = None
        return Args()

    def test_raises_when_fasterq_dump_generates_no_fastq_files(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(FileNotFoundError, match='did not generate FASTQ files'):
            run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

    def test_raises_when_fasterq_dump_output_path_is_directory(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()
        (tmp_path / '{}.fastq'.format(sra_id)).mkdir()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(IsADirectoryError, match='fasterq-dump output path exists but is not a file'):
            run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)


class TestSequenceExtractionSecondRound:
    def test_raises_when_second_round_fastq_is_missing(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        (tmp_path / '{}.amalgkit.fastq.gz'.format(sra_id)).write_text('AAAA\n')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
            'spot_start_2nd': [2],
            'spot_end_2nd': [5],
            'bp_written': [100],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        args = SimpleNamespace()
        g = {}

        monkeypatch.setattr('amalgkit.getfastq.sequence_extraction', lambda *a, **k: metadata)

        with pytest.raises(FileNotFoundError, match='Dumped fastq not found'):
            sequence_extraction_2nd_round(args=args, sra_stat=sra_stat, metadata=metadata, g=g)

    def test_raises_when_second_round_fastq_path_is_directory(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        (tmp_path / '{}.amalgkit.fastq.gz'.format(sra_id)).write_text('AAAA\n')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
            'spot_start_2nd': [2],
            'spot_end_2nd': [5],
            'bp_written': [100],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        args = SimpleNamespace()
        g = {}

        def fake_sequence_extraction(*_args, **_kwargs):
            tmp_first = tmp_path / '{}.amalgkit_1st.fastq.gz'.format(sra_id)
            if tmp_first.exists():
                tmp_first.unlink()
            tmp_first.mkdir()
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.sequence_extraction', fake_sequence_extraction)

        with pytest.raises(IsADirectoryError, match='Dumped fastq path exists but is not a file'):
            sequence_extraction_2nd_round(args=args, sra_stat=sra_stat, metadata=metadata, g=g)

    def test_paired_2nd_round_merges_subexts_with_parallel_workers(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        (tmp_path / '{}_1.amalgkit.fastq.gz'.format(sra_id)).write_bytes(b'FIRST1')
        (tmp_path / '{}_2.amalgkit.fastq.gz'.format(sra_id)).write_bytes(b'FIRST2')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['paired'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
            'spot_start_2nd': [2],
            'spot_end_2nd': [5],
            'bp_written': [100],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        args = SimpleNamespace(threads=2)
        g = {}
        observed = {'max_workers': None, 'task_items': None}

        def fake_sequence_extraction(*_args, **_kwargs):
            (tmp_path / '{}_1.amalgkit.fastq.gz'.format(sra_id)).write_bytes(b'SECOND1')
            (tmp_path / '{}_2.amalgkit.fastq.gz'.format(sra_id)).write_bytes(b'SECOND2')
            return metadata

        def fake_run_tasks(task_items, task_fn, max_workers):
            observed['max_workers'] = max_workers
            observed['task_items'] = list(task_items)
            results = {}
            failures = []
            for task_item in task_items:
                try:
                    results[task_item] = task_fn(task_item)
                except Exception as exc:
                    failures.append((task_item, exc))
            return results, failures

        monkeypatch.setattr('amalgkit.getfastq.sequence_extraction', fake_sequence_extraction)
        monkeypatch.setattr('amalgkit.getfastq.run_tasks_with_optional_threads', fake_run_tasks)

        sequence_extraction_2nd_round(args=args, sra_stat=sra_stat, metadata=metadata, g=g)

        assert observed['max_workers'] == 2
        assert set(observed['task_items']) == {'_1', '_2'}
        assert (tmp_path / '{}_1.amalgkit.fastq.gz'.format(sra_id)).read_bytes() == b'FIRST1SECOND1'
        assert (tmp_path / '{}_2.amalgkit.fastq.gz'.format(sra_id)).read_bytes() == b'FIRST2SECOND2'


class TestDownloadSraUrlSchemes:
    @staticmethod
    def _make_metadata(sra_id, aws_link, gcp_link, ncbi_link):
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'AWS_Link': [aws_link],
            'GCP_Link': [gcp_link],
            'NCBI_Link': [ncbi_link],
        }))

    @staticmethod
    def _make_args(gcp_project='', sra_download_method='urllib'):
        args = type('Args', (), {})()
        args.aws = True
        args.gcp = True
        args.ncbi = True
        args.gcp_project = gcp_project
        args.sra_download_method = sra_download_method
        return args

    def test_raises_when_existing_sra_path_is_directory(self, tmp_path):
        sra_id = 'SRR_DIR'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='',
            ncbi_link='https://example.invalid/path/to.sra',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args()
        (tmp_path / '{}.sra'.format(sra_id)).mkdir()

        with pytest.raises(IsADirectoryError, match='SRA path exists but is not a file'):
            download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

    def test_skips_non_http_schemes_before_urllib(self, tmp_path, monkeypatch, capsys):
        sra_id = 'SRR_SCHEME'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='s3://bucket/path/to.sra',
            gcp_link='gs://bucket/path/to.sra',
            ncbi_link='https://example.invalid/path/to.sra',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args()
        called_urls = []

        def fake_urlretrieve(url, _path):
            called_urls.append(url)
            raise urllib.error.URLError('network down')

        monkeypatch.setattr('amalgkit.getfastq.urllib.request.urlretrieve', fake_urlretrieve)

        with pytest.raises(FileNotFoundError):
            download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert called_urls == [
            'https://storage.googleapis.com/bucket/path/to.sra',
            'https://example.invalid/path/to.sra',
        ]
        stderr = capsys.readouterr().err
        assert 'unsupported URL scheme for urllib: s3' in stderr

    def test_downloads_when_http_source_succeeds(self, tmp_path, monkeypatch):
        sra_id = 'SRR_OK'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='',
            ncbi_link='https://example.invalid/path/to.sra',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args()
        out_path = tmp_path / '{}.sra'.format(sra_id)

        def fake_urlretrieve(_url, path):
            with open(path, 'w') as f:
                f.write('ok')
            return (path, None)

        monkeypatch.setattr('amalgkit.getfastq.urllib.request.urlretrieve', fake_urlretrieve)
        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)
        assert out_path.exists()

    def test_appends_user_project_for_gcp_requester_pays(self, tmp_path, monkeypatch):
        sra_id = 'SRR_GCP_PROJECT'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='gs://bucket/path/to.sra',
            ncbi_link='',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args(gcp_project='test-project')
        called_urls = []

        def fake_urlretrieve(url, _path):
            called_urls.append(url)
            raise urllib.error.URLError('network down')

        monkeypatch.setattr('amalgkit.getfastq.urllib.request.urlretrieve', fake_urlretrieve)

        with pytest.raises(FileNotFoundError):
            download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert called_urls == ['https://storage.googleapis.com/bucket/path/to.sra?userProject=test-project']

    def test_uses_curl_when_requested(self, tmp_path, monkeypatch):
        sra_id = 'SRR_CURL'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='',
            ncbi_link='https://example.invalid/path/to.sra',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args(sra_download_method='curl')
        observed = {'cmd': None}

        monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda exe: '/usr/bin/curl' if exe == 'curl' else None)

        def fail_urlretrieve(*_args, **_kwargs):
            raise AssertionError('urllib.request.urlretrieve should not be called when curl succeeds.')

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            out_path = cmd[cmd.index('-o') + 1]
            with open(out_path, 'w') as fh:
                fh.write('ok')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.urllib.request.urlretrieve', fail_urlretrieve)
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert observed['cmd'][0] == '/usr/bin/curl'
        assert (tmp_path / '{}.sra'.format(sra_id)).exists()


class TestSequenceExtractionPrivate:
    def test_handles_missing_read2_path_without_type_error(self, tmp_path, monkeypatch):
        read1_path = tmp_path / 'input_R1.fastq.gz'
        read1_path.write_text('dummy')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'read1_path': [str(read1_path)],
            'read2_path': [numpy.nan],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'work'
        sra_dir.mkdir()
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(sra_dir),
        }
        args = SimpleNamespace(fastp=False)

        monkeypatch.setattr('amalgkit.getfastq.set_current_intermediate_extension', lambda *_args, **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.get_or_detect_intermediate_extension', lambda *_args, **_kwargs: '.fastq.gz')
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *_args, **_kwargs: None)

        sequence_extraction_private(metadata=metadata, sra_stat=sra_stat, args=args)

        assert (sra_dir / 'input_R1.fastq.gz').exists()

    def test_warns_when_private_path_is_directory(self, tmp_path, monkeypatch, capsys):
        read1_dir = tmp_path / 'read1_dir'
        read1_dir.mkdir()
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'read1_path': [str(read1_dir)],
            'read2_path': ['unavailable'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'work'
        sra_dir.mkdir()
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(sra_dir),
        }
        args = SimpleNamespace(fastp=False)

        monkeypatch.setattr('amalgkit.getfastq.set_current_intermediate_extension', lambda *_args, **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.get_or_detect_intermediate_extension', lambda *_args, **_kwargs: '.fastq.gz')
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *_args, **_kwargs: None)

        sequence_extraction_private(metadata=metadata, sra_stat=sra_stat, args=args)
        err = capsys.readouterr().err

        assert 'exists but is not a file' in err

    def test_raises_when_private_output_path_is_directory(self, tmp_path, monkeypatch):
        read1_path = tmp_path / 'input_R1.fastq.gz'
        read1_path.write_text('dummy')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'read1_path': [str(read1_path)],
            'read2_path': ['unavailable'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'work'
        sra_dir.mkdir()
        (sra_dir / 'input_R1.fastq.gz').mkdir()
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(sra_dir),
        }
        args = SimpleNamespace(fastp=False)

        monkeypatch.setattr('amalgkit.getfastq.set_current_intermediate_extension', lambda *_args, **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.get_or_detect_intermediate_extension', lambda *_args, **_kwargs: '.fastq.gz')
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *_args, **_kwargs: None)

        with pytest.raises(IsADirectoryError, match='Private output path exists but is not a file/symlink'):
            sequence_extraction_private(metadata=metadata, sra_stat=sra_stat, args=args)


class TestGetfastqMainJobs:
    def test_rejects_nonpositive_jobs(self):
        args = SimpleNamespace(internal_jobs=0)
        with pytest.raises(ValueError, match='--internal_jobs must be > 0'):
            getfastq_main(args)

    def test_parallel_jobs_process_all_runs(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'lib_layout': ['single', 'single'],
            'total_spots': [10, 10],
            'total_bases': [1000, 1000],
            'spot_length': [100, 100],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        args = SimpleNamespace(
            internal_jobs=2,
            redo=False,
            out_dir=str(tmp_path),
            remove_sra=False,
            tol=1,
            fastp=False,
            read_name='default',
        )
        processed = []

        def fake_process_getfastq_run(args, row_index, sra_id, run_row_df, g):
            processed.append(sra_id)
            return {
                'row_index': row_index,
                'sra_id': sra_id,
                'row': run_row_df.iloc[0].copy(),
                'flag_any_output_file_present': True,
                'flag_private_file': False,
                'getfastq_sra_dir': os.path.join(args.out_dir, 'getfastq', sra_id),
            }

        monkeypatch.setattr('amalgkit.getfastq.check_getfastq_dependency', lambda _args: None)
        monkeypatch.setattr('amalgkit.getfastq.getfastq_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.getfastq.remove_experiment_without_run', lambda m: m)
        monkeypatch.setattr('amalgkit.getfastq.check_metadata_validity', lambda m: m)
        monkeypatch.setattr(
            'amalgkit.getfastq.initialize_global_params',
            lambda _args, _metadata: {
                'start_time': 0.0,
                'max_bp': 2000,
                'num_sra': 2,
                'num_bp_per_sra': 1000,
                'total_sra_bp': 2000,
            },
        )
        monkeypatch.setattr('amalgkit.getfastq.initialize_columns', lambda m, _g: m)
        monkeypatch.setattr('amalgkit.getfastq.process_getfastq_run', fake_process_getfastq_run)

        getfastq_main(args)

        assert set(processed) == {'SRR001', 'SRR002'}

    def test_cpu_budget_caps_jobs_to_serial(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'lib_layout': ['single', 'single'],
            'total_spots': [10, 10],
            'total_bases': [1000, 1000],
            'spot_length': [100, 100],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        args = SimpleNamespace(
            internal_jobs=4,
            threads=4,
            internal_cpu_budget=1,
            redo=False,
            out_dir=str(tmp_path),
            remove_sra=False,
            tol=1,
            fastp=False,
            read_name='default',
        )
        processed = []

        def fake_process_getfastq_run(args, row_index, sra_id, run_row_df, g):
            processed.append(sra_id)
            return {
                'row_index': row_index,
                'sra_id': sra_id,
                'row': run_row_df.iloc[0].copy(),
                'flag_any_output_file_present': True,
                'flag_private_file': False,
                'getfastq_sra_dir': os.path.join(args.out_dir, 'getfastq', sra_id),
            }

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('run_tasks_with_optional_threads should not be used when --internal_cpu_budget caps internal_jobs to 1.')

        monkeypatch.setattr('amalgkit.getfastq.check_getfastq_dependency', lambda _args: None)
        monkeypatch.setattr('amalgkit.getfastq.getfastq_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.getfastq.remove_experiment_without_run', lambda m: m)
        monkeypatch.setattr('amalgkit.getfastq.check_metadata_validity', lambda m: m)
        monkeypatch.setattr(
            'amalgkit.getfastq.initialize_global_params',
            lambda _args, _metadata: {
                'start_time': 0.0,
                'max_bp': 2000,
                'num_sra': 2,
                'num_bp_per_sra': 1000,
                'total_sra_bp': 2000,
            },
        )
        monkeypatch.setattr('amalgkit.getfastq.initialize_columns', lambda m, _g: m)
        monkeypatch.setattr('amalgkit.getfastq.process_getfastq_run', fake_process_getfastq_run)
        monkeypatch.setattr('amalgkit.getfastq.run_tasks_with_optional_threads', fail_if_called)

        getfastq_main(args)

        assert set(processed) == {'SRR001', 'SRR002'}
        assert args.threads == 1

    def test_private_file_in_any_run_skips_second_round(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'lib_layout': ['single', 'single'],
            'total_spots': [10, 10],
            'total_bases': [1000, 1000],
            'spot_length': [100, 100],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
            'bp_amalgkit': [500, 500],
        }))
        args = SimpleNamespace(
            internal_jobs=2,
            redo=False,
            out_dir=str(tmp_path),
            remove_sra=False,
            tol=1,
            fastp=False,
            read_name='default',
        )

        def fake_process_getfastq_run(args, row_index, sra_id, run_row_df, g):
            return {
                'row_index': row_index,
                'sra_id': sra_id,
                'row': run_row_df.iloc[0].copy(),
                'flag_any_output_file_present': False,
                'flag_private_file': (sra_id == 'SRR001'),
                'getfastq_sra_dir': os.path.join(args.out_dir, 'getfastq', sra_id),
            }

        monkeypatch.setattr('amalgkit.getfastq.check_getfastq_dependency', lambda _args: None)
        monkeypatch.setattr('amalgkit.getfastq.getfastq_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.getfastq.remove_experiment_without_run', lambda m: m)
        monkeypatch.setattr('amalgkit.getfastq.check_metadata_validity', lambda m: m)
        monkeypatch.setattr(
            'amalgkit.getfastq.initialize_global_params',
            lambda _args, _metadata: {
                'start_time': 0.0,
                'max_bp': 2000,
                'num_sra': 2,
                'num_bp_per_sra': 1000,
                'total_sra_bp': 2000,
            },
        )
        monkeypatch.setattr('amalgkit.getfastq.initialize_columns', lambda m, _g: m)
        monkeypatch.setattr('amalgkit.getfastq.process_getfastq_run', fake_process_getfastq_run)
        monkeypatch.setattr('amalgkit.getfastq.print_read_stats', lambda *_args, **_kwargs: None)

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('Second-round extraction should be skipped when any run is private.')

        monkeypatch.setattr('amalgkit.getfastq.is_2nd_round_needed', fail_if_called)

        getfastq_main(args)


class TestGetfastqDependencyChecks:
    def test_handles_missing_fastp_attributes(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'] == [['fasterq-dump', '-h'], ['seqkit', '--help']]

    def test_detects_missing_fasterq_spot_range_flags(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'

        def fake_run(cmd, stdout=None, stderr=None):
            if cmd[0] == 'fasterq-dump':
                return subprocess.CompletedProcess(
                    cmd,
                    0,
                    stdout=b'Usage:\\n  -x|--details\\n  -s|--split-spot\\n',
                    stderr=b'',
                )
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        args = Args()
        check_getfastq_dependency(args)
        assert args._fasterq_supports_spot_range is False

    def test_detects_fasterq_spot_range_flags_when_present(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'

        def fake_run(cmd, stdout=None, stderr=None):
            if cmd[0] == 'fasterq-dump':
                return subprocess.CompletedProcess(
                    cmd,
                    0,
                    stdout=b'Usage:\\n  -N|--minSpotId\\n  -X|--maxSpotId\\n',
                    stderr=b'',
                )
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        args = Args()
        check_getfastq_dependency(args)
        assert args._fasterq_supports_spot_range is True

    def test_uses_default_fastp_exe_when_fastp_enabled_without_override(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'
            fastp = True

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'][0] == ['fasterq-dump', '-h']
        assert called['cmds'][1] == ['seqkit', '--help']
        assert called['cmds'][2] == ['fastp', '--help']

    def test_uses_fasterq_dump_dependency(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'default'

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'][0][0] == 'fasterq-dump'

    def test_obsolete_flags_are_ignored_and_still_uses_fasterq_dump(self, monkeypatch):
        class Args:
            obsolete_pfd = True
            obsolete_pfd_exe = '/tmp/legacy_pfd_exe'
            obsolete_fastq_dump_exe = '/tmp/legacy_fastq_dump_exe'
            obsolete_prefetch_exe = '/tmp/legacy_prefetch_exe'
            fasterq_dump_exe = 'fasterq-dump'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'default'

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'][0][0] == 'fasterq-dump'

    def test_trinity_mode_uses_same_dependencies(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'trinity'

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'][0][0] == 'fasterq-dump'

    def test_raises_clear_error_when_fasterq_dump_missing(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'missing-fasterq-dump'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'default'

        def fake_run(cmd, stdout=None, stderr=None):
            raise FileNotFoundError(cmd[0])

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        with pytest.raises(FileNotFoundError, match='fasterq-dump executable not found'):
            check_getfastq_dependency(Args())

    def test_raises_clear_error_when_seqkit_missing(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'
            seqkit_exe = 'missing-seqkit'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'default'

        def fake_run(cmd, stdout=None, stderr=None):
            if cmd[0] == 'missing-seqkit':
                raise FileNotFoundError(cmd[0])
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        with pytest.raises(FileNotFoundError, match='seqkit executable not found'):
            check_getfastq_dependency(Args())

    def test_raises_clear_error_when_fastp_probe_fails(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'
            fastp = True
            fastp_exe = 'fastp'
            read_name = 'default'

        def fake_run(cmd, stdout=None, stderr=None):
            if cmd[0] == 'fastp':
                return subprocess.CompletedProcess(cmd, 127, stdout=b'', stderr=b'')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        with pytest.raises(RuntimeError, match='fastp dependency probe failed'):
            check_getfastq_dependency(Args())

    def test_rrna_filter_probes_sortmerna_and_prepares_download_dir(self, tmp_path, monkeypatch):
        download_dir = tmp_path / 'downloads'
        args = SimpleNamespace(
            fasterq_dump_exe='fasterq-dump',
            fastp=False,
            rrna_filter=True,
            download_dir=str(download_dir),
            sortmerna_exe='sortmerna',
        )

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(args)
        assert called['cmds'][0] == ['fasterq-dump', '-h']
        assert called['cmds'][1] == ['seqkit', '--help']
        assert called['cmds'][2] == ['sortmerna', '--version']
        assert download_dir.exists()
        assert download_dir.is_dir()

    def test_rrna_filter_raises_when_download_path_is_file(self, tmp_path, monkeypatch):
        bad_path = tmp_path / 'downloads_as_file'
        bad_path.write_text('not a directory')

        args = SimpleNamespace(
            fasterq_dump_exe='fasterq-dump',
            fastp=False,
            rrna_filter=True,
            download_dir=str(bad_path),
            sortmerna_exe='sortmerna',
        )

        monkeypatch.setattr(
            'amalgkit.getfastq.subprocess.run',
            lambda cmd, stdout=None, stderr=None: subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b''),
        )
        with pytest.raises(NotADirectoryError, match='SortMeRNA download path exists but is not a directory'):
            check_getfastq_dependency(args)


class TestSortmernaReferenceDownload:
    def test_downloads_silva_refs_to_default_out_dir_downloads(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        calls = {'n': 0}

        def fake_urlretrieve(_url, out_path):
            calls['n'] += 1
            with gzip.open(out_path, 'wt') as fout:
                fout.write('>rRNA\nACGT\n')
            return (out_path, None)

        monkeypatch.setattr('amalgkit.getfastq.urllib.request.urlretrieve', fake_urlretrieve)

        args = SimpleNamespace(
            out_dir=str(out_dir),
            download_dir='inferred',
        )

        refs_first = resolve_sortmerna_refs(args)
        assert len(refs_first) == 2
        for ref_path in refs_first:
            assert ref_path.startswith(os.path.join(str(out_dir), 'downloads'))
            assert os.path.exists(ref_path)
        assert calls['n'] == 2

        refs_second = resolve_sortmerna_refs(args)
        assert refs_second == refs_first
        assert calls['n'] == 2

    def test_downloads_silva_refs_to_custom_download_dir(self, tmp_path, monkeypatch):
        custom_dir = tmp_path / 'custom_downloads'

        def fake_urlretrieve(_url, out_path):
            with gzip.open(out_path, 'wt') as fout:
                fout.write('>rRNA\nACGT\n')
            return (out_path, None)

        monkeypatch.setattr('amalgkit.getfastq.urllib.request.urlretrieve', fake_urlretrieve)

        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(custom_dir),
        )

        refs = resolve_sortmerna_refs(args)
        assert len(refs) == 2
        for ref_path in refs:
            assert ref_path.startswith(str(custom_dir))
            assert os.path.exists(ref_path)


class TestFasterqSeqkitCompression:
    def _write_fastq(self, path, reads):
        with open(path, 'wt') as out:
            for i, seq in enumerate(reads):
                out.write('@r{}\n'.format(i))
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_uses_seqkit_for_compression(self, tmp_path, monkeypatch):
        sra_id = 'SRR999'
        fastq_path = tmp_path / '{}.fastq'.format(sra_id)
        self._write_fastq(str(fastq_path), ['ACGT', 'TGCA'])

        class Args:
            threads = 2
            dump_print = False
            seqkit_exe = 'seqkit'

        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
        }

        def fake_seqkit_run(cmd, stdout=None, stderr=None):
            assert cmd[0] == 'seqkit'
            assert cmd[1] == 'seq'
            out_index = cmd.index('-o')
            out_path = cmd[out_index + 1]
            in_path = cmd[out_index + 2]
            with open(in_path, 'rb') as fin, gzip.open(out_path, 'wb') as fout:
                fout.write(fin.read())
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_seqkit_run)
        compress_fasterq_output_files(sra_stat=sra_stat, args=Args())

        gz_path = tmp_path / '{}.fastq.gz'.format(sra_id)
        assert not fastq_path.exists()
        assert gz_path.exists()
        with gzip.open(str(gz_path), 'rt') as f:
            content = f.read()
        assert '@r0' in content
        assert '@r1' in content


class TestRenameReadsSeqkitCompression:
    def _write_fastq_gz(self, path, headers_and_seqs):
        with gzip.open(path, 'wt') as out:
            for header, seq in headers_and_seqs:
                out.write(header + '\n')
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_single_end_trinity_rename(self, tmp_path, monkeypatch):
        sra_id = 'SRR777'
        in_path = tmp_path / '{}.fastq.gz'.format(sra_id)
        self._write_fastq_gz(str(in_path), [('@r0 comment', 'ACGT'), ('@r1 x', 'TGCA')])
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            threads = 1
            remove_tmp = False
            dump_print = False
            seqkit_exe = 'seqkit'

        def fake_seqkit_run(cmd, stdout=None, stderr=None):
            assert cmd[0] == 'seqkit'
            assert cmd[1] == 'replace'
            assert '-p' in cmd
            assert '-r' in cmd
            out_index = cmd.index('-o')
            out_path = cmd[out_index + 1]
            in_path = cmd[-1]
            replacement = cmd[cmd.index('-r') + 1]
            suffix = replacement[len('$1'):] if replacement.startswith('$1') else ''
            in_open = gzip.open if str(in_path).endswith('.gz') else open
            out_open = gzip.open if str(out_path).endswith('.gz') else open
            with in_open(in_path, 'rt') as fin, out_open(out_path, 'wt') as fout:
                while True:
                    line1 = fin.readline()
                    if line1 == '':
                        break
                    line2 = fin.readline()
                    line3 = fin.readline()
                    line4 = fin.readline()
                    header = line1.rstrip('\n').split()[0]
                    fout.write(header + suffix + '\n')
                    fout.write(line2)
                    fout.write(line3)
                    fout.write(line4)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_seqkit_run)

        rename_reads(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path))
        out_path = tmp_path / '{}.rename.fastq.gz'.format(sra_id)
        assert out_path.exists()
        with gzip.open(str(out_path), 'rt') as f:
            lines = [f.readline().rstrip('\n') for _ in range(4)]
        assert lines[0] == '@r0/1'

    def test_paired_end_trinity_rename_uses_parallel_seqkit_workers(self, tmp_path, monkeypatch):
        sra_id = 'SRR778'
        in_path1 = tmp_path / '{}_1.fastq.gz'.format(sra_id)
        in_path2 = tmp_path / '{}_2.fastq.gz'.format(sra_id)
        self._write_fastq_gz(str(in_path1), [('@r0 comment1', 'ACGT')])
        self._write_fastq_gz(str(in_path2), [('@r0 comment2', 'TGCA')])
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            threads = 2
            remove_tmp = False
            dump_print = False
            seqkit_exe = 'seqkit'

        observed_threads = []

        def fake_seqkit_run(cmd, stdout=None, stderr=None):
            assert cmd[0] == 'seqkit'
            assert cmd[1] == 'replace'
            observed_threads.append(cmd[cmd.index('-j') + 1])
            out_index = cmd.index('-o')
            out_path = cmd[out_index + 1]
            in_path = cmd[-1]
            replacement = cmd[cmd.index('-r') + 1]
            suffix = replacement[len('$1'):] if replacement.startswith('$1') else ''
            with gzip.open(in_path, 'rt') as fin, gzip.open(out_path, 'wt') as fout:
                while True:
                    line1 = fin.readline()
                    if line1 == '':
                        break
                    line2 = fin.readline()
                    line3 = fin.readline()
                    line4 = fin.readline()
                    header = line1.rstrip('\n').split()[0]
                    fout.write(header + suffix + '\n')
                    fout.write(line2)
                    fout.write(line3)
                    fout.write(line4)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_seqkit_run)

        rename_reads(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path))
        out_path1 = tmp_path / '{}_1.rename.fastq.gz'.format(sra_id)
        out_path2 = tmp_path / '{}_2.rename.fastq.gz'.format(sra_id)
        assert out_path1.exists()
        assert out_path2.exists()
        with gzip.open(str(out_path1), 'rt') as f:
            first_header1 = f.readline().rstrip('\n')
        with gzip.open(str(out_path2), 'rt') as f:
            first_header2 = f.readline().rstrip('\n')
        assert first_header1 == '@r0/1'
        assert first_header2 == '@r0/2'
        assert sorted(observed_threads) == ['1', '1']
