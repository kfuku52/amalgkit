import re
import os
import pytest
import pandas

from amalgkit.curate import get_sample_group, run_curate_r_script, curate_main
from amalgkit.util import Metadata


# ---------------------------------------------------------------------------
# get_sample_group (wiki: curate extracts sample groups from args or metadata)
# ---------------------------------------------------------------------------

class TestGetSampleGroup:
    def test_from_args(self):
        """When --sample_group is specified, parse it."""
        class Args:
            sample_group = 'brain,liver,heart'
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'], 'exclusion': ['no'],
            'sample_group': ['brain'],
        }))
        result = get_sample_group(Args(), m)
        assert result == 'brain|liver|heart'

    def test_from_metadata(self):
        """When --sample_group is None, extract from metadata."""
        class Args:
            sample_group = None
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2', 'R3'],
            'exclusion': ['no', 'no', 'no'],
            'sample_group': ['brain', 'liver', 'brain'],
        }))
        result = get_sample_group(Args(), m)
        # Should contain brain and liver separated by pipe
        groups = result.split('|')
        assert 'brain' in groups
        assert 'liver' in groups

    def test_empty_sample_group_exits(self):
        """Wiki: curate exits with error if sample_group is empty."""
        class Args:
            sample_group = None
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'exclusion': ['no'],
            'sample_group': [float('nan')],
        }))
        with pytest.raises(SystemExit):
            get_sample_group(Args(), m)


class TestRunCurateRScript:
    @pytest.mark.parametrize('one_outlier_per_iter, expected_flag', [(False, '0'), (True, '1')])
    def test_passes_one_outlier_flag_to_rscript(self, tmp_path, monkeypatch, one_outlier_per_iter, expected_flag):
        class Args:
            pass
        args = Args()
        args.dist_method = 'pearson'
        args.mapping_rate = 0.2
        args.correlation_threshold = 0.3
        args.plot_intermediate = False
        args.sample_group = None
        args.sample_group_color = 'DEFAULT'
        args.norm = 'log2p1-fpkm'
        args.one_outlier_per_iter = one_outlier_per_iter
        args.batch_effect_alg = 'no'
        args.clip_negative = True
        args.maintain_zero = True
        args.skip_curation = False
        args.out_dir = str(tmp_path / 'out')

        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'sample_group': ['groupA'],
            'exclusion': ['no'],
            'scientific_name': ['sp'],
        }))
        input_dir = tmp_path / 'input_cstmm'
        (input_dir / 'SpA').mkdir(parents=True)
        (input_dir / 'metadata.tsv').write_text('run\tsample_group\texclusion\nR1\tgroupA\tno\n')
        (input_dir / 'SpA' / 'SpA_cstmm_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        (input_dir / 'SpA' / 'SpA_eff_length.tsv').write_text('target_id\tR1\nG1\t100\n')
        captured = {}

        def fake_call(cmd):
            captured['cmd'] = cmd
            return 0

        monkeypatch.setattr('amalgkit.curate.subprocess.call', fake_call)
        code = run_curate_r_script(
            args=args,
            metadata=metadata,
            sp='SpA',
            input_dir=str(input_dir),
        )
        assert code == 0
        assert captured['cmd'][13] == expected_flag


class TestCurateMain:
    @staticmethod
    def _args(out_dir):
        class Args:
            pass
        args = Args()
        args.input_dir = 'custom_input'
        args.out_dir = str(out_dir)
        args.norm = 'log2p1-fpkm'
        args.redo = False
        args.sample_group = None
        args.metadata = 'metadata.tsv'
        args.dist_method = 'pearson'
        args.mapping_rate = 0.2
        args.correlation_threshold = 0.3
        args.plot_intermediate = False
        args.sample_group_color = 'DEFAULT'
        args.one_outlier_per_iter = False
        args.batch_effect_alg = 'no'
        args.clip_negative = True
        args.maintain_zero = True
        args.skip_curation = False
        args.species_jobs = 1
        return args

    def test_creates_nested_out_dir_and_completion_flags(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A', 'Species B'],
            'run': ['R1', 'R2'],
            'sample_group': ['g1', 'g2'],
            'exclusion': ['no', 'no'],
        }))

        monkeypatch.setattr('amalgkit.curate.check_rscript', lambda: None)
        monkeypatch.setattr('amalgkit.curate.load_metadata', lambda *_args, **_kwargs: metadata)

        def fake_run_curate(args, metadata, sp, input_dir):
            os.makedirs(os.path.join(args.out_dir, 'curate', sp), exist_ok=True)
            return 0

        monkeypatch.setattr('amalgkit.curate.run_curate_r_script', fake_run_curate)
        curate_main(args)
        assert (out_dir / 'curate' / 'Species_A' / 'curate_completion_flag.txt').exists()
        assert (out_dir / 'curate' / 'Species_B' / 'curate_completion_flag.txt').exists()

    def test_exits_nonzero_when_any_species_fails(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A', 'Species B'],
            'run': ['R1', 'R2'],
            'sample_group': ['g1', 'g2'],
            'exclusion': ['no', 'no'],
        }))

        monkeypatch.setattr('amalgkit.curate.check_rscript', lambda: None)
        monkeypatch.setattr('amalgkit.curate.load_metadata', lambda *_args, **_kwargs: metadata)

        def fake_run_curate(args, metadata, sp, input_dir):
            os.makedirs(os.path.join(args.out_dir, 'curate', sp), exist_ok=True)
            return 1 if sp == 'Species_A' else 0

        monkeypatch.setattr('amalgkit.curate.run_curate_r_script', fake_run_curate)
        with pytest.raises(SystemExit) as e:
            curate_main(args)
        assert e.value.code == 1
        assert not (out_dir / 'curate' / 'Species_A' / 'curate_completion_flag.txt').exists()
        assert (out_dir / 'curate' / 'Species_B' / 'curate_completion_flag.txt').exists()

    def test_parallel_species_jobs_creates_completion_flags(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        args.species_jobs = 2
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A', 'Species B'],
            'run': ['R1', 'R2'],
            'sample_group': ['g1', 'g2'],
            'exclusion': ['no', 'no'],
        }))

        monkeypatch.setattr('amalgkit.curate.check_rscript', lambda: None)
        monkeypatch.setattr('amalgkit.curate.load_metadata', lambda *_args, **_kwargs: metadata)

        def fake_run_curate(args, metadata, sp, input_dir):
            os.makedirs(os.path.join(args.out_dir, 'curate', sp), exist_ok=True)
            return 0

        monkeypatch.setattr('amalgkit.curate.run_curate_r_script', fake_run_curate)
        curate_main(args)

        assert (out_dir / 'curate' / 'Species_A' / 'curate_completion_flag.txt').exists()
        assert (out_dir / 'curate' / 'Species_B' / 'curate_completion_flag.txt').exists()

    def test_rejects_nonpositive_species_jobs(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        args.species_jobs = 0
        monkeypatch.setattr('amalgkit.curate.check_rscript', lambda: None)
        with pytest.raises(ValueError, match='--species_jobs must be > 0'):
            curate_main(args)
