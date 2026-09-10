import json
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import pandas
import pytest

from amalgkit.fragment_length import (
    FRAGMENT_METADATA_COLUMNS, PROVENANCE_KEY, fragment_file_records,
    load_fragment_length_file, resolve_fragment_distribution, validate_fragment_args,
)
from amalgkit.main import build_main_parser
from amalgkit.command_context import QuantRuntimeContext
from amalgkit.getfastq import initialize_columns
from amalgkit.merge import collect_quant_models
from amalgkit.metadata_utils import Metadata, get_sra_stat
from amalgkit.output_contracts import validate_quant_run_info_json
from amalgkit.quant import call_kallisto, check_fragment_model_reuse, resolve_run_fragment_model, run_quant, quant_main
from amalgkit.rerun import rerun_quant_check
from amalgkit.select import apply_select_aggregate_rules, read_select_rules
from amalgkit.table_io import read_annotation_tsv


def metadata(**values):
    row = dict(run='R1', scientific_name='Species A', lib_layout='single',
               total_spots=10, total_bases=1000, spot_length=100)
    row.update(values)
    return Metadata.from_DataFrame(pandas.DataFrame([row]))


def args(**values):
    return SimpleNamespace(threads=1, **values)


@pytest.mark.parametrize('mean', [150, 199.5, 200, 350])
@pytest.mark.parametrize('use_file', [False, True])
def test_known_mean_and_independent_sd_reach_kallisto(tmp_path, monkeypatch, capsys, mean, use_file):
    observed = []

    def execute(command, **kwargs):
        observed.extend(command)
        (tmp_path / 'run_info.json').write_text('{"p_pseudoaligned":100,"kallisto_version":"test"}')
        (tmp_path / 'abundance.tsv').write_text('target_id\tlength\teff_length\test_counts\ttpm\nt1\t500\t350\t10\t1000000\n')
        return SimpleNamespace(returncode=0, stdout=b'', stderr=b'')

    monkeypatch.setattr('amalgkit.quant.subprocess.run', execute)
    data = metadata(nominal_length=mean, nominal_sdev=7)
    runtime = args()
    if use_file:
        path = write_run_file(tmp_path, [['R1', mean, 7, 'measured', 'insert assay L1']])
        runtime.fragment_length_file = str(path)
        data.df['nominal_length'] = mean + 50
        data.df['nominal_sdev'] = 9
    call_kallisto(runtime, ['reads.fq'], data, get_sra_stat('R1', data), str(tmp_path), 'index.idx')
    assert float(observed[observed.index('-l') + 1]) == mean
    assert float(observed[observed.index('-s') + 1]) == 7
    assert 'WARNING' not in capsys.readouterr().err
    info = json.loads((tmp_path / 'R1_run_info.json').read_text())
    assert info['kallisto_version'] == 'test'
    assert info[PROVENANCE_KEY]['command'] == observed
    assert info[PROVENANCE_KEY]['sd']['source'] == ('run_file:measured' if use_file else 'metadata_nominal')
    if use_file:
        assert info[PROVENANCE_KEY]['source_file'] == str(path)
    assert validate_quant_run_info_json(str(tmp_path / 'R1_run_info.json')) == ''


@pytest.mark.parametrize('row,expected,assumed', [
    ({}, (200, 20), ['mean', 'sd']),
    ({'nominal_length': 150}, (150, 20), ['sd']),
    ({'nominal_sdev': 7}, (200, 7), ['mean']),
    ({'nominal_length': 350}, (350, 20), ['sd']),
    ({'spot_length': 80, 'total_bases': 800, 'total_spots': 10, 'fastp_insert_size_peak': 123}, (200, 20), ['mean', 'sd']),
])
def test_assume_fills_only_missing_values_and_warns(row, expected, assumed, capsys):
    model = resolve_fragment_distribution(args(), row, 'R1', {})
    assert (model['mean']['value'], model['sd']['value']) == expected
    assert [field for field in ('mean', 'sd') if model[field]['source'] == 'assumed'] == assumed
    warning = capsys.readouterr().err
    assert 'WARNING: Run R1:' in warning and 'not measured' in warning
    for field in assumed:
        assert field + '=' in warning
    with pytest.raises(ValueError, match='missing fragment length'):
        resolve_fragment_distribution(args(fragment_length_policy='error'), row, 'R1', {})


@pytest.mark.parametrize('bad', ['broken', '150; 200', '100-200', 'NaN', 'inf', '-inf', 0, -1, True])
@pytest.mark.parametrize('field', ['nominal_length', 'nominal_sdev', 'fragment_length_mean', 'fragment_length_sd'])
def test_invalid_metadata_is_not_treated_as_missing(bad, field):
    with pytest.raises(ValueError, match='finite number > 0'):
        resolve_fragment_distribution(args(), {field: bad}, 'R1', {})


@pytest.mark.parametrize('values', [
    {'fragment_length_mean': 150}, {'fragment_length_sd': 7},
    {'fragment_length_mean': float('nan'), 'fragment_length_sd': 7},
    {'fragment_length_mean': 150, 'fragment_length_sd': float('inf')},
    {'fragment_length_policy': 'typo'},
])
def test_invalid_cli_fails_even_when_metadata_would_override_it(values):
    with pytest.raises(ValueError):
        resolve_fragment_distribution(args(**values), {'fragment_length_mean': 150, 'fragment_length_sd': 7}, 'R1', {})


def write_run_file(tmp_path, rows):
    path = tmp_path / 'fragments.tsv'
    pandas.DataFrame(rows, columns=['run', 'fragment_length_mean', 'fragment_length_sd', 'source', 'source_detail']).to_csv(path, sep='\t', index=False)
    return path


def test_file_metadata_cli_nominal_priority_and_lexical_run_ids(tmp_path):
    path = write_run_file(tmp_path, [['0001', 140, 8, 'measured', 'library L1, adapter-corrected'],
                                     ['NA', 175, 13, 'user', 'protocol P2']])
    records = load_fragment_length_file(str(path), {'0001', '1', 'NA'})
    assert set(records) == {'0001', 'NA'}
    row = {'fragment_length_mean': 150, 'fragment_length_sd': 7,
           'nominal_length': 250, 'nominal_sdev': 40}
    common = args(fragment_length_mean=180, fragment_length_sd=30)
    assert resolve_fragment_distribution(common, row, '0001', records)['mean']['value'] == 140
    assert resolve_fragment_distribution(common, row, '1', records)['mean']['value'] == 150
    nominal = {key: value for key, value in row.items() if key.startswith('nominal')}
    assert resolve_fragment_distribution(common, nominal, '1', records)['mean']['value'] == 180
    assert resolve_fragment_distribution(args(), nominal, '1', records)['mean']['value'] == 250
    assert resolve_fragment_distribution(common, row, 'NA', records)['sd']['value'] == 13
    # A batch worker receives the already validated full table, even if only one run remains.
    common._fragment_length_by_run = records
    assert fragment_file_records(common, {'0001'}) == records


@pytest.mark.parametrize('rows,match', [
    ([['R2', 150, 7, 'user', 'manual']], 'Unknown'),
    ([['R1', 150, 7, 'user', 'manual']] * 2, 'Duplicate run'),
    ([['R1', 150, '', 'user', 'manual']], 'both mean and SD'),
    ([['R1', 150, 7, '', '']], 'source and source_detail'),
    ([['R1', 'inf', 7, 'user', 'manual']], 'finite number'),
])
def test_run_file_errors(tmp_path, rows, match):
    path = write_run_file(tmp_path, rows)
    with pytest.raises(ValueError, match=match):
        load_fragment_length_file(str(path), {'R1'})


def test_measured_metadata_requires_a_method_and_provenance_survives_tsv(tmp_path):
    data = metadata(fragment_length_mean=150, fragment_length_sd=7, fragment_length_source='measured')
    with pytest.raises(ValueError, match='source_detail is required'):
        resolve_run_fragment_model(args(), data, 'R1', 'single')
    data.df['fragment_length_source_detail'] = 'L1 insert sizes after adapter subtraction; assay A'
    data.reorder(omit_misc=True)
    path = tmp_path / 'metadata.tsv'
    data.df.to_csv(path, sep='\t', index=False)
    restored = Metadata.from_DataFrame(read_annotation_tsv(path))
    model = resolve_run_fragment_model(args(), restored, 'R1', 'single')
    assert model['sd']['source'] == 'metadata:measured'
    assert 'assay A' in model['mean']['source_detail']
    assert all(field in restored.df for field in FRAGMENT_METADATA_COLUMNS)


def test_sra_nominal_fields_are_reported_insert_statistics():
    root = ET.fromstring('''<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE>
      <EXPERIMENT><IDENTIFIERS><PRIMARY_ID>E1</PRIMARY_ID></IDENTIFIERS><DESIGN><LIBRARY_DESCRIPTOR>
      <LIBRARY_LAYOUT><PAIRED NOMINAL_LENGTH="150" NOMINAL_SDEV="7"/></LIBRARY_LAYOUT>
      </LIBRARY_DESCRIPTOR></DESIGN></EXPERIMENT>
      <SAMPLE><SAMPLE_ATTRIBUTES>
      <SAMPLE_ATTRIBUTE><TAG>nominal_sdev</TAG><VALUE>9</VALUE></SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE><TAG>fragment_length_mean</TAG><VALUE>999</VALUE></SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE><TAG>fragment_length_sd</TAG><VALUE>50</VALUE></SAMPLE_ATTRIBUTE>
      </SAMPLE_ATTRIBUTES></SAMPLE>
      <RUN_SET><RUN><IDENTIFIERS><PRIMARY_ID>R1</PRIMARY_ID></IDENTIFIERS></RUN></RUN_SET>
      </EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>''')
    data = Metadata.from_xml(root)
    row = data.df.iloc[0]
    assert row['nominal_sdev'] == '7'
    assert row['sample_attribute_nominal_sdev'] == '9'
    assert row['sample_attribute_fragment_length_mean'] == '999'
    assert row['fragment_length_mean'] == ''
    model = resolve_run_fragment_model(args(), data, 'R1', 'single')
    assert model['mean']['value'] == 150
    assert model['sd']['source'] == 'sra_experiment'
    assert 'E1' in model['sd']['source_detail']
    assert model['metadata_layout'] == 'paired'


def test_selection_preserves_distinct_nominal_candidates_even_with_old_rules():
    frame = pandas.DataFrame([{'nominal_length': '150', 'mean_insert_size': '200', 'nominal_sdev': '7'}])
    old_rule = dict(stage='aggregate', columns=['mean_insert_size'], target_column='nominal_length')
    result = apply_select_aggregate_rules(frame, [old_rule])
    assert result.iloc[0]['nominal_length'] == '150'
    assert result.iloc[0]['mean_insert_size'] == '200'
    with pytest.raises(ValueError, match='conflicting'):
        resolve_fragment_distribution(args(), result.iloc[0], 'R1', {})
    for name in ('plantae', 'vertebrate', 'test'):
        rules = read_select_rules(Path(__file__).parents[1] / 'amalgkit' / 'select_rule_sets' / name / 'select_rules.tsv')
        assert not any(rule['stage'] == 'aggregate' and rule['target_column'] == 'nominal_length' for rule in rules)


def test_paired_and_oarfish_do_not_use_single_end_values(capsys):
    data = metadata(lib_layout='paired', nominal_length='unusable')
    common = args(fragment_length_mean=150, fragment_length_sd=7)
    model = resolve_run_fragment_model(common, data, 'R1', 'paired')
    assert model['source'] == 'kallisto_paired_estimation'
    assert 'ignoring common fragment length' in capsys.readouterr().err
    assert resolve_run_fragment_model(args(fragment_length_policy='error'), data, 'R1', 'single', backend='oarfish') is None
    for layout, backend in [('paired', 'kallisto'), ('single', 'oarfish')]:
        with pytest.raises(ValueError, match='only to single-end kallisto'):
            resolve_run_fragment_model(args(_fragment_length_by_run={'R1': {}}), data, 'R1', layout, backend=backend)


def write_completed(tmp_path, model=None):
    directory = tmp_path / 'quant' / 'R1'
    directory.mkdir(parents=True, exist_ok=True)
    info = {'p_pseudoaligned': 100}
    if model is not None:
        info[PROVENANCE_KEY] = model
    (directory / 'R1_run_info.json').write_text(json.dumps(info))
    (directory / 'R1_abundance.tsv').write_text('target_id\tlength\teff_length\test_counts\ttpm\nt1\t500\t350\t10\t1000000\n')
    return directory


def test_reuse_checks_distribution_and_provenance_without_fastqs(tmp_path):
    data = metadata(nominal_length=150, nominal_sdev=7)
    model = resolve_run_fragment_model(args(), data, 'R1', 'single')
    directory = write_completed(tmp_path, model)
    runtime = args(out_dir=str(tmp_path), redo=False, clean_fastq=True)
    run_quant(runtime, data, 'R1', 'unused.idx')
    original = (directory / 'R1_run_info.json').read_bytes()
    data.df['nominal_sdev'] = 9
    with pytest.raises(ValueError, match='settings differ.*--redo yes'):
        run_quant(runtime, data, 'R1', 'unused.idx')
    assert (directory / 'R1_run_info.json').read_bytes() == original


def test_legacy_reuse_warns_but_cannot_certify_new_explicit_settings(tmp_path, capsys):
    directory = write_completed(tmp_path)
    data = metadata(nominal_length=150)
    check_fragment_model_reuse(args(), data, 'R1', str(directory))
    assert 'unknown fragment length provenance' in capsys.readouterr().err
    for runtime in [args(fragment_length_policy='error'), args(fragment_length_mean=150, fragment_length_sd=7)]:
        with pytest.raises(ValueError, match='unknown.*--redo yes'):
            check_fragment_model_reuse(runtime, data, 'R1', str(directory))
    data.df['nominal_length'] = 'malformed'
    with pytest.raises(ValueError, match='finite number'):
        check_fragment_model_reuse(args(), data, 'R1', str(directory))


@pytest.mark.parametrize('change', ['invalid_sd', 'no_source', 'strict_assumption', 'schema'])
def test_consumers_reject_malformed_fragment_provenance(tmp_path, change):
    model = resolve_run_fragment_model(args(), metadata(), 'R1', 'single')
    if change == 'invalid_sd':
        model['sd']['value'] = 0
    elif change == 'no_source':
        model['mean']['source'] = ''
    elif change == 'strict_assumption':
        model['policy'] = 'error'
    else:
        model['schema_version'] = 2
    directory = write_completed(tmp_path, model)
    assert validate_quant_run_info_json(str(directory / 'R1_run_info.json'))
    with pytest.raises(ValueError):
        collect_quant_models(['R1'], [str(directory / 'R1_abundance.tsv')])


def test_cli_defaults_and_strict_mode():
    parser = build_main_parser()
    default = parser.parse_args(['quant', '--out_dir', 'out'])
    assert default.fragment_length_policy == 'assume'
    assert default.fragment_length_mean is default.fragment_length_sd is default.fragment_length_file is None
    strict = parser.parse_args(['quant', '--out_dir', 'out', '--fragment_length_policy', 'error',
                               '--fragment_length_mean', '150', '--fragment_length_sd', '7'])
    validate_fragment_args(strict)


@pytest.mark.parametrize('batch,expected_run,expected_mean', [(1, 'R1', 140), (2, 'R3', 160)])
def test_real_batch_loading_validates_full_fragment_file(tmp_path, monkeypatch, batch, expected_run, expected_mean):
    rows = [metadata(run=run, is_sampled=sampled).df.iloc[0].to_dict()
            for run, sampled in [('R1', 'yes'), ('R2', 'no'), ('R3', 'yes')]]
    metadata_path = tmp_path / 'input.tsv'
    pandas.DataFrame(rows).to_csv(metadata_path, sep='\t', index=False)
    path = write_run_file(tmp_path, [[run, mean, 7, 'user', 'test library']
                                    for run, mean in [('R1', 140), ('R2', 150), ('R3', 160)]])
    runtime = build_main_parser().parse_args(['quant', '--out_dir', str(tmp_path), '--metadata', str(metadata_path),
        '--fragment_length_file', str(path), '--batch', str(batch), '--threads', '1'])
    seen = []

    def dispatch(loaded_args, loaded_metadata, run, species, runtime_context=None):
        assert loaded_metadata.df['run'].tolist() == [expected_run]
        model = resolve_run_fragment_model(loaded_args, loaded_metadata, run, 'single')
        seen.append((run, model['mean']['value']))

    monkeypatch.setattr('amalgkit.quant.check_quant_dependencies', lambda *a, **k: None)
    monkeypatch.setattr('amalgkit.quant.prepare_quant_runtime_context', lambda *a, **k: QuantRuntimeContext())
    monkeypatch.setattr('amalgkit.quant.run_quant_for_sra', dispatch)
    quant_main(runtime)
    assert seen == [(expected_run, expected_mean)]
    assert runtime.batch == batch and not hasattr(runtime, '_fragment_length_by_run')
    # A typo outside this batch must still be caught against the full input.
    with path.open('a') as handle:
        handle.write('unknown\t180\t7\tuser\ttypo\n')
    with pytest.raises(ValueError, match='Unknown'):
        quant_main(runtime)


def test_getfastq_preserves_invalid_fragment_text_for_quant_validation():
    data = metadata(nominal_length='NaN', nominal_sdev='7')
    initialize_columns(data, {'num_bp_per_sra': 1000})
    assert data.df.loc[0, 'nominal_length'] == 'NaN'
    assert data.df.loc[0, 'nominal_sdev'] == '7'
    with pytest.raises(ValueError, match='finite number'):
        resolve_run_fragment_model(args(), data, 'R1', 'single')


def test_reuse_detects_actual_input_layout_change_and_retired_corrections(tmp_path):
    data = metadata(lib_layout='paired', nominal_length=150, nominal_sdev=7)
    model = resolve_run_fragment_model(args(), data, 'R1', 'single')
    directory = write_completed(tmp_path, model)
    runtime = args(out_dir=str(tmp_path))
    # Cleanup may leave no reads or a single retirement marker; preserve the
    # original paired-metadata -> single-input correction in either case.
    check_fragment_model_reuse(runtime, data, 'R1', str(directory))
    inputs = tmp_path / 'getfastq' / 'R1'
    inputs.mkdir(parents=True)
    marker = inputs / 'R1.fastq.gz.safely_removed'
    marker.touch()
    check_fragment_model_reuse(runtime, data, 'R1', str(directory))
    marker.unlink()
    (inputs / 'R1_1.fastq.gz').touch()
    (inputs / 'R1_2.fastq.gz').touch()
    with pytest.raises(ValueError, match='settings differ'):
        check_fragment_model_reuse(runtime, data, 'R1', str(directory))


@pytest.mark.parametrize('old_backend,new_backend', [('kallisto', 'oarfish'), ('oarfish', 'kallisto')])
def test_reuse_never_returns_output_from_a_different_backend(tmp_path, old_backend, new_backend):
    data = metadata(nominal_length=150, nominal_sdev=7)
    directory = write_completed(tmp_path, resolve_run_fragment_model(args(), data, 'R1', 'single'))
    info_path = directory / 'R1_run_info.json'
    info = json.loads(info_path.read_text())
    info['quant_backend'] = old_backend
    if old_backend == 'oarfish':
        info.pop(PROVENANCE_KEY)
    info_path.write_text(json.dumps(info))
    with pytest.raises(ValueError, match='backend differs'):
        check_fragment_model_reuse(args(quant_backend=new_backend), data, 'R1', str(directory))


def test_single_overhang_change_requires_redo(tmp_path):
    data = metadata(nominal_length=150, nominal_sdev=7)
    directory = write_completed(tmp_path, resolve_run_fragment_model(args(), data, 'R1', 'single'))
    with pytest.raises(ValueError, match='settings differ'):
        check_fragment_model_reuse(args(kallisto_options='--single-overhang'), data, 'R1', str(directory))


def test_legacy_paired_output_ignores_single_end_strict_policy(tmp_path):
    directory = write_completed(tmp_path)
    check_fragment_model_reuse(args(fragment_length_policy='error'), metadata(lib_layout='paired'), 'R1', str(directory))


@pytest.mark.parametrize('corruption', ['policy_list', 'schema_bool', 'value_list', 'value_dict', 'different_supplied', 'false_assumption'])
def test_corrupt_provenance_is_an_error_not_a_validator_crash(tmp_path, corruption):
    model = resolve_run_fragment_model(args(), metadata(), 'R1', 'single')
    if corruption == 'policy_list':
        model['policy'] = []
    elif corruption == 'schema_bool':
        model['schema_version'] = True
    elif corruption == 'value_list':
        model['mean']['value'] = [200]
    elif corruption == 'value_dict':
        model['sd']['value'] = {'value': 20}
    elif corruption == 'different_supplied':
        model['mean'].update(source='cli', supplied_value='250')
    else:
        model['sd']['value'] = 21
    directory = write_completed(tmp_path, model)
    assert validate_quant_run_info_json(str(directory / 'R1_run_info.json'))
    with pytest.raises(ValueError):
        collect_quant_models(['R1'], [str(directory / 'R1_abundance.tsv')])


def test_generic_rerun_refuses_to_replace_prior_cli_values(tmp_path, monkeypatch):
    data = metadata(nominal_length=200, nominal_sdev=20)
    old_args = args(fragment_length_mean=150, fragment_length_sd=7)
    directory = write_completed(tmp_path, resolve_run_fragment_model(old_args, data, 'R1', 'single'))
    # Abundance is damaged, which is why the generic repair was selected.
    (directory / 'R1_abundance.tsv').write_text('damaged output')
    before = (directory / 'R1_run_info.json').read_bytes()

    def dispatch(runtime):
        assert runtime._preserve_quant_settings
        run_quant(runtime, data, 'R1', 'unused.idx')

    monkeypatch.setattr('amalgkit.rerun.quant_main', dispatch)
    with pytest.raises(ValueError, match='rerun cannot safely recover'):
        rerun_quant_check(args(out_dir=str(tmp_path), redo=True), data, ['R1'])
    assert (directory / 'R1_run_info.json').read_bytes() == before
    assert (directory / 'R1_abundance.tsv').read_text() == 'damaged output'
