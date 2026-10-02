import json
import os
from types import SimpleNamespace

import pandas
import pytest

from amalgkit.fragment_length import PROVENANCE_KEY as FRAGMENT_KEY
from amalgkit.metadata_utils import Metadata
from amalgkit.quant import resolve_run_fragment_model, run_quant
from amalgkit.quant_provenance import PROVENANCE_KEY, validate_quant_provenance


@pytest.fixture(params=['kallisto', 'oarfish'])
def completed_quant(tmp_path, monkeypatch, request):
    backend = request.param
    inputs = tmp_path / 'getfastq' / 'R1'
    inputs.mkdir(parents=True)
    reads = inputs / 'R1.fastq'
    reads.write_text('@r1\nAAAA\n+\nIIII\n')
    index = tmp_path / 'reference.idx'
    index.write_bytes(b'reference A')
    data = Metadata.from_DataFrame(pandas.DataFrame([dict(
        run='R1', scientific_name='Species A', lib_layout='single',
        total_spots=1, total_bases=4, spot_length=4, nominal_length=200, nominal_sdev=20,
    )]))
    args = SimpleNamespace(out_dir=str(tmp_path), redo=False, clean_fastq=False, threads=1,
                           quant_backend=backend, oarfish_seq_tech='ont-cdna')
    calls = []

    def quantify(args, in_files, metadata, stat, output, index, *extra):
        calls.append(in_files)
        info = {'p_pseudoaligned': 100, 'quant_backend': backend}
        if backend == 'kallisto':
            info[FRAGMENT_KEY] = resolve_run_fragment_model(args, metadata, 'R1', 'single')
        else:
            info.update(oarfish_seq_tech='ont-cdna', oarfish_options=[])
        with open(os.path.join(output, 'R1_run_info.json'), 'w') as handle:
            json.dump(info, handle)
        with open(os.path.join(output, 'R1_abundance.tsv'), 'w') as handle:
            handle.write('target_id\tlength\teff_length\test_counts\ttpm\ng1\t100\t90\t1\t1000000\n')

    monkeypatch.setattr('amalgkit.quant.call_' + backend, quantify)
    run_quant(args, data, 'R1', str(index), backend=backend, oarfish_seq_tech='ont-cdna')
    return SimpleNamespace(args=args, data=data, index=index, reads=reads, calls=calls,
                           output=tmp_path / 'quant' / 'R1', backend=backend, quantify=quantify)


def rerun(fixture):
    run_quant(fixture.args, fixture.data, 'R1', str(fixture.index), backend=fixture.backend,
              oarfish_seq_tech='ont-cdna')


def test_unchanged_inputs_are_reused(completed_quant):
    rerun(completed_quant)
    assert len(completed_quant.calls) == 1
    info = json.loads((completed_quant.output / 'R1_run_info.json').read_text())
    assert validate_quant_provenance(info[PROVENANCE_KEY]) == ''


@pytest.mark.parametrize('changed', ['reads', 'index'])
def test_changed_content_requires_explicit_redo(completed_quant, changed):
    fixture = completed_quant
    path = getattr(fixture, changed)
    original_stat = path.stat()
    content = path.read_bytes()
    path.write_bytes(content.replace(b'A', b'C'))
    os.utime(path, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
    original = (fixture.output / 'R1_run_info.json').read_bytes()
    with pytest.raises(ValueError, match='differs? from existing output.*--redo yes'):
        rerun(fixture)
    assert len(fixture.calls) == 1
    assert (fixture.output / 'R1_run_info.json').read_bytes() == original
    fixture.args.redo = True
    rerun(fixture)
    assert len(fixture.calls) == 2


def test_cleanup_retains_verifiable_reuse(completed_quant):
    fixture = completed_quant
    fixture.args.redo = True
    fixture.args.clean_fastq = True
    rerun(fixture)
    assert not fixture.reads.exists()
    assert fixture.reads.with_name(fixture.reads.name + '.safely_removed').is_file()
    fixture.args.redo = False
    rerun(fixture)
    assert len(fixture.calls) == 2
    fixture.index.write_bytes(b'changed reference')
    with pytest.raises(ValueError, match='reference index differs'):
        rerun(fixture)


def test_missing_fastqs_without_cleanup_markers_cannot_certify_reuse(completed_quant):
    completed_quant.reads.unlink()
    with pytest.raises(ValueError, match='provenance cannot be verified'):
        rerun(completed_quant)


def test_legacy_output_remains_readable_but_requires_redo(completed_quant):
    fixture = completed_quant
    path = fixture.output / 'R1_run_info.json'
    info = json.loads(path.read_text())
    del info[PROVENANCE_KEY]
    path.write_text(json.dumps(info))
    from amalgkit.output_contracts import validate_quant_output_files
    assert validate_quant_output_files('R1', str(fixture.output))[0]
    with pytest.raises(ValueError, match='unknown or invalid input/reference provenance.*--redo yes'):
        rerun(fixture)


@pytest.mark.parametrize('changed', ['reads', 'index'])
def test_mutation_during_quant_rolls_back(completed_quant, monkeypatch, changed):
    fixture = completed_quant
    original = (fixture.output / 'R1_run_info.json').read_bytes()
    fixture.args.redo = True

    def mutate(*args):
        fixture.quantify(*args)
        getattr(fixture, changed).write_bytes(b'changed during quant')

    monkeypatch.setattr('amalgkit.quant.call_' + fixture.backend, mutate)
    with pytest.raises(ValueError, match='changed during quantification'):
        rerun(fixture)
    assert (fixture.output / 'R1_run_info.json').read_bytes() == original


def test_run_identity_is_lexical(completed_quant):
    path = completed_quant.output / 'R1_run_info.json'
    info = json.loads(path.read_text())
    info[PROVENANCE_KEY]['run'] = 'r1'
    info[PROVENANCE_KEY]['inputs'][0]['name'] = 'r1.fastq'
    path.write_text(json.dumps(info))
    with pytest.raises(ValueError, match='belongs to run r1'):
        rerun(completed_quant)


@pytest.mark.parametrize('field,value', [('run', ''), ('run', '../R1'), ('inputs', []),
                                       ('inputs', [{'name': '../reads.fastq'}]), ('index', None),
                                       ('index', {'size': -1, 'sha256': '0' * 64}),
                                       ('index', {'size': 10, 'sha256': 'wrong'})])
def test_invalid_provenance_is_rejected(completed_quant, field, value):
    info = json.loads((completed_quant.output / 'R1_run_info.json').read_text())[PROVENANCE_KEY]
    info[field] = value
    assert validate_quant_provenance(info)
