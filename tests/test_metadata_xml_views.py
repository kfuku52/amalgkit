import copy
import xml.etree.ElementTree as ET

import pandas
import pytest

from amalgkit import metadata_utils


def _legacy_views(root):
    """Original full-copy implementation, retained as an independent oracle."""
    for package in root.iter(tag='EXPERIMENT_PACKAGE'):
        runs = package.findall('./RUN_SET/RUN')
        if len(runs) <= 1:
            yield package
            continue
        for run_index in range(len(runs)):
            package_copy = copy.deepcopy(package)
            run_set_copy = package_copy.find('./RUN_SET')
            for copied_index, copied_run in enumerate(run_set_copy.findall('./RUN')):
                if copied_index != run_index:
                    run_set_copy.remove(copied_run)
            yield package_copy


def _xml(run_sets):
    return ET.fromstring(
        '<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE accession="ERP1">package text'
        '<EXPERIMENT><IDENTIFIERS><PRIMARY_ID>SRX000001</PRIMARY_ID></IDENTIFIERS>'
        '<DESIGN><LIBRARY_DESCRIPTOR><LIBRARY_LAYOUT><PAIRED NOMINAL_LENGTH="300" '
        'NOMINAL_SDEV="20"/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN></EXPERIMENT>experiment tail'
        '<SAMPLE><SAMPLE_NAME><SCIENTIFIC_NAME>Species alpha</SCIENTIFIC_NAME></SAMPLE_NAME>'
        '<SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE><TAG>tissue</TAG><VALUE>leaf</VALUE></SAMPLE_ATTRIBUTE>'
        '<SAMPLE_ATTRIBUTE><TAG>tissue</TAG><VALUE>young leaf</VALUE></SAMPLE_ATTRIBUTE>'
        '<SAMPLE_ATTRIBUTE><TAG>run</TAG><VALUE>submitted run</VALUE></SAMPLE_ATTRIBUTE>'
        '<SAMPLE_ATTRIBUTE><TAG>nominal_length</TAG><VALUE>350</VALUE></SAMPLE_ATTRIBUTE>'
        '</SAMPLE_ATTRIBUTES></SAMPLE>sample tail'
        + run_sets + '<LABEL>after run sets</LABEL></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>'
    )


def _run(accession):
    return (
        '<RUN accession="' + accession + '" total_spots="12" total_bases="2400">run text'
        '<IDENTIFIERS><PRIMARY_ID>' + accession + '</PRIMARY_ID>'
        '<EXTERNAL_ID namespace="BioProject">PRJNA' + accession + '</EXTERNAL_ID>'
        '<EXTERNAL_ID namespace="BioSample">SAMN' + accession + '</EXTERNAL_ID></IDENTIFIERS>'
        '<SRAFiles><SRAFile supertype="Primary ETL">'
        '<Alternatives org="NCBI" url="https://example.org/' + accession + '.sra"/>'
        '</SRAFile></SRAFiles></RUN>run tail'
    )


@pytest.mark.parametrize('run_sets', [
    '',
    '<RUN_SET source="SRA"><LABEL>empty</LABEL></RUN_SET>',
    '<RUN_SET source="SRA">' + _run('SRR000001') + '</RUN_SET>',
    '<RUN_SET source="SRA">set text<LABEL>before</LABEL>' + _run('SRR000001')
    + '<LABEL>between</LABEL>' + _run('SRR000002') + '<LABEL>after</LABEL></RUN_SET>set tail',
    '<RUN_SET>' + _run('SRR000001') + _run('SRR000002') + '</RUN_SET>'
    + '<RUN_SET>' + _run('SRR000003') + '</RUN_SET>',
])
def test_run_views_match_full_copy_oracle(run_sets):
    root = _xml(run_sets)
    original = ET.tostring(root)
    expected = [ET.tostring(package) for package in _legacy_views(root)]
    actual = [ET.tostring(package) for package in metadata_utils.iter_experiment_packages_per_run(root)]
    assert actual == expected
    assert ET.tostring(root) == original
    if len(root.findall('./EXPERIMENT_PACKAGE/RUN_SET/RUN')) <= 1:
        assert next(metadata_utils.iter_experiment_packages_per_run(root)) is root[0]


def test_multi_run_views_are_independent_deep_copies():
    root = _xml('<RUN_SET source="SRA">' + _run('SRR000001') + _run('SRR000002') + '</RUN_SET>')
    original = ET.tostring(root)
    first, second = metadata_utils.iter_experiment_packages_per_run(root)
    second_original = ET.tostring(second)
    for element in first.iter():
        element.set('mutated', 'yes')
        element.text = 'changed text'
        element.tail = 'changed tail'
    first.find('./RUN_SET').clear()
    first.remove(first.find('./SAMPLE'))
    assert ET.tostring(second) == second_original
    assert ET.tostring(root) == original


def test_metadata_run_views_preserve_complete_output_and_collision_info(monkeypatch):
    root = _xml('<RUN_SET source="SRA">' + _run('SRR000001') + '<LABEL>between</LABEL>'
                + _run('SRR000002') + _run('SRR000003') + '</RUN_SET>')
    original = ET.tostring(root)
    actual = metadata_utils.Metadata.from_xml(root)
    monkeypatch.setattr(metadata_utils, 'iter_experiment_packages_per_run', _legacy_views)
    expected = metadata_utils.Metadata.from_xml(root)
    pandas.testing.assert_frame_equal(actual.df, expected.df)
    assert actual.sample_attribute_collision_count == expected.sample_attribute_collision_count
    assert actual.sample_attribute_collision_examples == expected.sample_attribute_collision_examples
    assert actual.sample_attribute_collision_count > 0
    assert actual.df['run'].tolist() == [
        'SRR000001 | submitted run', 'SRR000002 | submitted run', 'SRR000003 | submitted run',
    ]
    assert actual.df['bioproject'].tolist() == ['PRJNASRR000001', 'PRJNASRR000002', 'PRJNASRR000003']
    assert actual.df['biosample'].tolist() == ['SAMNSRR000001', 'SAMNSRR000002', 'SAMNSRR000003']
    assert actual.df['NCBI_Link'].tolist() == [
        'https://example.org/SRR000001.sra', 'https://example.org/SRR000002.sra', 'https://example.org/SRR000003.sra',
    ]
    assert ET.tostring(root) == original
