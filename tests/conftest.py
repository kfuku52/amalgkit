import pytest
import pandas
from pathlib import Path

from amalgkit.util import Metadata


@pytest.fixture
def sample_metadata_df():
    """A small realistic metadata DataFrame."""
    data = {
        'scientific_name': ['Homo sapiens', 'Homo sapiens', 'Homo sapiens',
                            'Mus musculus', 'Mus musculus'],
        'sample_group': ['brain', 'brain', 'liver', 'brain', 'liver'],
        'tissue': ['brain', 'brain', 'liver', 'brain', 'liver'],
        'genotype': ['', '', '', '', ''],
        'sex': ['male', 'female', 'male', 'female', 'male'],
        'age': ['', '', '', '', ''],
        'treatment': ['', '', '', '', ''],
        'source_name': ['', '', '', '', ''],
        'is_sampled': ['yes', 'yes', 'yes', 'yes', 'yes'],
        'is_qualified': ['yes', 'yes', 'yes', 'yes', 'yes'],
        'exclusion': ['no', 'no', 'no', 'no', 'no'],
        'protocol': ['', '', '', '', ''],
        'bioproject': ['PRJNA1', 'PRJNA1', 'PRJNA2', 'PRJNA3', 'PRJNA3'],
        'biosample': ['SAMN1', 'SAMN2', 'SAMN3', 'SAMN4', 'SAMN5'],
        'experiment': ['SRX1', 'SRX2', 'SRX3', 'SRX4', 'SRX5'],
        'run': ['SRR001', 'SRR002', 'SRR003', 'SRR004', 'SRR005'],
        'sra_primary': ['SRA1', 'SRA2', 'SRA3', 'SRA4', 'SRA5'],
        'sra_sample': ['SRS1', 'SRS2', 'SRS3', 'SRS4', 'SRS5'],
        'sra_study': ['SRP1', 'SRP1', 'SRP2', 'SRP3', 'SRP3'],
        'study_title': ['Study1', 'Study1', 'Study2', 'Study3', 'Study3'],
        'exp_title': ['', '', '', '', ''],
        'design': ['', '', '', '', ''],
        'sample_title': ['', '', '', '', ''],
        'sample_description': ['', '', '', '', ''],
        'lib_name': ['', '', '', '', ''],
        'lib_layout': ['paired', 'paired', 'single', 'paired', 'single'],
        'lib_strategy': ['RNA-Seq', 'RNA-Seq', 'RNA-Seq', 'RNA-Seq', 'RNA-Seq'],
        'lib_source': ['TRANSCRIPTOMIC', 'TRANSCRIPTOMIC', 'TRANSCRIPTOMIC',
                        'TRANSCRIPTOMIC', 'TRANSCRIPTOMIC'],
        'lib_selection': ['cDNA', 'cDNA', 'cDNA', 'cDNA', 'cDNA'],
        'instrument': ['Illumina HiSeq 2500'] * 5,
        'total_spots': [10000000, 5000000, 100, 8000000, 200],
        'total_bases': [2000000000, 1000000000, 20000, 1600000000, 40000],
        'size': [500000000, 250000000, 5000, 400000000, 10000],
        'nominal_length': [200, 200, 0, 200, 0],
        'nominal_sdev': [0, 0, 0, 0, 0],
        'spot_length': [200, 200, 100, 200, 100],
        'read_index': ['', '', '', '', ''],
        'read_class': ['', '', '', '', ''],
        'read_type': ['', '', '', '', ''],
        'base_coord': ['', '', '', '', ''],
        'center': ['', '', '', '', ''],
        'submitter_id': ['', '', '', '', ''],
        'pubmed_id': ['', '', '', '', ''],
        'taxid': [9606, 9606, 9606, 10090, 10090],
        'published_date': ['2020-01-01', '2020-02-01', '2020-03-01',
                           '2020-04-01', '2020-05-01'],
        'NCBI_Link': ['', '', '', '', ''],
        'AWS_Link': ['', '', '', '', ''],
        'GCP_Link': ['', '', '', '', ''],
        'ENA_SRA_Link': ['', '', '', '', ''],
        'DDBJ_SRA_Link': ['', '', '', '', ''],
    }
    return pandas.DataFrame(data)


@pytest.fixture
def sample_metadata(sample_metadata_df):
    """A Metadata object built from sample data."""
    return Metadata.from_DataFrame(sample_metadata_df)


@pytest.fixture
def stub_pdf_rendering(monkeypatch):
    """Keep integration tests focused on plot orchestration, not PDF encoding."""
    from matplotlib.figure import Figure
    from amalgkit import cross_species_filter
    from amalgkit import per_species_finalize_python
    from amalgkit import per_species_python

    def find_output_path(args, kwargs):
        output_path = kwargs.get('out_pdf_path')
        if output_path is not None:
            return output_path
        for value in reversed(args):
            if isinstance(value, (str, Path)) and str(value).lower().endswith('.pdf'):
                return value
        raise AssertionError('Plot helper was called without a PDF output path.')

    def write_plot_placeholder(*args, **kwargs):
        output_path = find_output_path(args, kwargs)
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'%PDF-1.4\n% amalgkit test placeholder\n')
        return str(path)

    def write_pdf_placeholder(_figure, output_path, *args, **kwargs):
        _ = (args, kwargs)
        if hasattr(output_path, 'write'):
            output_path.write(b'%PDF-1.4\n% amalgkit test placeholder\n')
            return
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'%PDF-1.4\n% amalgkit test placeholder\n')

    monkeypatch.setattr(Figure, 'savefig', write_pdf_placeholder)
    for module, helper_names in (
        (
            cross_species_filter,
            (
                '_save_sample_number_heatmap_pdf',
                '_save_group_cor_scatter_plot',
                '_save_overview_pdf',
                '_save_csfilter_scatter_plot',
                '_save_heatmap_pdf',
                '_save_within_group_histogram',
                '_save_pca_pdf',
                '_save_unaveraged_tsne_pdf',
                '_save_averaged_heatmap_pdf',
                '_save_averaged_dendrogram_pdf',
                '_save_averaged_summary_pdf',
                '_save_averaged_boxplot_pdf',
                '_save_averaged_tsne_pdf',
                '_save_delta_pcc_plot',
                'save_exclusion_plot_pdf',
            ),
        ),
        (
            per_species_python,
            ('_save_ws_scatter_plot', 'save_state_overview_pdf', 'save_tau_histogram_pdf'),
        ),
        (
            per_species_finalize_python,
            ('save_quick_state_comparison_plot', 'save_tau_histogram_pdf'),
        ),
    ):
        for helper_name in helper_names:
            monkeypatch.setattr(module, helper_name, write_plot_placeholder)


@pytest.fixture(autouse=True)
def stub_pdf_rendering_for_integration_tests(request):
    """Use cheap PDF placeholders except in explicitly slow rendering tests."""
    is_integration = request.node.get_closest_marker('integration') is not None
    is_real_rendering = request.node.get_closest_marker('slow') is not None
    if is_integration and not is_real_rendering:
        request.getfixturevalue('stub_pdf_rendering')
