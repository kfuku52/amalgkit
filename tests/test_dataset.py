import os

import pytest

from amalgkit import dataset as dataset_module
from amalgkit.dataset import (
    extract_dataset,
    validate_dataset_name,
    resolve_dataset_source_dir,
)


def _prepare_fake_dataset(tmp_path, monkeypatch):
    dataset_root = tmp_path / 'datasets'
    source_dir = dataset_root / 'toy'
    source_dir.mkdir(parents=True)
    (source_dir / 'toy.fa.gz').write_text('FA')
    (source_dir / 'toy_busco.tsv').write_text('BUSCO')
    (source_dir / 'toy.config').write_text('CONF')
    monkeypatch.setattr(dataset_module, 'DATASETS', {
        'toy': {
            'description': 'toy dataset',
            'species': ['Toy_species'],
            'files': {
                'fasta': ['toy.fa.gz'],
                'busco': ['toy_busco.tsv'],
                'config': ['toy.config'],
            },
        },
    })
    monkeypatch.setattr(dataset_module, 'get_dataset_dir', lambda: str(dataset_root))
    return source_dir


class TestDatasetExtract:
    def test_extract_dataset_copies_files(self, tmp_path, monkeypatch):
        _prepare_fake_dataset(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'

        paths = extract_dataset(name='toy', out_dir=str(out_dir), overwrite=False)

        assert os.path.exists(os.path.join(paths['fasta'], 'toy.fa.gz'))
        assert os.path.exists(os.path.join(paths['busco'], 'toy_busco.tsv'))
        assert os.path.exists(os.path.join(paths['config'], 'toy.config'))

    def test_extract_dataset_skips_existing_when_not_overwrite(self, tmp_path, monkeypatch):
        _prepare_fake_dataset(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'
        existing_fasta = out_dir / 'fasta' / 'toy.fa.gz'
        existing_fasta.parent.mkdir(parents=True)
        existing_fasta.write_text('OLD')

        extract_dataset(name='toy', out_dir=str(out_dir), overwrite=False)

        assert existing_fasta.read_text() == 'OLD'

    def test_extract_dataset_overwrites_when_requested(self, tmp_path, monkeypatch):
        _prepare_fake_dataset(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'
        existing_fasta = out_dir / 'fasta' / 'toy.fa.gz'
        existing_fasta.parent.mkdir(parents=True)
        existing_fasta.write_text('OLD')

        extract_dataset(name='toy', out_dir=str(out_dir), overwrite=True)

        assert existing_fasta.read_text() == 'FA'


class TestDatasetValidation:
    def test_validate_dataset_name_unknown_exits(self):
        with pytest.raises(SystemExit):
            validate_dataset_name('unknown_dataset')

    def test_resolve_dataset_source_dir_missing_exits(self, tmp_path, monkeypatch):
        monkeypatch.setattr(dataset_module, 'get_dataset_dir', lambda: str(tmp_path / 'missing_root'))
        with pytest.raises(SystemExit):
            resolve_dataset_source_dir('toy')
