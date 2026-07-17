import io
import hashlib
import os
import shutil
import sqlite3
import stat
import tarfile
import time
from types import SimpleNamespace

import pytest

import amalgkit.download_utils as download_utils
from amalgkit.ncbi_taxonomy import (
    LEGACY_ETE4_DB_FORMAT,
    NATIVE_DB_FORMAT,
    NcbiTaxonomy,
    build_taxonomy_database,
    cleanup_stale_taxonomy_build_files,
    detect_taxonomy_database_format,
    is_taxonomy_database_compatible,
    native_database_matches_taxdump,
)


def _dmp_line(*fields):
    return "\t|\t".join(str(field) for field in fields) + "\t|\n"


def _write_tar_member(tar, name, text):
    data = text.encode("utf-8")
    info = tarfile.TarInfo(name=name)
    info.size = len(data)
    tar.addfile(info, io.BytesIO(data))


def write_tiny_taxdump(path, include_names=True, human_scientific_name="Homo sapiens"):
    nodes = "".join(
        [
            _dmp_line(1, 1, "no rank"),
            _dmp_line(2, 1, "domain"),
            _dmp_line(9605, 2, "genus"),
            _dmp_line(9606, 9605, "species"),
            _dmp_line(10088, 2, "genus"),
            _dmp_line(10090, 10088, "species"),
        ]
    )
    names = "".join(
        [
            _dmp_line(1, "root", "", "scientific name"),
            _dmp_line(2, "Bacteria", "", "scientific name"),
            _dmp_line(9605, "Homo", "", "scientific name"),
            _dmp_line(9606, human_scientific_name, "", "scientific name"),
            _dmp_line(10088, "Mus", "", "scientific name"),
            _dmp_line(10090, "Mus musculus", "", "scientific name"),
            _dmp_line(9606, "human", "", "synonym"),
            _dmp_line(9606, "shared alias", "", "synonym"),
            _dmp_line(10090, "shared alias", "", "equivalent name"),
            _dmp_line(9606, "Bacteria", "", "synonym"),
            _dmp_line(10090, "mouse", "", "genbank common name"),
        ]
    )
    merged = "".join(
        [
            _dmp_line(998, 999),
            _dmp_line(999, 9606),
        ]
    )
    with tarfile.open(path, "w:gz") as tar:
        _write_tar_member(tar, "nodes.dmp", nodes)
        if include_names:
            _write_tar_member(tar, "names.dmp", names)
        _write_tar_member(tar, "merged.dmp", merged)
    return path


def write_legacy_ete_database(path):
    db = sqlite3.connect(path)
    db.executescript(
        """
        CREATE TABLE stats (version INT PRIMARY KEY);
        CREATE TABLE species (
            taxid INT PRIMARY KEY,
            parent INT,
            spname TEXT COLLATE NOCASE,
            common TEXT COLLATE NOCASE,
            rank TEXT,
            track TEXT
        );
        CREATE TABLE synonym (
            taxid INT,
            spname TEXT COLLATE NOCASE,
            PRIMARY KEY (spname, taxid)
        );
        CREATE TABLE merged (taxid_old INT, taxid_new INT);
        INSERT INTO stats VALUES (2);
        INSERT INTO species VALUES (1, 0, 'root', '', 'no rank', '1');
        INSERT INTO species VALUES (2, 1, 'Bacteria', '', 'domain', '2,1');
        INSERT INTO species VALUES (9605, 2, 'Homo', '', 'genus', '9605,2,1');
        INSERT INTO species VALUES (9606, 9605, 'Homo sapiens', '', 'species', '9606,9605,2,1');
        INSERT INTO synonym VALUES (9606, 'human');
        INSERT INTO merged VALUES (999, 9606);
        """
    )
    db.commit()
    db.close()
    return path


def test_build_and_query_native_taxonomy_database(tmp_path):
    taxdump = write_tiny_taxdump(tmp_path / "taxdump.tar.gz")
    dbfile = tmp_path / "taxa.sqlite"

    result = build_taxonomy_database(str(dbfile), str(taxdump))

    assert result == str(dbfile)
    assert detect_taxonomy_database_format(str(dbfile)) == NATIVE_DB_FORMAT
    assert is_taxonomy_database_compatible(str(dbfile)) is True

    taxonomy = NcbiTaxonomy(str(dbfile), update=False)
    try:
        assert taxonomy.get_rank([1, 9605, 9606, None, ""]) == {
            1: "no rank",
            9605: "genus",
            9606: "species",
        }
        assert taxonomy.get_lineage(9606) == [1, 2, 9605, 9606]
        assert taxonomy.get_lineage_translator([9606, 10090]) == {
            9606: [1, 2, 9605, 9606],
            10090: [1, 2, 10088, 10090],
        }
        assert taxonomy.get_taxid_translator([9606, 10090]) == {
            9606: "Homo sapiens",
            10090: "Mus musculus",
        }
        assert taxonomy.get_name_translator(
            [
                "hOmO SaPiEnS",
                "human",
                "shared alias",
                "Bacteria",
                "mouse",
            ]
        ) == {
            "hOmO SaPiEnS": [9606],
            "human": [9606],
            "shared alias": [9606, 10090],
            "Bacteria": [2],
        }
        assert taxonomy.get_name_translator(["x') OR 1=1 --"]) == {}
    finally:
        taxonomy.close()


def test_native_taxonomy_resolves_merged_taxids(tmp_path):
    taxdump = write_tiny_taxdump(tmp_path / "taxdump.tar.gz")
    dbfile = tmp_path / "taxa.sqlite"
    build_taxonomy_database(str(dbfile), str(taxdump))
    taxonomy = NcbiTaxonomy(str(dbfile), update=False)
    try:
        assert taxonomy.get_lineage_translator([998, 999]) == {}
        with pytest.warns(UserWarning, match="taxid 999 was translated into 9606"):
            assert taxonomy.get_lineage(999) == [1, 2, 9605, 9606]
        with pytest.warns(UserWarning, match="taxid 998 was translated into 9606"):
            assert taxonomy.get_lineage(998) == [1, 2, 9605, 9606]
        assert taxonomy.get_taxid_translator([998, 999]) == {
            998: "Homo sapiens",
            999: "Homo sapiens",
        }
    finally:
        taxonomy.close()


def test_reads_existing_ete4_taxonomy_database(tmp_path):
    dbfile = write_legacy_ete_database(tmp_path / "taxa.sqlite")

    assert detect_taxonomy_database_format(str(dbfile)) == LEGACY_ETE4_DB_FORMAT
    taxonomy = NcbiTaxonomy(str(dbfile), update=False)
    try:
        assert taxonomy.get_lineage(9606) == [1, 2, 9605, 9606]
        assert taxonomy.get_rank([2, 9606]) == {2: "domain", 9606: "species"}
        assert taxonomy.get_taxid_translator([9606, 999]) == {
            9606: "Homo sapiens",
            999: "Homo sapiens",
        }
        assert taxonomy.get_name_translator(["HOMO SAPIENS", "human"]) == {
            "HOMO SAPIENS": [9606],
            "human": [9606],
        }
    finally:
        taxonomy.close()


def test_failed_rebuild_keeps_existing_database(tmp_path):
    valid_taxdump = write_tiny_taxdump(tmp_path / "valid.tar.gz")
    invalid_taxdump = write_tiny_taxdump(
        tmp_path / "invalid.tar.gz", include_names=False
    )
    dbfile = tmp_path / "taxa.sqlite"
    build_taxonomy_database(str(dbfile), str(valid_taxdump))

    with pytest.raises(ValueError, match="names.dmp"):
        build_taxonomy_database(str(dbfile), str(invalid_taxdump))

    taxonomy = NcbiTaxonomy(str(dbfile), update=False)
    try:
        assert taxonomy.get_taxid_translator([9606]) == {9606: "Homo sapiens"}
    finally:
        taxonomy.close()
    assert list(tmp_path.glob("taxa.sqlite.*.tmp")) == []


def test_rejects_corrupt_or_unsupported_database(tmp_path):
    dbfile = tmp_path / "taxa.sqlite"
    dbfile.write_bytes(b"not sqlite")

    assert detect_taxonomy_database_format(str(dbfile)) is None
    assert is_taxonomy_database_compatible(str(dbfile)) is False
    with pytest.raises(ValueError, match="Cannot open a supported"):
        NcbiTaxonomy(str(dbfile), update=False)


def test_builder_rejects_symlink_inputs_and_destination(tmp_path):
    taxdump = write_tiny_taxdump(tmp_path / "taxdump.tar.gz")
    taxdump_link = tmp_path / "taxdump-link.tar.gz"
    os.symlink(taxdump, taxdump_link)

    with pytest.raises(FileNotFoundError, match="must not be a symlink"):
        build_taxonomy_database(str(tmp_path / "taxa.sqlite"), str(taxdump_link))

    db_target = tmp_path / "target.sqlite"
    db_target.write_bytes(b"existing")
    db_link = tmp_path / "taxa-link.sqlite"
    os.symlink(db_target, db_link)
    with pytest.raises(IsADirectoryError, match="must not be a symlink"):
        build_taxonomy_database(str(db_link), str(taxdump))
    assert db_target.read_bytes() == b"existing"


def test_download_wrapper_builds_and_caches_native_database(tmp_path):
    source_taxdump = write_tiny_taxdump(tmp_path / "source-taxdump.tar.gz")
    args = SimpleNamespace(
        out_dir=str(tmp_path / "out"),
        download_dir=str(tmp_path / "downloads"),
    )
    cache = download_utils._get_thread_local_ncbi_taxonomy_cache()
    cache.clear()

    def copy_taxdump(_url, destination):
        shutil.copyfile(source_taxdump, destination)

    taxonomy = download_utils.get_ncbi_taxonomy(
        args=args,
        urlretrieve_fn=copy_taxdump,
    )
    try:
        assert taxonomy.get_taxid_translator([9606]) == {9606: "Homo sapiens"}
        assert download_utils.get_ncbi_taxonomy(args=args) is taxonomy
        data_dir = tmp_path / "downloads" / "ete_taxonomy"
        assert (
            detect_taxonomy_database_format(str(data_dir / "taxa.sqlite"))
            == NATIVE_DB_FORMAT
        )
        assert (data_dir / "taxdump.tar.gz").is_file()
    finally:
        taxonomy.close()
        cache.clear()


def test_default_data_dir_reuses_existing_ete_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path))
    legacy_dir = tmp_path / "ete"
    legacy_dir.mkdir()
    write_legacy_ete_database(legacy_dir / "taxa.sqlite")

    assert download_utils.resolve_default_ncbi_taxonomy_data_dir() == str(legacy_dir)

    preferred_dir = tmp_path / "amalgkit" / "ncbi_taxonomy"
    preferred_dir.mkdir(parents=True)
    assert download_utils.resolve_default_ncbi_taxonomy_data_dir() == str(legacy_dir)

    preferred_taxdump = write_tiny_taxdump(preferred_dir / "taxdump.tar.gz")
    build_taxonomy_database(preferred_dir / "taxa.sqlite", preferred_taxdump)
    assert download_utils.resolve_default_ncbi_taxonomy_data_dir() == str(preferred_dir)


@pytest.mark.parametrize("legacy", [False, True])
def test_format_detection_rejects_tables_with_missing_required_columns(
    tmp_path, legacy
):
    dbfile = tmp_path / "taxa.sqlite"
    db = sqlite3.connect(dbfile)
    if legacy:
        db.executescript(
            """
            CREATE TABLE stats (version INT);
            CREATE TABLE species (taxid INT, spname TEXT);
            CREATE TABLE synonym (taxid INT, spname TEXT);
            CREATE TABLE merged (taxid_old INT, taxid_new INT);
            INSERT INTO stats VALUES (2);
            INSERT INTO species VALUES (1, 'root');
            """
        )
    else:
        db.executescript(
            """
            CREATE TABLE metadata (key TEXT, value TEXT);
            CREATE TABLE nodes (wrong_taxid INT, parent_taxid INT, rank TEXT);
            CREATE TABLE names (taxid INT, name TEXT, is_scientific INT);
            CREATE TABLE merged (old_taxid INT, new_taxid INT);
            INSERT INTO metadata VALUES ('schema_version', '1');
            """
        )
    db.commit()
    db.close()

    assert detect_taxonomy_database_format(dbfile) is None
    assert is_taxonomy_database_compatible(dbfile) is False


@pytest.mark.skipif(os.name == "nt", reason="POSIX permission modes are required")
def test_builder_uses_shared_source_permissions_and_preserves_existing_mode(tmp_path):
    taxdump = write_tiny_taxdump(tmp_path / "taxdump.tar.gz")
    os.chmod(taxdump, 0o664)
    dbfile = tmp_path / "taxa.sqlite"

    build_taxonomy_database(dbfile, taxdump)
    assert stat.S_IMODE(os.stat(dbfile).st_mode) == 0o664

    os.chmod(dbfile, 0o640)
    build_taxonomy_database(dbfile, taxdump)
    assert stat.S_IMODE(os.stat(dbfile).st_mode) == 0o640


def test_native_database_freshness_uses_taxdump_fingerprint(tmp_path):
    taxdump = write_tiny_taxdump(tmp_path / "taxdump.tar.gz")
    dbfile = tmp_path / "taxa.sqlite"
    build_taxonomy_database(dbfile, taxdump)

    assert native_database_matches_taxdump(dbfile, taxdump) is True
    assert (
        download_utils.should_build_ncbi_taxonomy_db(dbfile, taxdump_file=taxdump)
        is False
    )

    write_tiny_taxdump(taxdump, human_scientific_name="Homo sapiens updated")
    assert native_database_matches_taxdump(dbfile, taxdump) is False
    assert (
        download_utils.should_build_ncbi_taxonomy_db(dbfile, taxdump_file=taxdump)
        is True
    )

    build_taxonomy_database(dbfile, taxdump)
    taxonomy = NcbiTaxonomy(dbfile, update=False)
    try:
        assert taxonomy.get_taxid_translator([9606]) == {9606: "Homo sapiens updated"}
    finally:
        taxonomy.close()


def test_builder_cleans_stale_but_not_recent_temp_files(tmp_path):
    taxdump = write_tiny_taxdump(tmp_path / "taxdump.tar.gz")
    dbfile = tmp_path / "taxa.sqlite"
    stale = tmp_path / "taxa.sqlite.abandoned.tmp"
    recent = tmp_path / "taxa.sqlite.active.tmp"
    stale.write_bytes(b"stale")
    recent.write_bytes(b"recent")
    old_time = time.time() - 7200
    os.utime(stale, (old_time, old_time))

    build_taxonomy_database(dbfile, taxdump)

    assert not stale.exists()
    assert recent.exists()
    assert cleanup_stale_taxonomy_build_files(dbfile, stale_seconds=0) == [str(recent)]


def test_download_wrapper_rejects_db_and_taxdump_symlinks(tmp_path):
    args = SimpleNamespace(
        out_dir=str(tmp_path / "out"),
        download_dir=str(tmp_path / "downloads"),
    )
    data_dir = tmp_path / "downloads" / "ete_taxonomy"
    data_dir.mkdir(parents=True)
    cache = download_utils._get_thread_local_ncbi_taxonomy_cache()
    cache.clear()

    db_victim = tmp_path / "db-victim.sqlite"
    os.symlink(db_victim, data_dir / "taxa.sqlite")
    with pytest.raises(IsADirectoryError, match="DB path exists but is not a file"):
        download_utils.get_ncbi_taxonomy(args=args)
    assert not db_victim.exists()

    os.unlink(data_dir / "taxa.sqlite")
    taxdump_victim = tmp_path / "taxdump-victim.tar.gz"
    os.symlink(taxdump_victim, data_dir / "taxdump.tar.gz")
    with pytest.raises(
        IsADirectoryError, match="taxdump path exists but is not a file"
    ):
        download_utils.get_ncbi_taxonomy(args=args)
    assert not taxdump_victim.exists()


def test_download_wrapper_rejects_symlinked_taxonomy_directory(tmp_path):
    download_dir = tmp_path / "downloads"
    download_dir.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    os.symlink(outside, download_dir / "ete_taxonomy")
    args = SimpleNamespace(
        out_dir=str(tmp_path / "out"), download_dir=str(download_dir)
    )

    with pytest.raises(NotADirectoryError, match="data directory"):
        download_utils.get_ncbi_taxonomy(args=args)
    assert list(outside.iterdir()) == []


def test_invalid_existing_taxdump_is_replaced_only_after_validation(tmp_path):
    source = write_tiny_taxdump(tmp_path / "source.tar.gz")
    destination = tmp_path / "taxdump.tar.gz"
    destination.write_bytes(b"broken")
    os.chmod(destination, 0o640)
    calls = []

    def copy_valid_dump(url, out_path):
        calls.append(url)
        shutil.copyfile(source, out_path)

    result = download_utils.ensure_ncbi_taxdump_file(
        destination, urlretrieve_fn=copy_valid_dump
    )

    assert result == str(destination)
    assert len(calls) == 1
    assert stat.S_IMODE(os.stat(destination).st_mode) == 0o640
    assert list(tmp_path.glob("taxdump.tar.gz.*.tmp")) == []


def test_invalid_download_is_not_installed(tmp_path):
    destination = tmp_path / "taxdump.tar.gz"

    def write_invalid_dump(_url, out_path):
        with open(out_path, "wb") as handle:
            handle.write(b"broken")

    with pytest.raises(ValueError, match="not a readable tar archive"):
        download_utils.ensure_ncbi_taxdump_file(
            destination, urlretrieve_fn=write_invalid_dump
        )

    assert not destination.exists()
    assert list(tmp_path.glob("taxdump.tar.gz.*.tmp")) == []


def test_default_download_validates_published_checksum(tmp_path, monkeypatch):
    source = write_tiny_taxdump(tmp_path / "source.tar.gz")
    destination = tmp_path / "taxdump.tar.gz"
    expected_md5 = hashlib.md5(  # noqa: S324 - mirrors the upstream integrity checksum
        source.read_bytes(),
        usedforsecurity=False,
    ).hexdigest()
    requested = []

    def copy_valid_dump(url, out_path):
        requested.append(url)
        shutil.copyfile(source, out_path)

    def serve_checksum(url, timeout):
        requested.append((url, timeout))
        return io.BytesIO("{}  taxdump.tar.gz\n".format(expected_md5).encode("ascii"))

    monkeypatch.setattr(download_utils.urllib.request, "urlretrieve", copy_valid_dump)
    monkeypatch.setattr(download_utils.urllib.request, "urlopen", serve_checksum)

    result = download_utils.ensure_ncbi_taxdump_file(destination)

    assert result == str(destination)
    assert requested == [
        download_utils.NCBI_TAXDUMP_URL,
        (download_utils.NCBI_TAXDUMP_MD5_URL, 30),
    ]


def test_checksum_mismatch_is_not_installed(tmp_path):
    source = write_tiny_taxdump(tmp_path / "source.tar.gz")
    destination = tmp_path / "taxdump.tar.gz"

    def copy_valid_dump(_url, out_path):
        shutil.copyfile(source, out_path)

    def serve_wrong_checksum(_url, timeout):
        assert timeout == 30
        return io.BytesIO(b"00000000000000000000000000000000  taxdump.tar.gz\n")

    with pytest.raises(ValueError, match="checksum mismatch"):
        download_utils.ensure_ncbi_taxdump_file(
            destination,
            urlretrieve_fn=copy_valid_dump,
            checksum_urlopen_fn=serve_wrong_checksum,
        )

    assert not destination.exists()
    assert list(tmp_path.glob("taxdump.tar.gz.*.tmp")) == []


def test_closed_cached_taxonomy_is_reopened(tmp_path):
    source_taxdump = write_tiny_taxdump(tmp_path / "source-taxdump.tar.gz")
    args = SimpleNamespace(
        out_dir=str(tmp_path / "out"),
        download_dir=str(tmp_path / "downloads"),
    )
    cache = download_utils._get_thread_local_ncbi_taxonomy_cache()
    cache.clear()

    def copy_taxdump(_url, destination):
        shutil.copyfile(source_taxdump, destination)

    first = download_utils.get_ncbi_taxonomy(args=args, urlretrieve_fn=copy_taxdump)
    first.close()
    second = download_utils.get_ncbi_taxonomy(args=args)
    try:
        assert second is not first
        assert second.get_taxid_translator([9606]) == {9606: "Homo sapiens"}
    finally:
        second.close()
        cache.clear()


def test_wrapper_removes_crash_temp_even_when_existing_database_is_usable(tmp_path):
    source_taxdump = write_tiny_taxdump(tmp_path / "source-taxdump.tar.gz")
    args = SimpleNamespace(
        out_dir=str(tmp_path / "out"),
        download_dir=str(tmp_path / "downloads"),
    )
    data_dir = tmp_path / "downloads" / "ete_taxonomy"
    data_dir.mkdir(parents=True)
    shutil.copyfile(source_taxdump, data_dir / "taxdump.tar.gz")
    build_taxonomy_database(data_dir / "taxa.sqlite", data_dir / "taxdump.tar.gz")
    abandoned = data_dir / "taxa.sqlite.crashed.tmp"
    abandoned.write_bytes(b"partial database")
    cache = download_utils._get_thread_local_ncbi_taxonomy_cache()
    cache.clear()

    taxonomy = download_utils.get_ncbi_taxonomy(args=args)
    try:
        assert taxonomy.get_taxid_translator([9606]) == {9606: "Homo sapiens"}
        assert not abandoned.exists()
    finally:
        taxonomy.close()
        cache.clear()


def test_invalid_api_inputs_are_ignored_or_reported_consistently(tmp_path):
    taxdump = write_tiny_taxdump(tmp_path / "taxdump.tar.gz")
    dbfile = tmp_path / "taxa.sqlite"
    build_taxonomy_database(dbfile, taxdump)
    taxonomy = NcbiTaxonomy(dbfile, update=False)
    try:
        assert taxonomy.get_rank(["not-a-taxid", 2**70, None]) == {}
        assert taxonomy.get_name_translator([None, 9606, "human"]) == {"human": [9606]}
        with pytest.raises(ValueError, match="Could not find taxid"):
            taxonomy.get_lineage("not-a-taxid")
        with pytest.raises(ValueError, match="Could not find taxid"):
            taxonomy.get_lineage(2**70)
    finally:
        taxonomy.close()
