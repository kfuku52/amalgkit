import datetime
import errno
import hashlib
import os
import sqlite3
import stat
import tarfile
import tempfile
import time
import warnings
from pathlib import Path


NATIVE_SCHEMA_VERSION = 1
NATIVE_DB_FORMAT = "amalgkit"
LEGACY_ETE4_DB_FORMAT = "ete4"
SQLITE_QUERY_CHUNK_SIZE = 500
MAX_LINEAGE_DEPTH = 1000
STALE_BUILD_TEMP_SECONDS = 3600

REQUIRED_TAXDUMP_MEMBERS = {"nodes.dmp", "names.dmp", "merged.dmp"}
NATIVE_REQUIRED_COLUMNS = {
    "metadata": {"key", "value"},
    "nodes": {"taxid", "parent_taxid", "rank"},
    "names": {"taxid", "name", "is_scientific"},
    "merged": {"old_taxid", "new_taxid"},
}
LEGACY_ETE4_REQUIRED_COLUMNS = {
    "stats": {"version"},
    "species": {"taxid", "parent", "spname", "rank", "track"},
    "synonym": {"taxid", "spname"},
    "merged": {"taxid_old", "taxid_new"},
}

SCIENTIFIC_NAME_CLASS = "scientific name"
SUPPORTED_SYNONYM_CLASSES = {
    "synonym",
    "equivalent name",
    "genbank equivalent name",
    "anamorph",
    "genbank synonym",
    "genbank anamorph",
    "teleomorph",
}


def _sqlite_read_only_uri(path):
    return Path(os.path.realpath(path)).as_uri() + "?mode=ro"


def _open_read_only_database(path):
    return sqlite3.connect(_sqlite_read_only_uri(path), uri=True)


def _iter_chunks(values, chunk_size=SQLITE_QUERY_CHUNK_SIZE):
    for start in range(0, len(values), chunk_size):
        yield values[start : start + chunk_size]


def _sql_placeholders(values):
    return ",".join(["?"] * len(values))


def _normalize_taxids(taxids):
    normalized = set()
    for taxid in taxids:
        if taxid is None or taxid == "":
            continue
        try:
            taxid = int(taxid)
        except (TypeError, ValueError, OverflowError):
            continue
        if -(2**63) <= taxid <= (2**63 - 1):
            normalized.add(taxid)
    return sorted(normalized)


def _split_taxdump_line(raw_line):
    return [field.strip() for field in raw_line.decode("utf-8").split("|")]


def _iter_parsed_batches(member, parse_line, batch_size=50000):
    batch = []
    for raw_line in member:
        row = parse_line(raw_line)
        if row is None:
            continue
        batch.append(row)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def _extract_required_member(tar, member_name):
    try:
        member = tar.extractfile(member_name)
    except KeyError as exc:
        raise ValueError(
            "NCBI taxonomy dump does not contain {}".format(member_name)
        ) from exc
    if member is None:
        raise ValueError("NCBI taxonomy dump does not contain {}".format(member_name))
    return member


def _regular_file_stat(path):
    raw_path = os.path.abspath(os.fspath(path))
    if os.path.islink(raw_path):
        raise OSError("Path must not be a symlink: {}".format(raw_path))
    result = os.stat(raw_path, follow_symlinks=False)
    if not stat.S_ISREG(result.st_mode):
        raise OSError("Path is not a regular file: {}".format(raw_path))
    return result


def _file_identity(stat_result):
    return (
        stat_result.st_dev,
        stat_result.st_ino,
        stat_result.st_size,
        stat_result.st_mtime_ns,
        stat_result.st_ctime_ns,
    )


def _sha256_file(path):
    raw_path = os.path.abspath(os.fspath(path))
    digest = hashlib.sha256()
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(raw_path, flags)
    with os.fdopen(fd, "rb") as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise OSError("Path is not a regular file: {}".format(raw_path))
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
        after = os.fstat(handle.fileno())
    if _file_identity(before) != _file_identity(after):
        raise RuntimeError("File changed while it was being read: {}".format(raw_path))
    return digest.hexdigest(), after


def validate_taxonomy_dump(taxdump_file):
    """Validate that a regular tar archive contains the files needed to build."""
    raw_path = os.path.abspath(os.fspath(taxdump_file))
    _regular_file_stat(raw_path)
    try:
        with tarfile.open(raw_path, "r:*") as tar:
            members = tar.getmembers()
    except (OSError, EOFError, tarfile.TarError) as exc:
        raise ValueError(
            "NCBI taxonomy dump is not a readable tar archive: {}".format(raw_path)
        ) from exc
    counts = {}
    for member in members:
        if member.name in REQUIRED_TAXDUMP_MEMBERS and member.isfile():
            counts[member.name] = counts.get(member.name, 0) + 1
    invalid = sorted(
        name for name in REQUIRED_TAXDUMP_MEMBERS if counts.get(name, 0) != 1
    )
    if invalid:
        raise ValueError(
            "NCBI taxonomy dump must contain one regular copy of each required "
            "member; invalid: {}".format(", ".join(invalid))
        )
    return raw_path


def is_taxonomy_dump_usable(taxdump_file):
    try:
        validate_taxonomy_dump(taxdump_file)
    except (OSError, ValueError):
        return False
    return True


def _database_has_required_columns(db, required_columns):
    for table, required in required_columns.items():
        actual = {
            row[1]
            for row in db.execute("PRAGMA table_info({})".format(table)).fetchall()
        }
        if not required.issubset(actual):
            return False
    return True


def detect_taxonomy_database_format(dbfile):
    """Return the supported taxonomy DB format, or ``None`` when unusable."""
    dbfile = os.path.abspath(os.fspath(dbfile))
    if os.path.islink(dbfile) or not os.path.isfile(dbfile):
        return None
    try:
        db = _open_read_only_database(dbfile)
    except (OSError, sqlite3.Error):
        return None
    try:
        tables = {
            row[0]
            for row in db.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()
        }
        if set(NATIVE_REQUIRED_COLUMNS).issubset(
            tables
        ) and _database_has_required_columns(db, NATIVE_REQUIRED_COLUMNS):
            row = db.execute(
                "SELECT value FROM metadata WHERE key='schema_version'"
            ).fetchone()
            root = db.execute("SELECT parent_taxid FROM nodes WHERE taxid=1").fetchone()
            root_name = db.execute(
                "SELECT name FROM names WHERE taxid=1 AND is_scientific=1"
            ).fetchone()
            if (
                row is not None
                and int(row[0]) == NATIVE_SCHEMA_VERSION
                and root is not None
                and int(root[0]) == 1
                and root_name is not None
            ):
                return NATIVE_DB_FORMAT
        if set(LEGACY_ETE4_REQUIRED_COLUMNS).issubset(
            tables
        ) and _database_has_required_columns(db, LEGACY_ETE4_REQUIRED_COLUMNS):
            row = db.execute("SELECT version FROM stats").fetchone()
            root = db.execute("SELECT spname FROM species WHERE taxid=1").fetchone()
            if row is not None and int(row[0]) == 2 and root is not None:
                return LEGACY_ETE4_DB_FORMAT
    except (ValueError, TypeError, sqlite3.Error):
        return None
    finally:
        db.close()
    return None


def is_taxonomy_database_compatible(dbfile):
    return detect_taxonomy_database_format(dbfile) is not None


def _read_native_database_metadata(dbfile):
    if detect_taxonomy_database_format(dbfile) != NATIVE_DB_FORMAT:
        return None
    db = _open_read_only_database(dbfile)
    try:
        return dict(db.execute("SELECT key, value FROM metadata").fetchall())
    except sqlite3.Error:
        return None
    finally:
        db.close()


def native_database_matches_taxdump(dbfile, taxdump_file):
    """Return whether a native DB was built from the current local taxdump."""
    metadata = _read_native_database_metadata(dbfile)
    if metadata is None:
        return None
    try:
        taxdump_stat = _regular_file_stat(taxdump_file)
    except OSError:
        return False
    stored_size = metadata.get("taxdump_size")
    stored_mtime_ns = metadata.get("taxdump_mtime_ns")
    stored_ctime_ns = metadata.get("taxdump_ctime_ns")
    if (
        stored_size is not None
        and stored_mtime_ns is not None
        and stored_ctime_ns is not None
    ):
        try:
            if (
                int(stored_size) == taxdump_stat.st_size
                and int(stored_mtime_ns) == taxdump_stat.st_mtime_ns
                and int(stored_ctime_ns) == taxdump_stat.st_ctime_ns
            ):
                return True
        except (TypeError, ValueError):
            pass
    stored_hash = metadata.get("taxdump_sha256")
    if not stored_hash:
        return False
    try:
        current_hash, _ = _sha256_file(taxdump_file)
    except (OSError, RuntimeError):
        return False
    return current_hash == stored_hash


def _create_native_schema(db):
    db.executescript(
        """
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;
        PRAGMA temp_store=MEMORY;
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        ) WITHOUT ROWID;
        CREATE TABLE nodes (
            taxid INTEGER PRIMARY KEY,
            parent_taxid INTEGER NOT NULL,
            rank TEXT NOT NULL
        );
        CREATE TABLE names (
            taxid INTEGER NOT NULL,
            name TEXT COLLATE NOCASE NOT NULL,
            is_scientific INTEGER NOT NULL CHECK (is_scientific IN (0, 1)),
            PRIMARY KEY (name, taxid, is_scientific)
        ) WITHOUT ROWID;
        CREATE TABLE merged (
            old_taxid INTEGER PRIMARY KEY,
            new_taxid INTEGER NOT NULL
        );
        """
    )


def _load_nodes(db, tar):
    def parse_node(raw_line):
        fields = _split_taxdump_line(raw_line)
        return int(fields[0]), int(fields[1]), fields[2]

    member = _extract_required_member(tar, "nodes.dmp")
    for batch in _iter_parsed_batches(member, parse_node):
        db.executemany(
            "INSERT INTO nodes (taxid, parent_taxid, rank) VALUES (?, ?, ?)",
            batch,
        )


def _load_names(db, tar):
    accepted_classes = SUPPORTED_SYNONYM_CLASSES | {SCIENTIFIC_NAME_CLASS}

    def parse_name(raw_line):
        fields = _split_taxdump_line(raw_line)
        name_class = fields[3].lower()
        if name_class not in accepted_classes:
            return None
        name = fields[1].rstrip('"').lstrip('"')
        return int(fields[0]), name, int(name_class == SCIENTIFIC_NAME_CLASS)

    member = _extract_required_member(tar, "names.dmp")
    for batch in _iter_parsed_batches(member, parse_name):
        db.executemany(
            "INSERT OR IGNORE INTO names (taxid, name, is_scientific) VALUES (?, ?, ?)",
            batch,
        )


def _load_merged_taxids(db, tar):
    def parse_merged(raw_line):
        fields = _split_taxdump_line(raw_line)
        return int(fields[0]), int(fields[1])

    member = _extract_required_member(tar, "merged.dmp")
    for batch in _iter_parsed_batches(member, parse_merged):
        db.executemany(
            "INSERT INTO merged (old_taxid, new_taxid) VALUES (?, ?)",
            batch,
        )


def _validate_native_database(dbfile):
    if detect_taxonomy_database_format(dbfile) != NATIVE_DB_FORMAT:
        raise ValueError("Generated NCBI taxonomy database has an unsupported schema")
    db = _open_read_only_database(dbfile)
    try:
        quick_check = db.execute("PRAGMA quick_check").fetchone()
        if quick_check is None or quick_check[0] != "ok":
            raise ValueError(
                "Generated NCBI taxonomy database failed SQLite quick_check"
            )
        root = db.execute("SELECT parent_taxid FROM nodes WHERE taxid=1").fetchone()
        if root is None or int(root[0]) != 1:
            raise ValueError(
                "Generated NCBI taxonomy database does not contain the root taxid"
            )
        root_name = db.execute(
            "SELECT name FROM names WHERE taxid=1 AND is_scientific=1"
        ).fetchone()
        if root_name is None:
            raise ValueError(
                "Generated NCBI taxonomy database does not contain the root name"
            )
    finally:
        db.close()


def cleanup_stale_taxonomy_build_files(dbfile, stale_seconds=STALE_BUILD_TEMP_SECONDS):
    """Remove abandoned builder temp files without touching active recent builds."""
    dbfile = os.path.abspath(os.fspath(dbfile))
    db_dir = os.path.dirname(dbfile) or "."
    if not os.path.isdir(db_dir):
        return []
    prefix = os.path.basename(dbfile) + "."
    now = time.time()
    removed = []
    with os.scandir(db_dir) as entries:
        for entry in entries:
            if not entry.name.startswith(prefix) or not entry.name.endswith(".tmp"):
                continue
            try:
                entry_stat = entry.stat(follow_symlinks=False)
            except FileNotFoundError:
                continue
            if stale_seconds > 0 and now - entry_stat.st_mtime < stale_seconds:
                continue
            if not (
                stat.S_ISREG(entry_stat.st_mode) or stat.S_ISLNK(entry_stat.st_mode)
            ):
                continue
            try:
                os.unlink(entry.path)
            except FileNotFoundError:
                continue
            removed.append(entry.path)
    return removed


def _database_output_mode(dbfile, taxdump_file):
    source_path = dbfile if os.path.exists(dbfile) else taxdump_file
    source_mode = stat.S_IMODE(_regular_file_stat(source_path).st_mode)
    return (source_mode & 0o666) | stat.S_IRUSR | stat.S_IWUSR


def _prepare_regular_file_for_replace(path, mode=None):
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(path, flags)
    try:
        opened_stat = os.fstat(fd)
        if not stat.S_ISREG(opened_stat.st_mode):
            raise OSError("Path is not a regular file: {}".format(path))
        if mode is not None:
            if hasattr(os, "fchmod"):
                os.fchmod(fd, mode)
            else:
                os.chmod(path, mode)
        os.fsync(fd)
        path_stat = os.stat(path, follow_symlinks=False)
        if (
            not stat.S_ISREG(path_stat.st_mode)
            or path_stat.st_dev != opened_stat.st_dev
            or path_stat.st_ino != opened_stat.st_ino
        ):
            raise RuntimeError(
                "File changed before atomic replacement: {}".format(path)
            )
    finally:
        os.close(fd)


def _fsync_parent_directory(path):
    directory = os.path.dirname(os.path.abspath(path)) or "."
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    try:
        fd = os.open(directory, flags)
    except OSError as exc:
        unsupported = {
            errno.EACCES,
            errno.EBADF,
            errno.EINVAL,
            getattr(errno, "ENOTSUP", errno.EINVAL),
        }
        if exc.errno in unsupported:
            return
        raise
    try:
        os.fsync(fd)
    except OSError as exc:
        unsupported = {
            errno.EBADF,
            errno.EINVAL,
            getattr(errno, "ENOTSUP", errno.EINVAL),
        }
        if exc.errno not in unsupported:
            raise
    finally:
        os.close(fd)


def build_taxonomy_database(dbfile, taxdump_file):
    """Build an AMALGKIT taxonomy SQLite DB and replace ``dbfile`` atomically."""
    raw_dbfile = os.path.abspath(os.fspath(dbfile))
    raw_taxdump_file = os.path.abspath(os.fspath(taxdump_file))
    if os.path.islink(raw_dbfile):
        raise IsADirectoryError(
            "NCBI taxonomy DB path must not be a symlink: {}".format(raw_dbfile)
        )
    if os.path.islink(raw_taxdump_file):
        raise FileNotFoundError(
            "NCBI taxonomy dump must not be a symlink: {}".format(raw_taxdump_file)
        )
    dbfile = os.path.realpath(dbfile)
    taxdump_file = os.path.realpath(taxdump_file)
    if not os.path.isfile(taxdump_file):
        raise FileNotFoundError(
            "NCBI taxonomy dump is not a regular file: {}".format(taxdump_file)
        )
    if os.path.lexists(dbfile) and not os.path.isfile(dbfile):
        raise IsADirectoryError(
            "NCBI taxonomy DB path exists but is not a file: {}".format(dbfile)
        )
    db_dir = os.path.dirname(dbfile)
    if db_dir:
        if os.path.exists(db_dir) and not os.path.isdir(db_dir):
            raise NotADirectoryError(
                "NCBI taxonomy DB parent is not a directory: {}".format(db_dir)
            )
        os.makedirs(db_dir, exist_ok=True)

    cleanup_stale_taxonomy_build_files(dbfile)
    source_stat = _regular_file_stat(taxdump_file)
    output_mode = _database_output_mode(dbfile, taxdump_file)

    temp_fd, temp_path = tempfile.mkstemp(
        prefix=os.path.basename(dbfile) + ".",
        suffix=".tmp",
        dir=db_dir or ".",
    )
    os.close(temp_fd)
    db = None
    try:
        db = sqlite3.connect(temp_path)
        _create_native_schema(db)
        with tarfile.open(taxdump_file, "r:*") as tar:
            _load_nodes(db, tar)
            _load_names(db, tar)
            _load_merged_taxids(db, tar)
        taxdump_hash, final_source_stat = _sha256_file(taxdump_file)
        if _file_identity(source_stat) != _file_identity(final_source_stat):
            raise RuntimeError(
                "NCBI taxonomy dump changed while the database was being built"
            )
        metadata_rows = [
            ("schema_version", str(NATIVE_SCHEMA_VERSION)),
            ("taxdump_sha256", taxdump_hash),
            ("taxdump_size", str(final_source_stat.st_size)),
            ("taxdump_mtime_ns", str(final_source_stat.st_mtime_ns)),
            ("taxdump_ctime_ns", str(final_source_stat.st_ctime_ns)),
            (
                "built_at_utc",
                datetime.datetime.now(datetime.timezone.utc).isoformat(),
            ),
        ]
        db.executemany("INSERT INTO metadata (key, value) VALUES (?, ?)", metadata_rows)
        db.executescript(
            """
            CREATE INDEX names_by_taxid ON names (taxid, is_scientific);
            ANALYZE;
            """
        )
        db.commit()
        db.close()
        db = None
        _validate_native_database(temp_path)
        _prepare_regular_file_for_replace(temp_path, mode=output_mode)
        os.replace(temp_path, dbfile)
        _fsync_parent_directory(dbfile)
    except BaseException:
        if db is not None:
            db.close()
        if os.path.exists(temp_path):
            os.remove(temp_path)
        raise
    return dbfile


class NcbiTaxonomy:
    """Small local connector for the NCBI taxonomy data used by AMALGKIT."""

    def __init__(self, dbfile, taxdump_file=None, memory=False, update=True):
        raw_dbfile = os.path.abspath(os.fspath(dbfile))
        if os.path.islink(raw_dbfile):
            raise ValueError(
                "NCBI taxonomy DB path must not be a symlink: {}".format(raw_dbfile)
            )
        self.dbfile = os.path.realpath(raw_dbfile)
        if taxdump_file is not None:
            build_taxonomy_database(self.dbfile, taxdump_file)
        db_format = detect_taxonomy_database_format(self.dbfile)
        if db_format is None:
            update_hint = ""
            if update and taxdump_file is None:
                update_hint = " A valid taxdump_file is required to build it."
            raise ValueError(
                "Cannot open a supported NCBI taxonomy database: {}.{}".format(
                    self.dbfile,
                    update_hint,
                )
            )
        self.db_format = db_format
        file_db = _open_read_only_database(self.dbfile)
        if memory:
            self.db = sqlite3.connect(":memory:")
            file_db.backup(self.db)
            file_db.close()
        else:
            self.db = file_db

    def close(self):
        db = getattr(self, "db", None)
        if db is not None:
            db.close()
            self.db = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def _get_merged_mapping(self, taxids):
        pending = _normalize_taxids(taxids)
        mapping = {}
        if not pending:
            return mapping
        if self.db_format == NATIVE_DB_FORMAT:
            table = "merged"
            old_column = "old_taxid"
            new_column = "new_taxid"
        else:
            table = "merged"
            old_column = "taxid_old"
            new_column = "taxid_new"
        for chunk in _iter_chunks(pending):
            query = "SELECT {old}, {new} FROM {table} WHERE {old} IN ({values})".format(
                old=old_column,
                new=new_column,
                table=table,
                values=_sql_placeholders(chunk),
            )
            for old_taxid, new_taxid in self.db.execute(query, chunk).fetchall():
                mapping[int(old_taxid)] = int(new_taxid)
        return mapping

    def _resolve_merged_taxids(self, taxids):
        resolved = {taxid: taxid for taxid in _normalize_taxids(taxids)}
        for _ in range(20):
            current_to_original = {}
            for original, current in resolved.items():
                current_to_original.setdefault(current, []).append(original)
            merged = self._get_merged_mapping(list(current_to_original))
            if not merged:
                break
            changed = False
            for current, new_taxid in merged.items():
                for original in current_to_original[current]:
                    if resolved[original] != new_taxid:
                        resolved[original] = new_taxid
                        changed = True
            if not changed:
                break
        return resolved

    def get_rank(self, taxids):
        taxids = _normalize_taxids(taxids)
        ranks = {}
        if not taxids:
            return ranks
        table = "nodes" if self.db_format == NATIVE_DB_FORMAT else "species"
        for chunk in _iter_chunks(taxids):
            query = "SELECT taxid, rank FROM {} WHERE taxid IN ({})".format(
                table,
                _sql_placeholders(chunk),
            )
            for taxid, rank in self.db.execute(query, chunk).fetchall():
                ranks[int(taxid)] = rank
        return ranks

    def get_lineage_translator(self, taxids):
        taxids = _normalize_taxids(taxids)
        lineages = {}
        if not taxids:
            return lineages
        if self.db_format == LEGACY_ETE4_DB_FORMAT:
            for chunk in _iter_chunks(taxids):
                query = "SELECT taxid, track FROM species WHERE taxid IN ({})".format(
                    _sql_placeholders(chunk)
                )
                for taxid, track in self.db.execute(query, chunk).fetchall():
                    lineages[int(taxid)] = [
                        int(value) for value in reversed(track.split(","))
                    ]
            return lineages

        for chunk in _iter_chunks(taxids):
            query = """
                WITH RECURSIVE lineage(query_taxid, taxid, parent_taxid, depth) AS (
                    SELECT taxid, taxid, parent_taxid, 0
                    FROM nodes
                    WHERE taxid IN ({values})
                    UNION ALL
                    SELECT lineage.query_taxid, nodes.taxid, nodes.parent_taxid, lineage.depth + 1
                    FROM lineage
                    JOIN nodes ON nodes.taxid = lineage.parent_taxid
                    WHERE lineage.taxid != lineage.parent_taxid
                      AND lineage.depth < {max_depth}
                )
                SELECT query_taxid, taxid, depth
                FROM lineage
                ORDER BY query_taxid, depth DESC
            """.format(
                values=_sql_placeholders(chunk),
                max_depth=MAX_LINEAGE_DEPTH,
            )
            for query_taxid, lineage_taxid, _depth in self.db.execute(
                query, chunk
            ).fetchall():
                lineages.setdefault(int(query_taxid), []).append(int(lineage_taxid))
        return lineages

    def get_lineage(self, taxid):
        if not taxid:
            return None
        try:
            taxid = int(taxid)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("Could not find taxid: {}".format(taxid)) from exc
        if not -(2**63) <= taxid <= (2**63 - 1):
            raise ValueError("Could not find taxid: {}".format(taxid))
        lineage = self.get_lineage_translator([taxid]).get(taxid)
        if lineage is not None:
            return lineage
        resolved_taxid = self._resolve_merged_taxids([taxid]).get(taxid, taxid)
        if resolved_taxid != taxid:
            lineage = self.get_lineage_translator([resolved_taxid]).get(resolved_taxid)
            if lineage is not None:
                warnings.warn(
                    "taxid {} was translated into {}".format(taxid, resolved_taxid),
                    UserWarning,
                )
                return lineage
        raise ValueError("Could not find taxid: {}".format(taxid))

    def _get_scientific_names(self, taxids):
        taxids = _normalize_taxids(taxids)
        names = {}
        if not taxids:
            return names
        for chunk in _iter_chunks(taxids):
            if self.db_format == NATIVE_DB_FORMAT:
                query = """
                    SELECT taxid, name FROM names
                    WHERE is_scientific=1 AND taxid IN ({})
                """.format(_sql_placeholders(chunk))
            else:
                query = "SELECT taxid, spname FROM species WHERE taxid IN ({})".format(
                    _sql_placeholders(chunk)
                )
            for taxid, name in self.db.execute(query, chunk).fetchall():
                names[int(taxid)] = name
        return names

    def get_taxid_translator(self, taxids, try_synonyms=True):
        taxids = _normalize_taxids(taxids)
        names = self._get_scientific_names(taxids)
        if not try_synonyms:
            return names
        missing = [taxid for taxid in taxids if taxid not in names]
        resolved = self._resolve_merged_taxids(missing)
        translated_taxids = sorted(
            {current for original, current in resolved.items() if current != original}
        )
        translated_names = self._get_scientific_names(translated_taxids)
        for original, current in resolved.items():
            if current in translated_names:
                names[original] = translated_names[current]
        return names

    def _query_name_to_taxids(self, lowercase_names, scientific):
        found = {}
        if not lowercase_names:
            return found
        for chunk in _iter_chunks(sorted(lowercase_names)):
            if self.db_format == NATIVE_DB_FORMAT:
                query = """
                    SELECT name, taxid FROM names
                    WHERE is_scientific=? AND name IN ({})
                    ORDER BY taxid
                """.format(_sql_placeholders(chunk))
                params = [int(scientific)] + chunk
            else:
                table = "species" if scientific else "synonym"
                query = "SELECT spname, taxid FROM {} WHERE spname IN ({}) ORDER BY taxid".format(
                    table,
                    _sql_placeholders(chunk),
                )
                params = chunk
            for name, taxid in self.db.execute(query, params).fetchall():
                key = name.lower()
                taxid = int(taxid)
                bucket = found.setdefault(key, [])
                if taxid not in bucket:
                    bucket.append(taxid)
        return found

    def get_name_translator(self, names):
        original_by_lowercase = {}
        for name in names:
            if not isinstance(name, str):
                continue
            original_by_lowercase[name.lower()] = name
        scientific = self._query_name_to_taxids(
            set(original_by_lowercase), scientific=True
        )
        missing = set(original_by_lowercase) - set(scientific)
        synonyms = self._query_name_to_taxids(missing, scientific=False)
        scientific.update(synonyms)
        return {
            original_by_lowercase[lowercase_name]: taxids
            for lowercase_name, taxids in scientific.items()
        }
