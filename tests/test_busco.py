import os
import pytest

from types import SimpleNamespace

from amalgkit.busco import (
    normalize_busco_table,
    find_full_table,
    resolve_species_fasta,
    select_tool,
)


def write_busco_table(path, rows):
    header = "# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n"
    with open(path, 'w') as f:
        f.write("# BUSCO version is: 6.0.0\n")
        f.write(header)
        for row in rows:
            f.write("\t".join(row) + "\n")


def test_normalize_busco_table(tmp_path):
    src = tmp_path / "full_table.tsv"
    rows = [
        ["BUSCO1", "Complete", "seq1", "100", "200", "url1", "desc1"],
        ["BUSCO2", "Missing", "-", "0", "0", "url2", "desc2"],
    ]
    write_busco_table(src, rows)
    dest = tmp_path / "normalized.tsv"
    normalize_busco_table(str(src), str(dest))
    content = dest.read_text()
    assert content.startswith("# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription")
    assert "BUSCO1" in content
    assert "BUSCO2" in content


def test_find_full_table_single(tmp_path):
    out_dir = tmp_path / "busco_out"
    out_dir.mkdir()
    table = out_dir / "full_table.tsv"
    table.write_text("Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
    found = find_full_table(str(out_dir))
    assert os.path.realpath(found) == os.path.realpath(str(table))


def test_find_full_table_multiple_raises(tmp_path):
    out_dir = tmp_path / "busco_out"
    out_dir.mkdir()
    (out_dir / "full_table.tsv").write_text("Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
    (out_dir / "full_table_copy.tsv").write_text("Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n")
    with pytest.raises(ValueError, match="Multiple BUSCO full_table files"):
        find_full_table(str(out_dir))


def test_resolve_species_fasta(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    fasta_path = fasta_dir / "Homo_sapiens.fa"
    fasta_path.write_text(">seq\nATGC\n")
    result = resolve_species_fasta("Homo sapiens", str(fasta_dir))
    assert os.path.realpath(result) == os.path.realpath(str(fasta_path))


def test_select_tool_auto_prefers_compleasm(monkeypatch):
    args = SimpleNamespace(tool='auto', compleasm_exe='compleasm', busco_exe='busco')

    def fake_which(name):
        if name == 'compleasm':
            return '/usr/bin/compleasm'
        return None

    monkeypatch.setattr('shutil.which', fake_which)
    assert select_tool(args) == 'compleasm'


def test_select_tool_auto_falls_back_to_busco(monkeypatch):
    args = SimpleNamespace(tool='auto', compleasm_exe='compleasm', busco_exe='busco')

    def fake_which(name):
        if name == 'busco':
            return '/usr/bin/busco'
        return None

    monkeypatch.setattr('shutil.which', fake_which)
    assert select_tool(args) == 'busco'
