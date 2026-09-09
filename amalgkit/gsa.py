"""Public CNCB GSA discovery and normalization.

BIG Search JSON supplies paginated experiment discovery. GSA currently exposes
run files and BioSample attributes in HTML; keep that provider-specific parser
here and fail closed when its required relationships disappear.
"""

import datetime
from html.parser import HTMLParser
import json
import re
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

import pandas

from amalgkit import __version__
from amalgkit.download_utils import maybe_acquire_download_semaphore
from contextlib import nullcontext
from amalgkit.metadata_utils import Metadata
from amalgkit.sra_sources import read_bounded_response

GSA_ROOT = "https://ngdc.cncb.ac.cn/gsa/"
SEARCH_API = "https://ngdc.cncb.ac.cn/search/api/specific"
GSA_ACCESSION = re.compile(r"(?:CRA|CRX|CRR|PRJCA)\d+\Z")
_REQUEST_LOCK = threading.Lock()


def is_gsa_accession(value):
    return GSA_ACCESSION.fullmatch(str(value).strip()) is not None


class GsaPage(HTMLParser):
    """Read links and table cells, including nested BioSample tables."""

    def __init__(self, text):
        super().__init__(convert_charrefs=True)
        self.links = []
        self.rows = []
        self._rows = []
        self._cells = []
        self.feed(text)
        self.close()

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if tag == "a" and attributes.get("href"):
            self.links.append(attributes["href"])
        if tag == "tr":
            self._rows.append([])
        if tag in {"td", "th"}:
            self._cells.append((tag, []))
        if tag == "br":
            self.handle_data(" ")

    def handle_data(self, data):
        if self._cells:
            self._cells[-1][1].append(data)

    def handle_endtag(self, tag):
        if tag in {"td", "th"} and self._cells:
            cell_tag, pieces = self._cells.pop()
            if self._rows:
                self._rows[-1].append((cell_tag, " ".join("".join(pieces).split())))
        if tag == "tr" and self._rows:
            self.rows.append(self._rows.pop())

    def fields(self):
        return {
            row[0][1].rstrip(":：").strip().lower(): row[1][1]
            for row in self.rows
            if len(row) == 2 and row[0][0] == "th" and row[1][0] == "td"
        }

    def accession_links(self, prefix):
        result = {}
        for link in self.links:
            absolute = urllib.parse.urljoin(GSA_ROOT, link)
            path = urllib.parse.urlparse(absolute).path
            match = re.fullmatch(r"/gsa/browse/(CRA\d+)/(" + prefix + r"\d+)", path)
            if match:
                result[match[2]] = absolute
        return result


def _field(fields, *names):
    return next((fields[name.lower()] for name in names if fields.get(name.lower())), "")


def _mate_marker(filename):
    # A sample/lane number must not shadow an explicit R1/R2 marker.
    for tokens in (r"read[12]|[rf][12]", r"[12]"):
        matches = list(re.finditer(r"(?:^|(?<=[._-]))(" + tokens + r")(?=[._-]|$)", filename, re.I))
        if len(matches) > 1:
            raise ValueError("Ambiguous GSA FASTQ mate markers: {}".format(filename))
        if matches:
            return matches[0]
    return None


def assign_file_mates(files, layout):
    """Require unambiguous matched lane names; never infer pairs from order."""
    if layout == "single":
        return [dict(entry, mate=0, group=i) for i, entry in enumerate(sorted(files, key=lambda x: x["filename"]))]
    if layout != "paired":
        raise ValueError("GSA FASTQ requires SINGLE or PAIRED library layout.")
    groups: dict[str, dict[int, dict[str, Any]]] = {}
    for entry in files:
        name = entry["filename"]
        marker = _mate_marker(name)
        if marker is None:
            raise ValueError("Cannot establish GSA FASTQ mate from filename: {}".format(name))
        mate = int(marker[1][-1])
        key = name[: marker.start(1)] + "{mate}" + name[marker.end(1) :]
        group = groups.setdefault(key, {})
        if mate in group:
            raise ValueError("Duplicate GSA FASTQ mate: {}".format(name))
        group[mate] = entry
    result: list[dict[str, Any]] = []
    for order, key in enumerate(sorted(groups)):
        pair = groups[key]
        if set(pair) != {1, 2}:
            raise ValueError("Missing GSA FASTQ mate for file group: {}".format(key))
        result.extend(dict(pair[mate], mate=mate, group=order) for mate in (1, 2))
    return result


def parse_run_files(page, run_id, layout):
    files: dict[str, dict[str, Any]] = {}
    for link in page.links:
        parsed = urllib.parse.urlparse(link)
        if parsed.scheme != "https" or parsed.hostname not in {
            "download.cncb.ac.cn",
            "download.big.ac.cn",
        }:
            continue
        name = urllib.parse.unquote(parsed.path.rsplit("/", 1)[-1])
        if not re.search(r"\.(fastq|fq)(\.(gz|bz2))?$", name, re.I):
            continue
        if "/" in name or "\\" in name or name in {".", ".."}:
            raise ValueError("Unsafe GSA FASTQ filename.")
        if run_id not in parsed.path.split("/"):
            raise ValueError("GSA FASTQ link does not belong to requested run: {}".format(run_id))
        entry = files.setdefault(name, {"filename": name, "sources": []})
        source = {"source_name": "GSA", "url": link}
        if source not in entry["sources"]:
            entry["sources"].append(source)
    if not files:
        raise ValueError(
            "No public supported FASTQ files found for {}; unavailable data or changed GSA page.".format(run_id)
        )
    # MB values on GSA pages are rounded; they cannot be used as exact byte counts.
    for entry in files.values():
        for row in page.rows:
            texts = [cell[1] for cell in row]
            if entry["filename"] in texts:
                checksums = {text.lower() for text in texts if re.fullmatch(r"[0-9a-fA-F]{32}", text)}
                if len(checksums) > 1:
                    raise ValueError("Conflicting GSA file checksums: {}".format(entry["filename"]))
                if checksums:
                    entry["expected_md5"] = checksums.pop()
        entry["sources"].sort(key=lambda source: not source["url"].startswith("https://"))
    return assign_file_mates(list(files.values()), layout)


def _validate_metadata_url(url):
    parsed = urllib.parse.urlparse(url)
    if (
        parsed.scheme != "https"
        or parsed.hostname != "ngdc.cncb.ac.cn"
        or parsed.username
        or parsed.password
        or parsed.port not in {None, 443}
    ):
        raise ValueError("Unexpected GSA metadata URL: {}".format(url))


class _GsaRedirectHandler(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        _validate_metadata_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _open_metadata(request, timeout):
    return urllib.request.build_opener(_GsaRedirectHandler()).open(request, timeout=timeout)


class GsaClient:
    def __init__(self, args=None):
        self.args = args
        self.timeout = float(getattr(args, "gsa_metadata_timeout_seconds", 30))
        if not 0 < self.timeout < float("inf"):
            raise ValueError("--gsa_metadata_timeout_seconds must be finite and > 0.")
        self.cache = {}

    def read(self, url):
        if url in self.cache:
            return self.cache[url]
        _validate_metadata_url(url)
        for attempt in range(3):
            try:
                request = urllib.request.Request(  # noqa: S310 (HTTPS and host validated above)
                    url,
                    headers={
                        "User-Agent": "amalgkit/{} (+https://github.com/kfuku52/amalgkit)".format(__version__),
                        "Accept-Language": "en-US,en;q=0.9",
                    },
                )
                # Limit metadata traffic across species workers in this process.
                slot = (
                    maybe_acquire_download_semaphore(
                        self.args,
                        "gsa_metadata_max_concurrency",
                        "gsa_metadata",
                        "GSA metadata",
                    )
                    if self.args is not None
                    else nullcontext()
                )
                with _REQUEST_LOCK, slot:
                    with _open_metadata(request, timeout=self.timeout) as response:
                        text = read_bounded_response(response, timeout=self.timeout).decode("utf-8")
                    time.sleep(0.25)
                self.cache[url] = text
                return text
            except (OSError, ValueError) as exc:
                if isinstance(exc, urllib.error.HTTPError) and exc.code not in {408, 429, 500, 502, 503, 504}:
                    raise RuntimeError("GSA metadata request failed: {}".format(url)) from exc
                if attempt == 2:
                    raise RuntimeError("GSA metadata request failed after retries: {}".format(url)) from exc
                time.sleep(2**attempt)
        raise AssertionError("Unreachable GSA retry state")

    def search(self, query):
        start = 0
        total = None
        records = []
        seen = set()
        while total is None or start < total:
            url = SEARCH_API + "?" + urllib.parse.urlencode({"db": "gsa", "q": query, "start": start, "length": 100})
            try:
                response = json.loads(self.read(url))
                if str(response["code"]) != "200":
                    raise ValueError("Non-success GSA search response")
                page = response["result"]["data"]
                if page.get("error"):
                    raise ValueError(str(page["error"]))
                count = int(page["recordsFiltered"])
                entries = page["data"]
                if not isinstance(entries, list) or count < 0:
                    raise ValueError("Invalid GSA search page")
                if total is not None and total != count:
                    raise ValueError("GSA search changed during pagination; rerun the query")
                total = count
                if not entries and start < total:
                    raise ValueError("GSA search ended before its reported total")
                for entry in entries:
                    if not isinstance(entry, dict) or entry["id"] in seen:
                        raise ValueError("Invalid or repeated GSA search record")
                    seen.add(entry["id"])
                    records.append(entry)
                start += len(entries)
                if start > total:
                    raise ValueError("GSA search returned more entries than its reported total")
            except (KeyError, TypeError, ValueError) as exc:
                raise RuntimeError("Invalid/incomplete GSA search response for {!r}: {}".format(query, exc)) from exc
        return records

    def experiment_urls(self, query):
        if re.fullmatch(r"CRR\d+", query):
            # BIG Search does not index CRR. The GSA search endpoint resolves it.
            page = GsaPage(self.read(GSA_ROOT + "search?" + urllib.parse.urlencode({"searchTerm": query})))
            if query not in page.accession_links("CRR"):
                raise ValueError("Public GSA Run accession not found: {}".format(query))
            urls = page.accession_links("CRX")
        else:
            entries = self.search(query)
            urls = {r["id"]: r["url"] for r in entries if re.fullmatch(r"CRX\d+", str(r.get("id", "")))}
            if is_gsa_accession(query):
                if query.startswith("CRX"):
                    urls = {key: value for key, value in urls.items() if key == query}
                elif query.startswith("CRA"):
                    urls = {key: value for key, value in urls.items() if "/browse/{}/".format(query) in value}
                else:
                    allowed = {r["id"] for r in entries if r.get("attrs", {}).get("BioProject") == query}
                    urls = {key: value for key, value in urls.items() if key in allowed}
        if not urls and is_gsa_accession(query):
            raise ValueError("No public GSA experiments found for {}".format(query))
        return urls

    def experiment_rows(self, experiment, url, requested_run=None):
        parsed_url = urllib.parse.urlparse(url)
        if not re.fullmatch(r"/gsa/browse/CRA\d+/" + re.escape(experiment), parsed_url.path):
            raise ValueError("Invalid GSA experiment URL: {}".format(url))
        page = GsaPage(self.read(url))
        fields = page.fields()
        if _field(fields, "实验编号", "Experiment accession", "Accession") != experiment:
            raise ValueError("GSA experiment identity missing or mismatched: {}".format(experiment))
        scientific_name = _field(fields, "物种名称", "Organism", "Organism name")
        taxids = [
            match[1]
            for link in page.links
            if "taxonomy" in link.lower() and (match := re.search(r"[?&]id=(\d+)", link))
        ]
        samples = [match[1] for link in page.links if (match := re.search(r"/(SAMC\d+)$", link))]
        projects = [match[1] for link in page.links if (match := re.search(r"/(PRJCA\d+)$", link))]
        runs = page.accession_links("CRR")
        if requested_run:
            runs = {key: value for key, value in runs.items() if key == requested_run}
        if not scientific_name or not samples or not projects or not runs:
            raise ValueError("Required GSA experiment relationships missing: {}".format(experiment))
        sample_id = samples[0]
        sample_page = GsaPage(self.read("https://ngdc.cncb.ac.cn/biosample/browse/" + sample_id))
        sample_fields = sample_page.fields()
        if _field(sample_fields, "样本编号", "Accession") != sample_id:
            raise ValueError("GSA BioSample identity missing or mismatched: {}".format(sample_id))
        sample_attrs = {
            "tissue": _field(sample_fields, "组织器官", "组织", "Tissue", "Tissue/organ", "Organ"),
            "genotype": _field(sample_fields, "基因型", "Genotype"),
            "sex": _field(sample_fields, "性别", "Sex"),
            "age": _field(sample_fields, "年龄", "Age"),
            "treatment": _field(sample_fields, "处理方法", "Treatment"),
            "sample_title": _field(
                sample_fields, "样本名称", "样品名称", "样本标题", "Sample name", "Sample title", "Title", "标题"
            ),
            "sample_description": _field(sample_fields, "描述信息", "Description"),
        }
        for run, run_url in sorted(runs.items()):
            if run_url.rsplit("/", 1)[0] != url.rsplit("/", 1)[0]:
                raise ValueError("GSA Run belongs to a different archive: {}".format(run))
            run_page = GsaPage(self.read(run_url))
            # The experiment summary in each run has a stable seven-column shape.
            summaries = [row for row in run_page.rows if len(row) == 7 and row[0][1] == experiment]
            if len(summaries) != 1:
                raise ValueError("GSA run experiment summary missing: {}".format(run))
            values = [cell[1] for cell in summaries[0]]
            run_rows = [row for row in run_page.rows if len(row) == 5 and row[0][1] == run]
            if len(run_rows) != 1:
                raise ValueError("GSA run identity missing: {}".format(run))
            run_values = [cell[1] for cell in run_rows[0]]
            layout = values[6].lower()
            file_format = run_values[2].lower()
            platform = values[2]
            exclusion = "no"
            if file_format != "fastq":
                exclusion = "unsupported_gsa_format"
            elif re.search(r"pacbio|pacific biosciences|nanopore|\bONT\b", platform, re.I):
                exclusion = "unsupported_gsa_long_read"
            elif layout not in {"single", "paired"}:
                exclusion = "unsupported_gsa_layout"
            manifest = parse_run_files(run_page, run, layout) if exclusion == "no" else []
            library_rows = [
                entry for entry in page.rows if len(entry) == 6 and entry[0][0] == "td" and entry[2][1] == values[3]
            ]
            design = library_rows[0][1][1] if len(library_rows) == 1 else ""
            yield {
                "run": run,
                "experiment": experiment,
                "biosample": sample_id,
                "bioproject": projects[0],
                "scientific_name": scientific_name,
                "taxid": taxids[0] if taxids else "",
                "sra_primary": url.rstrip("/").split("/")[-2],
                "lib_name": values[1],
                "lib_layout": layout,
                "instrument": platform,
                "platform": "ILLUMINA" if "illumina" in platform.lower() else platform,
                "lib_strategy": values[3],
                "design": design,
                "center": _field(fields, "Organization", "所属单位"),
                "lib_source": values[4],
                "lib_selection": values[5],
                "exp_title": _field(fields, "标题", "Title"),
                "published_date": run_values[4],
                "total_spots": "",
                "total_bases": "",
                "spot_length": "",
                "read_count_status": "unknown",
                "private_file": "no",
                "data_source": "gsa",
                "data_format": file_format,
                "exclusion": exclusion,
                "gsa_metadata_url": run_url,
                "gsa_fastq_files": json.dumps(manifest, sort_keys=True),
                "gsa_retrieved_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                **sample_attrs,
            }


def fetch_gsa_metadata(query, args=None, species_name=None, title_terms=None):
    query = str(query).strip()
    if not query:
        raise ValueError("A GSA accession or search query is required.")
    if "[" in query or "]" in query:
        raise ValueError("GSA search does not accept Entrez field tags; use BIG Search syntax or --species_tsv.")
    if query.upper().startswith("HRA"):
        raise ValueError("Controlled-access GSA-Human is not supported.")
    client = GsaClient(args)
    rows = []
    for experiment, url in client.experiment_urls(query).items():
        for row in client.experiment_rows(experiment, url, query if query.startswith("CRR") else None):
            if species_name and row["scientific_name"].casefold() != species_name.casefold():
                continue
            if species_name and row["lib_strategy"].lower() != "rna-seq":
                continue
            if title_terms and not any(term.casefold() in row["exp_title"].casefold() for term in title_terms):
                continue
            rows.append(row)
    metadata = Metadata.from_DataFrame(pandas.DataFrame(rows)) if rows else Metadata()
    if not metadata.df.empty and metadata.df["run"].duplicated().any():
        raise ValueError("GSA search resolved duplicate Run accessions.")
    metadata.df = metadata.df.sort_values("run", kind="stable").reset_index(drop=True)
    return metadata
