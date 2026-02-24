import json
import numpy
import pandas
import xml.etree.ElementTree as ET
import datetime
import ete4

import inspect
import os
import re
import shutil
import subprocess
import sys
import warnings
from bisect import bisect_left
from concurrent.futures import ThreadPoolExecutor, as_completed

def strtobool(val):
    val = val.lower()
    if val in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif val in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"invalid truth value {val!r}")


def parse_bool_flags(values, column_name='value', default=''):
    raw_values = list(values)
    normalized = (
        pandas.Series(raw_values)
        .fillna(default)
        .astype(str)
        .str.strip()
        .str.lower()
    )
    if default != '':
        normalized = normalized.replace('', default)
    parsed = []
    invalid_values = []
    for raw_value, normalized_value in zip(raw_values, normalized.tolist()):
        try:
            parsed.append(strtobool(normalized_value))
        except ValueError:
            invalid_values.append(raw_value)
    if invalid_values:
        unique_invalid = sorted({str(v) for v in invalid_values})
        raise ValueError(
            'Column "{}" contains invalid boolean flag(s): {}'.format(
                column_name,
                ', '.join(unique_invalid),
            )
        )
    return numpy.array(parsed, dtype=bool)

def run_tasks_with_optional_threads(task_items, task_fn, max_workers=1):
    tasks = list(task_items)
    results = dict()
    failures = list()
    if len(tasks) == 0:
        return results, failures
    if max_workers is None:
        worker_limit = 1
    elif is_auto_parallel_option(max_workers):
        worker_limit = 1
    else:
        try:
            worker_limit = int(max_workers)
        except (TypeError, ValueError) as exc:
            raise ValueError('max_workers must be an integer >= 1, None, or "auto".') from exc
        if worker_limit <= 1:
            worker_limit = 1
    if (worker_limit <= 1) or (len(tasks) <= 1):
        for task in tasks:
            try:
                results[task] = task_fn(task)
            except SystemExit as exc:
                failures.append((task, RuntimeError('Task requested exit with code {}.'.format(exc.code))))
            except Exception as exc:
                failures.append((task, exc))
        return results, failures
    worker_count = min(worker_limit, len(tasks))
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {
            executor.submit(task_fn, task): task
            for task in tasks
        }
        for future in as_completed(futures):
            task = futures[future]
            try:
                results[task] = future.result()
            except SystemExit as exc:
                failures.append((task, RuntimeError('Task requested exit with code {}.'.format(exc.code))))
            except Exception as exc:
                failures.append((task, exc))
    return results, failures

def validate_positive_int_option(value, option_name):
    if option_name in ('jobs', 'species_jobs'):
        option_name = 'internal_jobs'
    if is_auto_parallel_option(value):
        raise ValueError('--{} must be > 0.'.format(option_name))
    int_value = int(value)
    if int_value <= 0:
        raise ValueError('--{} must be > 0.'.format(option_name))
    return int_value


def is_auto_parallel_option(value):
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip().lower() in ('', 'auto')
    return False


def resolve_detected_cpu_count():
    detected = os.cpu_count()
    if (detected is None) or (detected <= 0):
        return 1
    return int(detected)


def resolve_cpu_budget(internal_cpu_budget='auto'):
    if is_auto_parallel_option(internal_cpu_budget):
        return resolve_detected_cpu_count()
    internal_cpu_budget = int(internal_cpu_budget)
    if internal_cpu_budget < 0:
        raise ValueError('--internal_cpu_budget must be >= 0.')
    if internal_cpu_budget == 0:
        return resolve_detected_cpu_count()
    return internal_cpu_budget


def resolve_total_core_budget(
    requested_threads='auto',
    internal_cpu_budget='auto',
    context='',
):
    if is_auto_parallel_option(requested_threads):
        requested_total = resolve_detected_cpu_count()
    else:
        requested_total = validate_positive_int_option(requested_threads, 'threads')
    budget_cap = resolve_cpu_budget(internal_cpu_budget=internal_cpu_budget)
    budget = min(requested_total, budget_cap)
    if budget < requested_total:
        print(
            '{} reducing total cores from {} to {} to fit --internal_cpu_budget {}.'.format(
                context if context else 'CPU budget:',
                requested_total,
                budget,
                budget_cap,
            ),
            flush=True,
        )
    return budget


def resolve_thread_worker_allocation(
    requested_threads='auto',
    requested_workers='auto',
    internal_cpu_budget='auto',
    worker_option_name='internal_jobs',
    context='',
    disable_workers=False,
):
    budget = resolve_total_core_budget(
        requested_threads=requested_threads,
        internal_cpu_budget=internal_cpu_budget,
        context=context,
    )
    if disable_workers:
        if not is_auto_parallel_option(requested_workers):
            workers = validate_positive_int_option(requested_workers, worker_option_name)
            if workers != 1:
                print(
                    '{} --batch is set. Forcing --{} to 1.'.format(
                        context if context else 'Parallel:',
                        worker_option_name,
                    ),
                    flush=True,
                )
        effective_workers = 1
        effective_threads = budget
    else:
        if is_auto_parallel_option(requested_workers):
            workers = budget
        else:
            workers = validate_positive_int_option(requested_workers, worker_option_name)
        effective_workers = min(workers, budget)
        if effective_workers < workers:
            print(
                '{} reducing --{} from {} to {} to fit total core budget {}.'.format(
                    context if context else 'CPU budget:',
                    worker_option_name,
                    workers,
                    effective_workers,
                    budget,
                ),
                flush=True,
            )
        effective_threads = max(1, budget // effective_workers)
    print(
        '{} effective parallelism: {} x {} = {} core(s) max.'.format(
            context if context else 'CPU budget:',
            effective_workers,
            effective_threads,
            effective_workers * effective_threads,
        ),
        flush=True,
    )
    return effective_threads, effective_workers, budget


def resolve_worker_allocation(
    requested_workers='auto',
    requested_threads='auto',
    internal_cpu_budget='auto',
    worker_option_name='internal_jobs',
    context='',
    disable_workers=False,
):
    budget = resolve_total_core_budget(
        requested_threads=requested_threads,
        internal_cpu_budget=internal_cpu_budget,
        context=context,
    )
    if disable_workers:
        if not is_auto_parallel_option(requested_workers):
            workers = validate_positive_int_option(requested_workers, worker_option_name)
            if workers != 1:
                print(
                    '{} --batch is set. Forcing --{} to 1.'.format(
                        context if context else 'CPU budget:',
                        worker_option_name,
                    ),
                    flush=True,
                )
        effective_workers = 1
    else:
        if is_auto_parallel_option(requested_workers):
            workers = budget
        else:
            workers = validate_positive_int_option(requested_workers, worker_option_name)
        effective_workers = min(workers, budget)
        if effective_workers < workers:
            print(
                '{} reducing --{} from {} to {} to fit total core budget {}.'.format(
                    context if context else 'CPU budget:',
                    worker_option_name,
                    workers,
                    effective_workers,
                    budget,
                ),
                flush=True,
            )
    print(
        '{} effective parallel workers: {} (total core budget {}).'.format(
            context if context else 'CPU budget:',
            effective_workers,
            budget,
        ),
        flush=True,
    )
    return effective_workers, budget


def find_prefixed_entries(entries, prefix, entries_sorted=False):
    if entries_sorted:
        left = bisect_left(entries, prefix)
        right = bisect_left(entries, prefix + '\uffff')
        matched = []
        for i in range(left, right):
            entry = entries[i]
            if entry.startswith(prefix):
                matched.append(entry)
        return matched
    return sorted([
        entry for entry in entries
        if entry.startswith(prefix)
    ])


def find_species_prefixed_entries(entries, species_prefix, entries_sorted=False):
    # Match exact species prefix and only delimiter-based variants (e.g., ".idx", "_k31.idx", "-v1.fa").
    matched = find_prefixed_entries(entries, species_prefix, entries_sorted=entries_sorted)
    allowed = []
    for entry in matched:
        if entry == species_prefix:
            allowed.append(entry)
            continue
        if entry.startswith(species_prefix + '.'):
            allowed.append(entry)
            continue
        if entry.startswith(species_prefix + '_'):
            allowed.append(entry)
            continue
        if entry.startswith(species_prefix + '-'):
            allowed.append(entry)
            continue
    return allowed


def find_run_prefixed_entries(entries, run_id, entries_sorted=False):
    # Match exact run ID and only known output delimiters (e.g., ".fastq.gz", "_1.fastq.gz").
    matched = find_prefixed_entries(entries, run_id, entries_sorted=entries_sorted)
    allowed = []
    for entry in matched:
        if entry == run_id:
            allowed.append(entry)
            continue
        if entry.startswith(run_id + '.'):
            allowed.append(entry)
            continue
        if entry.startswith(run_id + '_'):
            allowed.append(entry)
            continue
    return allowed

class Metadata:
    column_names = ['scientific_name', 'tissue', 'sample_group', 'genotype', 'sex', 'age',
                    'treatment', 'source_name',
                    'is_sampled', 'is_qualified', 'exclusion', 'protocol', 'bioproject', 'biosample',
                    'experiment', 'run', 'sra_primary', 'sra_sample', 'sra_study', 'study_title', 'exp_title', 'design',
                    'sample_title', 'sample_description', 'lib_name', 'lib_layout', 'lib_strategy', 'lib_source',
                    'lib_selection', 'instrument', 'total_spots', 'total_bases', 'size', 'nominal_length',
                    'nominal_sdev',
                    'spot_length', 'read_index', 'read_class', 'read_type', 'base_coord', 'center',
                    'submitter_id',
                    'pubmed_id', 'taxid', 'published_date', 'NCBI_Link', 'AWS_Link', 'GCP_Link', ]
    removed_metadata_columns = ['lab', 'biomaterial_provider', 'cell', 'location', 'antibody', 'batch', 'misc']
    id_cols = ['bioproject', 'biosample', 'experiment', 'run', 'sra_primary', 'sra_sample', 'sra_study']

    def __init__(self, column_names=column_names):
        self.config_dir = ''
        self.df = pandas.DataFrame(index=[], columns=column_names)

    def reorder(self, omit_misc=False, column_names=column_names):
        is_empty = (self.df.shape[0] == 0)
        legacy_columns = [col for col in self.removed_metadata_columns if col in self.df.columns]
        if len(legacy_columns) > 0:
            self.df = self.df.drop(columns=legacy_columns)
        self.df.loc[:, [col for col in column_names if col not in self.df.columns]] = ''
        if omit_misc:
            self.df = self.df.loc[:, column_names]
        else:
            misc_columns = [col for col in self.df.columns if col not in column_names]
            self.df = self.df.loc[:, column_names + misc_columns]
        exclusion_series = self.df.loc[:, 'exclusion'].fillna('').astype(str).str.strip()
        exclusion_series = exclusion_series.replace('', 'no')
        self.df['exclusion'] = exclusion_series
        # reorder sample_group to the front
        if 'sample_group' in self.df.columns:
            cols = list(self.df)
            cols.insert(1, cols.pop(cols.index('sample_group')))
            self.df = self.df.loc[:, cols]
        self.df = self.df.reset_index(drop=True)
        if is_empty:
            return None

    def from_DataFrame(df):
        metadata = Metadata()
        metadata.df = df.copy(deep=True)
        metadata.reorder(omit_misc=False)
        return metadata

    def from_xml(xml_root):
        if isinstance(xml_root, ET.ElementTree):
            root = xml_root.getroot()
        elif isinstance(xml_root, ET.Element):
            root = xml_root
        elif hasattr(xml_root, 'getroot'):
            root = xml_root.getroot()
        else:
            raise TypeError("Unknown input type.")

        def get_first_text(entry, path):
            element = entry.find(path)
            if element is None:
                return ""
            text = element.text
            if text is None:
                return ""
            return str(text)

        def get_first_attr(entry, path, attr_name):
            element = entry.find(path)
            if element is None:
                return ""
            value = element.get(attr_name)
            if value is None:
                return ""
            return str(value)

        def get_external_id_map(entry):
            external_ids = {}
            for element in entry.findall(".//EXTERNAL_ID"):
                namespace = element.attrib.get('namespace')
                if (namespace is None) or (namespace in external_ids):
                    continue
                text = element.text
                external_ids[namespace] = '' if text is None else str(text)
            return external_ids

        row_list = list()
        counter = 0
        metadata = Metadata()
        for entry in root.iter(tag="EXPERIMENT_PACKAGE"):
            if counter % 1000 == 0:
                now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                print('{}: Converting {:,}th sample from XML to DataFrame'.format(now, counter), flush=True)
            external_ids = get_external_id_map(entry)
            bioproject = external_ids.get('BioProject', '')
            if bioproject == "":
                labels = entry.findall('.//LABEL')
                for label in labels:
                    text = label.text
                    if (text is not None) and text.startswith("PRJ"):
                        bioproject = text
                        break
            is_single = entry.find('.//LIBRARY_LAYOUT/SINGLE') is not None
            is_paired = entry.find('.//LIBRARY_LAYOUT/PAIRED') is not None
            if is_single:
                library_layout = "single"
            elif is_paired:
                library_layout = "paired"
            else:
                library_layout = ""
            row = {
                "bioproject": bioproject,
                "scientific_name": get_first_text(entry, './SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME'),
                "biosample": external_ids.get('BioSample', ''),
                "experiment": get_first_text(entry, './EXPERIMENT/IDENTIFIERS/PRIMARY_ID'),
                "run": get_first_text(entry, './RUN_SET/RUN/IDENTIFIERS/PRIMARY_ID'),
                "sra_primary": get_first_text(entry, './SUBMISSION/IDENTIFIERS/PRIMARY_ID'),
                "sra_sample": get_first_text(entry, './SAMPLE/IDENTIFIERS/PRIMARY_ID'),
                "sra_study": get_first_text(entry, './EXPERIMENT/STUDY_REF/IDENTIFIERS/PRIMARY_ID'),
                "published_date": get_first_attr(entry, './RUN_SET/RUN', 'published'),
                "exp_title": get_first_text(entry, './EXPERIMENT/TITLE'),
                "design": get_first_text(entry, './EXPERIMENT/DESIGN/DESIGN_DESCRIPTION'),
                "lib_name": get_first_text(entry, './EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME'),
                "lib_strategy": get_first_text(entry, './EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY'),
                "lib_source": get_first_text(entry, './EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE'),
                "lib_selection": get_first_text(entry, './EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION'),
                "lib_layout": library_layout,
                "nominal_length": get_first_attr(
                    entry,
                    "./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED",
                    "NOMINAL_LENGTH",
                ),
                "nominal_sdev": get_first_attr(
                    entry,
                    "./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED",
                    "NOMINAL_SDEV",
                ),
                "spot_length": get_first_text(entry, './EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/SPOT_LENGTH'),
                "read_index": get_first_text(entry, './EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/READ_INDEX'),
                "read_class": get_first_text(entry, './EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/READ_CLASS'),
                "read_type": get_first_text(entry, './EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/READ_TYPE'),
                "base_coord": get_first_text(entry, './EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/BASE_COORD'),
                "instrument": get_first_text(entry, './EXPERIMENT/PLATFORM/ILLUMINA/INSTRUMENT_MODEL'),
                "center": get_first_attr(entry, './SUBMISSION', 'center_name'),
                "submitter_id": get_first_text(entry, './SUBMISSION/IDENTIFIERS/SUBMITTER_ID'),
                "study_title": get_first_text(entry, './STUDY/DESCRIPTOR/STUDY_TITLE'),
                "pubmed_id": get_first_text(entry, './STUDY/STUDY_LINKS/STUDY_LINK/XREF_LINK/ID'),
                "sample_title": get_first_text(entry, './SAMPLE/TITLE'),
                "taxid": get_first_text(entry, './SAMPLE/SAMPLE_NAME/TAXON_ID'),
                "sample_description": get_first_text(entry, './SAMPLE/DESCRIPTION'),
                "total_spots": get_first_attr(entry, './RUN_SET/RUN', 'total_spots'),
                "total_bases": get_first_attr(entry, './RUN_SET/RUN', 'total_bases'),
                "size": get_first_attr(entry, './RUN_SET/RUN', 'size'),
                "NCBI_Link": get_first_attr(
                    entry,
                    "./RUN_SET/RUN/SRAFiles/SRAFile[@supertype='Primary ETL']/Alternatives[@org='NCBI']",
                    "url",
                ),
                "AWS_Link": get_first_attr(
                    entry,
                    "./RUN_SET/RUN/SRAFiles/SRAFile[@supertype='Primary ETL']/Alternatives[@org='AWS']",
                    "url",
                ),
                "GCP_Link": get_first_attr(
                    entry,
                    "./RUN_SET/RUN/SRAFiles/SRAFile[@supertype='Primary ETL']/Alternatives[@org='GCP']",
                    "url",
                ),
            }
            blocked_tags = set(metadata.removed_metadata_columns)
            sas = entry.findall('./SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
            for sa in sas:
                tag = get_first_text(sa, './TAG')
                if tag != "":
                    tag = tag.lower()
                    tag = re.sub(r" \(.*", "", tag)
                    tag = re.sub(r" ", "_", tag)
                    if (tag not in row) and (tag not in blocked_tags):
                        value = get_first_text(sa, './VALUE')
                        if value != "":
                            row[tag] = value
            row_list.append(row)
            counter += 1
        if len(row_list)==0:
            return metadata
        df = pandas.DataFrame.from_records(row_list)
        now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print('{}: Finished converting {:,} samples'.format(now, counter), flush=True)
        metadata.df = df
        metadata.reorder(omit_misc=False)
        return metadata

    def add_standard_rank_taxids(self):
        self._require_nullable_int_taxid()
        standard_ranks = ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        ncbi = ete4.NCBITaxa()
        lineage_taxid_row_list = []
        lineage_failures = []
        for taxid in self.df['taxid'].dropna().unique().tolist():
            lineage_taxid_row = {'taxid': taxid}
            for rank in standard_ranks:
                lineage_taxid_row['taxid_' + rank] = numpy.nan
            try:
                lineage = ncbi.get_lineage(taxid)
                rank_dict = ncbi.get_rank(lineage)
                for lineage_taxid, rank in rank_dict.items():
                    if rank in standard_ranks:
                        lineage_taxid_row['taxid_' + rank] = lineage_taxid
            except Exception as exc:
                lineage_failures.append((taxid, exc))
            lineage_taxid_row_list.append(lineage_taxid_row)
        if len(lineage_failures) > 0:
            preview = ', '.join(
                ['{} ({})'.format(taxid, exc.__class__.__name__) for taxid, exc in lineage_failures[:5]]
            )
            if len(lineage_failures) > 5:
                preview += ', ...'
            warnings.warn(
                'Failed to resolve NCBI lineage for {} taxid(s): {}'.format(
                    len(lineage_failures),
                    preview,
                )
            )
        lineage_taxid_df = pandas.DataFrame(lineage_taxid_row_list)
        lineage_taxid_df = lineage_taxid_df.reindex(columns=['taxid'] + ['taxid_' + rank for rank in standard_ranks], fill_value=numpy.nan)
        lineage_taxid_df = lineage_taxid_df.astype('Int64')
        self.df = self.df.merge(lineage_taxid_df, on='taxid', how='left')

    def resolve_scientific_names(self):
        self._require_nullable_int_taxid()
        self.df['scientific_name_original'] = self.df['scientific_name']
        ncbi = ete4.NCBITaxa()
        taxid2sciname = ncbi.get_taxid_translator(self.df['taxid'].dropna().unique().tolist())
        self.df['scientific_name'] = self.df['taxid'].map(taxid2sciname).fillna(self.df['scientific_name_original'])

    def _require_nullable_int_taxid(self):
        if 'taxid' not in self.df.columns:
            raise KeyError('taxid column not found in metadata.')
        if str(self.df['taxid'].dtype) != 'Int64':
            raise TypeError('taxid column must be Int64 dtype')

    def _load_tab_config(self, dir_config, config_filename):
        config_path = os.path.join(dir_config, config_filename)
        if os.path.exists(config_path) and (not os.path.isfile(config_path)):
            raise IsADirectoryError('Config path exists but is not a file: {}'.format(config_path))
        try:
            config = pandas.read_csv(
                config_path,
                parse_dates=False,
                quotechar='"',
                sep='\t',
                header=None,
                index_col=None,
                skip_blank_lines=True,
                comment='#',
            )
        except (FileNotFoundError, pandas.errors.EmptyDataError):
            config = pandas.DataFrame()
        return config.replace(numpy.nan, '')

    def _build_text_series_getter(self):
        text_series_cache = {}
        def get_text_series(col):
            if col not in text_series_cache:
                text_series_cache[col] = self.df.loc[:, col].astype(str)
            return text_series_cache[col]
        return get_text_series

    def group_attributes(self, dir_config):
        config = self._load_tab_config(dir_config=dir_config, config_filename='group_attribute.config')
        if (config.shape[0] > 0) and (config.shape[1] < 2):
            raise ValueError('group_attribute.config must contain at least 2 tab-separated columns.')
        for i in config.index:
            aggregate_to = config.iloc[i, 0]
            aggregate_from = config.iloc[i, 1]
            if (aggregate_from in self.df.columns) & (aggregate_from != ''):
                print('{}: Aggregating column "{}" to column "{}"'.format(datetime.datetime.now(), aggregate_from, aggregate_to), flush=True)
                if not aggregate_to in self.df.columns:
                    self.df[aggregate_to] = ''
                else:
                    # Prevent assigning strings into float columns (pandas FutureWarning).
                    self.df[aggregate_to] = self.df[aggregate_to].fillna('').astype(str)
                is_from_empty = (self.df.loc[:, aggregate_from].isnull()) | (
                            self.df.loc[:, aggregate_from].astype(str) == '')
                is_to_empty = (self.df.loc[:, aggregate_to].isnull()) | (self.df.loc[:, aggregate_to].astype(str) == '')
                new_annotations = self.df.loc[(~is_from_empty) & (is_to_empty), aggregate_from].astype(
                    str) + '[' + aggregate_from + ']'
                self.df.loc[(~is_from_empty) & (is_to_empty), aggregate_to] = new_annotations
                new_annotations = self.df.loc[(~is_from_empty) & (~is_to_empty), aggregate_to].astype(str) + "; " + \
                                  self.df.loc[(~is_from_empty) & (~is_to_empty), aggregate_from].astype(
                                      str) + '[' + aggregate_from + ']'
                self.df.loc[(~is_from_empty) & (~is_to_empty), aggregate_to] = new_annotations
                self.df = self.df.drop(labels=aggregate_from, axis=1)
        self.reorder(omit_misc=False)

    def mark_exclude_keywords(self, dir_config):
        config = self._load_tab_config(dir_config=dir_config, config_filename='exclude_keyword.config')
        if (config.shape[0] > 0) and (config.shape[1] < 3):
            raise ValueError('exclude_keyword.config must contain at least 3 tab-separated columns.')
        if config.shape[0]>0:
            print('{}: Marking SRAs with bad keywords'.format(datetime.datetime.now()), flush=True)
        get_text_series = self._build_text_series_getter()

        for i in config.index:
            cols = [col.strip() for col in str(config.iloc[i, 0]).split(',') if col.strip() != '']
            if len(cols) == 0:
                continue
            missing_cols = [col for col in cols if col not in self.df.columns]
            if len(missing_cols) > 0:
                raise ValueError(
                    'Column(s) in exclude_keyword.config were not found in metadata: {}'.format(
                        ', '.join(missing_cols)
                    )
                )
            reason = str(config.iloc[i, 1]).strip()
            exclude_keyword = str(config.iloc[i, 2]).strip()
            if reason == '':
                warnings.warn(
                    'Skipping exclude_keyword.config row {} because exclusion reason is empty.'.format(i + 1)
                )
                continue
            if exclude_keyword == '':
                warnings.warn(
                    'Skipping exclude_keyword.config row {} because keyword is empty.'.format(i + 1)
                )
                continue
            num_detected = 0
            for col in cols:
                text_col = get_text_series(col)
                try:
                    has_bad_keyword = text_col.str.contains(exclude_keyword, regex=True, case=False).fillna(False)
                except re.error as exc:
                    raise ValueError(
                        'Invalid regex pattern in exclude_keyword.config row {}: {}'.format(i + 1, exclude_keyword)
                    ) from exc
                self.df.loc[has_bad_keyword, 'exclusion'] = reason
                num_detected += has_bad_keyword.sum()
            txt = '{}: Marking {:,} SRAs with keyword "{}"'
            print(txt.format(datetime.datetime.now(), num_detected, exclude_keyword), flush=True)
    def mark_treatment_terms(self, dir_config):
        config = self._load_tab_config(dir_config=dir_config, config_filename='control_term.config')
        if (config.shape[0] > 0) and (config.shape[1] < 2):
            raise ValueError('control_term.config must contain at least 2 tab-separated columns.')
        if config.shape[0]>0:
            print('{}: Marking SRAs with non-control terms'.format(datetime.datetime.now()), flush=True)
        get_text_series = self._build_text_series_getter()

        for i in config.index:
            cols = [col.strip() for col in str(config.iloc[i, 0]).split(',') if col.strip() != '']
            if len(cols) == 0:
                continue
            missing_cols = [col for col in cols if col not in self.df.columns]
            if len(missing_cols) > 0:
                raise ValueError(
                    'Column(s) in control_term.config were not found in metadata: {}'.format(
                        ', '.join(missing_cols)
                    )
                )
            control_term = str(config.iloc[i, 1]).strip()
            if control_term == '':
                continue
            num_control = 0
            num_treatment = 0
            for col in cols:
                text_col = get_text_series(col)
                try:
                    is_control = text_col.str.contains(control_term, regex=True, case=False).fillna(False)
                except re.error as exc:
                    raise ValueError(
                        'Invalid regex pattern in control_term.config row {}: {}'.format(i + 1, control_term)
                    ) from exc
                if not any(is_control):
                    continue
                control_bioprojects = self.df.loc[is_control, 'bioproject'].unique()
                if len(control_bioprojects) == 0:
                    continue
                is_control_bioproject = self.df.loc[:, 'bioproject'].isin(control_bioprojects)
                is_treatment = is_control_bioproject & (~is_control)
                self.df.loc[is_treatment, 'exclusion'] = 'non_control'
                num_control += int((is_control_bioproject & is_control).sum())
                num_treatment += int(is_treatment.sum())
            txt = '{}: Applying control term "{}" to "{}": Detected control and treatment SRAs: {:,} and {:,}'
            print(txt.format(datetime.datetime.now(), control_term, ','.join(cols), num_control, num_treatment), flush=True)

    def nspot_cutoff(self, min_nspots):
        print('{}: Marking SRAs with less than {:,} reads'.format(datetime.datetime.now(), min_nspots), flush=True)
        self.df['total_spots'] = pandas.to_numeric(self.df['total_spots'], errors='coerce').fillna(0).astype(int)
        self.df.loc[-(self.df.loc[:, 'total_spots'] == 0) & (
                    self.df.loc[:, 'total_spots'] < min_nspots), 'exclusion'] = 'low_nspots'

    def mark_missing_rank(self, rank_name):
        if rank_name=='none':
            return
        print('{}: Marking SRAs with missing taxid at the {} level'.format(datetime.datetime.now(), rank_name), flush=True)
        rank_col = 'taxid_' + rank_name
        if rank_col not in self.df.columns:
            raise ValueError(
                'Column "{}" is required in metadata for missing-rank filtering.'.format(rank_col)
            )
        is_empty = (self.df[rank_col].isna())
        self.df.loc[is_empty, 'exclusion'] = 'missing_taxid'

    def mark_redundant_biosample(self, exe_flag):
        if exe_flag:
            print('{}: Marking SRAs with redundant BioSample IDs'.format(datetime.datetime.now()), flush=True)
            redundant_bool = self.df.duplicated(subset=['bioproject', 'biosample'], keep='first')
            self.df.loc[redundant_bool, 'exclusion'] = 'redundant_biosample'

    def _maximize_bioproject_sampling(self, df, target_n=10):
        if 'exclusion' not in df.columns:
            raise ValueError('Column "exclusion" is required for sample selection.')
        exclusion_series = (
            df.loc[:, 'exclusion']
            .fillna('')
            .astype(str)
            .str.strip()
            .str.lower()
        )
        is_eligible = (exclusion_series == 'no')
        if int(is_eligible.sum()) <= target_n:
            df.loc[is_eligible, 'is_sampled'] = 'yes'
            return df

        is_selected = (df.loc[:, 'is_sampled'] == 'yes') & is_eligible
        selected_n = int(is_selected.sum())
        if selected_n >= target_n:
            return df

        is_unselected_eligible = (df.loc[:, 'is_sampled'] == 'no') & is_eligible
        if not bool(is_unselected_eligible.any()):
            return df

        grouped = df.loc[is_unselected_eligible, :].groupby('bioproject', sort=False).groups
        index_pools = {}
        for bioproject, indices in grouped.items():
            shuffled = list(numpy.random.permutation(indices.to_numpy()))
            if len(shuffled) > 0:
                index_pools[bioproject] = shuffled

        active_bioprojects = list(index_pools.keys())
        while (selected_n < target_n) and (len(active_bioprojects) > 0):
            exhausted = set()
            bioproject_order = numpy.random.permutation(active_bioprojects)
            for bioproject in bioproject_order:
                if selected_n >= target_n:
                    break
                pool = index_pools[bioproject]
                if len(pool) == 0:
                    exhausted.add(bioproject)
                    continue
                selected_index = pool.pop()
                df.at[selected_index, 'is_sampled'] = 'yes'
                selected_n += 1
                if len(pool) == 0:
                    exhausted.add(bioproject)
            if len(exhausted) > 0:
                active_bioprojects = [bp for bp in active_bioprojects if bp not in exhausted]
        return df

    def label_sampled_data(self, max_sample=10):
        previous_mode = pandas.get_option('mode.chained_assignment')
        pandas.set_option('mode.chained_assignment', None)
        try:
            txt = '{}: Selecting subsets of SRA IDs for >{:,} samples per sample_group per species'
            print(txt.format(datetime.datetime.now(), max_sample), flush=True)
            self.df['sample_group'] = self.df['sample_group'].fillna('').astype(str)
            self.df['exclusion'] = self.df['exclusion'].fillna('no').astype(str).str.strip().str.lower()
            self.df['bioproject'] = self.df['bioproject'].fillna('unknown').astype(str).values
            self.df['is_sampled'] = 'no'
            self.df['is_qualified'] = 'no'
            is_empty = (self.df['sample_group'] == '')
            self.df.loc[is_empty,'exclusion'] = 'no_tissue_label'
            self.df.loc[(self.df.loc[:, 'exclusion'] == 'no'), 'is_qualified'] = 'yes'
            grouped = self.df.groupby(['scientific_name', 'sample_group'], sort=False, dropna=False)
            for (_, _), sp_sample_group in grouped:
                sampled_group = self._maximize_bioproject_sampling(
                    df=sp_sample_group.copy(),
                    target_n=max_sample,
                )
                self.df.loc[sampled_group.index, 'is_sampled'] = sampled_group['is_sampled'].to_numpy()
            self.reorder(omit_misc=False)
        finally:
            pandas.set_option('mode.chained_assignment', previous_mode)

    def remove_specialchars(self):
        for col, dtype in zip(self.df.dtypes.index, self.df.dtypes.values):
            if any([key in str(dtype) for key in ['str', 'object']]):
                self.df.loc[:, col] = self.df[col].replace(r"[\r\n'\"|]", '', regex=True)

    def pivot(self, n_sp_cutoff=0, qualified_only=True, sampled_only=False):
        df = self.df
        if qualified_only:
            df = df.loc[(df.loc[:, 'is_qualified'] == 'yes'), :]
        if sampled_only:
            df = df.loc[(df.loc[:, 'is_sampled'] == 'yes'), :]
        df_reduced = df.loc[:, ['scientific_name', 'biosample', 'sample_group']]
        pivot = df_reduced.pivot_table(columns='sample_group', index='scientific_name', aggfunc='count')
        pivot.columns = pivot.columns.get_level_values(1)
        column_sort = pivot.count(axis='index').sort_values(ascending=False).index
        index_sort = pivot.count(axis='columns').sort_values(ascending=False).index
        pivot = pivot.loc[index_sort, column_sort]
        pivot_reduced = pivot.loc[:, pivot.count(axis='index') >= n_sp_cutoff]
        column_sort = pivot_reduced.count(axis='index').sort_values(ascending=False).index
        index_sort = pivot_reduced.count(axis='columns').sort_values(ascending=False).index
        pivot_reduced = pivot_reduced.loc[index_sort, column_sort]
        return pivot_reduced

def read_config_file(file_name, dir_path):
    config_path = os.path.join(dir_path, file_name)
    if os.path.exists(config_path) and (not os.path.isfile(config_path)):
        raise IsADirectoryError('Config path exists but is not a file: {}'.format(config_path))
    try:
        df = pandas.read_csv(config_path,
                             parse_dates=False, quotechar='"', sep='\t',
                             header=None, index_col=None, skip_blank_lines=True, comment='#')
    except (FileNotFoundError, pandas.errors.EmptyDataError):
        df = pandas.DataFrame([])
    if df.shape[1]==1:
        df = df.iloc[:,0]
    return df

def load_metadata(args, dir_subcommand='metadata'):
    if args.metadata=='inferred':
        relative_path = os.path.join(args.out_dir, dir_subcommand, 'metadata.tsv')
        real_path = os.path.realpath(relative_path)
    else:
        real_path = os.path.realpath(args.metadata)
    if not os.path.exists(real_path):
        raise FileNotFoundError('Metadata file not found: {}'.format(real_path))
    if not os.path.isfile(real_path):
        raise IsADirectoryError('Metadata path exists but is not a file: {}'.format(real_path))
    print('{}: Loading metadata from: {}'.format(datetime.datetime.now(), real_path), flush=True)
    df = pandas.read_csv(real_path, sep='\t', header=0, low_memory=False)
    metadata = Metadata.from_DataFrame(df)
    if 'batch' not in dir(args):
        return metadata
    if args.batch is None:
        return metadata
    batch = int(args.batch)
    if batch <= 0:
        raise ValueError('--batch must be >= 1.')
    # --batch must be handled species-wise in curate.py
    # so we need to find out where the call came from
    frm = inspect.stack()[1]
    mod = inspect.getmodule(frm[0])
    mod_name = getattr(mod, '__name__', '')
    if mod_name == 'amalgkit.curate':
        print('Entering --batch mode for amalgkit curate. processing 1 species', flush=True)
        txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table.'
        species_series = metadata.df.loc[:, 'scientific_name'].fillna('').astype(str).str.strip()
        spp = species_series.loc[species_series != ''].drop_duplicates().sort_values().values
        if len(spp) == 0:
            raise ValueError('No valid scientific_name was found in metadata for curate --batch mode.')
        if batch > len(spp):
            sys.stderr.write('--batch {} is too large. Exiting.\n'.format(args.batch))
            sys.exit(0)
        print(txt.format(batch, len(spp)), flush=True)
        sp = spp[batch - 1]
        print('Processing species: {}'.format(sp), flush=True)
        is_sp = (species_series == sp)
        metadata.df = metadata.df.loc[is_sp,:].reset_index(drop=True)
        return metadata
    else:
        print('--batch is specified. Processing one SRA per job.', flush=True)
        if 'is_sampled' not in df.columns:
            raise ValueError('Column "is_sampled" is required when --batch is specified.')
        is_sampled = parse_bool_flags(
            df.loc[:, 'is_sampled'],
            column_name='is_sampled',
            default='no',
        )
        num_sampled = int(is_sampled.sum())
        if num_sampled == 0:
            print('No sample is "sampled". Please check the "is_sampled" column in the metadata. Exiting.')
            sys.exit(1)
        txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table. {:,} '
        txt += 'SRAs were excluded from the table (is_sampled==no).'
        print(txt.format(batch, num_sampled, len(numpy.where(is_sampled == False)[0])), flush=True)
        if batch > num_sampled:
            sys.stderr.write('--batch {} is too large. Exiting.\n'.format(args.batch))
            sys.exit(0)
        metadata.df = metadata.df.loc[is_sampled,:]
        metadata.df = metadata.df.reset_index(drop=True)
        metadata.df = metadata.df.loc[[batch - 1,],:]
        return metadata

def _build_run_index_cache(metadata):
    if 'run' not in metadata.df.columns:
        raise ValueError('Column "run" is required in metadata.')
    run_to_idx = dict()
    duplicate_runs = set()
    for idx, run_id in zip(metadata.df.index, metadata.df['run'].values):
        if run_id in run_to_idx:
            duplicate_runs.add(run_id)
        else:
            run_to_idx[run_id] = idx
    metadata._run_to_idx = run_to_idx
    metadata._run_duplicate_ids = duplicate_runs
    metadata._run_cache_df_id = id(metadata.df)
    metadata._run_cache_len = len(metadata.df)
    return run_to_idx, duplicate_runs

def get_metadata_row_index_by_run(metadata, sra_id):
    run_to_idx = getattr(metadata, '_run_to_idx', None)
    duplicate_runs = getattr(metadata, '_run_duplicate_ids', None)
    cache_df_id = getattr(metadata, '_run_cache_df_id', None)
    cache_len = getattr(metadata, '_run_cache_len', None)
    is_cache_valid = (
        isinstance(run_to_idx, dict)
        and isinstance(duplicate_runs, set)
        and (cache_df_id == id(metadata.df))
        and (cache_len == len(metadata.df))
    )
    if not is_cache_valid:
        run_to_idx, duplicate_runs = _build_run_index_cache(metadata)
    idx = run_to_idx.get(sra_id, None)
    if (idx is None) or (idx not in metadata.df.index) or (metadata.df.at[idx, 'run'] != sra_id):
        run_to_idx, duplicate_runs = _build_run_index_cache(metadata)
    if sra_id in duplicate_runs:
        raise AssertionError('There are multiple metadata rows with the same SRA ID: ' + sra_id)
    idx = run_to_idx.get(sra_id, None)
    if idx is None:
        raise AssertionError('SRA ID not found in metadata: ' + sra_id)
    return idx

def get_sra_stat(sra_id, metadata, num_bp_per_sra=None):
    required_columns = ['run', 'lib_layout', 'total_spots', 'spot_length', 'total_bases']
    missing_columns = [col for col in required_columns if col not in metadata.df.columns]
    if len(missing_columns) > 0:
        raise ValueError('Missing required metadata column(s) for get_sra_stat: {}'.format(', '.join(missing_columns)))
    sra_stat = dict()
    sra_stat['sra_id'] = sra_id
    idx = get_metadata_row_index_by_run(metadata, sra_id)
    sra_stat['metadata_idx'] = idx
    layout = str(metadata.df.at[idx, 'lib_layout']).strip().lower()
    if layout in ('nan', 'none'):
        layout = ''
    if 'layout_amalgkit' in metadata.df.columns:
        layout_amalgkit = str(metadata.df.at[idx, 'layout_amalgkit']).strip().lower()
        if layout_amalgkit in ['single', 'paired']:
            layout = layout_amalgkit
    if layout not in ['single', 'paired']:
        raise ValueError(
            'Unsupported lib_layout "{}" for SRA ID {}. Expected "single" or "paired".'.format(
                metadata.df.at[idx, 'lib_layout'],
                sra_id,
            )
        )
    sra_stat['layout'] = layout
    total_spot = pandas.to_numeric(metadata.df.at[idx, 'total_spots'], errors='coerce')
    if pandas.isna(total_spot) or (float(total_spot) <= 0):
        raise ValueError('total_spots must be > 0 for SRA ID: {}'.format(sra_id))
    sra_stat['total_spot'] = int(total_spot)
    original_spot_len = pandas.to_numeric(metadata.df.at[idx, 'spot_length'], errors='coerce')
    if pandas.isna(original_spot_len) or (float(original_spot_len) <= 0):
        total_bases = pandas.to_numeric(metadata.df.at[idx, 'total_bases'], errors='coerce')
        if pandas.isna(total_bases) or (float(total_bases) <= 0):
            raise ValueError('spot_length cannot be inferred because total_bases is missing or <= 0 for SRA ID: {}'.format(sra_id))
        inferred_spot_len = float(total_bases) / float(sra_stat['total_spot'])
        if inferred_spot_len <= 0:
            raise ValueError('spot_length inference failed for SRA ID: {}'.format(sra_id))
        sra_stat['spot_length'] = int(inferred_spot_len)
        txt = 'spot_length cannot be obtained directly from metadata. Using total_bases/total_spots instead: {:,}'
        print(txt.format(sra_stat['spot_length']))
    else:
        sra_stat['spot_length'] = int(original_spot_len)
    if 'nominal_length' in metadata.df.columns:
        nominal_length_value = metadata.df.at[idx, 'nominal_length']
    else:
        nominal_length_value = numpy.nan
    sra_stat['nominal_length'] = pandas.to_numeric(nominal_length_value, errors='coerce')
    if num_bp_per_sra is not None:
        sra_stat['num_read_per_sra'] = int(num_bp_per_sra/sra_stat['spot_length'])
    return sra_stat

def get_newest_intermediate_file_extension(sra_stat, work_dir, files=None):
    ext_out = 'no_extension_found'
    # Order is important in this list. More downstream should come first.
    extensions = ['.amalgkit.fastq.gz','.rename.fastq.gz','.fastp.fastq.gz','.fastq.gz','.fastq']
    if 'getfastq_sra_dir' not in sra_stat:
        sra_stat['getfastq_sra_dir'] = work_dir
    if files is None:
        try:
            with os.scandir(work_dir) as entries:
                files = {
                    entry.name
                    for entry in entries
                    if entry.is_file()
                }
        except FileNotFoundError:
            files = set()
    else:
        files = set(files)
    sra_stat = detect_layout_from_file(sra_stat, files=files)
    sra_id = sra_stat['sra_id']
    if sra_stat['layout'] == 'single':
        for ext in extensions:
            if (sra_id + ext) in files:
                ext_out = ext
                break
    elif sra_stat['layout'] == 'paired':
        for ext in extensions:
            pair1 = sra_id + '_1' + ext
            pair2 = sra_id + '_2' + ext
            if (pair1 in files) and (pair2 in files):
                ext_out = ext
                break
    if ext_out == 'no_extension_found':
        safe_delete_names = []
        for ext in extensions:
            safe_single = sra_id + ext + '.safely_removed'
            safe_r1 = sra_id + '_1' + ext + '.safely_removed'
            safe_r2 = sra_id + '_2' + ext + '.safely_removed'
            if safe_single in files:
                safe_delete_names.append(safe_single)
            if (safe_r1 in files) and (safe_r2 in files):
                safe_delete_names.extend([safe_r1, safe_r2])
        safe_delete_files = sorted([
            os.path.join(work_dir, f)
            for f in set(safe_delete_names)
        ])
        if len(safe_delete_files):
            txt = 'getfastq safely_removed flag was detected. `amalgkit quant` has been completed in this sample: {}\n'
            sys.stdout.write(txt.format(work_dir))
            for safe_delete_file in safe_delete_files:
                sys.stdout.write('{}\n'.format(safe_delete_file))
            return '.safely_removed'
    return ext_out

def is_there_unpaired_file(sra_stat, extensions, files=None):
    if files is None:
        try:
            with os.scandir(sra_stat['getfastq_sra_dir']) as entries:
                files = {
                    entry.name
                    for entry in entries
                    if entry.is_file()
                }
        except FileNotFoundError:
            files = set()
    is_unpaired_file = False
    for ext in extensions:
        candidates = [sra_stat['sra_id'] + ext]
        if not ext.endswith('.safely_removed'):
            candidates.append(sra_stat['sra_id'] + ext + '.safely_removed')
        if any([candidate in files for candidate in candidates]):
            is_unpaired_file = True
            break
    return is_unpaired_file

def detect_layout_from_file(sra_stat, files=None):
    # Order is important in this list. More downstream should come first.
    extensions = ['.amalgkit.fastq.gz','.rename.fastq.gz','.fastp.fastq.gz','.fastq.gz','.fastq']
    if files is None:
        try:
            with os.scandir(sra_stat['getfastq_sra_dir']) as entries:
                files = {
                    entry.name
                    for entry in entries
                    if entry.is_file()
                }
        except FileNotFoundError:
            files = set()
    is_paired_end = False
    for ext in extensions:
        pair1_candidates = [sra_stat['sra_id'] + '_1' + ext]
        pair2_candidates = [sra_stat['sra_id'] + '_2' + ext]
        if not ext.endswith('.safely_removed'):
            pair1_candidates.append(sra_stat['sra_id'] + '_1' + ext + '.safely_removed')
            pair2_candidates.append(sra_stat['sra_id'] + '_2' + ext + '.safely_removed')
        if any([f in files for f in pair1_candidates]) and any([f in files for f in pair2_candidates]):
            is_paired_end = True
            break
    is_unpaired_file = is_there_unpaired_file(sra_stat, extensions, files=files)
    if (not is_paired_end) & is_unpaired_file:
        is_single_end = True
    else:
        is_single_end = False
    if is_single_end & (sra_stat['layout'] == 'paired'):
        txt = 'Single-end fastq was generated even though layout in the metadata = {}. '
        txt += 'This sample will be treated as single-end reads: {}\n'
        txt = txt.format(sra_stat['layout'], sra_stat['sra_id'])
        sys.stderr.write(txt)
        sra_stat['layout'] = 'single'
    if is_paired_end & (sra_stat['layout'] == 'single'):
        txt = 'Paired-end fastq was generated even though layout in the metadata = {}. '
        txt += 'This sample will be treated as paired-end reads: {}\n'
        txt = txt.format(sra_stat['layout'], sra_stat['sra_id'])
        sys.stderr.write(txt)
        sra_stat['layout'] = 'paired'
    return sra_stat

def write_updated_metadata(metadata, outpath, args, max_workers='auto'):
    if os.path.exists(outpath):
        print('Updated metadata file was detected. Will be overwritten: {}'.format(outpath), flush=True)
    quant_dir = os.path.join(args.out_dir, 'quant')
    metadata = get_mapping_rate(metadata, quant_dir, max_workers=max_workers)
    print('Writing curate metadata containing mapping rate: {}'.format(outpath))
    metadata.df.to_csv(outpath, sep='\t', index=False)

def get_mapping_rate(metadata, quant_dir, max_workers='auto'):
    if 'run' not in metadata.df.columns:
        raise ValueError('Column "run" is required in metadata to compute mapping_rate.')
    if os.path.exists(quant_dir):
        if not os.path.isdir(quant_dir):
            txt = 'quant directory not found. Mapping rate cutoff will not be applied: {}\n'
            sys.stderr.write(txt.format(quant_dir))
            sys.stderr.write('Path exists but is not a directory: {}\n'.format(quant_dir))
            return metadata
        print('quant directory found: {}'.format(quant_dir))
        metadata.df.loc[:, 'mapping_rate'] = numpy.nan
        normalized_runs = metadata.df.loc[:, 'run'].fillna('').astype(str).str.strip()
        sra_dirs = []
        for sra_id in dict.fromkeys(normalized_runs.tolist()):
            if pandas.isna(sra_id):
                continue
            sra_id = str(sra_id)
            if sra_id == '':
                continue
            if os.path.isdir(os.path.join(quant_dir, sra_id)):
                sra_dirs.append(sra_id)
        print('Number of quant sub-directories that matched to metadata: {:,}'.format(len(sra_dirs)))

        def load_mapping_rate(sra_id):
            run_info_path = os.path.join(quant_dir, sra_id, sra_id + '_run_info.json')
            if not os.path.exists(run_info_path):
                return sra_id, None, 'missing'
            try:
                with open(run_info_path) as f:
                    run_info = json.load(f)
            except Exception as e:
                return sra_id, None, 'Failed to read run_info.json for {}: {}'.format(sra_id, e)
            if 'p_pseudoaligned' not in run_info:
                return sra_id, None, 'p_pseudoaligned missing in run_info.json for {}.'.format(sra_id)
            mapping_rate = pandas.to_numeric(run_info['p_pseudoaligned'], errors='coerce')
            if pandas.isna(mapping_rate) or (not numpy.isfinite(float(mapping_rate))):
                return sra_id, None, 'Invalid p_pseudoaligned value in run_info.json for {}: {}'.format(
                    sra_id,
                    run_info['p_pseudoaligned'],
                )
            return sra_id, float(mapping_rate), None

        if is_auto_parallel_option(max_workers):
            worker_cap = 8
        else:
            worker_cap = validate_positive_int_option(max_workers, 'threads')
        max_workers = min(worker_cap, len(sra_dirs))
        results_by_sra, failures = run_tasks_with_optional_threads(
            task_items=sra_dirs,
            task_fn=load_mapping_rate,
            max_workers=max_workers,
        )
        for sra_id, exc in failures:
            results_by_sra[sra_id] = (
                sra_id,
                None,
                'Failed to read run_info.json for {}: {}'.format(sra_id, exc),
            )
        results = [results_by_sra[sra_id] for sra_id in sra_dirs if sra_id in results_by_sra]

        for sra_id, mapping_rate, error in results:
            if mapping_rate is None:
                if error == 'missing':
                    sys.stderr.write('run_info.json not found. Skipping {}.\n'.format(sra_id))
                else:
                    sys.stderr.write('{}\n'.format(error))
                continue
            is_run = (normalized_runs == sra_id)
            if not bool(is_run.any()):
                sys.stderr.write('Run ID from quant output was not found in metadata. Skipping {}.\n'.format(sra_id))
                continue
            metadata.df.loc[is_run, 'mapping_rate'] = mapping_rate
    else:
        txt = 'quant directory not found. Mapping rate cutoff will not be applied: {}\n'
        sys.stderr.write(txt.format(quant_dir))
    return metadata

def check_rscript():
    try:
        probe = subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        print("R (Rscript) is not installed. Exiting.")
        sys.exit(1)
    if probe.returncode != 0:
        sys.stderr.write('Rscript dependency probe failed with exit code {}: Rscript --help\n'.format(probe.returncode))
        sys.exit(1)

def orthogroup2genecount(file_orthogroup, file_genecount, spp):
    df = pandas.read_csv(file_orthogroup, sep='\t', header=0, low_memory=False)
    if 'busco_id' not in df.columns:
        raise ValueError('Column "busco_id" is required in orthogroup table: {}'.format(file_orthogroup))
    missing_species = [sp for sp in spp if sp not in df.columns]
    if len(missing_species) > 0:
        raise ValueError(
            'Species column(s) not found in orthogroup table ({}): {}'.format(
                file_orthogroup,
                ', '.join(missing_species),
            )
        )
    orthogroup_df = pandas.DataFrame({'orthogroup_id': df['busco_id'].to_numpy()})
    is_spp = df.columns.isin(spp)
    df = df.loc[:,is_spp].fillna('').replace('-', '').astype(str)
    values = df.to_numpy(dtype=str)
    non_empty = (values != '')
    comma_counts = numpy.char.count(values, ',').astype(int)
    gc_values = non_empty.astype(int) + (comma_counts * non_empty)
    gc = pandas.DataFrame(gc_values, index=df.index, columns=df.columns)
    gc = pandas.concat([orthogroup_df, gc], axis=1)
    col_order = ['orthogroup_id'] + [col for col in gc.columns if col != 'orthogroup_id']
    gc = gc[col_order]
    gc.to_csv(file_genecount, index=False, sep='\t')

def check_ortholog_parameter_compatibility(args):
    orthogroup_table = getattr(args, 'orthogroup_table', None)
    dir_busco = getattr(args, 'dir_busco', None)
    if isinstance(orthogroup_table, str):
        orthogroup_table = orthogroup_table.strip()
        if orthogroup_table == '':
            orthogroup_table = None
    if isinstance(dir_busco, str):
        dir_busco = dir_busco.strip()
        if dir_busco == '':
            dir_busco = None
    setattr(args, 'orthogroup_table', orthogroup_table)
    setattr(args, 'dir_busco', dir_busco)
    if (orthogroup_table is None) and (dir_busco is None):
        raise ValueError('One of --orthogroup_table and --dir_busco should be specified.')
    if (orthogroup_table is not None) and (dir_busco is not None):
        raise ValueError('Only one of --orthogroup_table and --dir_busco should be specified.')

BUSCO_TABLE_COLUMNS = ['busco_id', 'status', 'sequence', 'score', 'length', 'orthodb_url', 'description']
BUSCO_TABLE_USE_COLUMNS = ['busco_id', 'sequence', 'orthodb_url', 'description']
BUSCO_SPECIES_SUFFIX_PATTERN = re.compile(r'\.tsv(?:\.gz)?$', re.IGNORECASE)
BUSCO_SPECIES_CLEANUP_PATTERN = re.compile(r'(_busco|_full_table.*)$', re.IGNORECASE)
BUSCO_SPECIES_MATCH_PATTERN = re.compile(r'^([^_]+_[^_]+)')


def parse_busco_species_name(species_infile):
    species_colname = BUSCO_SPECIES_SUFFIX_PATTERN.sub('', species_infile)
    species_colname = BUSCO_SPECIES_CLEANUP_PATTERN.sub('', species_colname)
    matched = BUSCO_SPECIES_MATCH_PATTERN.match(species_colname)
    if matched is not None:
        species_colname = matched.group(1)
    return species_colname


def read_busco_species_table(path_to_table):
    tmp_table = pandas.read_table(
        path_to_table,
        sep='\t',
        header=None,
        comment='#',
        names=BUSCO_TABLE_COLUMNS,
        usecols=BUSCO_TABLE_USE_COLUMNS,
    )
    # Some BUSCO tables include an uncommented header row. Exclude it to avoid a bogus BUSCO ID entry.
    busco_id_key = (
        tmp_table.loc[:, 'busco_id']
        .fillna('')
        .astype(str)
        .str.lower()
        .str.replace(r'[^a-z0-9]', '', regex=True)
    )
    tmp_table = tmp_table.loc[busco_id_key != 'buscoid', :].copy()
    tmp_table.loc[:, 'sequence'] = tmp_table.loc[:, 'sequence'].str.replace(r':[-\.0-9]*$', '', regex=True)
    for col in ['sequence', 'orthodb_url', 'description']:
        tmp_table[col] = tmp_table[col].fillna('').astype(str)
        tmp_table.loc[(tmp_table[col] == ''), col] = '-'
    return tmp_table


def parse_busco_species_table(dir_busco, species_infile):
    path_to_table = os.path.join(dir_busco, species_infile)
    if not os.path.isfile(path_to_table):
        if os.path.exists(path_to_table):
            warnings.warn('full_table.tsv path exists but is not a file. Skipping: {}'.format(species_infile))
        else:
            warnings.warn('full_table.tsv does not exist. Skipping: {}'.format(species_infile))
        return None
    tmp_table = read_busco_species_table(path_to_table)
    meta_rows = tmp_table.loc[:, ['busco_id', 'orthodb_url', 'description']].drop_duplicates(
        subset=['busco_id'],
        keep='first',
        inplace=False,
    )
    grouped = tmp_table.loc[:, ['busco_id', 'sequence']].groupby('busco_id', sort=False)['sequence'].agg(','.join)
    species_colname = parse_busco_species_name(species_infile)
    return species_colname, meta_rows, grouped


def append_unique_busco_ids(busco_id_seen, busco_id_order, busco_ids):
    for busco_id in busco_ids:
        if busco_id in busco_id_seen:
            continue
        busco_id_seen.add(busco_id)
        busco_id_order.append(busco_id)


def update_busco_meta(busco_meta, busco_ids, urls, descriptions):
    for busco_id, orthodb_url, description in zip(busco_ids, urls, descriptions):
        if busco_id not in busco_meta:
            busco_meta[busco_id] = {
                'orthodb_url': orthodb_url,
                'description': description,
            }
            continue
        if (busco_meta[busco_id]['orthodb_url'] == '-') and (orthodb_url != '-'):
            busco_meta[busco_id]['orthodb_url'] = orthodb_url
        if (busco_meta[busco_id]['description'] == '-') and (description != '-'):
            busco_meta[busco_id]['description'] = description


def generate_multisp_busco_table(dir_busco, outfile):
    dir_busco = os.path.realpath(dir_busco)
    if not os.path.exists(dir_busco):
        raise FileNotFoundError('BUSCO directory not found: {}'.format(dir_busco))
    if not os.path.isdir(dir_busco):
        raise NotADirectoryError('BUSCO path exists but is not a directory: {}'.format(dir_busco))
    print('Generating multi-species BUSCO table.', flush=True)
    species_infiles = [
        f for f in os.listdir(path=dir_busco)
        if BUSCO_SPECIES_SUFFIX_PATTERN.search(f)
        and os.path.isfile(os.path.join(dir_busco, f))
    ]
    species_infiles = sorted(species_infiles)
    print('BUSCO full tables for {} species were detected at: {}'.format(len(species_infiles), dir_busco), flush=True)
    if len(species_infiles) == 0:
        raise FileNotFoundError('No BUSCO full table file (.tsv) was detected in: {}'.format(dir_busco))

    busco_id_seen = set()
    busco_id_order = []
    busco_meta = dict()
    species_series = dict()
    species_order = []

    max_workers = min(8, len(species_infiles))
    parsed_by_file, parse_failures = run_tasks_with_optional_threads(
        task_items=species_infiles,
        task_fn=lambda species_infile: parse_busco_species_table(dir_busco, species_infile),
        max_workers=max_workers,
    )
    for species_infile, exc in parse_failures:
        warnings.warn('Failed to parse BUSCO table {}: {}'.format(species_infile, exc))
        parsed_by_file[species_infile] = None
    parsed_results = [parsed_by_file.get(species_infile) for species_infile in species_infiles]

    species_source_file = dict()
    for species_infile, parsed in zip(species_infiles, parsed_results):
        if parsed is None:
            continue
        species_colname, meta_rows, grouped = parsed
        existing_infile = species_source_file.get(species_colname)
        if (existing_infile is not None) and (existing_infile != species_infile):
            raise ValueError(
                'Duplicate species label was detected across BUSCO tables: {} ({} vs {}). '
                'Rename input BUSCO table files to avoid species label collision.'.format(
                    species_colname,
                    existing_infile,
                    species_infile,
                )
            )
        species_source_file[species_colname] = species_infile
        busco_ids = meta_rows['busco_id'].to_numpy()
        urls = meta_rows['orthodb_url'].to_numpy()
        descriptions = meta_rows['description'].to_numpy()
        append_unique_busco_ids(busco_id_seen, busco_id_order, busco_ids)
        update_busco_meta(busco_meta, busco_ids, urls, descriptions)
        species_order.append(species_colname)
        species_series[species_colname] = grouped

    if len(species_series) == 0:
        raise ValueError('Failed to parse any BUSCO table under: {}'.format(dir_busco))

    merged_table = pandas.DataFrame({'busco_id': busco_id_order})
    merged_table['orthodb_url'] = merged_table['busco_id'].map(
        lambda bid: busco_meta.get(bid, {}).get('orthodb_url', '-')
    )
    merged_table['description'] = merged_table['busco_id'].map(
        lambda bid: busco_meta.get(bid, {}).get('description', '-')
    )
    for species_colname in species_order:
        merged_table[species_colname] = merged_table['busco_id'].map(species_series[species_colname])
    merged_table.to_csv(outfile, sep='\t', index=None, doublequote=False)

def check_config_dir(dir_path, mode):
    mode_to_files = {
        'select': [
            'group_attribute.config',
            'exclude_keyword.config',
            'control_term.config',
        ],
    }
    asserted_files = mode_to_files.get(mode)
    if asserted_files is None:
        raise ValueError('Unsupported config check mode: {}'.format(mode))
    if not os.path.exists(dir_path):
        raise FileNotFoundError('Config directory not found: {}'.format(dir_path))
    if not os.path.isdir(dir_path):
        raise NotADirectoryError('Config path exists but is not a directory: {}'.format(dir_path))
    missing_count = 0
    for af in asserted_files:
        af_path = os.path.join(dir_path, af)
        if os.path.isfile(af_path):
            print('Config file found: {}'.format(af))
        elif os.path.exists(af_path):
            sys.stderr.write('Config entry exists but is not a file: {}\n'.format(af))
            missing_count += 1
        else:
            sys.stderr.write('Config file not found: {}\n'.format(af))
            missing_count += 1
    if (missing_count>0):
        txt = 'Please refer to the AMALGKIT Wiki for more info: https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata\n'
        sys.stderr.write(txt)

def cleanup_tmp_amalgkit_files(work_dir='.'):
    try:
        with os.scandir(work_dir) as entries:
            for entry in entries:
                if not entry.name.startswith('tmp.amalgkit.'):
                    continue
                path = entry.path
                try:
                    if entry.is_dir(follow_symlinks=False):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)
                except FileNotFoundError:
                    # Another process may remove the same temp path concurrently.
                    continue
    except FileNotFoundError:
        return

def get_getfastq_run_dir(args, sra_id):
    amalgkit_out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(amalgkit_out_dir) and (not os.path.isdir(amalgkit_out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(amalgkit_out_dir))
    getfastq_dir = os.path.join(amalgkit_out_dir, 'getfastq')
    if os.path.exists(getfastq_dir) and (not os.path.isdir(getfastq_dir)):
        raise NotADirectoryError('getfastq path exists but is not a directory: {}'.format(getfastq_dir))
    run_output_dir = os.path.join(getfastq_dir, sra_id)
    if os.path.exists(run_output_dir) and (not os.path.isdir(run_output_dir)):
        raise NotADirectoryError('getfastq run path exists but is not a directory: {}'.format(run_output_dir))
    os.makedirs(run_output_dir, exist_ok=True)
    return run_output_dir
