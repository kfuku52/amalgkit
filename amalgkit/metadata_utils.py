import datetime
import json
import os
import re
import sys
import threading
import warnings
import xml.etree.ElementTree as ET

import numpy
import pandas

from amalgkit.download_utils import get_ete_ncbitaxa
from amalgkit.exceptions import AmalgkitExit
from amalgkit.output_utils import atomic_write_dataframe
from amalgkit.parallel_utils import (
    is_auto_parallel_option,
    run_tasks_with_optional_threads,
    validate_positive_int_option,
)

PRIVATE_FASTQ_SCIENTIFIC_NAME_PLACEHOLDER = 'Please add in format: Genus species'
_TAXONOMY_LOOKUP_CACHE_LOCK = threading.Lock()
_LINEAGE_BY_TAXID_CACHE = {}
_RANK_BY_TAXID_CACHE = {}
_SCIENTIFIC_NAME_BY_TAXID_CACHE = {}


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


def normalize_scientific_name_text(value):
    if pandas.isna(value):
        return ''
    text = str(value).strip()
    if text.lower() in ['', 'nan', 'none']:
        return ''
    return re.sub(r'\s+', ' ', text)


def is_private_fastq_scientific_name_placeholder(value):
    normalized = normalize_scientific_name_text(value)
    if normalized == '':
        return False
    return normalized.lower() == PRIVATE_FASTQ_SCIENTIFIC_NAME_PLACEHOLDER.lower()


def _taxonomy_cache_namespace(ncbi):
    return getattr(ncbi, '_amalgkit_cache_key', ('object', id(ncbi)))


def _get_cached_taxonomy_entries(cache_store, namespace, taxids):
    with _TAXONOMY_LOOKUP_CACHE_LOCK:
        namespace_cache = cache_store.setdefault(namespace, {})
        found = {taxid: namespace_cache[taxid] for taxid in taxids if taxid in namespace_cache}
    missing = [taxid for taxid in taxids if taxid not in found]
    return found, missing


def _update_cached_taxonomy_entries(cache_store, namespace, mapping):
    if len(mapping) == 0:
        return
    with _TAXONOMY_LOOKUP_CACHE_LOCK:
        cache_store.setdefault(namespace, {}).update(mapping)


class Metadata:
    column_names = ['scientific_name', 'tissue', 'sample_group', 'genotype', 'sex', 'age',
                    'treatment', 'source_name',
                    'is_sampled', 'is_qualified', 'exclusion', 'protocol', 'bioproject', 'biosample',
                    'experiment', 'run', 'sra_primary', 'sra_sample', 'sra_study', 'study_title', 'exp_title', 'design',
                    'sample_title', 'sample_description', 'lib_name', 'lib_layout', 'lib_strategy', 'lib_source',
                    'lib_selection', 'platform', 'instrument', 'total_spots', 'total_bases', 'size', 'nominal_length',
                    'nominal_sdev',
                    'spot_length', 'read_index', 'read_class', 'read_type', 'base_coord', 'center',
                    'submitter_id',
                    'pubmed_id', 'taxid', 'published_date', 'NCBI_Link', 'AWS_Link', 'GCP_Link', ]
    removed_metadata_columns = ['lab', 'biomaterial_provider', 'cell', 'location', 'antibody', 'batch', 'misc']
    id_cols = ['bioproject', 'biosample', 'experiment', 'run', 'sra_primary', 'sra_sample', 'sra_study']

    def __init__(self, column_names=column_names):
        self.config_dir = ''
        self.df = pandas.DataFrame(index=[], columns=column_names)
        self.sample_attribute_collision_count = 0
        self.sample_attribute_collision_examples = []

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

    def _normalize_xml_root(xml_root):
        if isinstance(xml_root, ET.ElementTree):
            return xml_root.getroot()
        if isinstance(xml_root, ET.Element):
            return xml_root
        if hasattr(xml_root, 'getroot'):
            return xml_root.getroot()
        raise TypeError("Unknown input type.")

    def from_xml_roots(xml_roots):
        metadata = Metadata()

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

        def get_platform_info(entry):
            platform = ''
            instrument = ''
            platform_root = entry.find('./EXPERIMENT/PLATFORM')
            if platform_root is None:
                return platform, instrument
            for child in list(platform_root):
                tag = str(getattr(child, 'tag', '')).strip()
                if tag == '':
                    continue
                platform = tag
                instrument = get_first_text(child, './INSTRUMENT_MODEL')
                if instrument == '':
                    instrument = str(child.attrib.get('instrument_model', '')).strip()
                break
            return platform, instrument

        def normalize_sample_attribute_tag(tag):
            normalized = str(tag).strip().lower()
            normalized = re.sub(r" \(.*", "", normalized)
            normalized = re.sub(r"\s+", "_", normalized)
            normalized = re.sub(r"[^0-9a-zA-Z_]+", "_", normalized)
            normalized = re.sub(r"_+", "_", normalized)
            return normalized.strip('_')

        def append_unique_text(existing_value, new_value):
            if existing_value == "":
                return new_value
            existing_parts = [part.strip() for part in str(existing_value).split(' | ') if part.strip() != ""]
            if new_value in existing_parts:
                return existing_value
            return existing_value + ' | ' + new_value

        def record_sample_attribute_collision(normalized_tag, target_tag, existing_value, new_value):
            metadata.sample_attribute_collision_count += 1
            if len(metadata.sample_attribute_collision_examples) >= 10:
                return
            metadata.sample_attribute_collision_examples.append(
                {
                    'normalized_tag': normalized_tag,
                    'target_tag': target_tag,
                    'existing_value': existing_value,
                    'new_value': new_value,
                }
            )

        blocked_tags = set(metadata.removed_metadata_columns)
        core_column_tags = set(metadata.column_names)
        row_list = list()
        counter = 0
        for xml_root in xml_roots:
            root = Metadata._normalize_xml_root(xml_root)
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
                platform_name, instrument_model = get_platform_info(entry)
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
                    "platform": platform_name,
                    "instrument": instrument_model,
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
                sas = entry.findall('./SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
                for sa in sas:
                    tag = get_first_text(sa, './TAG')
                    if tag == "":
                        continue
                    normalized_tag = normalize_sample_attribute_tag(tag)
                    if normalized_tag == "":
                        continue
                    value = get_first_text(sa, './VALUE')
                    if value == "":
                        continue
                    if normalized_tag in blocked_tags:
                        target_tag = 'sample_attribute_' + normalized_tag
                    else:
                        target_tag = normalized_tag

                    def store_sample_attribute(target_column):
                        existing_value = str(row.get(target_column, ''))
                        if existing_value != "":
                            updated_value = append_unique_text(existing_value, value)
                            if updated_value != existing_value:
                                record_sample_attribute_collision(
                                    normalized_tag=normalized_tag,
                                    target_tag=target_column,
                                    existing_value=existing_value,
                                    new_value=value,
                                )
                            row[target_column] = updated_value
                            return
                        if target_column != normalized_tag:
                            record_sample_attribute_collision(
                                normalized_tag=normalized_tag,
                                target_tag=target_column,
                                existing_value=str(row.get(normalized_tag, '')),
                                new_value=value,
                            )
                        row[target_column] = value

                    store_sample_attribute(target_tag)
                    if normalized_tag in core_column_tags:
                        preserved_tag = 'sample_attribute_' + normalized_tag
                        if preserved_tag != target_tag:
                            store_sample_attribute(preserved_tag)
                row_list.append(row)
                counter += 1
        if len(row_list) == 0:
            return metadata
        df = pandas.DataFrame.from_records(row_list)
        now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print('{}: Finished converting {:,} samples'.format(now, counter), flush=True)
        metadata.df = df
        metadata.reorder(omit_misc=False)
        return metadata

    def from_xml(xml_root):
        return Metadata.from_xml_roots([xml_root])

    def add_standard_rank_taxids(self, args=None):
        self._require_nullable_int_taxid()
        standard_ranks = ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        ncbi = get_ete_ncbitaxa(args=args)
        cache_namespace = _taxonomy_cache_namespace(ncbi)
        lineage_columns = ['taxid_' + rank for rank in standard_ranks]
        unique_taxids = [int(taxid) for taxid in self.df['taxid'].dropna().unique().tolist()]
        if len(unique_taxids) == 0:
            for col in lineage_columns:
                self.df.loc[:, col] = pandas.Series([pandas.NA] * len(self.df), dtype='Int64')
            return

        row_map = {
            taxid: dict({'taxid': taxid}, **{col: numpy.nan for col in lineage_columns})
            for taxid in unique_taxids
        }
        lineage_failures = {}
        lineage_map, missing_taxids = _get_cached_taxonomy_entries(
            cache_store=_LINEAGE_BY_TAXID_CACHE,
            namespace=cache_namespace,
            taxids=unique_taxids,
        )

        translated_lineages = {}
        if (len(missing_taxids) > 0) and hasattr(ncbi, 'get_lineage_translator'):
            try:
                translated_lineages = {
                    int(taxid): [int(lineage_taxid) for lineage_taxid in lineage]
                    for taxid, lineage in ncbi.get_lineage_translator(missing_taxids).items()
                }
                if len(translated_lineages) > 0:
                    lineage_map.update(translated_lineages)
                    _update_cached_taxonomy_entries(
                        cache_store=_LINEAGE_BY_TAXID_CACHE,
                        namespace=cache_namespace,
                        mapping=translated_lineages,
                    )
                missing_taxids = [taxid for taxid in missing_taxids if taxid not in translated_lineages]
            except KeyboardInterrupt:
                raise
            except Exception:
                translated_lineages = {}

        for taxid in missing_taxids:
            try:
                lineage_map[taxid] = [int(lineage_taxid) for lineage_taxid in ncbi.get_lineage(taxid)]
                _update_cached_taxonomy_entries(
                    cache_store=_LINEAGE_BY_TAXID_CACHE,
                    namespace=cache_namespace,
                    mapping={taxid: lineage_map[taxid]},
                )
            except KeyboardInterrupt:
                raise
            except Exception as exc:
                lineage_failures.setdefault(taxid, exc)

        resolved_lineage_taxids = sorted(
            {
                lineage_taxid
                for lineage in lineage_map.values()
                for lineage_taxid in lineage
            }
        )
        rank_dict = None
        if len(resolved_lineage_taxids) > 0:
            rank_dict, missing_rank_taxids = _get_cached_taxonomy_entries(
                cache_store=_RANK_BY_TAXID_CACHE,
                namespace=cache_namespace,
                taxids=resolved_lineage_taxids,
            )
            if len(missing_rank_taxids) > 0:
                try:
                    fetched_rank_dict = ncbi.get_rank(missing_rank_taxids)
                    if fetched_rank_dict is None:
                        fetched_rank_dict = {}
                    fetched_rank_dict = {
                        int(lineage_taxid): rank
                        for lineage_taxid, rank in fetched_rank_dict.items()
                    }
                    _update_cached_taxonomy_entries(
                        cache_store=_RANK_BY_TAXID_CACHE,
                        namespace=cache_namespace,
                        mapping=fetched_rank_dict,
                    )
                    rank_dict.update(fetched_rank_dict)
                except KeyboardInterrupt:
                    raise
                except Exception:
                    rank_dict = None

        if rank_dict is not None:
            for taxid, lineage in lineage_map.items():
                lineage_taxid_row = row_map[taxid]
                for lineage_taxid in lineage:
                    rank = rank_dict.get(lineage_taxid)
                    if rank in standard_ranks:
                        lineage_taxid_row['taxid_' + rank] = lineage_taxid
        else:
            for taxid, lineage in lineage_map.items():
                lineage_taxid_row = row_map[taxid]
                try:
                    per_taxid_rank_dict = ncbi.get_rank(lineage)
                except KeyboardInterrupt:
                    raise
                except Exception as exc:
                    lineage_failures.setdefault(taxid, exc)
                    continue
                for lineage_taxid, rank in per_taxid_rank_dict.items():
                    if rank in standard_ranks:
                        lineage_taxid_row['taxid_' + rank] = lineage_taxid

        lineage_taxid_row_list = [row_map[taxid] for taxid in unique_taxids]
        if len(lineage_failures) > 0:
            preview = ', '.join(
                ['{} ({})'.format(taxid, exc.__class__.__name__) for taxid, exc in list(lineage_failures.items())[:5]]
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
        self.df = self.df.drop(columns=lineage_columns, errors='ignore')
        self.df = self.df.merge(lineage_taxid_df, on='taxid', how='left')

    def resolve_scientific_names(self, args=None):
        self._require_nullable_int_taxid()
        self.df['scientific_name_original'] = self.df['scientific_name']
        ncbi = get_ete_ncbitaxa(args=args)
        cache_namespace = _taxonomy_cache_namespace(ncbi)
        unique_taxids = [int(taxid) for taxid in self.df['taxid'].dropna().unique().tolist()]
        taxid2sciname, missing_taxids = _get_cached_taxonomy_entries(
            cache_store=_SCIENTIFIC_NAME_BY_TAXID_CACHE,
            namespace=cache_namespace,
            taxids=unique_taxids,
        )
        if len(missing_taxids) > 0:
            fetched = ncbi.get_taxid_translator(missing_taxids)
            fetched = {int(taxid): name for taxid, name in fetched.items()}
            _update_cached_taxonomy_entries(
                cache_store=_SCIENTIFIC_NAME_BY_TAXID_CACHE,
                namespace=cache_namespace,
                mapping=fetched,
            )
            taxid2sciname.update(fetched)
        self.df['scientific_name'] = self.df['taxid'].map(taxid2sciname).fillna(self.df['scientific_name_original'])

    def _require_nullable_int_taxid(self):
        if 'taxid' not in self.df.columns:
            raise KeyError('taxid column not found in metadata.')
        if str(self.df['taxid'].dtype) != 'Int64':
            raise TypeError('taxid column must be Int64 dtype')

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
            self.df.loc[(self.df.loc[:, 'exclusion'] == 'no') & (~is_empty), 'is_qualified'] = 'yes'
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


def load_metadata(args, dir_subcommand='metadata', batch_scope='run'):
    if args.metadata == 'inferred':
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
    normalized_scope = str(batch_scope).strip().lower()
    if normalized_scope == 'species':
        print('Entering species batch mode for per-species table generation. processing 1 species', flush=True)
        txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table.'
        species_series = metadata.df.loc[:, 'scientific_name'].fillna('').astype(str).str.strip()
        spp = species_series.loc[species_series != ''].drop_duplicates().sort_values().values
        if len(spp) == 0:
            raise ValueError('No valid scientific_name was found in metadata for species batch mode.')
        if batch > len(spp):
            raise AmalgkitExit('--batch {} is too large. Exiting.'.format(args.batch), exit_code=0)
        print(txt.format(batch, len(spp)), flush=True)
        sp = spp[batch - 1]
        print('Processing species: {}'.format(sp), flush=True)
        is_sp = (species_series == sp)
        metadata.df = metadata.df.loc[is_sp, :].reset_index(drop=True)
        return metadata
    if normalized_scope != 'run':
        raise ValueError('Unknown batch_scope "{}". Expected "run" or "species".'.format(batch_scope))
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
        raise ValueError('No sample is "sampled". Please check the "is_sampled" column in the metadata.')
    txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table. {:,} '
    txt += 'SRAs were excluded from the table (is_sampled==no).'
    print(txt.format(batch, num_sampled, len(numpy.where(~is_sampled)[0])), flush=True)
    if batch > num_sampled:
        raise AmalgkitExit('--batch {} is too large. Exiting.'.format(args.batch), exit_code=0)
    metadata.df = metadata.df.loc[is_sampled, :]
    metadata.df = metadata.df.reset_index(drop=True)
    metadata.df = metadata.df.loc[[batch - 1, ], :]
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
        sra_stat['num_read_per_sra'] = int(num_bp_per_sra / sra_stat['spot_length'])
    return sra_stat


def get_newest_intermediate_file_extension(sra_stat, work_dir, files=None):
    ext_out = 'no_extension_found'
    extensions = ['.amalgkit.fastq.gz', '.rename.fastq.gz', '.contam-filtered.fastq.gz', '.rrna-filtered.fastq.gz', '.fastp.fastq.gz', '.fastq.gz', '.fastq']
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
    extensions = ['.amalgkit.fastq.gz', '.rename.fastq.gz', '.contam-filtered.fastq.gz', '.rrna-filtered.fastq.gz', '.fastp.fastq.gz', '.fastq.gz', '.fastq']
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


def get_mapping_rate(
    metadata,
    quant_dir,
    max_workers='auto',
    run_tasks_with_optional_threads_fn=run_tasks_with_optional_threads,
    is_auto_parallel_option_fn=is_auto_parallel_option,
    validate_positive_int_option_fn=validate_positive_int_option,
):
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

        if is_auto_parallel_option_fn(max_workers):
            worker_cap = 8
        else:
            worker_cap = validate_positive_int_option_fn(max_workers, 'threads')
        max_workers = min(worker_cap, len(sra_dirs))
        results_by_sra, failures = run_tasks_with_optional_threads_fn(
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


def write_updated_metadata(
    metadata,
    outpath,
    args,
    max_workers='auto',
    get_mapping_rate_fn=get_mapping_rate,
):
    if os.path.exists(outpath):
        print('Updated metadata file was detected. Will be overwritten: {}'.format(outpath), flush=True)
    quant_dir = os.path.join(args.out_dir, 'quant')
    metadata = get_mapping_rate_fn(metadata, quant_dir, max_workers=max_workers)
    print('Writing per-species metadata containing mapping rate: {}'.format(outpath))
    atomic_write_dataframe(metadata.df, outpath, sep='\t', index=False)
