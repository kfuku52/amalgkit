import json
import numpy
import pandas
import xml.etree.ElementTree as ET
import datetime
import ete4

import glob
import inspect
import os
import re
import subprocess
import sys
import warnings

def strtobool(val):
    val = val.lower()
    if val in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif val in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"invalid truth value {val!r}")

class Metadata:
    column_names = ['scientific_name', 'tissue', 'sample_group', 'genotype', 'sex', 'age',
                    'treatment', 'source_name',
                    'is_sampled', 'is_qualified', 'exclusion', 'protocol', 'bioproject', 'biosample',
                    'experiment', 'run', 'sra_primary', 'sra_sample', 'sra_study', 'study_title', 'exp_title', 'design',
                    'sample_title', 'sample_description', 'lib_name', 'lib_layout', 'lib_strategy', 'lib_source',
                    'lib_selection', 'instrument', 'total_spots', 'total_bases', 'size', 'nominal_length',
                    'nominal_sdev',
                    'spot_length', 'read_index', 'read_class', 'read_type', 'base_coord', 'lab', 'center',
                    'submitter_id',
                    'pubmed_id', 'taxid', 'published_date', 'biomaterial_provider', 'cell', 'location', 'antibody',
                    'batch',
                    'misc', 'NCBI_Link', 'AWS_Link', 'GCP_Link', ]
    id_cols = ['bioproject', 'biosample', 'experiment', 'run', 'sra_primary', 'sra_sample', 'sra_study']

    def __init__(self, column_names=column_names):
        self.config_dir = ''
        self.df = pandas.DataFrame(index=[], columns=column_names)

    def reorder(self, omit_misc=False, column_names=column_names):
        if (self.df.shape[0] == 0):
            return None
        self.df.loc[:, [col for col in column_names if col not in self.df.columns]] = ''
        if omit_misc:
            self.df = self.df.loc[:, column_names]
        else:
            misc_columns = [col for col in self.df.columns if col not in column_names]
            self.df = self.df.loc[:, column_names + misc_columns]
        self.df.loc[:, 'exclusion'] = self.df.loc[:, 'exclusion'].replace('', 'no')
        # reorder sample_group to the front
        if 'sample_group' in self.df.columns:
            cols = list(self.df)
            cols.insert(1, cols.pop(cols.index('sample_group')))
            self.df = self.df.loc[:, cols]
        self.df = self.df.reset_index(drop=True)

    def from_DataFrame(df):
        metadata = Metadata()
        metadata.df = df
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

        def get_first_external_id(entry, namespace):
            matches = entry.findall(".//EXTERNAL_ID[@namespace='{}']".format(namespace))
            if len(matches) == 0:
                return ""
            text = matches[0].text
            if text is None:
                return ""
            return str(text)

        df_list = list()
        counter = 0
        metadata = Metadata()
        for entry in root.iter(tag="EXPERIMENT_PACKAGE"):
            if counter % 1000 == 0:
                now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                print('{}: Converting {:,}th sample from XML to DataFrame'.format(now, counter), flush=True)
            bioproject = get_first_external_id(entry, 'BioProject')
            if bioproject == "":
                labels = entry.findall('.//LABEL')
                for label in labels:
                    text = label.text
                    if (text is not None) and text.startswith("PRJ"):
                        bioproject = text
                        break
            is_single = len(entry.findall('.//LIBRARY_LAYOUT/SINGLE'))
            is_paired = len(entry.findall('.//LIBRARY_LAYOUT/PAIRED'))
            if is_single:
                library_layout = "single"
            elif is_paired:
                library_layout = "paired"
            else:
                library_layout = ""
            row_df = pandas.DataFrame([{
                "bioproject": bioproject,
                "scientific_name": get_first_text(entry, './SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME'),
                "biosample": get_first_external_id(entry, 'BioSample'),
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
                "lab": get_first_attr(entry, './SUBMISSION', 'lab_name'),
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
            }])
            sas = entry.findall('./SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
            for sa in sas:
                tag = get_first_text(sa, './TAG')
                if tag != "":
                    tag = tag.lower()
                    tag = re.sub(r" \(.*", "", tag)
                    tag = re.sub(r" ", "_", tag)
                    if tag not in row_df.columns:
                        value = get_first_text(sa, './VALUE')
                        if value != "":
                            row_df[tag] = value
            df_list.append(row_df)
            counter += 1
        if len(df_list)==0:
            return metadata
        if len(df_list) <= 1000:
            df = pandas.concat(df_list, ignore_index=True)
        else:
            chunked = [pandas.concat(df_list[i:i+1000], ignore_index=True) for i in range(0, len(df_list), 1000)]
            df = pandas.concat(chunked, ignore_index=True)
        now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print('{}: Finished converting {:,} samples'.format(now, counter), flush=True)
        metadata.df = df
        metadata.reorder(omit_misc=False)
        return metadata

    def add_standard_rank_taxids(self):
        standard_ranks = ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        ncbi = ete4.NCBITaxa()
        assert self.df['taxid'].dtype == 'Int64', 'taxid column must be Int64 dtype'
        lineage_taxid_row_list = []
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
            except:
                pass
            lineage_taxid_row_list.append(lineage_taxid_row)
        lineage_taxid_df = pandas.DataFrame(lineage_taxid_row_list)
        lineage_taxid_df = lineage_taxid_df.reindex(columns=['taxid'] + ['taxid_' + rank for rank in standard_ranks], fill_value=numpy.nan)
        lineage_taxid_df = lineage_taxid_df.astype('Int64')
        self.df = self.df.merge(lineage_taxid_df, on='taxid', how='left')

    def resolve_scientific_names(self):
        self.df['scientific_name_original'] = self.df['scientific_name']
        ncbi = ete4.NCBITaxa()
        assert self.df['taxid'].dtype == 'Int64', 'taxid column must be Int64 dtype'
        taxid2sciname = ncbi.get_taxid_translator(self.df['taxid'].dropna().unique().tolist())
        self.df['scientific_name'] = self.df['taxid'].map(taxid2sciname).fillna(self.df['scientific_name_original'])

    def group_attributes(self, dir_config):
        try:
            config = pandas.read_csv(os.path.join(dir_config, 'group_attribute.config'),
                                     parse_dates=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
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
        try:
            config = pandas.read_csv(os.path.join(dir_config, 'exclude_keyword.config'),
                                     parse_dates=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        if config.shape[0]>0:
            print('{}: Marking SRAs with bad keywords'.format(datetime.datetime.now()), flush=True)
        for i in config.index:
            cols = config.iloc[i, 0].split(',')
            reason = config.iloc[i, 1]
            exclude_keyword = config.iloc[i, 2]
            num_detected = 0
            for col in cols:
                has_bad_keyword = self.df.loc[:, col].astype(str).str.contains(exclude_keyword, regex=True, case=False).fillna(False)
                self.df.loc[has_bad_keyword, 'exclusion'] = reason
                num_detected += has_bad_keyword.sum()
            txt = '{}: Marking {:,} SRAs with keyword "{}"'
            print(txt.format(datetime.datetime.now(), num_detected, exclude_keyword), flush=True)
        self.df.loc[~((self.df.loc[:, 'antibody'].isnull()) | (
                    self.df.loc[:, 'antibody'] == '')), 'exclusion'] = 'immunoprecipitation'
        self.df.loc[~((self.df.loc[:, 'cell'].isnull()) | (self.df.loc[:, 'cell'] == '')), 'exclusion'] = 'cell_culture'

    def mark_treatment_terms(self, dir_config):
        try:
            config = pandas.read_csv(os.path.join(dir_config, 'control_term.config'),
                                     parse_dates=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        if config.shape[0]>0:
            print('{}: Marking SRAs with non-control terms'.format(datetime.datetime.now()), flush=True)
        for i in config.index:
            cols = config.iloc[i, 0].split(',')
            control_term = config.iloc[i, 1]
            if control_term == '':
                continue
            num_control = 0
            num_treatment = 0
            for col in cols:
                is_control = self.df.loc[:, col].astype(str).str.contains(control_term, regex=True, case=False).fillna(False)
                if not any(is_control):
                    continue
                bioprojects = self.df.loc[is_control, 'bioproject'].unique()
                for bioproject in bioprojects:
                    is_bioproject = (self.df.loc[:, 'bioproject'] == bioproject)
                    self.df.loc[(is_bioproject & -is_control), 'exclusion'] = 'non_control'
                    num_control += (is_bioproject & is_control).sum()
                    num_treatment += (is_bioproject & -is_control).sum()
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
        is_empty = (self.df['taxid_' + rank_name].isna())
        self.df.loc[is_empty, 'exclusion'] = 'missing_taxid'

    def mark_redundant_biosample(self, exe_flag):
        if exe_flag:
            print('{}: Marking SRAs with redundant BioSample IDs'.format(datetime.datetime.now()), flush=True)
            redundant_bool = self.df.duplicated(subset=['bioproject', 'biosample'], keep='first')
            self.df.loc[redundant_bool, 'exclusion'] = 'redundant_biosample'

    def _maximize_bioproject_sampling(self, df, target_n=10):
        while len(df.loc[(df.loc[:, 'is_sampled'] == 'yes') & (df.loc[:, 'exclusion'] == 'no'), :]) < target_n:
            if len(df) <= target_n:
                df.loc[(df.loc[:, 'exclusion'] == 'no'), 'is_sampled'] = 'yes'
                break
            else:
                df_unselected = df.loc[(df.loc[:, 'is_sampled'] == 'no') & (df.loc[:, 'exclusion'] == 'no'), :]
                bioprojects = df_unselected.loc[:, 'bioproject'].unique()
                if len(bioprojects) == 0:
                    break
                remaining_n = target_n - (df.loc[:, 'is_sampled'] == 'yes').sum()
                select_n = min([len(bioprojects), remaining_n])
                selected_bioprojects = numpy.random.choice(bioprojects, size=select_n, replace=False)
                selected_index = []
                for bioproject in selected_bioprojects:
                    is_bp = (df_unselected.loc[:, 'bioproject'] == bioproject)
                    index = numpy.random.choice(df_unselected.index[is_bp], size=1, replace=False)
                    selected_index.append(int(index.item()))
                df.loc[selected_index, 'is_sampled'] = 'yes'
        return df

    def label_sampled_data(self, max_sample=10):
        pandas.set_option('mode.chained_assignment', None)
        txt = '{}: Selecting subsets of SRA IDs for >{:,} samples per sample_group per species'
        print(txt.format(datetime.datetime.now(), max_sample), flush=True)
        self.df['sample_group'] = self.df['sample_group'].fillna('').astype(str)
        self.df['exclusion'] = self.df['exclusion'].fillna('no').astype(str)
        self.df['bioproject'] = self.df['bioproject'].fillna('unknown').astype(str).values
        self.df['is_sampled'] = 'no'
        self.df['is_qualified'] = 'no'
        is_empty = (self.df['sample_group'] == '')
        self.df.loc[is_empty,'exclusion'] = 'no_tissue_label'
        self.df.loc[(self.df.loc[:, 'exclusion'] == 'no'), 'is_qualified'] = 'yes'
        df_list = list()
        species = self.df.loc[:, 'scientific_name'].unique()
        for sp in species:
            sp_table = self.df.loc[(self.df.loc[:, 'scientific_name'] == sp), :]
            sample_groups = sp_table.loc[:, 'sample_group'].unique()
            for sample_group in sample_groups:
                sp_sample_group = sp_table.loc[(sp_table.loc[:, 'sample_group'] == sample_group), :]
                if sp_sample_group.shape[0] == 0:
                    continue
                sp_sample_group = self._maximize_bioproject_sampling(df=sp_sample_group, target_n=max_sample)
                df_list.append(sp_sample_group)
        if len(df_list) <= 100:
            self.df = pandas.concat(df_list, ignore_index=True)
        else:
            chunked = [pandas.concat(df_list[i:i+100], ignore_index=True) for i in range(0, len(df_list), 100)]
            self.df = pandas.concat(chunked, ignore_index=True)
        self.reorder(omit_misc=False)
        pandas.set_option('mode.chained_assignment', 'warn')

    def remove_specialchars(self):
        for col, dtype in zip(self.df.dtypes.index, self.df.dtypes.values):
            if any([key in str(dtype) for key in ['str', 'object']]):
                self.df.loc[:, col] = self.df[col].replace(r'\r', '', regex=True)
                self.df.loc[:, col] = self.df[col].replace(r'\n', '', regex=True)
                self.df.loc[:, col] = self.df[col].replace(r'\'', '', regex=True)
                self.df.loc[:, col] = self.df[col].replace(r'\"', '', regex=True)
                self.df.loc[:, col] = self.df[col].replace(r'\|', '', regex=True)

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
    try:
        df = pandas.read_csv(os.path.join(dir_path, file_name),
                             parse_dates=False, quotechar='"', sep='\t',
                             header=None, index_col=None, skip_blank_lines=True, comment='#')
    except:
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
    print('{}: Loading metadata from: {}'.format(datetime.datetime.now(), real_path), flush=True)
    df = pandas.read_csv(real_path, sep='\t', header=0, low_memory=False)
    metadata = Metadata.from_DataFrame(df)
    if 'batch' not in dir(args):
        return metadata
    if args.batch is None:
        return metadata
    # --batch must be handled species-wise in curate.py
    # so we need to find out where the call came from
    frm = inspect.stack()[1]
    mod = inspect.getmodule(frm[0])
    if mod.__name__ == 'amalgkit.curate':
        print('Entering --batch mode for amalgkit curate. processing 1 species', flush=True)
        txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table.'
        spp = metadata.df.loc[:, 'scientific_name'].drop_duplicates().sort_values().values
        print(txt.format(args.batch, len(spp)), flush=True)
        sp = spp[args.batch - 1]
        print('Processing species: {}'.format(sp), flush=True)
        is_sp = (metadata.df['scientific_name'] == sp)
        metadata.df = metadata.df.loc[is_sp,:].reset_index(drop=True)
        return metadata
    else:
        print('--batch is specified. Processing one SRA per job.', flush=True)
        is_sampled = numpy.array([strtobool(yn) for yn in df.loc[:, 'is_sampled']], dtype=bool)
        txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table. {:,} '
        txt += 'SRAs were excluded from the table (is_sampled==no).'
        print(txt.format(args.batch, sum(is_sampled), len(numpy.where(is_sampled == False)[0])), flush=True)
        if args.batch>sum(is_sampled):
            sys.stderr.write('--batch {} is too large. Exiting.\n'.format(args.batch))
            sys.exit(0)
        if is_sampled.sum()==0:
            print('No sample is "sampled". Please check the "is_sampled" column in the metadata. Exiting.')
            sys.exit(1)
        metadata.df = metadata.df.loc[is_sampled,:]
        metadata.df = metadata.df.reset_index()
        metadata.df = metadata.df.loc[[args.batch-1,],:]
        return metadata

def get_sra_stat(sra_id, metadata, num_bp_per_sra=None):
    sra_stat = dict()
    sra_stat['sra_id'] = sra_id
    is_sra = (metadata.df.loc[:,'run']==sra_id)
    assert is_sra.sum()==1, 'There are multiple metadata rows with the same SRA ID: '+sra_id
    layout = metadata.df.loc[is_sra,'lib_layout'].values[0]
    if 'layout_amalgkit' in metadata.df.columns:
        layout_amalgkit = metadata.df.loc[is_sra,'layout_amalgkit'].values[0]
        if layout_amalgkit in ['single', 'paired']:
            layout = layout_amalgkit
    sra_stat['layout'] = layout
    sra_stat['total_spot'] = int(metadata.df.loc[is_sra,'total_spots'].values[0])
    original_spot_len = metadata.df.loc[is_sra,'spot_length'].values[0]
    if (numpy.isnan(original_spot_len) | (original_spot_len==0)):
        inferred_spot_len = int(metadata.df.loc[is_sra,'total_bases'].values[0]) / int(sra_stat['total_spot'])
        sra_stat['spot_length'] = int(inferred_spot_len)
        txt = 'spot_length cannot be obtained directly from metadata. Using total_bases/total_spots instead: {:,}'
        print(txt.format(sra_stat['spot_length']))
    else:
        sra_stat['spot_length'] = int(original_spot_len)
    if num_bp_per_sra is not None:
        sra_stat['num_read_per_sra'] = int(num_bp_per_sra/sra_stat['spot_length'])
    return sra_stat

def get_newest_intermediate_file_extension(sra_stat, work_dir):
    ext_out = 'no_extension_found'
    # Order is important in this list. More downstream should come first.
    extensions = ['.amalgkit.fastq.gz','.rename.fastq.gz','.fastp.fastq.gz','.fastq.gz']
    if 'getfastq_sra_dir' not in sra_stat:
        sra_stat['getfastq_sra_dir'] = work_dir
    sra_stat = detect_layout_from_file(sra_stat)
    if sra_stat['layout']=='single':
        subext = ''
    elif sra_stat['layout']=='paired':
        subext = '_1'
    files = os.listdir(work_dir)
    for ext in extensions:
        if any([ f==sra_stat['sra_id']+subext+ext for f in files ]):
            ext_out = ext
            break
    if ext_out == 'no_extension_found':
        safe_delete_files = glob.glob(os.path.join(work_dir, sra_stat['sra_id']+"*.safely_removed"))
        if len(safe_delete_files):
            txt = 'getfastq safely_removed flag was detected. `amalgkit quant` has been completed in this sample: {}\n'
            sys.stdout.write(txt.format(work_dir))
            for safe_delete_file in safe_delete_files:
                sys.stdout.write('{}\n'.format(safe_delete_file))
            return '.safely_removed'
    return ext_out

def is_there_unpaired_file(sra_stat, extensions):
    is_unpaired_file = False
    for ext in extensions:
        single_fastq_file = os.path.join(sra_stat['getfastq_sra_dir'], sra_stat['sra_id'] + ext)
        if os.path.exists(single_fastq_file):
            is_unpaired_file = True
            break
    return is_unpaired_file

def detect_layout_from_file(sra_stat):
    # Order is important in this list. More downstream should come first.
    extensions = ['.amalgkit.fastq.gz.safely_removed','.amalgkit.fastq.gz','.rename.fastq.gz','.fastp.fastq.gz','.fastq.gz']
    is_paired_end = False
    for ext in extensions:
        paired_fastq_files = [
            os.path.join(sra_stat['getfastq_sra_dir'], sra_stat['sra_id'] + '_1'+ext),
            os.path.join(sra_stat['getfastq_sra_dir'], sra_stat['sra_id'] + '_2'+ext),
        ]
        if all([os.path.exists(f) for f in paired_fastq_files]):
            is_paired_end = True
            break
    is_unpaired_file = is_there_unpaired_file(sra_stat, extensions)
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

def write_updated_metadata(metadata, outpath, args):
    if os.path.exists(outpath):
        print('Updated metadata file was detected. Will be overwritten: {}'.format(outpath), flush=True)
    quant_dir = os.path.join(args.out_dir, 'quant')
    metadata = get_mapping_rate(metadata, quant_dir)
    print('Writing curate metadata containing mapping rate: {}'.format(outpath))
    metadata.df.to_csv(outpath, sep='\t', index=False)

def get_mapping_rate(metadata, quant_dir):
    if os.path.exists(quant_dir):
        print('quant directory found: {}'.format(quant_dir))
        metadata.df.loc[:, 'mapping_rate'] = numpy.nan
        sra_ids = metadata.df.loc[:, 'run'].values
        sra_dirs = [d for d in os.listdir(quant_dir) if d in sra_ids]
        print('Number of quant sub-directories that matched to metadata: {:,}'.format(len(sra_dirs)))
        for sra_id in sra_dirs:
            run_info_path = os.path.join(quant_dir, sra_id, sra_id + '_run_info.json')
            if not os.path.exists(run_info_path):
                sys.stderr.write('run_info.json not found. Skipping {}.\n'.format(sra_id))
                continue
            is_sra = (metadata.df.loc[:, 'run'] == sra_id)
            with open(run_info_path) as f:
                run_info = json.load(f)
            metadata.df.loc[is_sra, 'mapping_rate'] = run_info['p_pseudoaligned']
    else:
        txt = 'quant directory not found. Mapping rate cutoff will not be applied: {}\n'
        sys.stderr.write(txt.format(quant_dir))
    return metadata

def check_rscript():
    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        print("R (Rscript) is not installed. Exiting.")
        sys.exit(1)

def orthogroup2genecount(file_orthogroup, file_genecount, spp):
    df = pandas.read_csv(file_orthogroup, sep='\t', header=0, low_memory=False)
    orthogroup_df = pandas.DataFrame({'orthogroup_id': df['busco_id'].to_numpy()})
    is_spp = df.columns.isin(spp)
    df = df.loc[:,is_spp]
    df[df.isnull()] = ''
    df[df=='-'] = ''
    gc = pandas.DataFrame(0, index=df.index, columns=df.columns)
    no_comma = (df != '') & (~df.apply(lambda x: x.str.contains(',')))
    gc[no_comma] = 1
    has_comma = df.apply(lambda x: x.str.contains(','))
    gc[has_comma] = df[has_comma].apply(lambda x: x.str.count(',') + 1)
    gc = pandas.concat([orthogroup_df, gc], axis=1)
    col_order = ['orthogroup_id'] + [col for col in gc.columns if col != 'orthogroup_id']
    gc = gc[col_order]
    gc.to_csv(file_genecount, index=False, sep='\t')

def check_ortholog_parameter_compatibility(args):
    if (args.orthogroup_table is None)&(args.dir_busco is None):
        raise Exception('One of --orthogroup_table and --dir_busco should be specified.')
    if (args.orthogroup_table is not None)&(args.dir_busco is not None):
        raise Exception('Only one of --orthogroup_table and --dir_busco should be specified.')

def generate_multisp_busco_table(dir_busco, outfile):
    print('Generating multi-species BUSCO table.', flush=True)
    col_names = ['busco_id', 'status', 'sequence', 'score', 'length', 'orthodb_url', 'description']
    species_infiles = [f for f in os.listdir(path=dir_busco) if f.endswith('.tsv')]
    species_infiles = sorted(species_infiles)
    print('BUSCO full tables for {} species were detected at: {}'.format(len(species_infiles), dir_busco), flush=True)
    for species_infile in species_infiles:
        path_to_table = os.path.join(dir_busco, species_infile)
        if not os.path.exists(path_to_table):
            warnings.warn('full_table.tsv does not exist. Skipping: {}'.format(species_infile))
            continue
        tmp_table = pandas.read_table(path_to_table, sep='\t', header=None, comment='#', names=col_names)
        tmp_table.loc[:, 'sequence'] = tmp_table.loc[:, 'sequence'].str.replace(r':[-\.0-9]*$', '', regex=True)
        for col in ['sequence', 'orthodb_url', 'description']:
            tmp_table[col] = tmp_table[col].fillna('').astype(str)
            tmp_table.loc[(tmp_table[col]==''), col] = '-'
        if species_infile == species_infiles[0]:
            merged_table = tmp_table.loc[:, ['busco_id', 'orthodb_url', 'description']]
            merged_table = merged_table.drop_duplicates(keep='first', inplace=False, ignore_index=True)
        else:
            is_mt_missing = (merged_table.loc[:, 'orthodb_url'] == '-')
            if is_mt_missing.sum() > 0:
                tmp_table2 = tmp_table.loc[:, ['busco_id', 'orthodb_url', 'description']]
                tmp_table2 = tmp_table2.drop_duplicates(keep='first', inplace=False, ignore_index=True)
                merged_table.loc[is_mt_missing, 'orthodb_url'] = tmp_table2.loc[is_mt_missing, 'orthodb_url']
                merged_table.loc[is_mt_missing, 'description'] = tmp_table2.loc[is_mt_missing, 'description']
        tmp_table = tmp_table.loc[:, ['busco_id', 'sequence']].groupby(['busco_id'])['sequence'].apply(
            lambda x: ','.join(x))
        tmp_table = tmp_table.reset_index()
        species_colname = species_infile
        species_colname = re.sub(r'\.tsv(?:\.gz)?$', '', species_colname, flags=re.IGNORECASE)
        species_colname = re.sub(r'(_busco|_full_table.*)$', '', species_colname, flags=re.IGNORECASE)
        matched = re.match(r'^([^_]+_[^_]+)', species_colname)
        if matched is not None:
            species_colname = matched.group(1)
        tmp_table = tmp_table.rename(columns={'sequence': species_colname})
        merged_table = merged_table.merge(tmp_table, on='busco_id', how='outer')
    merged_table.to_csv(outfile, sep='\t', index=None, doublequote=False)

def check_config_dir(dir_path, mode):
    files = os.listdir(dir_path)
    if mode=='select':
        asserted_files = [
            'group_attribute.config',
            'exclude_keyword.config',
            'control_term.config',
        ]
    missing_count = 0
    for af in asserted_files:
        if af in files:
            print('Config file found: {}'.format(af))
        else:
            sys.stderr.write('Config file not found: {}\n'.format(af))
            missing_count += 1
    if (missing_count>0):
        txt = 'Please refer to the AMALGKIT Wiki for more info: https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata\n'
        sys.stderr.write(txt)

def get_getfastq_run_dir(args, sra_id):
    amalgkit_out_dir = os.path.realpath(args.out_dir)
    run_output_dir = os.path.join(amalgkit_out_dir, 'getfastq', sra_id)
    if not os.path.exists(run_output_dir):
        os.makedirs(run_output_dir)
    return run_output_dir
