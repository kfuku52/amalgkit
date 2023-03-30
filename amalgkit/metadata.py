from Bio import Entrez
import lxml.etree
import numpy
import pandas

import datetime
import os
import re
import sys
import time
import warnings

from urllib.error import HTTPError

def check_config_dir(dir_path):
    files = os.listdir(dir_path)
    asserted_files = [
        'group_attribute.config',
        'group_tissue.config',
        'give_value.config',
        'replace_value.config',
        'exclude_keyword.config',
        'control_term.config',
        'rescue_id.config',
        'search_term_exclusion.config',
        'search_term_other.config',
        'search_term_species.config',
        'search_term_keyword.config',
        'orthographical_variant.config',
    ]
    missing_count = 0
    for af in asserted_files:
        if af in files:
            print('Config file found: {}'.format(af))
        else:
            sys.stderr.write('Config file not found: {}\n'.format(af))
            missing_count += 1
    if missing_count>0:
        txt = 'Please refer to the AMALGKIT Wiki for more info: https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata\n'
        sys.stderr.write(txt)

def get_search_term(species_name="", bioprojects=[], biosamples=[], keywords=[],
                    other_conditions=pandas.DataFrame(), excluded_conditions=pandas.DataFrame()):
    assert ((len(bioprojects)>0)+(len(biosamples)>0))!=2, "bioprojects and biosamples cannot be specified simultaneously."

    if species_name == '':
        species_term = ''
    else:
        species_term = '"'+species_name+'"'+"[Organism]"
    keyword_term = "(" + " OR ".join(keywords) + ")"

    other_terms = list()
    for i in numpy.arange(other_conditions.shape[0]):
        other_terms.append('\"'+other_conditions.loc[i,1]+'\"['+other_conditions.loc[i,0]+']')
    other_term = " AND ".join(other_terms)

    excluded_terms = list()
    for i in numpy.arange(excluded_conditions.shape[0]):
        excluded_terms.append('\"'+excluded_conditions.loc[i,1]+'\"['+excluded_conditions.loc[i,0]+']')
    excluded_term = '('+" OR ".join(excluded_terms)+')'

    if len(bioprojects):
        bioproject_term = "(" + " OR ".join(bioprojects) + ")"
        search_term = species_term + " AND " + bioproject_term + " AND " + other_term + " NOT " + excluded_term
    elif len(biosamples):
        biosample_term = "(" + " OR ".join(biosamples) + ")"
        search_term = biosample_term
    else:
        search_term = species_term + " AND " + keyword_term + " AND " + other_term + " NOT " + excluded_term
    search_term = re.sub(" AND  AND ", " AND ", search_term)
    search_term = re.sub(" AND  NOT ", " NOT ", search_term)
    search_term = re.sub("^ AND ", "", search_term)
    return search_term

def fetch_sra_xml(species_name, search_term, save_xml=True, read_from_existing_file=False, retmax=100):
    file_xml = "SRA_"+species_name.replace(" ", "_")+".xml"
    flag = True
    if (read_from_existing_file)&(os.path.exists(file_xml)):
        with open(file_xml) as f:
            if '<Error>' in f.read():
                print(species_name, ': <Error> found in the saved file. Deleting...')
                os.remove(file_xml)
            else:
                print(species_name, ': reading xml from file')
                root = lxml.etree.parse(file_xml, parser=lxml.etree.XMLParser())
                flag = False
    if flag:
        try:
            sra_handle = Entrez.esearch(db="sra", term=search_term, retmax=10000000)
        except HTTPError as e:
            print(e, '- Trying Entrez.esearch() again...')
            sra_handle = Entrez.esearch(db="sra", term=search_term, retmax=10000000)
        sra_record = Entrez.read(sra_handle)
        record_ids = sra_record["IdList"]
        num_record = len(record_ids)
        print('Number of SRA records:', num_record)
        start_time = time.time()
        query_search_time = datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')
        root = None
        for i in numpy.arange(numpy.ceil(num_record//retmax)+1):
            start = int(i*retmax)
            end = int(((i+1)*retmax)-1) if num_record >= int(((i+1)*retmax)-1) else num_record
            print('Processing SRA records:', start, '-', end, flush=True)
            max_retry = 10
            for i in range(max_retry):
                try:
                    handle = Entrez.efetch(db="sra", id=record_ids[start:end], rettype="full", retmode="xml", retmax=retmax)
                except HTTPError as e:
                    sleep_second = 60
                    print('{} - Trying Entrez.efetch() again after {} seconds...'.format(e, sleep_second), flush=True)
                    time.sleep(sleep_second)
                    continue
                try:
                    chunk = lxml.etree.parse(handle).getroot()
                except:
                    print('XML may be truncated. Retrying...', flush=True)
                    continue
                break
            if root is None:
                root = chunk
            else:
                root.append(chunk)
        elapsed_time = int(time.time() - start_time)
        xml_string = lxml.etree.tostring(root, pretty_print=True)
        for line in str(xml_string).split('\n'):
            if '<Error>' in line:
                print(line)
                if os.path.exists(file_xml):
                    os.remove(file_xml)
                raise Exception(species_name, ': <Error> found in the xml.')
        if save_xml:
            with open(file_xml, 'wb') as f:
                f.write(xml_string)
    return root

def create_run_dir(run_path, run_no=1):
        if not os.path.exists(os.path.join(run_path + '_' + str(run_no))):
            os.makedirs(os.path.join(run_path + '_' + str(run_no)))
            return(os.path.join(run_path + '_' + str(run_no)))
        else:
            run_no = run_no + 1
            return create_run_dir(run_path, run_no)

class Metadata:
    column_names = ['scientific_name','tissue','curate_group','tissue_original','genotype','sex','age','treatment','source_name',
                    'is_sampled','is_qualified','exclusion','protocol','bioproject','biosample',
                    'experiment','run','sra_primary','sra_sample','sra_study','study_title','exp_title','design',
                    'sample_title','sample_description','lib_name','lib_layout','lib_strategy','lib_source',
                    'lib_selection','instrument','total_spots','total_bases','size','nominal_length','nominal_sdev',
                    'spot_length','read_index','read_class','read_type','base_coord','lab','center','submitter_id',
                    'pubmed_id','taxid','published_date','biomaterial_provider','cell','location','antibody','batch',
                    'misc','NCBI_Link','AWS_Link','GCP_Link',]
    id_cols = ['bioproject','biosample','experiment','run','sra_primary','sra_sample','sra_study']

    def __init__(self, column_names=column_names):
        self.config_dir = ''
        self.df = pandas.DataFrame(index=[], columns=column_names)

    def reorder(self, omit_misc=False, column_names=column_names):
        if (self.df.shape[0]==0):
            return None
        for col in column_names:
            if not col in self.df.columns:
                self.df.loc[:,col] = ''
        if omit_misc:
            self.df = self.df.loc[:,column_names]
        else:
            misc_columns = [ col for col in self.df.columns if col not in column_names ]
            self.df = self.df.loc[:,column_names+misc_columns]
        self.df.loc[:,'exclusion'] = self.df.loc[:,'exclusion'].replace('', 'no')
        # reorder curate_group to the front
        if 'curate_group' in self.df.columns:
            cols = list(self.df)
            cols.insert(1, cols.pop(cols.index('curate_group')))
            self.df = self.df.loc[:, cols]

    def from_DataFrame(df):
        metadata = Metadata()
        metadata.df = df
        metadata.reorder(omit_misc=False)
        return metadata
        
    def from_xml(xml_root):
        if isinstance(xml_root, lxml.etree._Element):
            xml_root = lxml.etree.ElementTree(xml_root)
        root = xml_root
        assert isinstance(root, lxml.etree._ElementTree), "Unknown input type."
        df = pandas.DataFrame()
        for entry in root.iter(tag="EXPERIMENT_PACKAGE"):
            items = []
            bioproject = entry.findall('.//EXTERNAL_ID[@namespace="BioProject"]')
            if not len(bioproject):
                labels = entry.findall('.//LABEL')
                for label in labels:
                    text = label.text
                    if text.startswith("PRJ"):
                        bioproject = [label]
                        break
            is_single = len(entry.findall('.//LIBRARY_LAYOUT/SINGLE'))
            is_paired = len(entry.findall('.//LIBRARY_LAYOUT/PAIRED'))
            if is_single:
                library_layout = ["single"]
            elif is_paired:
                library_layout = ["paired"]
            else:
                library_layout = [""]
            values = entry.findall('.//VALUE')
            is_protected = ["No"]
            if len(values):
                for value in values:
                    text = value.text
                    if not text is None:
                        if text.endswith("PROTECTED"):
                            is_protected = ["Yes"]
                            break
            items.append(["bioproject", bioproject])
            items.append(["scientific_name", entry.xpath('./SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME')])
            items.append(["biosample", entry.findall('.//EXTERNAL_ID[@namespace="BioSample"]')])
            items.append(["experiment", entry.xpath('./EXPERIMENT/IDENTIFIERS/PRIMARY_ID')])
            items.append(["run", entry.xpath('./RUN_SET/RUN/IDENTIFIERS/PRIMARY_ID')])
            items.append(["sra_primary", entry.xpath('./SUBMISSION/IDENTIFIERS/PRIMARY_ID')])
            items.append(["sra_sample", entry.xpath('./SAMPLE/IDENTIFIERS/PRIMARY_ID')])
            items.append(["sra_study", entry.xpath('./EXPERIMENT/STUDY_REF/IDENTIFIERS/PRIMARY_ID')])
            items.append(["published_date", entry.xpath('./RUN_SET/RUN/@published')])
            items.append(["exp_title", entry.xpath('./EXPERIMENT/TITLE')])
            items.append(["design", entry.xpath('./EXPERIMENT/DESIGN/DESIGN_DESCRIPTION')])
            items.append(["lib_name", entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME')])
            items.append(["lib_strategy", entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY')])
            items.append(["lib_source", entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE')])
            items.append(["lib_selection", entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION')])
            items.append(["lib_layout", library_layout])
            items.append(["nominal_length", entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED/@NOMINAL_LENGTH')])
            items.append(["nominal_sdev", entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED/@NOMINAL_SDEV')])
            items.append(["spot_length", entry.xpath('./EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/SPOT_LENGTH')])
            items.append(["read_index", entry.xpath('./EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/READ_INDEX')])
            items.append(["read_class", entry.xpath('./EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/READ_CLASS')])
            items.append(["read_type", entry.xpath('./EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/READ_TYPE')])
            items.append(["base_coord", entry.xpath('./EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/BASE_COORD')])
            items.append(["instrument", entry.xpath('./EXPERIMENT/PLATFORM/ILLUMINA/INSTRUMENT_MODEL')])
            items.append(["lab", entry.xpath('./SUBMISSION/@lab_name')])
            items.append(["center", entry.xpath('./SUBMISSION/@center_name')])
            items.append(["submitter_id", entry.xpath('./SUBMISSION/IDENTIFIERS/SUBMITTER_ID')])
            items.append(["study_title", entry.xpath('./STUDY/DESCRIPTOR/STUDY_TITLE')])
            items.append(["pubmed_id", entry.xpath('./STUDY/STUDY_LINKS/STUDY_LINK/XREF_LINK/ID')])
            items.append(["sample_title", entry.xpath('./SAMPLE/TITLE')])
            items.append(["taxid", entry.xpath('./SAMPLE/SAMPLE_NAME/TAXON_ID')])
            items.append(["sample_description", entry.xpath('./SAMPLE/DESCRIPTION')])
            items.append(["total_spots", entry.xpath('./RUN_SET/RUN/@total_spots')])
            items.append(["total_bases", entry.xpath('./RUN_SET/RUN/@total_bases')])
            items.append(["size", entry.xpath('./RUN_SET/RUN/@size')])
            items.append(["NCBI_Link", entry.xpath('./RUN_SET/RUN/SRAFiles/SRAFile[@supertype="Primary ETL"]/Alternatives[@org="NCBI"]/@url')])
            items.append(["AWS_Link", entry.xpath('./RUN_SET/RUN/SRAFiles/SRAFile[@supertype="Primary ETL"]/Alternatives[@org="AWS"]/@url')])
            items.append(["GCP_Link", entry.xpath('./RUN_SET/RUN/SRAFiles/SRAFile[@supertype="Primary ETL"]/Alternatives[@org="GCP"]/@url')])

            row = []
            for item in items:
                try:
                    if isinstance(item[1][0], (lxml.etree._ElementUnicodeResult, int, str)):
                        row.append(str(item[1][0]))
                    else:
                        row.append(item[1][0].text)
                except:
                    row.append("")
            colnames = []
            for item in items:
                colnames.append(item[0])
            row_df = pandas.DataFrame(row).T
            row_df.columns = colnames
            sas = entry.xpath('./SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
            for sa in sas:
                tag = sa.xpath('./TAG')
                if not tag[0].text == None:
                    tag = tag[0].text.lower()
                    tag = re.sub(" \(.*", "", tag)
                    tag = re.sub(" ", "_", tag)
                    if not tag in row_df.columns:
                        value = sa.xpath('./VALUE')
                        if len(value):
                            value = value[0].text
                            if tag in colnames:
                                tag = tag+"_2"
                            sa_df = pandas.DataFrame([value])
                            sa_df.columns = [tag]
                            row_df = pandas.concat([row_df,sa_df], axis=1)
            df = pandas.concat([df, row_df], ignore_index=True, sort=False)
        if "scientific_name" in df.columns and len(df.loc[(df.loc[:,"scientific_name"]==""), "scientific_name"]):
            species_names = df.loc[~(df.loc[:,"scientific_name"]==""), "scientific_name"]
            if species_names.shape[0]>0:
                species_name = species_names.iloc[0]
                df.loc[(df.loc[:,"scientific_name"]==""), "scientific_name"] = species_name
        metadata = Metadata()
        metadata.df = df
        metadata.reorder(omit_misc=False)
        return metadata

    def mark_exclude_ids(self, id_cols=id_cols):
        try:
            config = pandas.read_csv(os.path.join(self.config_dir, 'exclude_id.config'),
                                 parse_dates=False, infer_datetime_format=False, quotechar='"', sep='\t',
                                 header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        for i in numpy.arange(config.shape[0]):
            reason = config.iloc[i,0]
            exclude_id = config.iloc[i,1]
            for col in id_cols:
                is_exclude_id = (self.df.loc[:,col]==exclude_id).fillna(False)
                if any(is_exclude_id):
                    self.df.loc[is_exclude_id,'exclusion'] = reason

    def unmark_rescue_ids(self, id_cols=id_cols):
        try:
            config = pandas.read_csv(os.path.join(self.config_dir, 'rescue_id.config'),
                                     parse_dates=False, infer_datetime_format=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        for i in numpy.arange(config.shape[0]):
            rescue_id = config.iloc[i,0]
            for col in id_cols:
                is_rescue_id = (self.df.loc[:,col]==rescue_id).fillna(False)
                if any(is_rescue_id):
                    self.df.loc[is_rescue_id,'exclusion'] = 'no'

    def group_attributes(self):
        try:
            config = pandas.read_csv(os.path.join(self.config_dir, 'group_attribute.config'),
                                     parse_dates=False, infer_datetime_format=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        for i in numpy.arange(config.shape[0]):
            aggregate_to = config.iloc[i,0]
            aggregate_from = config.iloc[i,1]
            if (aggregate_from in self.df.columns)&(aggregate_from!=''):
                if not aggregate_to in self.df.columns:
                    self.df.loc[:,aggregate_to] = ''
                is_from_empty = (self.df.loc[:,aggregate_from].isnull())|(self.df.loc[:,aggregate_from].astype(str)=='')
                is_to_empty = (self.df.loc[:,aggregate_to].isnull())|(self.df.loc[:,aggregate_to].astype(str)=='')
                new_annotations = self.df.loc[(~is_from_empty)&(is_to_empty), aggregate_from].astype(str)+'['+aggregate_from+']'
                self.df.loc[(~is_from_empty)&(is_to_empty), aggregate_to] = new_annotations
                new_annotations = self.df.loc[(~is_from_empty)&(~is_to_empty), aggregate_to].astype(str)+"; "+self.df.loc[(~is_from_empty)&(~is_to_empty), aggregate_from].astype(str)+'['+aggregate_from+']'
                self.df.loc[(~is_from_empty)&(~is_to_empty), aggregate_to] = new_annotations
                self.df = self.df.drop(labels=aggregate_from, axis=1)
        self.reorder(omit_misc=False)

    def correct_orthographical_variants(self):
        self.df.loc[:,"scientific_name"] = [re.sub(r'(.+)(\s)(.+)(\s)(.+)', r"\1\2\3", sp) for sp in self.df.loc[:,"scientific_name"]]
        try:
            config = pandas.read_csv(os.path.join(self.config_dir, 'orthographical_variant.config'),
                                     parse_dates=False, infer_datetime_format=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        for i in numpy.arange(config.shape[0]):
            cols = config.iloc[i,0].split(',')
            replace_from = config.iloc[i,2]
            replace_to = config.iloc[i,1]
            for col in cols:
                self.df.loc[:,col] = self.df.loc[:,col].str.replace(replace_from, replace_to, regex=True, case=False)

    def group_tissues_auto(self):
        import nltk
        from nltk.stem import WordNetLemmatizer
        import obonet
        nltk.download('wordnet')
        # retrieve metadata
        # retrieve tissue query from config
        tissues = pandas.read_csv(os.path.join(self.config_dir, 'search_term_keyword.config'),
                                    parse_dates=False, infer_datetime_format=False, quotechar='"', sep='\t',
                                    header=None, index_col=None, skip_blank_lines=True, comment='#').iloc[:,0]
        tissues = tissues.tolist()
        # init lemmatizer
        lemmatizer = WordNetLemmatizer()
        # isolate tissue from metadata, force lower case and remove tailing whitespaces
        content = self.df.loc[:,'tissue']
        content = ['' if pandas.isnull(x) else x.lower() for x in content]
        content = [x if pandas.isnull(x) else x.strip() for x in content]
        # lemmatize all words in tissue list
        lemm = [lemmatizer.lemmatize(i) for i in content]

        # load ontology
        graph = obonet.read_obo("http://purl.obolibrary.org/obo/po.obo")
        dat_list =[]

        # declare functions for later

        id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
        name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}

        def search(values, searchFor):
            for k in values:
                for v in values[k]:
                    if searchFor in v:
                       return k
            return None

        # for each item in the lemmatized tissue list, check if the query can be recovered
        for i in range(len(lemm)):

            sp = re.sub(r"[^a-zA-Z]+", ' ', lemm[i])
            sp = re.split('[_,:;" "/]', sp)
            sp = [lemmatizer.lemmatize(x) for x in sp]


            if len(sp) > 0 :

                list_set = set(sp)
                for tissue in tissues:
                    if tissue in list_set:
                        lemm[i] = tissue
                        break
                    else:
                        try:
                            lemm[i] = graph.node[name_to_id[tissue]]['name']
                        except Exception as e:
                            for id, dat in graph.nodes(data=True):
                                if search(dat, tissue) != None:
                                    lemm[i] = dat['name']


                    lemm[i] = ''.join(sp)

        self.df.loc[:,'tissue'] = lemm

    def group_tissues_by_config(self):
        try:
            config = pandas.read_csv(os.path.join(self.config_dir, 'group_tissue.config'),
                                 parse_dates=False, infer_datetime_format=False, quotechar='"', sep='\t',
                                 header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        for i in numpy.arange(config.shape[0]):
            replace_from = config.iloc[i,1]
            replace_to = config.iloc[i,0]
            is_matching = self.df.loc[:,'tissue'].str.match(replace_from, case=False, na=False)
            self.df.loc[is_matching,'tissue'] = replace_to
        self.df.loc[:,'tissue'] = self.df.loc[:,'tissue'].str.lower()
        self.df.loc[:, 'curate_group'] = self.df.loc[:, 'tissue']

    def mark_exclude_keywords(self):
        try:
            config = pandas.read_csv(os.path.join(self.config_dir, 'exclude_keyword.config'),
                                     parse_dates=False, infer_datetime_format=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        for i in numpy.arange(config.shape[0]):
            cols = config.iloc[i,0].split(',')
            reason = config.iloc[i,1]
            exclude_keyword = config.iloc[i,2]
            for col in cols:
                self.df.loc[self.df.loc[:,col].str.contains(exclude_keyword, regex=True, case=False).fillna(False), 'exclusion'] = reason
        self.df.loc[~((self.df.loc[:,'antibody'].isnull())|(self.df.loc[:,'antibody']=='')), 'exclusion'] = 'immunoprecipitation'
        self.df.loc[~((self.df.loc[:,'cell'].isnull())|(self.df.loc[:,'cell']=='')), 'exclusion'] = 'cell_culture'

    def replace_values(self):
        try:
            config = pandas.read_csv(os.path.join(self.config_dir, 'replace_value.config'),
                                     parse_dates=False, infer_datetime_format=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        for i in numpy.arange(config.shape[0]):
            col1 = config.iloc[i,0]
            col1_value = config.iloc[i,1]
            col2 = config.iloc[i,2]
            col3 = config.iloc[i,3]
            self.df.loc[(self.df.loc[:,col1]==col1_value), col2] = self.df.loc[(self.df.loc[:,col1]==col1_value), col3]

    def give_values(self, id_cols=id_cols):
        try:
            config = pandas.read_csv(os.path.join(self.config_dir, 'give_value.config'),
                                     parse_dates=False, infer_datetime_format=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        for i in numpy.arange(config.shape[0]):
            target_id = config.iloc[i,0]
            column = config.iloc[i,1]
            value = config.iloc[i,2]
            for col in id_cols:
                is_target_id = (self.df.loc[:,col]==target_id).fillna(False)
                if any(is_target_id):
                    self.df.loc[is_target_id,column] = value

    def mark_treatment_terms(self):
        try:
            config = pandas.read_csv(os.path.join(self.config_dir, 'control_term.config'),
                                     parse_dates=False, infer_datetime_format=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        for i in numpy.arange(config.shape[0]):
            col = config.iloc[i,0]
            control_term = config.iloc[i,1]
            if control_term=='':
                continue
            is_control = self.df.loc[:,col].str.contains(control_term, regex=True, case=False).fillna(False)
            if not any(is_control):
                continue
            bioprojects = self.df.loc[:,'bioproject'][is_control].unique()
            for bioproject in bioprojects:
                if bioproject == '':
                    continue
                is_bioproject = (self.df.loc[:,'bioproject']==bioproject)
                self.df.loc[(is_bioproject & -is_control), 'exclusion'] = 'treatment'

    def nspot_cutoff(self, min_nspots):
        self.df['total_spots'] = self.df.loc[:,'total_spots'].replace('', 0)
        self.df['total_spots'] = self.df.loc[:,'total_spots'].fillna(0).astype(int)
        self.df.loc[-(self.df.loc[:,'total_spots']==0) & (self.df.loc[:,'total_spots'] < min_nspots), 'exclusion'] = 'low_nspots'

    def remove_redundant_biosample(self):
        redundant_bool = self.df.duplicated(subset=['bioproject','biosample'], keep='first')
        self.df.loc[redundant_bool, 'exclusion'] = 'redundant_biosample'

    def _maximize_bioproject_sampling(self, df, target_n=10):
        pandas.set_option('mode.chained_assignment', None)
        df.loc[:, 'bioproject'] = df.loc[:, 'bioproject'].fillna('unknown').values
        df.loc[:, 'is_sampled'] = 'no'
        df.loc[:, 'is_qualified'] = 'no'
        df.loc[(df.loc[:, 'exclusion']=='no'), 'is_qualified'] = 'yes'
        while len(df.loc[(df.loc[:,'is_sampled']=='yes')&(df.loc[:,'exclusion']=='no'),:]) < target_n:
            if len(df) <= target_n:
                df.loc[(df.loc[:,'exclusion']=='no'), 'is_sampled'] = 'yes'
                break
            else:
                df_unselected = df.loc[(df.loc[:,'is_sampled']=='no')&(df.loc[:,'exclusion']=='no'),:]
                bioprojects = df_unselected.loc[:, 'bioproject'].unique()
                if len(bioprojects) == 0:
                    break
                remaining_n = target_n - (df.loc[:, 'is_sampled']=='yes').sum()
                select_n = min([len(bioprojects), remaining_n])
                selected_bioprojects = numpy.random.choice(bioprojects, size=select_n, replace=False)
                selected_index = []
                for bioproject in selected_bioprojects:
                    is_bp = (df_unselected.loc[:,'bioproject']==bioproject)
                    index = numpy.random.choice(df_unselected.index[is_bp], size=1, replace=False)
                    selected_index.append(int(index))
                df.loc[selected_index, 'is_sampled'] = 'yes'
        pandas.set_option('mode.chained_assignment', 'warn')
        return df

    def label_sampled_data(self, max_sample=10):
        df_labeled = pandas.DataFrame()
        species = self.df.loc[:,'scientific_name'].unique()
        for sp in species:
            sp_table = self.df.loc[(self.df.loc[:,'scientific_name']==sp), :]
            tissues = sp_table.loc[:,'tissue'].unique()
            for tissue in tissues:
                sp_tissue = sp_table.loc[(sp_table.loc[:,'tissue']==tissue), :]
                sp_tissue = self._maximize_bioproject_sampling(df=sp_tissue, target_n=max_sample)
                df_labeled = pandas.concat([df_labeled, sp_tissue], axis=0)
        self.df = df_labeled
        self.reorder(omit_misc=False)

    def remove_specialchars(self):
        for col,dtype in zip(self.df.dtypes.index, self.df.dtypes.values):
            if any([ key in str(dtype) for key in ['str','object'] ]):
                self.df.loc[:,col] = self.df[col].replace('\r','',regex=True)
                self.df.loc[:,col] = self.df[col].replace('\n','',regex=True)
                self.df.loc[:,col] = self.df[col].replace('\'','',regex=True)
                self.df.loc[:,col] = self.df[col].replace('\"','',regex=True)
                self.df.loc[:,col] = self.df[col].replace('\|','',regex=True)

    def pivot(self, n_sp_cutoff=0, qualified_only=True, sampled_only=False):
        df = self.df
        if qualified_only:
            df = df.loc[(df.loc[:,'is_qualified']=='yes'),:]
        if sampled_only:
            df = df.loc[(df.loc[:,'is_sampled']=='yes'),:]
        df_reduced = df.loc[:,['scientific_name', 'biosample', 'tissue']]
        pivot = df_reduced.pivot_table(columns='tissue',index='scientific_name', aggfunc='count')
        pivot.columns = pivot.columns.get_level_values(1)
        column_sort = pivot.count(axis='index').sort_values(ascending=False).index
        index_sort = pivot.count(axis='columns').sort_values(ascending=False).index
        pivot = pivot.loc[index_sort,column_sort]
        pivot_reduced = pivot.loc[:,pivot.count(axis='index') >= n_sp_cutoff]
        column_sort = pivot_reduced.count(axis='index').sort_values(ascending=False).index
        index_sort = pivot_reduced.count(axis='columns').sort_values(ascending=False).index
        pivot_reduced = pivot_reduced.loc[index_sort,column_sort]
        return pivot_reduced

def read_config_file(file_name, dir_path):
    try:
        df = pandas.read_csv(os.path.join(dir_path, file_name),
                             parse_dates=False, infer_datetime_format=False, quotechar='"', sep='\t',
                             header=None, index_col=None, skip_blank_lines=True, comment='#')
    except:
        df = pandas.DataFrame([])
    if df.shape[1]==1:
        df = df.iloc[:,0]
    return df

def metadata_main(args):
    metadata_dir = os.path.join(args.out_dir, 'metadata')
    metadata_tmp_dir = os.path.join(metadata_dir, 'tmp')
    metadata_results_dir = os.path.join(metadata_dir, 'metadata')
    pivot_table_dir = os.path.join(metadata_dir, 'pivot_tables')
    for path_dir in [args.out_dir, metadata_dir, metadata_tmp_dir, metadata_results_dir, pivot_table_dir]:
        if not os.path.exists(path_dir):
            print('Creating directory: {}'.format(path_dir))
            os.mkdir(path_dir)
    path_metadata_table = os.path.join(metadata_results_dir, 'metadata.tsv')
    if os.path.exists(path_metadata_table) & (args.overwrite==False):
        print('Exiting. --overwrite is set to "no". Output file exists at: {}'.format(path_metadata_table))
        sys.exit()

    if args.config_dir=='inferred':
        dir_config = os.path.join(args.out_dir, 'config')
    else:
        dir_config = args.config_dir
    check_config_dir(dir_path=dir_config)

    Entrez.email = args.entrez_email

    search_spp = read_config_file(file_name='search_term_species.config', dir_path=dir_config)
    print('Number of species for Entrez search: {:,}'.format(search_spp.shape[0]))
    if search_spp.shape[0]==0:
        search_spp = pandas.Series(['',])

    search_keywords = read_config_file(file_name='search_term_keyword.config', dir_path=dir_config)
    print('Number of free keywords for Entrez search: {:,}'.format(search_keywords.shape[0]))

    other_conditions = read_config_file(file_name='search_term_other.config', dir_path=dir_config)
    print('Number of other conditions for Entrez search: {:,}'.format(other_conditions.shape[0]))

    excluded_conditions = read_config_file(file_name='search_term_exclusion.config', dir_path=dir_config)
    print('Number of excluded conditions for Entrez search: {:,}'.format(excluded_conditions.shape[0]))

    metadata_species = dict()
    for sp in search_spp:
        print('Entrez search:', sp)
        sp_file_name = "tmp_"+sp.replace(' ', '_')+".tsv"
        bioprojects = []
        if (os.path.exists(sp_file_name))&(args.overwrite=='no'):
            if (os.path.getsize(sp_file_name)>1):
                print(sp, ': reading from tsv file.')
                metadata_species[sp] = pandas.read_csv(sp_file_name, sep='\t', header=0, low_memory=False)
            else:
                print(sp, ': empty tsv file. skipped.')
        else:
            search_term = get_search_term(species_name=sp, bioprojects=bioprojects, biosamples=[], keywords=search_keywords,
                                          other_conditions=other_conditions, excluded_conditions=excluded_conditions)
            print('Entrez search term:', search_term)
            root = fetch_sra_xml(species_name=sp, search_term=search_term, save_xml=False, read_from_existing_file=False)
            metadata_species[sp] = Metadata.from_xml(xml_root=root).df
            metadata_species[sp].to_csv(os.path.join(metadata_tmp_dir, sp_file_name), sep="\t", index=False)
        print('')
    metadata = Metadata.from_DataFrame(df=pandas.concat(metadata_species.values(), sort=False))
    metadata.config_dir = dir_config
    metadata.reorder(omit_misc=False)
    multisp_file_name = "metadata_01_raw.tsv"
    metadata.df.to_csv(os.path.join(metadata_tmp_dir, multisp_file_name), sep="\t", index=False)
    del metadata_species
    if metadata.df.shape[0]==0:
        print('No entry was found/survived in the metadata processing. Please check the config files.')
        return None
    metadata.df.loc[:, 'tissue_original'] = metadata.df.loc[:, 'tissue']
    if args.tissue_detect:
        metadata.remove_linebreak()
        # metadata.mark_exclude_ids() # TODO to Matthias, this should be activated even when --tissue_detect yes. Any conflicting feature?
        # metadata.group_attributes() # TODO to Matthias, this should be activated even when --tissue_detect yes. Any conflicting feature?
        metadata.correct_orthographical_variants()
        metadata.replace_values()
        metadata.give_values()
        #metadata.mark_exclude_keywords() # TODO to Matthias, this should be activated even when --tissue_detect yes. Any conflicting feature?
        metadata.group_tissues_auto()
    else:
        metadata.remove_specialchars()
        metadata.mark_exclude_ids()
        metadata.group_attributes()
        metadata.correct_orthographical_variants()
        metadata.replace_values()
        metadata.give_values()
        metadata.mark_exclude_keywords()
        metadata.group_tissues_by_config()

    metadata.reorder(omit_misc=False)
    metadata.df.to_csv(os.path.join(metadata_tmp_dir, 'metadata_02_grouped.tsv'), sep='\t', index=False)

    metadata.mark_treatment_terms()
    metadata.nspot_cutoff(args.min_nspots)
    metadata.remove_redundant_biosample()
    metadata.unmark_rescue_ids()
    metadata.label_sampled_data(args.max_sample)
    metadata.reorder(omit_misc=True)
    metadata.df.to_csv(path_metadata_table, sep='\t', index=False)

    sra_qualified_pivot = metadata.pivot(n_sp_cutoff=0, qualified_only=True, sampled_only=False)
    sra_qualified_pivot.to_csv(os.path.join(pivot_table_dir, 'pivot_qualified.tsv'), sep='\t')
    sra_selected_pivot = metadata.pivot(n_sp_cutoff=0, qualified_only=True, sampled_only=True)
    sra_selected_pivot.to_csv(os.path.join(pivot_table_dir, 'pivot_selected.tsv'), sep='\t')
