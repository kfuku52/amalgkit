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

from amalgkit.util import *

from urllib.error import HTTPError

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

def fetch_sra_xml(species_name, search_term, save_xml=True, read_from_existing_file=False, retmax=1000):
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
        print('Number of SRA records: {:,}'.format(num_record))
        start_time = time.time()
        print('{}: SRA XML retrieval started.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        root = None
        max_retry = 10
        for i in numpy.arange(numpy.ceil(num_record//retmax)+1):
            start = int(i*retmax)
            end = int(((i+1)*retmax)-1) if num_record >= int(((i+1)*retmax)-1) else num_record
            now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print('{}: Retrieving SRA XML: {:,}-{:,} of {:,} records'.format(now, start, end, num_record), flush=True)
            for j in range(max_retry):
                try:
                    handle = Entrez.efetch(db="sra", id=record_ids[start:end], rettype="full", retmode="xml", retmax=retmax)
                except HTTPError as e:
                    sleep_second = 60
                    print('{} - Trying Entrez.efetch() again after {:,} seconds...'.format(e, sleep_second), flush=True)
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
        print('{}: SRA XML retrieval ended.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        print('SRA XML retrieval time: {:,.1f} sec'.format(elapsed_time), flush=True)
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

def metadata_main(args):
    metadata_dir = os.path.join(args.out_dir, 'metadata')
    metadata_tmp_dir = os.path.join(metadata_dir, 'tmp')
    multisp_file_path = os.path.join(metadata_dir, "metadata.tsv")
    if os.path.exists(multisp_file_path):
        if args.redo:
            os.remove(multisp_file_path)
        else:
            print('Exiting. Output file already exists at: {}'.format(multisp_file_path))
            sys.exit(0)
    for path_dir in [args.out_dir, metadata_dir, metadata_tmp_dir]:
        if not os.path.exists(path_dir):
            print('Creating directory: {}'.format(path_dir))
            os.mkdir(path_dir)
    if args.config_dir=='inferred':
        dir_config = os.path.join(args.out_dir, 'config')
    else:
        dir_config = args.config_dir
    check_config_dir(dir_path=dir_config, mode='metadata')

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
        if sp=='':
            sp_print = 'All'
        else:
            sp_print = sp
        print('Entrez search species: {}'.format(sp_print), flush=True)
        sp_file_name = "tmp_"+sp_print.replace(' ', '_')+".tsv"
        bioprojects = []
        if (os.path.exists(sp_file_name))&(not args.redo):
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
            metadata_species[sp_print] = Metadata.from_xml(xml_root=root).df
            metadata_species[sp_print].to_csv(os.path.join(metadata_tmp_dir, sp_file_name), sep="\t", index=False)
        print('')
    metadata = Metadata.from_DataFrame(df=pandas.concat(metadata_species.values(), sort=False, ignore_index=True))
    metadata.config_dir = dir_config
    metadata.df['tissue'] = metadata.df['tissue'].astype(str)
    metadata.df.loc[(metadata.df['tissue']=='nan'), 'tissue'] = ''
    metadata.df.loc[:, 'curate_group'] = metadata.df.loc[:, 'tissue'].str.lower()
    metadata.reorder(omit_misc=False)
    metadata.df.to_csv(multisp_file_path, sep="\t", index=False)
    if metadata.df.shape[0]==0:
        sys.stderr.write('No entry was found/survived in the metadata processing. Please reconsider your config files.\n')
        return None
