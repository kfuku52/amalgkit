"""Small synthetic GSA pages with the structure observed on public provider pages."""

import gzip
import json

from amalgkit.gsa import GSA_ROOT

RUN = 'CRR0001'
EXPERIMENT = 'CRX0001'
EXP_URL = GSA_ROOT + 'browse/CRA0001/CRX0001'
RUN_URL = GSA_ROOT + 'browse/CRA0001/CRR0001'
SAMPLE_URL = 'https://ngdc.cncb.ac.cn/biosample/browse/SAMC0001'


def source_url(filename):
    return 'https://download.cncb.ac.cn/gsa2/CRA0001/CRR0001/' + filename


def fake_pages(language='en', paired=True):
    accession, organism, title, tissue = (('Accession', 'Organism', 'Title', 'Tissue') if language == 'en'
                                         else ('实验编号', '物种名称', '标题', '组织器官'))
    experiment = f'''<table><tr><th>{accession}</th><td>CRX0001</td></tr>
        <tr><th>{organism}</th><td><a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=3702">Arabidopsis thaliana</a></td></tr>
        <tr><th>{title}</th><td>Leaf RNA</td></tr>
        <tr><th>BioSample</th><td><a href="{SAMPLE_URL}">SAMC0001</a></td></tr>
        <tr><th>BioProject</th><td><a href="https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA0001">PRJCA0001</a></td></tr>
        <tr><th>Run</th><td><a href="browse/CRA0001/CRR0001">CRR0001</a></td></tr></table>'''
    sample_accession = 'Accession' if language == 'en' else '样本编号'
    sample = f'''<table><tr><th>{sample_accession}</th><td>SAMC0001</td></tr>
        <tr><th>Attributes</th><td><table><tr><th>{tissue}</th><td>leaf</td></tr>
        <tr><th>Genotype</th><td>Col-0</td></tr></table></td></tr></table>'''
    layout = 'PAIRED' if paired else 'SINGLE'
    names = ['CRR0001_f1.fastq.gz', 'CRR0001_r2.fastq.gz'] if paired else ['CRR0001.fastq.gz']
    links = ''.join(f'<a href="{source_url(name)}">{name}</a>' for name in names)
    run = f'''<table><tr><td>CRR0001</td><td>Leaf RNA</td><td>fastq</td><td>2024-01-01</td><td>2024-02-01</td></tr></table>
        <table><tr><td><a href="browse/CRA0001/CRX0001">CRX0001</a></td><td>RNA library</td>
        <td>Illumina NovaSeq 6000</td><td>RNA-Seq</td><td>TRANSCRIPTOMIC</td><td>PolyA</td><td>{layout}</td></tr></table>
        {links}'''
    search = '<a href="browse/CRA0001/CRX0001">CRX0001</a><a href="browse/CRA0001/CRR0001">CRR0001</a>'
    return {EXP_URL: experiment, SAMPLE_URL: sample, RUN_URL: run, GSA_ROOT + 'search?searchTerm=CRR0001': search}


def search_payload(entries, total=None):
    return json.dumps({'code': '200', 'result': {'data': {
        'recordsFiltered': len(entries) if total is None else total, 'data': entries, 'error': None,
    }}})


def manifest_row(paired=True, groups=1, run=RUN):
    files = []
    for group in range(groups):
        for mate in ([1, 2] if paired else [0]):
            name = f'{run}_L{group}' + (f'_r{mate}' if paired else '') + '.fastq.gz'
            files.append({'filename': name, 'mate': mate, 'group': group,
                          'sources': [{'source_name': 'GSA', 'url': source_url(name).replace('/CRR0001/', '/' + run + '/')} ]})
    return {'run': run, 'lib_layout': 'paired' if paired else 'single', 'data_source': 'gsa',
            'data_format': 'fastq', 'gsa_fastq_files': json.dumps(files), 'scientific_name': 'Arabidopsis thaliana',
            'private_file': 'no', 'instrument': 'Illumina NovaSeq 6000', 'total_spots': '', 'total_bases': '',
            'spot_length': '', 'read_count_status': 'unknown', 'exclusion': 'no', 'is_sampled': 'yes',
            'biosample': 'SAMC0001', 'bioproject': 'PRJCA0001', 'tissue': 'leaf', 'sample_group': 'leaf'}


def fastq_bytes(count=4, mate=0, offset=0, length=10):
    return b''.join(f'@read{i + offset}' .encode() + (f'/{mate}'.encode() if mate else b'') + b'\n'
                    + b'A' * length + b'\n+\n' + b'I' * length + b'\n' for i in range(count))


def payloads_for_row(row, count=4, length=10):
    return {entry['filename']: gzip.compress(fastq_bytes(count, entry['mate'], entry['group'] * count, length), mtime=0)
            for entry in json.loads(row['gsa_fastq_files'])}
