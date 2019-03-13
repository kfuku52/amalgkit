
# coding: utf-8

# In[1]:


import pandas, os, argparse, sys


# In[2]:


if sys.argv[1] == "-f":
    mode = "debug"
else:
    mode = "batch"
    
if mode=='debug':
    args = {
        'gtf':'/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/gtf/Homo_sapiens.GRCh38.91.gtf.gz',
        'search_attrs':'gene_biotype|transcript_name',
        'search_values':'lincRNA|macro_lncRNA|miRNA|misc_RNA|Mt_rRNA|Mt_tRNA|rRNA|scaRNA|scRNA|snoRNA|snRNA|sRNA|tRNA|mt-.*',
        'out_attr':'transcript_id',
    }
elif mode=='batch':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--search_attrs', metavar='STR', default='gene_biotype|transcript_name', type=str, help='')
    parser.add_argument('--search_values', metavar='STR', default='lincRNA|macro_lncRNA|miRNA|misc_RNA|Mt_rRNA|Mt_tRNA|rRNA|scaRNA|scRNA|snoRNA|snRNA|sRNA|tRNA|mt-.*', type=str, help='')
    parser.add_argument('--out_attr', metavar='STR', default='transcript_id', type=str, help='')
    args = parser.parse_args()  
    g = dict()
    for attr in [a for a in dir(args) if not a.startswith('_')]:
        g[attr] = getattr(args, attr)
    args = g


# In[3]:


df = pandas.read_csv(args['gtf'], sep='\t', comment='#', header=None, low_memory=False)
df.columns = ['seqname','source','feature','start','end','score','strand','frame','attribute']
df = df.loc[(df['feature']=='transcript'),:]

for sa in args['search_attrs'].split('|'):
    potential_values = df['attribute'].replace('.*'+sa+' \"','',regex=True).replace('\".*','',regex=True).unique()
    if len(potential_values) < 100:
        print('Potential search values for', sa, ':', potential_values)
    else:
        print('Potential search values for', sa, ':', 'too many to print,', len(potential_values))


# In[4]:


df2 = pandas.DataFrame()
for sv in args['search_values'].split('|'):
    for sa in args['search_attrs'].split('|'):
        search_term = sa+' \"'+sv+'\";'
        tmp = df.loc[(df.attribute.str.contains(search_term, regex=True)),:]
        df2 = pandas.concat([df2, tmp], ignore_index=True)
df2 = df2.drop_duplicates()


# In[5]:


df3 = df2.copy()
df3[args['out_attr']] = df3['attribute'].replace('.*'+args['out_attr']+' \"','',regex=True).replace('\".*','',regex=True)
for sa in args['search_attrs'].split('|'):
    df3.loc[:,sa] = df3['attribute'].replace('.*'+sa+' \"','',regex=True).replace('\".*','',regex=True)


# In[6]:


df2.to_csv('ensembl_gtf2geneid.gtf', sep='\t', index=False, header=False)
df3.loc[:,args['out_attr']].to_csv('ensembl_gtf2geneid.out_attr.tsv', sep='\t', index=False)

