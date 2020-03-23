from amalgkit.metadata import Metadata
import pandas

def load_metadata(args):
    df = pandas.read_csv(args.metadata, sep='\t', header=0)
    metadata = Metadata.from_DataFrame(df)
    if args.batch is not None:
        print('--batch is specified. Processing one SRA per job.')
        is_sampled = (metadata.df.loc[:,'is_sampled']=='Yes')
        txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table. {:,} SRAs were excluded.'
        print(txt.format(args.batch, is_sampled.sum(), (is_sampled==False).sum()))
        assert args.batch<is_sampled.sum(), '--batch ({}) is too large.'.format(args.batch)
        metadata.df = metadata.df.loc[is_sampled,:]
        metadata.df = metadata.df.reset_index()
        metadata.df = metadata.df.loc[[args.batch,],:]
    return metadata