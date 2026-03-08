import argparse

from amalgkit.cli_utils import (
    build_help_command_handler,
    int_or_auto,
    nonnegative_int_or_auto,
    strtobool,
)


def build_parser(command_handlers, command_names, version):
    parser = argparse.ArgumentParser(description='A toolkit for cross-species transcriptome amalgamation')
    parser.add_argument('--version', action='version', version='amalgkit version ' + version)
    subparsers = parser.add_subparsers()

    pp_meta = argparse.ArgumentParser(add_help=False)
    pp_meta.add_argument('--metadata', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: "inferred" = out_dir/metadata/metadata.tsv. '
                          'PATH to metadata table, the output file of `amalgkit metadata`. '
                          'In wsfilter/csfilter/finalize, if "inferred" and prior filter metadata exists, the newest of '
                          'out_dir/wsfilter/metadata.tsv and out_dir/csfilter/metadata.tsv is used automatically.')
    pp_out = argparse.ArgumentParser(add_help=False)
    pp_out.add_argument('--out_dir', metavar='PATH', default='./', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to the directory where intermediate and output files are generated.')
    pp_batch = argparse.ArgumentParser(add_help=False)
    pp_batch.add_argument('--batch', metavar='INT', default=None, type=int, required=False, action='store',
                     help='default=%(default)s: One-based index of metadata table (--metadata). '
                          'If set, process only one SRA record. This function is intended for array job processing. '
                          'If multiple batch jobs run at the same time, total CPU demand scales with the number of concurrent batch jobs.')
    pp_redo = argparse.ArgumentParser(add_help=False)
    pp_redo.add_argument('--redo', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Redo the analysis even if previous output files are detected.')
    pp_threads = argparse.ArgumentParser(add_help=False)
    pp_threads.add_argument('--threads', metavar='INT|auto', default='auto', type=int_or_auto, required=False, action='store',
                     help='default=%(default)s: Total CPU cores to use per amalgkit process. '
                          'This is a global core budget (not per-worker threads). '
                          'AMALGKIT auto-splits this budget into internal workers and per-worker threads.')
    pp_internal_jobs = argparse.ArgumentParser(add_help=False)
    pp_internal_jobs.add_argument('--internal_jobs', metavar='INT|auto', default='auto', type=int_or_auto, required=False, action='store',
                     help='default=%(default)s: Advanced override for internal parallel workers. '
                          'In auto mode, worker count is derived from --threads. '
                          'When --batch is set, this is forced to 1.')
    pp_cpu_budget = argparse.ArgumentParser(add_help=False)
    pp_cpu_budget.add_argument('--internal_cpu_budget', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto, required=False, action='store',
                     help='default=%(default)s: Advanced CPU budget cap for internal auto-parallelization. '
                          'The effective core budget is min(--threads, --internal_cpu_budget). '
                          '"auto" uses os.cpu_count().')
    pp_download = argparse.ArgumentParser(add_help=False)
    pp_download.add_argument('--download_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: Shared download/cache directory. "inferred" = out_dir/downloads')
    pp_sg = argparse.ArgumentParser(add_help=False)
    pp_sg.add_argument('--sample_group', metavar='tissueA,tissueB,tissueC,...', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: "Comma separated list of sample groups. '
                          'By default, all sample_group values in metadata.tsv are passed.')
    pp_sgc = argparse.ArgumentParser(add_help=False)
    pp_sgc.add_argument('--sample_group_color', metavar='#d95f02ff,#1b9e77ff,#7570b3ff,...', default='DEFAULT', type=str, required=False, action='store',
                     help='default=%(default)s: "Comma separated list of sample groups colors. '
                          'The order should be the same as --sample_group. '
                          'The number of colors should be the same as the number of selected sample groups. '
                          'By default, all colors are automatically assigned.')

    pme_help = 'NCBI SRA metadata retrieval and curation. See `amalgkit metadata -h`'
    pme = subparsers.add_parser('metadata', help=pme_help, parents=[pp_out, pp_redo, pp_download])
    pme.add_argument('--search_string', metavar='PATH', default=None, type=str, required=True, action='store',
                     help='default=%(default)s: Entrez search string. See https://www.ncbi.nlm.nih.gov/books/NBK25499/ for details. '
                          'The search string is used to identify SRA entries that can be found at https://www.ncbi.nlm.nih.gov/sra/ using the same string. '
                          'Example: "Cephalotus follicularis"[Organism] AND "Illumina"[Platform] AND "RNA-seq"[Strategy]')
    pme.add_argument('--entrez_email', metavar='aaa@bbb.com', default='', type=str, required=False, action='store',
                     help='default=%(default)s: Your email address. See https://www.ncbi.nlm.nih.gov/books/NBK25497/')
    pme.add_argument('--resolve_names', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Whether to resolve scientific names based on NCBI Taxonomy IDs.')
    pme.set_defaults(handler=command_handlers['metadata'])

    pse_help = 'Selecting SRA entries for analysis. See `amalgkit select -h`'
    pse = subparsers.add_parser('select', help=pse_help, parents=[pp_out, pp_meta, pp_sg])
    pse.add_argument('--config_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to the config directory. "inferred" = out_dir/config')
    pse.add_argument('--min_nspots', metavar='INT', default=5000000, type=int, required=False, action='store',
                     help='default=%(default)s: Minimum number of RNA-seq reads per sample.')
    pse.add_argument('--max_sample', metavar='INT', default=99999, type=int, required=False, action='store',
                     help='default=%(default)s: Maximum number of RNA-seq data to retain for one sample group in a species.')
    pse.add_argument('--mark_missing_rank', metavar='species|genus|family|order|class|phylum|kingdom|domain|none', default='species', type=str, required=False, action='store',
                     choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'domain', 'none'],
                     help='default=%(default)s: Mark samples lacking taxid information at the specified rank as unqualified.')
    pse.add_argument('--mark_redundant_biosamples', metavar='no|yes', default='no', type=strtobool,
                     required=False, action='store',
                     help='default=%(default)s: Whether to label SRAs with the same BioSample ID as unqualified.')
    pse.set_defaults(handler=command_handlers['select'])

    pge_help = 'Retrieving fastq files. See `amalgkit getfastq -h`'
    pge = subparsers.add_parser('getfastq', help=pge_help, parents=[pp_out, pp_meta, pp_threads, pp_internal_jobs, pp_cpu_budget, pp_redo, pp_batch, pp_download])
    pge.add_argument('--entrez_email', metavar='aaa@bbb.com', default='', type=str, required=False, action='store',
                     help='default=%(default)s: Your email address. See https://www.ncbi.nlm.nih.gov/books/NBK25497/')
    pge.add_argument('--id', metavar='XXXXX0000', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: BioProject/BioSample/SRR ID. This option can be used to directly specify '
                          'an ID to start FASTQ generation without running `amalgkit metadata` beforehand.')
    pge.add_argument('--id_list', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Location of file containing a list of SRA IDs. Otherwise works like --id')
    pge.add_argument('--layout', metavar='single|paired|auto', default='auto', type=str, required=False, action='store',
                     choices=['single', 'paired', 'auto'],
                     help='default=%(default)s: Library layout of RNA-seq data to be dumped. '
                          '"auto" prioritizes paird-end libraries if both types are available.')
    pge.add_argument('--max_bp', metavar='INT', default='999,999,999,999,999', type=str, required=False, action='store',
                     help='default=%(default)s: Target sequence size (bp) to be dumped.')
    pge.add_argument('--min_read_length', metavar='INT', default=25, type=int, required=False, action='store',
                     help='default=%(default)s: Minimum read length.')
    pge.add_argument('--pfd', dest='obsolete_pfd', metavar='yes|no', default=None, type=strtobool,
                     required=False, action='store', help=argparse.SUPPRESS)
    pge.add_argument('--pfd_exe', dest='obsolete_pfd_exe', metavar='PATH', default=None, type=str,
                     required=False, action='store', help=argparse.SUPPRESS)
    pge.add_argument('--fastq_dump_exe', dest='obsolete_fastq_dump_exe', metavar='PATH', default=None, type=str,
                     required=False, action='store', help=argparse.SUPPRESS)
    pge.add_argument('--fasterq_dump_exe', '--fasterq-dump_exe', '--fasterq-dump-exe',
                     dest='fasterq_dump_exe', metavar='PATH', default='fasterq-dump', type=str,
                     required=False, action='store',
                     help='default=%(default)s: PATH to fasterq-dump executable for public SRA extraction.')
    pge.add_argument('--seqkit_exe', metavar='PATH', default='seqkit', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to seqkit executable used for FASTQ .gz output.')
    pge.add_argument('--fasterq_size_check', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Toggle fasterq-dump pre-run disk-size checks. '
                          '"yes" forwards --size-check on, "no" forwards --size-check off.')
    pge.add_argument('--fasterq_disk_limit', metavar='STR', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Optional value forwarded to fasterq-dump --disk-limit '
                          '(for example, "200G").')
    pge.add_argument('--fasterq_disk_limit_tmp', metavar='STR', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Optional value forwarded to fasterq-dump --disk-limit-tmp '
                          '(for example, "200G").')
    pge.add_argument('--prefetch_exe', dest='obsolete_prefetch_exe', metavar='PATH', default=None, type=str,
                     required=False, action='store', help=argparse.SUPPRESS)
    pge.add_argument('--fastp', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Run fastp.')
    pge.add_argument('--fastp_exe', metavar='PATH', default='fastp', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to fastp executable.')
    pge.add_argument('--fastp_option', metavar='STR', default='-j /dev/null -h /dev/null', type=str, required=False,
                     action='store',
                     help='default=%(default)s: Options to be passed to fastp. Do not include --length_required option here. '
                          'It can be specified throught --min_read_length in amalgkit. ')
    pge.add_argument('--rrna_filter', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Remove rRNA reads using MMseqs2 before final FASTQ output. '
                          'Typical cost: minutes to tens of minutes per SRA, ~2-8 GB RAM; first run also builds '
                          'the SILVA DB.')
    pge.add_argument('--filter_order', metavar='ORDER',
                     default='fastp,rrna,contam', type=str, required=False, action='store',
                     help='default=%(default)s: Order of optional filters. Use comma or ">" separators, for example '
                          '"fastp,rrna,contam" or "rrna,contam,fastp". You may list all filters or only the enabled '
                          'ones.')
    pge.add_argument('--contam_filter', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Remove contaminant reads using MMseqs2 taxonomy (unclassified reads '
                          'are retained). Typical cost with UniRef90: tens of minutes to hours per SRA, ~32-128 GB '
                          'RAM; first run also downloads/builds the DB.')
    pge.add_argument('--contam_filter_rank', metavar='species|genus|family|order|class|phylum|kingdom|superkingdom',
                     choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom', 'domain'],
                     default='superkingdom', type=str, required=False, action='store',
                     help='default=%(default)s: Taxonomic rank used for contaminant filtering against metadata taxid lineage. '
                          '"domain" is accepted as an alias for "superkingdom".')
    pge.add_argument('--contam_filter_db_name', metavar='STR', default='UniRef90', type=str, required=False, action='store',
                     help='default=%(default)s: MMseqs2 downloadable DB name used for contaminant filtering (passed to `mmseqs databases`).')
    pge.add_argument('--contam_filter_db', metavar='PATH|inferred', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: MMseqs2 taxonomy DB prefix path. "inferred" = out_dir/downloads/mmseqs_<db_name>.')
    pge.add_argument('--mmseqs_exe', metavar='PATH', default='mmseqs', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to mmseqs executable used for contaminant filtering.')
    pge.add_argument('--remove_sra', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Remove downloaded SRA files after fastq extraction.')
    pge.add_argument('--remove_tmp', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Remove temporary files.')
    pge.add_argument('--dump_print', '--pfd_print', dest='dump_print', metavar='yes|no', default='no', type=strtobool,
                     required=False, action='store',
                     help='default=%(default)s: Show sequence extraction (fasterq-dump/compression) stdout and stderr. '
                          '--pfd_print is supported as an obsolete alias.')
    pge.add_argument('--fastp_print', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Show fastp stdout and stderr.')
    pge.add_argument('--sci_name', metavar='STR', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Species name in case the BioProject covers multiple species. Example: "Homo sapiens"')
    pge.add_argument('--ncbi', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Download SRA files using wget from NCBI cloud, if available.')
    pge.add_argument('--aws', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Download SRA files from Amazon Cloud (AWS), if available.')
    pge.add_argument('--gcp', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Download SRA files from Google Cloud (GCP), if available.')
    pge.add_argument('--gcp_project', metavar='STR', default='', type=str, required=False, action='store',
                     help='default=%(default)s: Google Cloud project for requester-pays GCP buckets (used when GCP_Link is gs://).')
    pge.add_argument('--sra_download_method', metavar='auto|urllib|curl', default='auto', type=str, required=False, action='store',
                     choices=['auto', 'urllib', 'curl'],
                     help='default=%(default)s: Method for downloading SRA objects from cloud URLs.')
    pge.add_argument('--read_name', metavar='default|trinity', default='default', type=str, required=False, action='store',
                     choices=['default', 'trinity'],
                     help='default=%(default)s: read name formatting for downstream analysis.')
    pge.add_argument('--entrez_additional_search_term', metavar='STR',
                     default=None,
                     type=str, required=False, action='store',
                     help='default=%(default)s: Entrez search terms in addition to --id option to further restrict the SRA entry.')
    pge.add_argument('--tol', metavar='FLOAT', default=1, type=float, required=False, action='store',
                     help='default=%(default)s: Acceptable percentage loss of reads relative to --max_bp. If the 1st-round sequence '
                          'generation could not produce enough reads, the 2nd-round sequence generation is activated to '
                          'compensate the loss.')
    pge.set_defaults(handler=command_handlers['getfastq'])

    pqu_help = 'Estimating transcript abundance with kallisto. See `amalgkit quant -h`'
    pqu = subparsers.add_parser('quant', help=pqu_help, parents=[pp_out, pp_meta, pp_threads, pp_internal_jobs, pp_cpu_budget, pp_redo, pp_batch])
    pqu.add_argument('--index_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to index directory. Only required if index directory is not '
                          'out_dir/index/')
    pqu.add_argument('--clean_fastq', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Remove getfastq-processed fastq files when quant is successfully completed.')
    pqu.add_argument('--fasta_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: "inferred" = out_dir/fasta. '
                          'PATH to directory containing reference transcriptome fasta files required for kallisto index building (see --build_index). '
                          'In this directory, file names of fasta files are expected to start with the string '
                          'in the "scientific_name" column of the metadata table, with a space replaced with an underbar. '
                          'Example: Arabidopsis_thaliana_v1.fasta for Arabidopsis thaliana.')
    pqu.add_argument('--build_index', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Allows AMALGKIT to build kallisto index from reference fasta files. '
                          'Will only do this for a species if an index file is not already present. One fasta file per species should be put in --fasta_dir. '
                          'AMALGKIT will read the species from the metadata, try to find the fasta file (.fa or .fasta) and build the index for further use.')
    pqu.add_argument('--index_lock_poll', metavar='INT', default=5, type=int, required=False, action='store',
                     help='default=%(default)s: Poll interval in seconds while waiting for another batch process to finish building a species index.')
    pqu.add_argument('--index_lock_timeout', metavar='INT', default=3600, type=int, required=False, action='store',
                     help='default=%(default)s: Maximum wait time in seconds for index build lock release before aborting.')
    pqu.set_defaults(handler=command_handlers['quant'])

    pmg_help = 'Generating transcript abundance tables. See `amalgkit merge -h`'
    pmg = subparsers.add_parser('merge', help=pmg_help, parents=[pp_out, pp_meta, pp_threads, pp_internal_jobs, pp_cpu_budget])
    pmg.set_defaults(handler=command_handlers['merge'])

    pbu_help = 'Generating BUSCO tables for amalgkit cstmm/csfilter. See `amalgkit busco -h`'
    pbu = subparsers.add_parser('busco', help=pbu_help, parents=[pp_out, pp_meta, pp_threads, pp_internal_jobs, pp_cpu_budget, pp_redo, pp_download])
    pbu.add_argument('--tool', metavar='auto|busco|compleasm', default='auto', type=str, required=False, action='store',
                     choices=['auto', 'busco', 'compleasm'],
                     help='default=%(default)s: Tool for BUSCO table generation.')
    pbu.add_argument('--lineage', metavar='STR', default=None, type=str, required=True, action='store',
                     help='Lineage dataset name (e.g., eukaryota_odb12).')
    pbu.add_argument('--fasta_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: "inferred" = out_dir/fasta. Directory with transcriptome fasta files.')
    pbu.add_argument('--fasta', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='Optional: single fasta path to process (requires --species).')
    pbu.add_argument('--species', metavar='STR', default=None, type=str, required=False, action='store',
                     help='Optional: species name used with --fasta (e.g., "Homo sapiens").')
    pbu.add_argument('--busco_exe', metavar='PATH', default='busco', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to busco executable.')
    pbu.add_argument('--compleasm_exe', metavar='PATH', default='compleasm', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to compleasm executable.')
    pbu.add_argument('--tool_args', metavar='STR', default=None, type=str, required=False, action='store',
                     help='Additional arguments passed to the selected tool.')
    pbu.set_defaults(handler=command_handlers['busco'])

    pcs_help = 'Applying cross-species TMM normalization using single-copy genes. See `amalgkit cstmm -h`'
    pcs = subparsers.add_parser('cstmm', help=pcs_help, parents=[pp_out, pp_meta, pp_redo])
    pcs.add_argument('--orthogroup_table', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to orthogroup table, which is, for example, Orthogroups.tsv and N0.tsv in OrthoFinder.'
                          'Specify `--orthogroup_table ""` for single-species TMM normalization.')
    pcs.add_argument('--dir_busco', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to the directory where per-species BUSCO full tables are stored. '
                          'File names in this directory are expected to be GENUS_SPECIES_busco.tsv: e.g., Arabidopsis_thaliana_busco.tsv')
    pcs.add_argument('--dir_count', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: AMALGKIT subfolder PATH to per-species transcript abundance data as produced by `amalgkit merge`. '
                          '"inferred" = out_dir/merge')
    pcs.set_defaults(handler=command_handlers['cstmm'])

    pws_help = 'Within-species outlier filtering. Outputs metadata.tsv + excluded.tsv + species PDFs (no plots/). See `amalgkit wsfilter -h`'
    pws = subparsers.add_parser('wsfilter', help=pws_help, parents=[pp_out, pp_meta, pp_batch, pp_threads, pp_internal_jobs, pp_cpu_budget, pp_sg, pp_sgc, pp_redo])
    pws.add_argument('--input_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to `amalgkit merge` or `amalgkit cstmm` output folder. '
                          '"inferred" = out_dir/cstmm if exist, else out_dir/merge.')
    pws.add_argument('--dist_method', metavar='STR', default='pearson', type=str, required=False, action='store',
                     help='default=%(default)s: Method for calculating distance.')
    pws.add_argument('--mapping_rate', metavar='FLOAT', default=0.20, type=float, required=False, action='store',
                     help='default=%(default)s: Cutoff for mapping rate.')
    pws.add_argument('--correlation_threshold', metavar='FLOAT', default=0.30, type=float, required=False, action='store',
                     help='default=%(default)s: Lower cutoff for pearson r during outlier removal (legacy fallback).')
    pws.add_argument('--plot_intermediate', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: If yes, writes intermediate plots during filtering.')
    pws.add_argument('--one_outlier_per_iter', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: If yes, removes at most one outlier per sample_group/BioProject per iteration.')
    pws.add_argument('--norm', metavar='(logn|log2|lognp1|log2p1|none)-(fpkm|tpm|none)',
                     default='log2p1-fpkm', type=str, required=False, action='store',
                     help='default=%(default)s: Expression transformation before filtering.')
    pws.add_argument('--margin_threshold', metavar='FLOAT', default=0.0, type=float, required=False, action='store',
                     help='default=%(default)s: Margin threshold for robust-margin outlier detection.')
    pws.add_argument('--robust_z_threshold', metavar='FLOAT', default=-2.5, type=float, required=False, action='store',
                     help='default=%(default)s: Robust z-score threshold for robust-margin outlier detection.')
    pws.set_defaults(handler=command_handlers['wsfilter'])

    pcsf_help = 'Cross-species outlier filtering. Outputs metadata.tsv + excluded.tsv + PDFs (no plots/). See `amalgkit csfilter -h`'
    pcsf = subparsers.add_parser('csfilter', help=pcsf_help, parents=[pp_out, pp_meta, pp_batch, pp_threads, pp_internal_jobs, pp_cpu_budget, pp_sg, pp_sgc, pp_redo])
    pcsf.add_argument('--input_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                      help='default=%(default)s: PATH to `amalgkit merge` or `amalgkit cstmm` output folder. '
                           '"inferred" = out_dir/cstmm if exist, else out_dir/merge.')
    pcsf.add_argument('--norm', metavar='(logn|log2|lognp1|log2p1|none)-(fpkm|tpm|none)',
                      default='log2p1-fpkm', type=str, required=False, action='store',
                      help='default=%(default)s: Expression transformation used during temporary table generation.')
    pcsf.add_argument('--orthogroup_table', metavar='PATH', default=None, type=str, required=False, action='store',
                      help='default=%(default)s: PATH to orthogroup table, for example Orthogroups.tsv or N0.tsv in OrthoFinder.')
    pcsf.add_argument('--dir_busco', metavar='PATH', default=None, type=str, required=False, action='store',
                      help='default=%(default)s: PATH to per-species BUSCO full tables.')
    pcsf.add_argument('--missing_strategy', metavar='em_pca|nipals|row_mean',
                      choices=['em_pca', 'nipals', 'row_mean'],
                      default='em_pca', type=str, required=False, action='store',
                      help='default=%(default)s: Missing-value handling strategy before dimensionality reduction.')
    pcsf.add_argument('--margin_threshold', metavar='FLOAT', default=0.0, type=float, required=False, action='store',
                      help='default=%(default)s: Margin threshold for robust-margin outlier detection.')
    pcsf.add_argument('--robust_z_threshold', metavar='FLOAT', default=-2.5, type=float, required=False, action='store',
                      help='default=%(default)s: Robust z-score threshold for robust-margin outlier detection.')
    pcsf.set_defaults(handler=command_handlers['csfilter'])

    pfi_help = 'Final table export from filtered metadata. See `amalgkit finalize -h`'
    pfi = subparsers.add_parser('finalize', help=pfi_help, parents=[pp_out, pp_meta, pp_batch, pp_threads, pp_internal_jobs, pp_cpu_budget, pp_sg, pp_sgc, pp_redo])
    pfi.add_argument('--input_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to `amalgkit merge` or `amalgkit cstmm` output folder. '
                          '"inferred" = out_dir/cstmm if exist, else out_dir/merge.')
    pfi.add_argument('--norm', metavar='(logn|log2|lognp1|log2p1|none)-(fpkm|tpm|none)',
                     default='log2p1-fpkm', type=str, required=False, action='store',
                     help='default=%(default)s: Expression transformation before optional batch correction.')
    pfi.add_argument('--batch_effect_alg', metavar='(no|sva|ruvseq|combatseq)',
                     choices=['no', 'sva', 'ruvseq', 'combatseq'],
                     default='no', type=str, required=False, action='store',
                     help='default=%(default)s: Batch-effect removal algorithm for finalized output tables.')
    pfi.add_argument('--clip_negative', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Clip negative values to zero after batch effect removal.')
    pfi.add_argument('--maintain_zero', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Preserve zero values from input after batch effect removal.')
    pfi.add_argument('--ruvseq_control_genes', metavar='auto|all', choices=['auto', 'all'],
                     default='auto', type=str, required=False, action='store',
                     help='default=%(default)s: Control-gene selection strategy when --batch_effect_alg ruvseq is used.')
    pfi.add_argument('--ruvseq_k', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto, required=False, action='store',
                     help='default=%(default)s: Number of unwanted factors for RUVSeq. "auto" selects k from PCA score tradeoff.')
    pfi.add_argument('--ruvseq_k_max', metavar='INT', default=5, type=int, required=False, action='store',
                     help='default=%(default)s: Maximum k considered when --ruvseq_k auto.')
    pfi.add_argument('--ruvseq_control_top_n', metavar='INT', default=1000, type=int, required=False, action='store',
                     help='default=%(default)s: Top non-DE genes considered as control candidates when --ruvseq_control_genes auto.')
    pfi.add_argument('--ruvseq_min_controls', metavar='INT', default=100, type=int, required=False, action='store',
                     help='default=%(default)s: Minimum number of control genes required for RUVSeq auto selection.')
    pfi.add_argument('--seed', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto, required=False, action='store',
                     help='default=%(default)s: Random seed for stochastic steps (SVA/RUVSeq/t-SNE). "auto" keeps default RNG behavior.')
    pfi.add_argument('--sva_nsv', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto, required=False, action='store',
                     help='default=%(default)s: Number of surrogate variables for SVA. "auto" lets sva estimate n.sv.')
    pfi.add_argument('--sva_B', metavar='INT|auto', default='auto', type=int_or_auto, required=False, action='store',
                     help='default=%(default)s: Number of permutation iterations used by SVA. "auto" chooses B from sample size.')
    pfi.add_argument('--sva_B_auto_max', metavar='INT', default=100, type=int, required=False, action='store',
                     help='default=%(default)s: Upper bound for auto-selected SVA permutation iterations.')
    pfi.set_defaults(handler=command_handlers['finalize'])

    psa_help = 'Checking the integrity of AMALGKIT input and output files. See `amalgkit sanity -h`'
    psa = subparsers.add_parser('sanity', help=psa_help, parents=[pp_out, pp_meta])
    psa.add_argument('--index', required=False, action='store_true',
                     help='set this option if you want to check for availability of index files '
                          'based on species name in metadata file.')
    psa.add_argument('--quant', required=False, action='store_true',
                     help='set this option if you want to check for availability quant output files '
                          'based on SRA IDs in metadata file.')
    psa.add_argument('--getfastq', required=False, action='store_true',
                     help='set this option if you want to check for availability of getfastq output files '
                          'based on SRA IDs in metadata file.')
    psa.add_argument('--all', required=False, action='store_true',
                     help='setting this option runs amalgkit sanity as if --index, --quant, --getfastq were set')
    psa.add_argument('--index_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to index directory. Only required if index directory is not '
                          'out_dir/index/')
    psa.add_argument('--quant_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to quant directory. Only required if quant directory is not '
                          'out_dir/quant/')
    psa.add_argument('--getfastq_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to index directory. Only required if getfastq directory is not '
                          'out_dir/getfastq/')
    psa.add_argument('--quiet', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Suppress per-run status lines and print only summary results.')
    psa.add_argument('--verbose_runs', metavar='INT', default=20, type=int, required=False, action='store',
                     help='default=%(default)s: Maximum number of runs for per-run logging in sanity checks. '
                          'Set a negative value to always print per-run logs.')
    psa.set_defaults(handler=command_handlers['sanity'])

    pin_help = 'Appending local fastq info to a metadata table. See `amalgkit integrate -h`'
    pin = subparsers.add_parser('integrate', help=pin_help, parents=[pp_out, pp_meta, pp_threads])
    pin.add_argument('--fastq_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to input directory where fastq files are stored. '
                          'Nested subdirectories are scanned recursively. The first subdirectory below --fastq_dir '
                          'is parsed as the species scientific_name (underscores are converted to spaces), so '
                          'files like PATH/Homo_sapiens/brain/sample1.fq.gz are assigned to "Homo sapiens". '
                          'If the same FASTQ basename appears under multiple species directories, the generated run '
                          'ID is prefixed with the species directory name, for example Homo_sapiens_sample1.')
    pin.add_argument('--seqkit_exe', metavar='PATH', default='seqkit', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to seqkit executable used for FASTQ statistics scanning.')
    pin.add_argument('--getfastq_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to index directory. Only required if getfastq directory is not '
                          'out_dir/getfastq/')
    pin.add_argument('--remove_tmp', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Remove temporary files.')
    pin.add_argument('--accurate_size', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: ONLY APPLIES TO .gz COMPRESSED FASTQ FILES. If no, scans the first 1,000 reads '
                          'to estimate average read length. If yes, scans the whole file for exact statistics.')
    pin.set_defaults(handler=command_handlers['integrate'])

    pco_help = 'Creating a series of config files for the metadata search. See `amalgkit config -h`'
    pco = subparsers.add_parser('config', help=pco_help, parents=[pp_out])
    pco.add_argument('--config', metavar='base|test|plantae|vertebrate', default='base', type=str, required=False, action='store',
                     help='default=%(default)s: Name of config dataset to be exported. Options: '
                          '"base": a minimal set of .config files for the purpose of creating custom config files. '
                          '"base_all": a complete set of near-empty .config files. '
                          '"test": short animal set for testing amalgkit metadtata. '
                          '"vertebrate" preconfigured set of config files for vertebrate animal data.'
                          '"plantae": preconfigured set of config files for plant data.')
    pco.add_argument('--overwrite', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: allow to overwrite config files in out_dir/config/config_name/ .')
    pco.set_defaults(handler=command_handlers['config'])

    pda_help = 'Extracting bundled test datasets. See `amalgkit dataset -h`'
    pda = subparsers.add_parser('dataset', help=pda_help, parents=[pp_out])
    pda.add_argument('--name', metavar='yeast|...', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Name of the dataset to extract. Use --list to see available datasets.')
    pda.add_argument('--list', default=False, required=False, action='store_true',
                     help='List available datasets and exit.')
    pda.add_argument('--overwrite', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Allow overwriting existing files.')
    pda.set_defaults(handler=command_handlers['dataset'])

    parser_help = subparsers.add_parser('help', help='Printing help messages')
    parser_help.add_argument(
        'topic',
        nargs='?',
        default=None,
        choices=command_names + ['help'],
        help='Optional command name to show detailed help for.',
    )
    parser_help.set_defaults(handler=build_help_command_handler(parser))
    return parser
