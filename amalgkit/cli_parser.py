import argparse
import math

from amalgkit.cli_utils import (
    EXPRESSION_NORMALIZATION_METHODS,
    build_help_command_handler,
    int_or_auto,
    nonnegative_int_or_auto,
    positive_float_or_auto,
    strtobool,
)
from amalgkit.subprocess_utils import DEPENDENCY_PROBE_TIMEOUT_SECONDS


def positive_int(value):
    int_value = int(value)
    if int_value <= 0:
        raise argparse.ArgumentTypeError('must be > 0')
    return int_value


def single_copy_threshold(value):
    try:
        threshold = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError('must be > 0 and <= 100') from exc
    if (not math.isfinite(threshold)) or (threshold <= 0.0) or (threshold > 100.0):
        raise argparse.ArgumentTypeError('must be > 0 and <= 100')
    return threshold


def finite_float(value):
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError('must be a finite number') from exc
    if not math.isfinite(parsed):
        raise argparse.ArgumentTypeError('must be a finite number')
    return parsed


def mapping_rate_threshold(value):
    parsed = finite_float(value)
    if parsed < 0.0 or parsed > 100.0:
        raise argparse.ArgumentTypeError('must be between 0 and 100')
    return parsed


def correlation_threshold(value):
    parsed = finite_float(value)
    if parsed < -1.0 or parsed > 1.0:
        raise argparse.ArgumentTypeError('must be between -1 and 1')
    return parsed


def correlation_margin(value):
    parsed = finite_float(value)
    if parsed < -2.0 or parsed > 2.0:
        raise argparse.ArgumentTypeError('must be between -2 and 2')
    return parsed


def build_parser(command_handlers, command_names, version, prog=None):
    parser = argparse.ArgumentParser(
        description='A toolkit for cross-species transcriptome amalgamation',
        prog=prog,
    )
    parser.add_argument('--version', action='version', version='amalgkit version ' + version)
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Show a full traceback for unexpected errors. This global option must precede the command name.',
    )
    parser.add_argument(
        '--log_level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default=None,
        help='default=disabled: Emit structured JSON logs to stderr at this level. This global option must precede the command name.',
    )
    parser.add_argument(
        '--log_file',
        default=None,
        metavar='PATH',
        help='Optional JSONL log file. This global option must precede the command name.',
    )
    subparsers = parser.add_subparsers()

    pp_meta = argparse.ArgumentParser(add_help=False)
    pp_meta.add_argument('--metadata', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: "inferred" = out_dir/metadata/metadata.tsv. '
                          'PATH to a metadata table. After integrate, explicitly pass its output metadata path.')
    pp_count_meta = argparse.ArgumentParser(add_help=False)
    pp_count_meta.add_argument('--metadata', metavar='PATH', default='inferred',
                              help='default=%(default)s: "inferred" = dir_count/metadata.tsv '
                                   '(normally out_dir/merge/metadata.tsv). Explicit metadata must describe every count column; '
                                   'exclusion and is_sampled control which runs are normalized.')
    pp_filter_meta = argparse.ArgumentParser(add_help=False)
    pp_filter_meta.add_argument('--metadata', metavar='PATH', default='inferred',
                               help='default=%(default)s: Use the last successful filter recorded in filter_metadata_state.json; '
                                    'legacy workspaces prefer csfilter over wsfilter with a warning. Otherwise use '
                                    'input_dir/metadata.tsv (inferred input prefers cstmm over merge). '
                                    'CSTMM-derived metadata is required for FPKM from CSTMM counts.')
    pp_out = argparse.ArgumentParser(add_help=False)
    pp_out.add_argument('--out_dir', metavar='PATH', default='./', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to the directory where intermediate and output files are generated.')
    pp_batch = argparse.ArgumentParser(add_help=False)
    pp_batch.add_argument('--batch', metavar='INT', default=None, type=int, required=False, action='store',
                     help='default=%(default)s: One-based run index after filtering metadata to is_sampled=yes, in table order. '
                          'Process only this run. This function is intended for array job processing. '
                          'If multiple batch jobs run at the same time, total CPU demand scales with the number of concurrent batch jobs.')
    pp_species_batch = argparse.ArgumentParser(add_help=False)
    pp_species_batch.add_argument('--batch', metavar='INT', default=None, type=int,
                                 help='default=%(default)s: One-based index into sorted unique scientific_name values '
                                      'in the input metadata; select all rows of that species. Not a run index. '
                                      'Species jobs are not equivalent to joint cross-species filtering; '
                                      'do not run concurrent filter/finalize jobs into the same output directory.')
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
                          'The effective worker count is capped by the number of tasks to process. '
                          'When --batch is set, this is forced to 1.')
    pp_cpu_budget = argparse.ArgumentParser(add_help=False)
    pp_cpu_budget.add_argument('--internal_cpu_budget', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto, required=False, action='store',
                     help='default=%(default)s: Advanced CPU budget cap for internal auto-parallelization. '
                          'The effective core budget is min(--threads, --internal_cpu_budget). '
                          '"auto" uses os.cpu_count().')
    pp_download = argparse.ArgumentParser(add_help=False)
    pp_download.add_argument('--download_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: Shared download/cache directory. "inferred" = out_dir/downloads')
    pp_download.add_argument('--download_lock_dir', metavar='PATH|inferred', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: Shared filesystem lock/semaphore directory for download throttling. '
                          '"inferred" = download_dir/locks')
    pp_sg = argparse.ArgumentParser(add_help=False)
    pp_sg.add_argument('--sample_group', metavar='tissueA,tissueB,tissueC,...', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: "Comma separated list of sample groups. '
                          'Use a backslash to escape a literal comma or pipe in a group name. '
                          'By default, all sample_group values in metadata.tsv are passed.')
    pp_sgc = argparse.ArgumentParser(add_help=False)
    pp_sgc.add_argument('--sample_group_color', metavar='#d95f02ff,#1b9e77ff,#7570b3ff,...', default='DEFAULT', type=str, required=False, action='store',
                     help='default=%(default)s: "Comma separated list of sample groups colors. '
                          'The order should be the same as --sample_group. '
                          'The number of colors should be the same as the number of selected sample groups. '
                          'By default, all colors are automatically assigned.')

    pme_help = 'NCBI SRA metadata retrieval and curation. See `amalgkit metadata -h`'
    pme = subparsers.add_parser('metadata', help=pme_help, parents=[pp_out, pp_redo, pp_download, pp_threads, pp_internal_jobs, pp_cpu_budget])
    pme.add_argument('--search_string', metavar='STR', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Entrez search string for one metadata query. '
                          'Required unless --species_tsv is used. See https://www.ncbi.nlm.nih.gov/books/NBK25499/ for details. '
                          'The search string is used to identify SRA entries that can be found at https://www.ncbi.nlm.nih.gov/sra/ using the same string. '
                          'Example: "Cephalotus follicularis"[Organism] AND "Illumina"[Platform] AND "RNA-seq"[Strategy]')
    pme.add_argument('--species_tsv', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: TSV with a scientific_name column. '
                          'When provided, amalgkit metadata runs species-wise batch queries and writes outputs to out_dir/metadata_specieswise/.')
    pme.add_argument('--mode', metavar='base|title_union|title_split', default='base', type=str, required=False, action='store',
                     choices=['base', 'title_union', 'title_split'],
                     help='default=%(default)s: Query construction mode for --species_tsv batch execution.')
    pme.add_argument('--title_terms', metavar='flower,leaf,root', default='flower,leaf,root', type=str, required=False, action='store',
                     help='default=%(default)s: Comma-separated title terms used by --mode title_union/title_split when --organ_terms_tsv is not set.')
    pme.add_argument('--organ_terms_tsv', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Optional TSV with columns sample_group and title_terms. '
                          'title_terms must be semicolon-separated. Used by --species_tsv with --mode title_union/title_split.')
    pme.add_argument('--species_limit', metavar='INT', default=None, type=int, required=False, action='store',
                     help='default=%(default)s: Optional maximum number of species to process from --species_tsv.')
    pme.add_argument('--merge', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: In --species_tsv mode, merge per-query metadata files into one species-level metadata table.')
    pme.add_argument('--entrez_email', metavar='aaa@bbb.com', default='', type=str, required=False, action='store',
                     help='default=%(default)s: Your email address. See https://www.ncbi.nlm.nih.gov/books/NBK25497/')
    pme.add_argument('--resolve_names', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Whether to resolve scientific names based on NCBI Taxonomy IDs.')
    pme.add_argument('--ncbi_metadata_max_concurrency', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: Maximum concurrent NCBI Entrez metadata requests across processes when '
                          'download_lock_dir is shared. Set to 0 or "auto" to disable throttling.')
    pme.add_argument('--ncbi_metadata_timeout_seconds', metavar='INT', default=300, type=int,
                     required=False, action='store',
                     help='default=%(default)s (5 min): Socket timeout for each NCBI Entrez metadata request. '
                          '0 disables the timeout.')
    pme.set_defaults(handler=command_handlers['metadata'])

    pse_help = 'Selecting SRA entries for analysis. See `amalgkit select -h`'
    pse = subparsers.add_parser('select', help=pse_help, parents=[pp_out, pp_meta, pp_threads, pp_internal_jobs, pp_cpu_budget])
    pse.add_argument('--species_tsv', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: TSV with a scientific_name column for native batch select mode. '
                          'When set, select builds per-species workspaces from metadata_specieswise and writes '
                          'batch summaries/queues/manifests under --out_dir.')
    pse.add_argument('--metadata_specieswise_dir', metavar='PATH|inferred', default='inferred', type=str,
                     required=False, action='store',
                     help='default=%(default)s: Root directory containing <species_token>/<species_token>.metadata.tsv '
                          'for --species_tsv batch mode. "inferred" = dirname(out_dir)/metadata_specieswise.')
    pse.add_argument('--summary_tsv', metavar='PATH|inferred', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: Batch select summary output path. "inferred" = out_dir/select_summary.tsv.')
    pse.add_argument('--queue_tsv', metavar='PATH|inferred', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: Batch select queue output path. "inferred" = out_dir/select_queue.tsv.')
    pse.add_argument('--manifest_tsv', metavar='PATH|inferred', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: Batch external-manifest output path. '
                          '"inferred" = out_dir/external_manifest.tsv. '
                          'Sidecars for all_tissues_ge30, all_tissues_ge3, all_tissues_ge1, and any_tissues_ge1 are written next to this path.')
    pse.add_argument('--batch_label', metavar='STR|inferred', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: Label written to batch manifests in --species_tsv mode. '
                          '"inferred" = basename(out_dir).')
    pse.add_argument('--select_rules_tsv', metavar='PATH|inferred', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to select_rules.tsv. "inferred" = out_dir/select_rules.tsv')
    pse.add_argument('--random_seed', metavar='INT', default=0, type=int, required=False, action='store',
                     help='default=%(default)s: Non-negative random seed used for deterministic sample selection.')
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
    pge.add_argument('--treat_identical_paired_as_single', metavar='yes|no', default='no', type=strtobool,
                     required=False, action='store',
                     help='default=%(default)s: Explicitly allow distributed sampling of up to 2,000 read pairs and '
                          'conversion to single-end when at least 99%% of sampled pairs are identical. Enabling this '
                          'removes read2 after sampling; the default preserves both mates.')
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
    pge.add_argument('--fastp', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Run fastp.')
    pge.add_argument('--fastp_exe', metavar='PATH', default='fastp', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to fastp executable.')
    pge.add_argument('--fastp_option', metavar='STR', default='-j /dev/null -h /dev/null', type=str, required=False,
                     action='store',
                     help='default=%(default)s: Additional options passed to fastp. Input/output, thread, and '
                          'length_required options are managed by amalgkit and are rejected here. Use '
                          '--min_read_length to set the minimum retained read length.')
    pge.add_argument('--rrna_filter', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Remove rRNA reads using MMseqs2 before final FASTQ output. '
                          'The first run builds both the SILVA sequence DB and its reusable MMseqs search index.')
    pge.add_argument('--rrna_filter_sensitivity', metavar='FLOAT|auto', default=1.0, type=positive_float_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: MMseqs2 sensitivity (-s) for rRNA filtering. '
                          '"auto" keeps the MMseqs2 default; lower values are faster but less sensitive.')
    pge.add_argument('--rrna_filter_max_seqs', metavar='INT|auto', default=20, type=int_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: MMseqs2 --max-seqs for rRNA filtering. '
                          '"auto" keeps the MMseqs2 default; lower values are faster but less sensitive.')
    pge.add_argument('--rrna_filter_chunk_spots', metavar='INT', default=5000000, type=int,
                     required=False, action='store',
                     help='default=%(default)s: Maximum FASTQ spots per synchronized MMseqs2 rRNA query chunk. '
                          'Paired mates stay in the same chunk and are removed together.')
    pge.add_argument('--rrna_filter_memory_limit', metavar='SIZE|auto', default='32G', type=str,
                     required=False, action='store',
                     help='default=%(default)s: Positive integer byte count with an optional uppercase B/K/M/G/T '
                          'unit, forwarded to MMseqs2 --split-memory-limit for SILVA index creation and each '
                          'rRNA search (for example, 32G). "auto" keeps the MMseqs2 default.')
    pge.add_argument('--rrna_filter_jobs', metavar='1|2|auto', default=1, type=int_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: Maximum concurrently processed SRA runs while rRNA filtering is '
                          'enabled. "auto" resolves conservatively to 1. Values are limited to 1 or 2; remaining '
                          '--threads are assigned to each MMseqs search.')
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
                     help='default=%(default)s: MMseqs2 taxonomy DB prefix path. "inferred" = out_dir/downloads/mmseqs2/<db_name>_DB.')
    pge.add_argument('--contam_filter_sensitivity', metavar='FLOAT|auto', default='auto', type=positive_float_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: MMseqs2 sensitivity (-s) for contaminant filtering. '
                          '"auto" keeps the MMseqs2 default; lower values are faster but less sensitive.')
    pge.add_argument('--contam_filter_max_seqs', metavar='INT|auto', default='auto', type=int_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: MMseqs2 --max-seqs for contaminant filtering. '
                          '"auto" keeps the MMseqs2 default; lower values are faster but less sensitive.')
    pge.add_argument('--contam_filter_chunk_spots', metavar='INT', default=5000000, type=positive_int,
                     required=False, action='store',
                     help='default=%(default)s: Maximum read pairs (or single-end reads) processed per MMseqs2 '
                          'taxonomy chunk. Smaller chunks reduce peak Python memory and temporary disk usage.')
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
    pge.add_argument('--ena', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Download SRA files from ENA, if derivable from the run accession. '
                          'Tried after AWS/GCP/NCBI and can absorb load when higher-priority sources are busy.')
    pge.add_argument('--ddbj', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Download SRA files from DDBJ for DRA runs, if derivable from DRR/DRX accessions. '
                          'Tried after AWS/GCP/NCBI/ENA and can absorb load when higher-priority sources are busy.')
    pge.add_argument('--gcp_project', metavar='STR', default='', type=str, required=False, action='store',
                     help='default=%(default)s: Google Cloud project for requester-pays GCP buckets (used when GCP_Link is gs://).')
    pge.add_argument('--ncbi_download_max_concurrency', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: Maximum concurrent NCBI cloud-object downloads across processes when '
                          'download_lock_dir is shared. Set to 0 or "auto" to disable throttling. '
                          'When NCBI is at its limit, getfastq tries the next enabled source before waiting.')
    pge.add_argument('--aws_download_max_concurrency', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: Maximum concurrent AWS cloud-object downloads across processes when '
                          'download_lock_dir is shared. Set to 0 or "auto" to disable throttling. '
                          'When AWS is at its limit, getfastq tries the next enabled source before waiting.')
    pge.add_argument('--gcp_download_max_concurrency', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: Maximum concurrent GCP cloud-object downloads across processes when '
                          'download_lock_dir is shared. Set to 0 or "auto" to disable throttling. '
                          'When GCP is at its limit, getfastq tries the next enabled source before waiting.')
    pge.add_argument('--ena_download_max_concurrency', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: Maximum concurrent ENA downloads across processes when '
                          'download_lock_dir is shared. Set to 0 or "auto" to disable throttling. '
                          'When ENA is at its limit, getfastq tries the next enabled source before waiting.')
    pge.add_argument('--ddbj_download_max_concurrency', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: Maximum concurrent DDBJ downloads across processes when '
                          'download_lock_dir is shared. Set to 0 or "auto" to disable throttling. '
                          'When DDBJ is at its limit, getfastq tries the next enabled source before waiting.')
    pge.add_argument('--ncbi_metadata_max_concurrency', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto,
                     required=False, action='store',
                     help='default=%(default)s: Maximum concurrent NCBI Entrez metadata requests for --id/--id_list across '
                          'processes when download_lock_dir is shared. Set to 0 or "auto" to disable throttling.')
    pge.add_argument('--ncbi_metadata_timeout_seconds', metavar='INT', default=300, type=int,
                     required=False, action='store',
                     help='default=%(default)s (5 min): Socket timeout for each NCBI Entrez metadata request. '
                          '0 disables the timeout.')
    pge.add_argument('--sra_download_method', metavar='auto|urllib|curl', default='auto', type=str, required=False, action='store',
                     choices=['auto', 'urllib', 'curl'],
                     help='default=%(default)s: Method for downloading SRA objects from cloud URLs.')
    pge.add_argument('--sra_download_wait_timeout_seconds', metavar='INT', default=86400, type=positive_int,
                     required=False, action='store',
                     help='default=%(default)s: Maximum seconds to wait for an available provider download slot.')
    pge.add_argument('--sra_download_transfer_timeout_seconds', metavar='INT', default=21600, type=positive_int,
                     required=False, action='store',
                     help='default=%(default)s: Maximum seconds allowed for one SRA or original FASTQ transfer.')
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
    pge.add_argument('--tool_timeout_seconds', metavar='INT', default=172800, type=int,
                     required=False, action='store',
                     help='default=%(default)s (48 h): Wall-clock timeout applied to each external\n'
                          ' tool invocation. 0 disables the timeout.')
    pge.add_argument('--dependency_probe_timeout_seconds', metavar='INT', default=DEPENDENCY_PROBE_TIMEOUT_SECONDS, type=int,
                     required=False, action='store',
                     help='default=%(default)s sec: Wall-clock timeout applied to each dependency\n'
                          ' version/help probe run at startup. 0 disables the timeout.')
    pge.set_defaults(handler=command_handlers['getfastq'])

    pqu_help = 'Estimating transcript abundance with auto-selected kallisto/oarfish backend. See `amalgkit quant -h`'
    pqu = subparsers.add_parser('quant', help=pqu_help, parents=[pp_out, pp_meta, pp_threads, pp_internal_jobs, pp_cpu_budget, pp_redo, pp_batch])
    pqu.add_argument('--quant_backend', metavar='auto|kallisto|oarfish', default='auto', type=str, required=False, action='store',
                     choices=['auto', 'kallisto', 'oarfish'],
                     help='default=%(default)s: Quantification backend. "auto" uses metadata to choose kallisto for short-read runs and oarfish for long-read runs.')
    pqu.add_argument('--oarfish_seq_tech', metavar='auto|ont-cdna|ont-drna|pac-bio|pac-bio-hifi', default='auto', type=str, required=False, action='store',
                     choices=['auto', 'ont-cdna', 'ont-drna', 'pac-bio', 'pac-bio-hifi'],
                     help='default=%(default)s: Override oarfish sequencing-technology preset. "auto" infers ONT/PacBio subtype from metadata where possible.')
    pqu.add_argument('--kallisto_options', metavar='STR', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Additional shell-style option string passed through to `kallisto quant`. Example: --kallisto_options "--bias --seed 42".')
    pqu.add_argument('--oarfish_options', metavar='STR', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Additional shell-style option string passed through to `oarfish`. Example: --oarfish_options "--filter-group no-filters --model-coverage".')
    pqu.add_argument('--index_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to index directory. Only required if index directory is not '
                          'out_dir/index/')
    pqu.add_argument('--clean_fastq', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Remove getfastq-processed fastq files when quant is successfully completed.')
    pqu.add_argument('--fasta_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: "inferred" = out_dir/fasta. '
                          'PATH to directory containing reference transcriptome fasta files used when building backend-specific quant indices '
                          '(see --build_index; kallisto .idx or oarfish .mmi). '
                          'In this directory, file names of fasta files are expected to start with the string '
                          'in the "scientific_name" column of the metadata table, with a space replaced with an underbar. '
                          'Accepted suffixes include .fa, .fasta, .fa.gz, and .fasta.gz. '
                          'Example: Arabidopsis_thaliana_v1.fasta for Arabidopsis thaliana.')
    pqu.add_argument('--build_index', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Allows AMALGKIT to build quant backend indices from reference fasta files '
                          '(kallisto .idx or oarfish .mmi, depending on the selected backend). '
                          'Will only do this for a species/backend target if an index file is not already present. '
                          'One fasta file per species should be put in --fasta_dir. '
                          'AMALGKIT will read species names from the metadata, try to find the matching fasta file, and build the index for further use.')
    pqu.add_argument('--index_lock_poll', metavar='INT', default=5, type=int, required=False, action='store',
                     help='default=%(default)s: Poll interval in seconds while waiting for another batch process to finish building a species index.')
    pqu.add_argument('--index_lock_timeout', metavar='INT', default=3600, type=int, required=False, action='store',
                     help='default=%(default)s: Maximum wait time in seconds for index build lock release before aborting.')
    pqu.add_argument('--tool_timeout_seconds', metavar='INT', default=43200, type=int,
                     required=False, action='store',
                     help='default=%(default)s (12 h): Wall-clock timeout applied to each external\n'
                          ' tool invocation. 0 disables the timeout.')
    pqu.add_argument('--dependency_probe_timeout_seconds', metavar='INT', default=DEPENDENCY_PROBE_TIMEOUT_SECONDS, type=int,
                     required=False, action='store',
                     help='default=%(default)s sec: Wall-clock timeout applied to each dependency\n'
                          ' version/help probe run at startup. 0 disables the timeout.')
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
                     help='Additional arguments passed to the selected tool. For BUSCO, --download is managed from --lineage and cannot be supplied here.')
    pbu.add_argument('--tool_timeout_seconds', metavar='INT', default=172800, type=int,
                     required=False, action='store',
                     help='default=%(default)s (48 h): Wall-clock timeout applied to each external\n'
                          ' tool invocation. 0 disables the timeout.')
    pbu.set_defaults(handler=command_handlers['busco'])

    pcs_help = 'Applying cross-species TMM normalization using single-copy genes. See `amalgkit cstmm -h`'
    pcs = subparsers.add_parser('cstmm', help=pcs_help, parents=[pp_out, pp_count_meta, pp_redo])
    pcs.add_argument('--orthogroup_table', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to an orthogroup table, such as Orthogroups.tsv or N0.tsv from OrthoFinder. '
                          'Not required when the count input contains only one species.')
    pcs.add_argument('--dir_busco', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to the directory where per-species BUSCO full tables are stored. '
                          'File names in this directory are expected to be GENUS_SPECIES_busco.tsv: e.g., Arabidopsis_thaliana_busco.tsv')
    pcs.add_argument('--dir_count', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: AMALGKIT subfolder PATH to per-species transcript abundance data as produced by `amalgkit merge`. '
                          '"inferred" = out_dir/merge')
    pcs.add_argument('--tmm_backend', metavar='python', default='python', choices=['python'], type=str, required=False, action='store',
                     help='default=%(default)s: Backend used for cstmm.')
    pcs.add_argument('--single_copy_threshold', metavar='PERCENT', default=50.0, type=single_copy_threshold,
                     required=False, action='store',
                     help='default=%(default)s: Minimum percentage of species in which an orthogroup must be single-copy.')
    pcs.set_defaults(handler=command_handlers['cstmm'])

    pws_help = 'Within-species outlier filtering. Outputs metadata.tsv + excluded.tsv + species PDFs (no plots/). See `amalgkit wsfilter -h`'
    pws = subparsers.add_parser('wsfilter', help=pws_help, parents=[pp_out, pp_filter_meta, pp_species_batch, pp_threads, pp_internal_jobs, pp_cpu_budget, pp_sg, pp_sgc, pp_redo])
    pws.add_argument('--input_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to `amalgkit merge` or `amalgkit cstmm` output folder. '
                          '"inferred" = out_dir/cstmm if exist, else out_dir/merge.')
    pws.add_argument('--dist_method', metavar='pearson|spearman|kendall', default='pearson', type=str,
                     choices=['pearson', 'spearman', 'kendall'], required=False, action='store',
                     help='default=%(default)s: Method for calculating distance.')
    pws.add_argument('--mapping_rate', metavar='FLOAT', default=0.20, type=mapping_rate_threshold, required=False, action='store',
                     help='default=%(default)s: Cutoff for mapping rate.')
    pws.add_argument('--correlation_threshold', metavar='FLOAT', default=0.30, type=correlation_threshold, required=False, action='store',
                     help='default=%(default)s: Lower cutoff for pearson r during outlier removal (legacy fallback).')
    pws.add_argument('--plot_intermediate', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: If yes, writes intermediate plots during filtering.')
    pws.add_argument('--one_outlier_per_iter', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: If yes, removes at most one outlier per sample_group/BioProject per iteration.')
    pws.add_argument('--norm', metavar='(logn|log2|lognp1|log2p1|none)-(fpkm|tpm|none)',
                     default='log2p1-fpkm', choices=EXPRESSION_NORMALIZATION_METHODS,
                     type=str, required=False, action='store',
                     help='default=%(default)s: Expression transformation before filtering.')
    pws.add_argument('--margin_threshold', metavar='FLOAT', default=0.0, type=correlation_margin, required=False, action='store',
                     help='default=%(default)s: Margin threshold for robust-margin outlier detection.')
    pws.add_argument('--robust_z_threshold', metavar='FLOAT', default=-2.5, type=finite_float, required=False, action='store',
                     help='default=%(default)s: Robust z-score threshold for robust-margin outlier detection.')
    pws.add_argument('--small_group_policy', metavar='margin_fallback|retain',
                     choices=['margin_fallback', 'retain'],
                     default='margin_fallback', type=str, required=False, action='store',
                     help='default=%(default)s: How to screen sample_groups too small for a robust z-score. '
                          '"margin_fallback" applies --margin_threshold alone; '
                          '"retain" keeps them without margin-based screening.')
    pws.set_defaults(handler=command_handlers['wsfilter'])

    pcsf_help = 'Cross-species outlier filtering. Outputs metadata.tsv + excluded.tsv + PDFs (no plots/). See `amalgkit csfilter -h`'
    pcsf = subparsers.add_parser('csfilter', help=pcsf_help, parents=[pp_out, pp_filter_meta, pp_species_batch, pp_threads, pp_internal_jobs, pp_cpu_budget, pp_sg, pp_sgc, pp_redo])
    pcsf.add_argument('--input_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                      help='default=%(default)s: PATH to `amalgkit merge` or `amalgkit cstmm` output folder. '
                           '"inferred" = out_dir/cstmm if exist, else out_dir/merge.')
    pcsf.add_argument('--norm', metavar='(logn|log2|lognp1|log2p1|none)-(fpkm|tpm|none)',
                      default='log2p1-fpkm', choices=EXPRESSION_NORMALIZATION_METHODS,
                      type=str, required=False, action='store',
                      help='default=%(default)s: Expression transformation used during temporary table generation.')
    pcsf.add_argument('--orthogroup_table', metavar='PATH', default=None, type=str, required=False, action='store',
                      help='default=%(default)s: PATH to orthogroup table, for example Orthogroups.tsv or N0.tsv in OrthoFinder.')
    pcsf.add_argument('--dir_busco', metavar='PATH', default=None, type=str, required=False, action='store',
                      help='default=%(default)s: PATH to per-species BUSCO full tables.')
    pcsf.add_argument('--missing_strategy', metavar='em_pca|nipals|row_mean',
                      choices=['em_pca', 'nipals', 'row_mean'],
                      default='em_pca', type=str, required=False, action='store',
                      help='default=%(default)s: Missing-value handling strategy before dimensionality reduction.')
    pcsf.add_argument('--margin_threshold', metavar='FLOAT', default=0.0, type=correlation_margin, required=False, action='store',
                      help='default=%(default)s: Margin threshold for robust-margin outlier detection.')
    pcsf.add_argument('--robust_z_threshold', metavar='FLOAT', default=-2.5, type=finite_float, required=False, action='store',
                      help='default=%(default)s: Robust z-score threshold for robust-margin outlier detection.')
    pcsf.add_argument('--small_group_policy', metavar='margin_fallback|retain',
                      choices=['margin_fallback', 'retain'],
                      default='margin_fallback', type=str, required=False, action='store',
                      help='default=%(default)s: How to screen sample_groups too small for a robust z-score. '
                           '"margin_fallback" applies --margin_threshold alone; '
                           '"retain" keeps them without margin-based screening.')
    pcsf.add_argument('--single_copy_threshold', metavar='PERCENT', default=None, type=single_copy_threshold,
                      required=False, action='store',
                      help='default=inferred: Minimum percentage of species in which an orthogroup must be single-copy. '
                           'Inherits the value recorded by cstmm, or uses 50 when no cstmm value is available.')
    pcsf.set_defaults(handler=command_handlers['csfilter'])

    pfi_help = 'Final table export from filtered metadata. See `amalgkit finalize -h`'
    pfi = subparsers.add_parser('finalize', help=pfi_help, parents=[pp_out, pp_filter_meta, pp_species_batch, pp_threads, pp_internal_jobs, pp_cpu_budget, pp_sg, pp_sgc, pp_redo])
    pfi.add_argument('--input_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to `amalgkit merge` or `amalgkit cstmm` output folder. '
                          '"inferred" = out_dir/cstmm if exist, else out_dir/merge.')
    pfi.add_argument('--norm', metavar='(logn|log2|lognp1|log2p1|none)-(fpkm|tpm|none)',
                     default='log2p1-fpkm', choices=EXPRESSION_NORMALIZATION_METHODS,
                     type=str, required=False, action='store',
                     help='default=%(default)s: Expression transformation before optional batch correction.')
    pfi.add_argument('--batch_effect_alg', metavar='(no|sva|ruvseq|combatseq|latent_glm)',
                     choices=['no', 'sva', 'ruvseq', 'combatseq', 'latent_glm'],
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
    pfi.add_argument('--seed', metavar='INT|auto', default=0, type=nonnegative_int_or_auto, required=False, action='store',
                     help='default=%(default)s: Random seed for stochastic steps (SVA/RUVSeq/t-SNE). '
                          '"auto" uses OS entropy and may produce non-reproducible results.')
    pfi.add_argument('--sva_nsv', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto, required=False, action='store',
                     help='default=%(default)s: Number of surrogate variables for SVA. "auto" lets sva estimate n.sv.')
    pfi.add_argument('--sva_B', metavar='INT|auto', default='auto', type=int_or_auto, required=False, action='store',
                     help='default=%(default)s: Number of permutation iterations used by SVA. "auto" chooses B from sample size.')
    pfi.add_argument('--sva_B_auto_max', metavar='INT', default=100, type=int, required=False, action='store',
                     help='default=%(default)s: Upper bound for auto-selected SVA permutation iterations.')
    pfi.add_argument('--sva_backend', metavar='python', choices=['python'],
                     default='python', type=str, required=False, action='store',
                     help='default=%(default)s: Backend used for --batch_effect_alg sva.')
    pfi.add_argument('--combatseq_backend', metavar='python', choices=['python'],
                     default='python', type=str, required=False, action='store',
                     help='default=%(default)s: Backend used for --batch_effect_alg combatseq.')
    pfi.add_argument('--ruvseq_backend', metavar='python', choices=['python'],
                     default='python', type=str, required=False, action='store',
                     help='default=%(default)s: Backend used for --batch_effect_alg ruvseq.')
    pfi.add_argument('--latent_family', metavar='poisson|nb', choices=['poisson', 'nb'],
                     default='nb', type=str, required=False, action='store',
                     help='default=%(default)s: Likelihood family for --batch_effect_alg latent_glm.')
    pfi.add_argument('--latent_k', metavar='INT|auto', default='auto', type=nonnegative_int_or_auto, required=False, action='store',
                     help='default=%(default)s: Number of latent factors for latent_glm. "auto" selects k from model diagnostics.')
    pfi.add_argument('--latent_k_max', metavar='INT', default=5, type=int, required=False, action='store',
                     help='default=%(default)s: Maximum k considered when --latent_k auto.')
    pfi.add_argument('--latent_max_iter', metavar='INT', default=200, type=int, required=False, action='store',
                     help='default=%(default)s: Maximum number of optimization iterations for latent_glm.')
    pfi.add_argument('--latent_tol', metavar='FLOAT', default=1e-5, type=float, required=False, action='store',
                     help='default=%(default)s: Convergence tolerance for latent_glm.')
    pfi.set_defaults(handler=command_handlers['finalize'])

    psa_help = 'Checking the integrity of AMALGKIT input and output files. See `amalgkit sanity -h`'
    psa = subparsers.add_parser('sanity', help=psa_help, parents=[pp_out, pp_meta, pp_threads])
    psa.add_argument('--index', required=False, action='store_true',
                     help='set this option if you want to check for availability of index files '
                          'based on species name in metadata file.')
    psa.add_argument('--quant', required=False, action='store_true',
                     help='set this option if you want to check for availability quant output files '
                          'based on SRA IDs in metadata file.')
    psa.add_argument('--getfastq', required=False, action='store_true',
                     help='set this option if you want to check for availability of getfastq output files '
                          'based on SRA IDs in metadata file.')
    psa.add_argument('--merge', required=False, action='store_true',
                     help='set this option if you want to check merge outputs based on species in metadata file.')
    psa.add_argument('--busco', required=False, action='store_true',
                     help='set this option if you want to check busco outputs based on species in metadata file.')
    psa.add_argument('--finalize', required=False, action='store_true',
                     help='set this option if you want to check finalize outputs based on species in metadata file.')
    psa.add_argument('--all', required=False, action='store_true',
                     help='setting this option runs amalgkit sanity as if --index, --quant, --getfastq, --merge, '
                          '--busco, --finalize were set. '
                          'If none of these target flags are specified, --all is assumed.')
    psa.add_argument('--check', metavar='LIST', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Comma-separated sanity checks to run. Accepted values: '
                          'getfastq,index,quant,merge,busco,finalize,all.')
    psa.add_argument('--run', metavar='RUN1,RUN2,...', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Optional comma-separated run IDs to limit sanity checks.')
    psa.add_argument('--species', metavar='SPECIES1,SPECIES2,...', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Optional comma-separated scientific_name values to limit sanity checks.')
    psa.add_argument('--index_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to index directory. Only required if index directory is not '
                          'out_dir/index/')
    psa.add_argument('--quant_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to quant directory. Only required if quant directory is not '
                          'out_dir/quant/')
    psa.add_argument('--getfastq_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to index directory. Only required if getfastq directory is not '
                          'out_dir/getfastq/')
    psa.add_argument('--merge_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to merge directory. Only required if merge directory is not '
                          'out_dir/merge/')
    psa.add_argument('--busco_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to busco directory. Only required if busco directory is not '
                          'out_dir/busco/')
    psa.add_argument('--finalize_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: PATH to finalize directory. Only required if finalize directory is not '
                          'out_dir/finalize/')
    psa.add_argument('--quiet', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Suppress per-run status lines and print only summary results.')
    psa.add_argument('--verbose_runs', metavar='INT', default=20, type=int, required=False, action='store',
                     help='default=%(default)s: Maximum number of runs for per-run logging in sanity checks. '
                          'Set a negative value to always print per-run logs.')
    psa.add_argument('--strict', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Exit with code 1 if any selected sanity check reports missing or '
                          'ambiguous inputs/outputs.')
    psa.add_argument('--strict_level', metavar='none|error|warning', default='none', type=str, required=False, action='store',
                     choices=['none', 'error', 'warning'],
                     help='default=%(default)s: Failure threshold for sanity issues. '
                          '"error" fails only on errors, "warning" fails on warnings or errors.')
    psa.set_defaults(handler=command_handlers['sanity'])

    pre_help = 'Rerunning failed outputs based on amalgkit sanity report. See `amalgkit rerun -h`'
    pre = subparsers.add_parser('rerun', help=pre_help, parents=[pp_threads, pp_internal_jobs, pp_cpu_budget, pp_download])
    pre.add_argument('--out_dir', metavar='PATH', default='./', type=str, required=False, action='store',
                     help='default=%(default)s: Workspace root. Used to infer --report when --report inferred.')
    pre.add_argument('--report', metavar='PATH|inferred', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: PATH to sanity_report.json. "inferred" = out_dir/sanity/sanity_report.json')
    pre.add_argument('--metadata', metavar='PATH|report|inferred', default='report', type=str, required=False, action='store',
                     help='default=%(default)s: Metadata table to use for rerun target selection. '
                          '"report" uses metadata_path stored in sanity_report.json. '
                          '"inferred" = out_dir/metadata/metadata.tsv.')
    pre.add_argument('--check', metavar='LIST', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Comma-separated rerun targets. Accepted values: '
                          'getfastq,index,quant,merge,busco,finalize,all. '
                          'If omitted, rerun uses the checks recorded in sanity_report.json.')
    pre.add_argument('--run', metavar='RUN1,RUN2,...', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Optional comma-separated run IDs to further limit run-based reruns.')
    pre.add_argument('--species', metavar='SPECIES1,SPECIES2,...', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Optional comma-separated scientific_name values to further limit species-based reruns.')
    pre.add_argument('--redo', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Overwrite broken outputs when rerunning targets.')
    pre.add_argument('--dry_run', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Print resolved rerun targets without executing commands.')
    pre.add_argument('--manifest', metavar='PATH|inferred|none', default='inferred', type=str, required=False, action='store',
                     help='default=%(default)s: JSON path for the resolved rerun execution plan. '
                          '"inferred" = out_dir/sanity/rerun_manifest.json. '
                          '"none" disables manifest output.')
    pre.add_argument('--include_warnings', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Include warning-level sanity issues in rerun target selection. '
                          'Useful for regenerating global summary outputs such as BUSCO/merge/finalize PDFs.')
    pre.set_defaults(handler=command_handlers['rerun'])

    pin_help = 'Appending local fastq info to a metadata table. See `amalgkit integrate -h`'
    pin = subparsers.add_parser('integrate', help=pin_help, parents=[pp_out, pp_meta, pp_threads, pp_download])
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
    pin.add_argument('--output_metadata', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Explicit output path for the generated metadata table. '
                          'Without this option, standalone integrate writes out_dir/metadata_private_fastq.tsv and '
                          'merge mode writes out_dir/metadata/metadata_updated_for_private_fastq.tsv.')
    pin.add_argument('--remove_tmp', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Remove temporary files.')
    pin.add_argument('--accurate_size', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: ONLY APPLIES TO .gz COMPRESSED FASTQ FILES. If no, scans the first 1,000 reads '
                          'to estimate average read length. If yes, scans the whole file for exact statistics.')
    pin.add_argument('--tool_timeout_seconds', metavar='INT', default=7200, type=int,
                     required=False, action='store',
                     help='default=%(default)s (2 h): Wall-clock timeout applied to each external\n'
                          ' tool invocation. 0 disables the timeout.')
    pin.set_defaults(handler=command_handlers['integrate'])

    pda_help = 'Extracting bundled datasets, initializing workspaces, and exporting select rule sets. See `amalgkit dataset -h`'
    pda = subparsers.add_parser('dataset', help=pda_help, parents=[pp_out])
    pda.add_argument('--name', metavar='init|yeast|...', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Name of the dataset to extract. '
                          '`init` creates an empty workspace scaffold with template TSVs and default select_rules.tsv. '
                          'Use --list to see available datasets.')
    pda.add_argument('--rule_set', metavar='base|test|plantae|vertebrate', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Bundled select rule set to export to out_dir/select_rules.tsv. '
                          'Use --list to see available rule sets.')
    pda.add_argument('--list', default=False, required=False, action='store_true',
                     help='List available datasets and rule sets and exit.')
    pda.add_argument('--overwrite', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Allow overwriting existing files, including out_dir/select_rules.tsv.')
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
