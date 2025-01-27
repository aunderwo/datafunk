import argparse
import sys

import datafunk
import datafunk.subcommands

def main(args=None):
    parser = argparse.ArgumentParser(
        prog="datafunk",
        usage="datafunk <subcommand> <options>",
        description="Miscellaneous data manipulation tools",
    )

    parser.add_argument("--version", action="version", version=datafunk.__version__)
    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )

    # _________________________________ repair_names ____________________________#
    subparser_repair_names = subparsers.add_parser("repair_names",
        usage="datafunk repair_names --fasta <fasta> --tree <tree> --out <outfile>",
        help="Returns iqtree taxa names to their former glory",
    )

    subparser_repair_names.add_argument("--fasta", action="store", type=str, dest="fasta")
    subparser_repair_names.add_argument("--tree", action="store", type=str, dest="tree")
    subparser_repair_names.add_argument("--out", action="store", type=str, dest="out")

    subparser_repair_names.set_defaults(func=datafunk.subcommands.repair_names.run)

    # _________________________________ remove_fasta ____________________________#
    subparser_remove_fasta = subparsers.add_parser(
        "remove_fasta",
        usage="datafunk remove_fasta -i <input_fasta_file> -f <filter_file> -o <output_file>",
        help="Removing sequences based on a filtering file",
    )

    subparser_remove_fasta.add_argument(
        "-i",
        "--input-fasta",
        dest="input_fasta",
        action="store",
        type=str,
        help="Input file: something about the input file format",
    )

    subparser_remove_fasta.add_argument(
        "-f",
        "--filter-file",
        dest="filter_file",
        action="store",
        type=str,
        help="Filter file for filtering based on filter file",
    )

    subparser_remove_fasta.add_argument(
        "-o",
        "--output-fasta",
        dest="output_file",
        action="store",
        type=str,
        default='removed.fasta',
        help="Output file name for resulting filtered fasta file",
    )

    subparser_remove_fasta.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_remove_fasta.set_defaults(func=datafunk.subcommands.remove_fasta.run)

    # _________________________________ clean_names ____________________________#
    subparser_clean_names = subparsers.add_parser(
        "clean_names",
        usage="datafunk clean_names -i <input_metadata> -t <trait> -o <output_file>",
        help="Standardizing country names based on trait given",
    )

    subparser_clean_names.add_argument(
        "-i",
        "--input-metadata",
        dest="input_metadata",
        action="store",
        type=str,
        help="Input file: metadata (csv) for location curation",
    )

    subparser_clean_names.add_argument(
        "-t",
        "--trait",
        dest="trait",
        action="store",
        type=str,
        help="Column name containing the locations",
    )

    subparser_clean_names.add_argument(
        "-o",
        "--output-metadata",
        dest="output_metadata",
        action="store",
        type=str,
        default='cleaned.csv',
        help="Output file name for resulting filtered metadata",
    )

    subparser_clean_names.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_clean_names.set_defaults(func=datafunk.subcommands.clean_names.run)

    # _________________________________ merge_fasta ____________________________#
    subparser_merge_fasta = subparsers.add_parser(
        "merge_fasta",
        usage="datafunk merge_fasta -f <folder> -i <input_metadata> -o <output_file>",
        help="Merge fasta file into one single file with removal of duplicates",
    )

    subparser_merge_fasta.add_argument(
        "-f",
        "--folder",
        dest="folder",
        action="store",
        type=str,
        help="Folder containing all fasta files needed to be merged",
    )

    subparser_merge_fasta.add_argument(
        "-i",
        "--input-metadata",
        dest="input_metadata",
        action="store",
        type=str,
        help="Input metadata (csv) for validating sequence information",
    )

    subparser_merge_fasta.add_argument(
        "-o",
        "--output-fasta",
        dest="output_fasta",
        action="store",
        type=str,
        default='filtered.fasta',
        help="Output for merged fasta file",
    )

    subparser_merge_fasta.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_merge_fasta.set_defaults(func=datafunk.subcommands.merge_fasta.run)

    # _________________________________ filter_fasta_by_covg_and_length ____________________________#
    subparser_filter_fasta_by_covg_and_length = subparsers.add_parser(
        "filter_fasta_by_covg_and_length",
        usage="datafunk filter_fasta_by_covg_and_length -i <input_fasta> -t <threshold> [-o <output_fasta>] [-f <failed_qc_output_fasta>]",
        help="Removes sequences where the fraction of non-N bases falls below the threshold",
    )

    subparser_filter_fasta_by_covg_and_length.add_argument(
        "-i",
        "--input-fasta",
        dest="input_fasta",
        action="store",
        required=True,
        type=str,
        help="Input FASTA",
    )

    subparser_filter_fasta_by_covg_and_length.add_argument(
        "--min-covg",
        dest="min_covg",
        action="store",
        required=False,
        type=int,
        help="Integer representing the minimum coverage percentage threshold. Sequences with coverage "
             "(strictly) less than this will be excluded from the filtered file.",
    )

    subparser_filter_fasta_by_covg_and_length.add_argument(
        "--min-length",
        dest="min_length",
        action="store",
        required=False,
        type=int,
        help="Integer representing the minimum length threshold. Sequences with length (strictly) "
             "less than this will be excluded from the filtered file. Default: ?",
    )

    subparser_filter_fasta_by_covg_and_length.add_argument(
        "-o",
        "--output-fasta",
        dest="output_fasta",
        action="store",
        default=None,
        type=str,
        help="Output file name for resulting filtered FASTA (default adds .filtered to input file name)",
    )
    
    subparser_filter_fasta_by_covg_and_length.add_argument(
        "-f",
        "--failed_qc_output_fasta",
        dest="failed_qc_output_fasta",
        action="store",
        default=None,
        type=str,
        help="Output file name for fasta files that do not pass filtering (default do not write)",
    )

    subparser_filter_fasta_by_covg_and_length.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_filter_fasta_by_covg_and_length.set_defaults(func=datafunk.subcommands.filter_fasta_by_covg_and_length.run)

    # _________________________________ process_gisaid_sequence_data ____________________________#
    subparser_process_gisaid_sequence_data = subparsers.add_parser(
        "process_gisaid_sequence_data",
        usage="datafunk process_gisaid_sequence_data -i <input.json OR input.fasta> [-o <output.fasta>] [-e file1 -e file2 ...] [--stdout]",
        help="Process raw sequence data in fasta or json format",
        description="Process raw sequence data in fasta or json format",
    )

    subparser_process_gisaid_sequence_data._action_groups.pop()
    required_process_gisaid_sequence_data = subparser_process_gisaid_sequence_data.add_argument_group('required arguments')
    optional_process_gisaid_sequence_data = subparser_process_gisaid_sequence_data.add_argument_group('optional arguments')

    required_process_gisaid_sequence_data.add_argument(
        '-i',
        '--input',
        required=True,
        metavar='GISAID.fasta OR GISAID.json',
        help='Sequence data in FASTA/json format'
    )

    optional_process_gisaid_sequence_data.add_argument(
        '-o',
        '--output-fasta',
        required=False,
        metavar = 'OUTPUT.fasta',
        help='FASTA format file to write, print to stdout if unspecified'
    )

    optional_process_gisaid_sequence_data.add_argument(
        '-e',
        '--exclude',
        action='append',
        required=False,
        metavar = 'FILE',
        help='A file that contains (anywhere) EPI_ISL_###### IDs to exclude (can provide many files, '
             'e.g. -e FILE1 -e FILE2 ...)'
    )

    optional_process_gisaid_sequence_data.add_argument(
        '--exclude-uk',
        required=False,
        action='store_true',
        help='Removes all GISAID entries with containing England, Ireland, Scotland or Wales',
    )

    optional_process_gisaid_sequence_data.add_argument(
        '--exclude-undated',
        required=False,
        action='store_true',
        help='Removes all GISAID entries with an incomplete date',
    )

    subparser_process_gisaid_sequence_data.set_defaults(func=datafunk.subcommands.process_gisaid_sequence_data.run)

    # _________________________________ sam_2_fasta _____________________________#
    subparser_sam_2_fasta = subparsers.add_parser(
        "sam_2_fasta",
        usage="datafunk sam_2_fasta -s <input.sam> -r <reference.fasta> [-o <output.fasta>] [-t [INT]:[INT]] [--prefix-ref] [--stdout]",
        help="Convert sam format alignment to fasta format multiple alignment, with optional trimming",
        description="aligned sam -> fasta (with optional trim to user-defined (reference) co-ordinates)",
    )

    subparser_sam_2_fasta._action_groups.pop()
    required_sam_2_fasta = subparser_sam_2_fasta.add_argument_group('required arguments')
    optional_sam_2_fasta = subparser_sam_2_fasta.add_argument_group('optional arguments')

    required_sam_2_fasta.add_argument(
        '-s', '--sam',
        help='samfile',
        required=True,
        metavar='in.sam'
                        )
    required_sam_2_fasta.add_argument(
        '-r', '--reference',
        help='reference',
        required=True,
        metavar='reference.fasta'
                        )
    optional_sam_2_fasta.add_argument(
        '-o', '--output-fasta',
        dest='output_fasta',
        help='FASTA format file to write. Prints to stdout if not specified',
        required=False,
        metavar='out.fasta'
                        )
    optional_sam_2_fasta.add_argument(
        '-t', '--trim',
        help='trim the alignment to these coordinates (0-based, half-open)',
        required=False,
        metavar='[[start]:[end]]'
                        )
    optional_sam_2_fasta.add_argument(
        '--pad',
        help='if --trim, pad trimmed ends with Ns, to retain reference length',
        required=False,
        action='store_true'
                        )
    optional_sam_2_fasta.add_argument(
        '--prefix-ref',
        dest='prefix_ref',
        help='write the reference sequence at the beginning of the file',
        required=False,
        action='store_true'
                        )
    optional_sam_2_fasta.add_argument(
        '--log-inserts',
        dest='log_inserts',
        help='log non-singleton insertions relative to the reference',
        required=False,
        action='store_true'
                        )
    optional_sam_2_fasta.add_argument(
        '--log-all-inserts',
        help='log all (including singleton) insertions relative to the reference',
        dest='log_all_inserts',
        required=False,
        action='store_true'
                        )
    optional_sam_2_fasta.add_argument(
        '--log-deletions',
        dest='log_dels',
        help='log non-singleton deletions relative to the reference',
        required=False,
        action='store_true'
                        )
    optional_sam_2_fasta.add_argument(
        '--log-all-deletions',
        help='log all (including singleton) deletions relative to the reference',
        dest='log_all_dels',
        required=False,
        action='store_true'
                        )
    optional_sam_2_fasta.add_argument(
        '--stdout',
        help='Overides -o/--output-fasta if present and prints output to stdout',
        required=False,
        action='store_true'
                        )

    subparser_sam_2_fasta.set_defaults(func=datafunk.subcommands.sam_2_fasta.run)

    # _________________________________ phylotype_consensus ____________________________#
    subparser_phylotype_consensus = subparsers.add_parser(
        "phylotype_consensus",
        usage="datafunk phylotype_consensus -i <input_fasta> -m <input_metadata> -c <clade_file> -o <output_folder>",
        help="Split sequences according to lineage and create an consensus",
    )

    subparser_phylotype_consensus.add_argument(
        "-i",
        "--input-fasta",
        dest="input_fasta",
        action="store",
        type=str,
        help="Fasta file for splitting into phylotypes",
    )

    subparser_phylotype_consensus.add_argument(
        "-m",
        "--input-metadata",
        dest="input_metadata",
        action="store",
        type=str,
        help="Input metadata (csv) with phylotype information",
    )


    subparser_phylotype_consensus.add_argument(
        "-c",
        "--clade-file",
        dest="clade_file",
        action="store",
        type=str,
        help="Clade file stating the phylotypes needed to be grouped",
    )

    subparser_phylotype_consensus.add_argument(
        "-o",
        "--output-folder",
        dest="output_folder",
        action="store",
        default="./",
        type=str,
        help="Output folder for the phylotype fasta files and consensus file",
    )

    subparser_phylotype_consensus.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )

    subparser_phylotype_consensus.set_defaults(func=datafunk.subcommands.phylotype_consensus.run)

    # _________________________________ gisaid_json_2_metadata ____________________________#
    subparser_gisaid_json_2_metadata = subparsers.add_parser(
        """gisaid_json_2_metadata""",
        usage="""datafunk gisaid_json_2_metadata [-h] -n gisaid.json -c <OLD_metadata.csv / False> -o NEW_metadata.csv -e omissions.txt""",
        description="""Add the info from a Gisaid json dump to an existing metadata table (or create a new one)""",
        help="""Add the info from a Gisaid json dump to an existing metadata table (or create a new one)""")

    subparser_gisaid_json_2_metadata._action_groups.pop()
    required_gisaid_json_2_metadata = subparser_gisaid_json_2_metadata.add_argument_group('required arguments')
    optional_gisaid_json_2_metadata = subparser_gisaid_json_2_metadata.add_argument_group('optional arguments')

    required_gisaid_json_2_metadata.add_argument('-n', '--new',
                        help='Most recent Gisaid json dump',
                        required=True,
                        metavar = 'gisaid.json')
    required_gisaid_json_2_metadata.add_argument('-c', '--csv',
                        help='Last metadata table (csv format), or \'False\' if you really want '
                             '(but you will lose date stamp information from previous dumps)',
                        required=True,
                        metavar = '<OLD_metadata.csv / False>')
    required_gisaid_json_2_metadata.add_argument('-o', '--output-metadata',
                        help='New csv file to write',
                        required=True,
                        metavar='NEW_metadata.csv')
    optional_gisaid_json_2_metadata.add_argument(
        '-l',
        '--lineages',
        required=False,
        metavar = 'lineages.csv',
        help='csv file of ineages to include'
    )
    required_gisaid_json_2_metadata.add_argument('-e', '--exclude',
                        action='append',
                        required=True,
                        metavar = 'omissions.txt',
                        help='A file that contains (anywhere) EPI_ISL_###### IDs to exclude (can provide more than one file, '
                             'e.g. -e FILE1 -e FILE2 ...)')

    subparser_gisaid_json_2_metadata.set_defaults(func=datafunk.subcommands.gisaid_json_2_metadata.run)

    # _________________________________ set_uniform_header ____________________________#
    subparser_set_uniform_header = subparsers.add_parser(
        "set_uniform_header",
        help="Adds a header column to metadata table and converts fasta to have a uniform header",
    )

    subparser_set_uniform_header.add_argument(
        "--input-fasta",
        dest="input_fasta",
        required=True,
        type=str,
        help="Input FASTA",
    )
    subparser_set_uniform_header.add_argument(
        "--input-metadata",
        dest="input_metadata",
        required=True,
        type=str,
        help="Input CSV or TSV",
    )
    subparser_set_uniform_header.add_argument(
        "--output-fasta",
        dest="output_fasta",
        required=True,
        type=str,
        help="Input FASTA",
    )
    subparser_set_uniform_header.add_argument(
        "--output-metadata",
        dest="output_metadata",
        required=True,
        type=str,
        help="Input CSV or TSV",
    )
    subparser_set_uniform_header.add_argument(
        "--gisaid",
        dest="gisaid",
        action="store_true",
        required=False,
        help="Input data is from GISAID",
    )
    subparser_set_uniform_header.add_argument(
        "--cog-uk",
        dest="cog_uk",
        action="store_true",
        required=False,
        help="Input data is from COG-UK",
    )
    subparser_set_uniform_header.add_argument(
        "--log",
        dest="log_file",
        required=False,
        help="Log file to use (otherwise uses stdout)"
    )
    subparser_set_uniform_header.add_argument(
        "--column-name",
        dest="column_name",
        required=False,
        default='sequence_name',
        help="Name of column in metadata corresponding to fasta header"
    )
    subparser_set_uniform_header.add_argument(
        "--index-column",
        dest="index_column",
        required=False,
        help="Name of column in metadata to parse for string matching with fasta header"
    )
    subparser_set_uniform_header.add_argument(
        "--extended",
        dest="extended",
        action="store_true",
        required=False,
        help="Longer fasta name"
    )

    subparser_set_uniform_header.set_defaults(func=datafunk.subcommands.set_uniform_header.run)

    # _________________________________ add_epi_week ____________________________#
    subparser_add_epi_week = subparsers.add_parser(
        "add_epi_week",
        usage="datafunk add_epi_week -i <input_metadata> -o <output_metadata> --date_column <column> [--epi-column-name <column>]",
        help="Adds a column for epi week to metadata table, and optionally a column for epi day",
    )

    subparser_add_epi_week.add_argument(
        "-i", "--input-metadata",
        dest="input_metadata",
        required=True,
        type=str,
        help="Input CSV or TSV",
    )
    subparser_add_epi_week.add_argument(
        "-o", "--output-metadata",
        dest="output_metadata",
        required=True,
        type=str,
        help="Input CSV or TSV",
    )
    subparser_add_epi_week.add_argument(
        "--date-column",
        dest="date_column",
        required=True,
        help="Column name to convert to epi week",
    )
    subparser_add_epi_week.add_argument(
        "--epi-week-column-name",
        dest="epi_week_column_name",
        default="edin_epi_week",
        required=False,
        help="Column name for epi week column",
    )
    subparser_add_epi_week.add_argument(
        "--epi-day-column-name",
        dest="epi_day_column_name",
        required=False,
        help="Column name for epi day column",
    )
    subparser_add_epi_week.set_defaults(func=datafunk.subcommands.add_epi_week.run)

    # _______________________________ process_gisaid_data ____________________________________________#
    subparser_process_gisaid_data = subparsers.add_parser(
        """process_gisaid_data""",
        usage="""datafunk process_gisaid_data --input-json <export.json> --input-metadata <in.csv>
                 --output-fasta <out.fa> --output-metadata <out.csv> --exclude-file <omissions.txt> --exclude-uk --exclude-undated""",
        description="""Gisaid json (+ metadata) -> (new) gisaid.fasta + metadata""",
        help="""Gisaid json (+ metadata) -> (new) gisaid.fasta + metadata""")

    subparser_process_gisaid_data._action_groups.pop()
    required_process_gisaid_data = subparser_process_gisaid_data.add_argument_group('required arguments')
    optional_process_gisaid_data = subparser_process_gisaid_data.add_argument_group('optional arguments')

    required_process_gisaid_data.add_argument('--input-json',
                        dest='json',
                        help='Gisaid json data',
                        required=True,
                        metavar = 'gisaid.json')
    required_process_gisaid_data.add_argument('--input-metadata',
                        help='previous metadata (can be \'False\')',
                        required=True,
                        dest='input_metadata',
                        metavar='metadata.in.csv')
    optional_process_gisaid_data.add_argument('--output-fasta',
                        help='fasta format file to write.',
                        dest='output_sequences',
                        required=False,
                        metavar = 'output.fasta')
    optional_process_gisaid_data.add_argument('--output-metadata',
                        help='metadata file to write.',
                        dest='output_metadata',
                        required=False,
                        metavar = 'metadata.out.csv')
    optional_process_gisaid_data.add_argument('--exclude-file',
                        action='append',
                        required=False,
                        dest='exclude',
                        metavar = 'FILE',
                        help='A file that contains (anywhere) EPI_ISL_###### IDs to exclude (can provide many files, e.g. -e FILE1 -e FILE2 ...)')
    optional_process_gisaid_data.add_argument('--exclude-uk',
                        action='store_true',
                        dest='exclude_uk',
                        required=False,
                        help='Excludes GISAID entries from England, Ireland, Scotland or Wales from being written to fasta (default is to include them)')
    optional_process_gisaid_data.add_argument('--exclude-undated',
                        action='store_true',
                        dest='exclude_undated',
                        required=False,
                        help='Excludes GISAID entries with an incomplete date from being written to fasta (default is to include them)')
    optional_process_gisaid_data.add_argument('--include-subsampled',
                        action='store_true',
                        dest='include_subsampled',
                        required=False,
                        help='Write GISAID entries previously flagged as duplicated to fasta (default is to exclude them)')
    optional_process_gisaid_data.add_argument('--include-omitted-file',
                        action='store_true',
                        dest='include_omitted_file',
                        required=False,
                        help='Write GISAID entries excluded in --exclude-file FILE to fasta (default is to exclude them)')

    subparser_process_gisaid_data.set_defaults(func=datafunk.subcommands.process_gisaid_data.run)

    # ___________________________________pad_alignment________________________________________#

    subparser_pad_alignment = subparsers.add_parser(
        """pad_alignment -i <input.fasta> -o <output.fasta> --left-pad <int> --right-pad <int> [--stdout]""",
        usage="""datafunk pad_alignment""",
        description="""pad alignment with leading and trailing Ns""",
        help="""pad alignment with leading and trailing Ns""")

    subparser_pad_alignment._action_groups.pop()
    required_pad_alignment = subparser_pad_alignment.add_argument_group('required arguments')
    optional_pad_alignment = subparser_pad_alignment.add_argument_group('optional arguments')


    required_pad_alignment.add_argument('-i', '--input-fasta',
                        help='Fasta file to pad',
                        required=True,
                        dest='fasta_in',
                        metavar='input.fasta')
    required_pad_alignment.add_argument('-l', '--left-pad',
                        help='Number of leading Ns to pad with',
                        dest='left_pad',
                        required=False,
                        metavar = 'int')
    required_pad_alignment.add_argument('-r', '--right-pad',
                        help='Number of trailing Ns to pad with',
                        dest='right_pad',
                        required=False,
                        metavar = 'int')
    optional_pad_alignment.add_argument('-o','--output-fasta',
                        help='Padded fasta file to write. Prints to stdout if not specified',
                        dest='fasta_out',
                        required=False,
                        metavar = 'output.fasta')
    optional_pad_alignment.add_argument('--stdout',
                        help='Overrides --output-fasta and writes to stdout',
                        action = 'store_true',
                        required=False)

    subparser_pad_alignment.set_defaults(func=datafunk.subcommands.pad_alignment.run)

    # ___________________________ exclude_uk_seqs _____________________________________#

    subparser_exclude_uk_seqs = subparsers.add_parser(
        """exclude_uk_seqs""",
        usage="""datafunk exclude_uk_seqs -i <input.fasta> -o <output.fasta>""",
        description="""exclude UK sequences from fasta""",
        help="""exclude UK sequences from fasta""")

    subparser_exclude_uk_seqs._action_groups.pop()
    required_exclude_uk_seqs = subparser_exclude_uk_seqs.add_argument_group('required arguments')


    required_exclude_uk_seqs.add_argument('-i', '--input-fasta',
                        help='Fasta file to read',
                        required=True,
                        dest='fasta_in',
                        metavar='input.fasta')
    required_exclude_uk_seqs.add_argument('-o','--output-fasta',
                        help='Fasta file to write',
                        dest='fasta_out',
                        required=False,
                        metavar = 'output.fasta')

    subparser_exclude_uk_seqs.set_defaults(func=datafunk.subcommands.exclude_uk_seqs.run)

    # ___________________________ get_CDS _____________________________________#

    subparser_get_CDS = subparsers.add_parser(
        """get_CDS""",
        usage="""datafunk get_CDS -i <input.fasta> -o <output.fasta> [--translate]""",
        description="""Extracts CDS from alignments in Wuhan-Hu-1 coordinates""",
        help="""Extracts CDS from alignments in Wuhan-Hu-1 coordinates""")

    subparser_get_CDS._action_groups.pop()
    required_get_CDS = subparser_get_CDS.add_argument_group('required arguments')
    optional_get_CDS = subparser_get_CDS.add_argument_group('optional arguments')


    required_get_CDS.add_argument('-i', '--input-fasta',
                        help='Fasta file to read',
                        required=True,
                        dest='fasta_in',
                        metavar='input.fasta')
    optional_get_CDS.add_argument('-o', '--output-fasta',
                        help='Fasta file to write. Prints to stdout if not specified',
                        required=False,
                        dest='fasta_out',
                        metavar='output.fasta')
    optional_get_CDS.add_argument('--translate',
                        help='output amino acid sequence (default is nucleotides)',
                        required=False,
                        action='store_true')

    subparser_get_CDS.set_defaults(func=datafunk.subcommands.get_CDS.run)

    # ___________________________ distance_to_root _____________________________________#

    subparser_distance_to_root = subparsers.add_parser(
        """distance_to_root""",
        usage="""datafunk distance_to_root --input-fasta <file> --input-metadata <file>""",
        description="""calculates per sample genetic distance to WH04 and writes it to 'distances.tsv'""",
        help="""calculates per sample genetic distance to WH04 and writes it to 'distances.tsv""")

    subparser_distance_to_root._action_groups.pop()
    required_distance_to_root = subparser_distance_to_root.add_argument_group('required arguments')
    optional_distance_to_root = subparser_distance_to_root.add_argument_group('optional arguments')

    required_distance_to_root.add_argument('--input-fasta',
                        help='Fasta file to read. Must be aligned to Wuhan-Hu-1',
                        required=True,
                        dest='fasta_in',
                        metavar='input.fasta')
    required_distance_to_root.add_argument('--input-metadata',
                        help='Metadata to read. Must contain epi week information',
                        required=True,
                        dest='metadata_in',
                        metavar='input.csv')


    subparser_distance_to_root.set_defaults(func=datafunk.subcommands.distance_to_root.run)

    # ___________________________ mask _____________________________________#

    subparser_mask = subparsers.add_parser(
        """mask""",
        usage="""datafunk mask -i <file.fasta> -m <mask.txt> -o <file.masked.fasta>""",
        description="""mask regions of a fasta file using information in an external file""",
        help="""mask regions of a fasta file using information in an external file""")

    subparser_mask._action_groups.pop()
    required_mask = subparser_mask.add_argument_group('required arguments')

    required_mask.add_argument('-i', '--input-fasta',
                        help='Fasta file to mask',
                        required=True,
                        dest='fasta_in',
                        metavar='input.fasta')
    required_mask.add_argument('-o', '--output-fasta',
                            help='Fasta file to write',
                            required=True,
                            dest='fasta_out',
                            metavar='output.fasta')
    required_mask.add_argument('-m,','--mask-file',
                        help='File with mask instructions to parse',
                        required=True,
                        dest='mask_file',
                        metavar='mask.txt')


    subparser_mask.set_defaults(func=datafunk.subcommands.mask.run)

    # _________________________________ curate_lineages ____________________________#

    subparser_curate_lineages = subparsers.add_parser(
        "curate_lineages",
        description="Find new lineages, merge ones that need merging, split ones that need splitting",
        help="Find new lineages, merge ones that need merging, split ones that need splitting",
        usage="datafunk curate_lineages -i <input_directory>",
    )

    subparser_curate_lineages.add_argument(
        "-i", "--input-directory",
        dest="input_directory",
        required=True,
        type=str,
        help="Path to input directory containing traits.csv files",
    )

    subparser_curate_lineages.add_argument(
        "-o", "--output_file",
        dest="output_file",
        required=False,
        default="updated_lineage_assignments.csv",
        type=str,
        help="Name of output CSV",
    )

    subparser_curate_lineages.set_defaults(func=datafunk.subcommands.curate_lineages.run)

    # _________________________________ snp_finder ____________________________#

    subparser_snp_finder = subparsers.add_parser(
        "snp_finder",
        description="Query an alignment position for informative SNP",
        help="Query an alignment position for informative SNP",
        usage="datafunk snp_finder -i <input_directory>",
    )

    subparser_snp_finder.add_argument("-a", action="store", type=str, dest="a")
    subparser_snp_finder.add_argument("--snp-csv", action="store", type=str, dest="snp")
    subparser_snp_finder.add_argument("-o", action="store", type=str, dest="o")

    subparser_snp_finder.set_defaults(func=datafunk.subcommands.snp_finder.run)

    # _________________________________ del_finder ____________________________#

    subparser_del_finder = subparsers.add_parser(
        "del_finder",
        description="Query an alignment position for deletions",
        help="Query an alignment position for deletions",
        usage="datafunk del_finder -i <input.fasta> --deletions-file <deletions.csv> --genotypes-table <results.csv>",
    )
    subparser_del_finder._action_groups.pop()
    required_del_finder = subparser_del_finder.add_argument_group('required arguments')
    optional_del_finder = subparser_del_finder.add_argument_group('optional arguments')

    required_del_finder.add_argument('-i', '--input-fasta',
                        help='Alignment (to Wuhan-Hu-1) in Fasta format to type',
                        required=True,
                        dest='fasta_in',
                        metavar='input.fasta')
    optional_del_finder.add_argument('-o', '--output-fasta',
                            help='Fasta file to write',
                            required=False,
                            dest='fasta_out',
                            metavar='output.fasta')
    required_del_finder.add_argument('--deletions-file',
                        help='Input CSV file with deletions type. Format is: 1-based start position of deletion,length of deletion (dont include a header line), eg: 1605,3',
                        required=True,
                        dest='del_file',
                        metavar='deletions.csv')
    required_del_finder.add_argument('--genotypes-table',
                        help='CSV file with deletion typing results to write. Returns the genotype for each deletion in --deletions-file for each sequence in --input-fasta: either "ref", "del" or "X" (for missing data)',
                        required=True,
                        dest='genotypes_file',
                        metavar='results.csv')
    optional_del_finder.add_argument('--append-as-SNP',
                        help='If invoked, then append the genotype of the deletion as a SNP on the end of the alignment',
                        required=False,
                        dest='append_snp',
                        action='store_true')

    subparser_del_finder.set_defaults(func=datafunk.subcommands.del_finder.run)

    # _________________________________ add_header_column ____________________________#

    subparser_add_header_column = subparsers.add_parser(
        "add_header_column",
        description="Add header column to metadata table corresponding to fasta record ids",
        help="Add header column to metadata table corresponding to fasta record ids"
    )

    subparser_add_header_column.add_argument(
        "--input-fasta",
        dest="input_fasta",
        required=True,
        type=str,
        help="Input FASTA",
    )
    subparser_add_header_column.add_argument(
        "--input-metadata",
        dest="input_metadata",
        required=True,
        type=str,
        help="Input CSV or TSV",
    )
    subparser_add_header_column.add_argument(
        "--output-metadata",
        dest="output_metadata",
        required=True,
        type=str,
        help="Output CSV",
    )
    subparser_add_header_column.add_argument(
        "--output-fasta",
        dest="output_fasta",
        required=True,
        type=str,
        help="Output FASTA",
    )
    subparser_add_header_column.add_argument(
        "--gisaid",
        dest="gisaid",
        action="store_true",
        required=False,
        help="Input data is from GISAID",
    )
    subparser_add_header_column.add_argument(
        "--cog-uk",
        dest="cog_uk",
        action="store_true",
        required=False,
        help="Input data is from COG-UK",
    )
    subparser_add_header_column.add_argument(
        "--log",
        dest="log_file",
        required=False,
        help="Log file to use (otherwise uses stdout)"
    )
    subparser_add_header_column.add_argument(
        "--column-name",
        dest="column_name",
        required=False,
        default='header',
        help="Name of column in metadata corresponding to fasta header"
    )
    subparser_add_header_column.add_argument(
        "--columns",
        dest="columns",
        nargs='+',
        required=False,
        help="List of columns in metadata to parse for string matching with fasta header"
    )

    subparser_add_header_column.set_defaults(func=datafunk.subcommands.add_header_column.run)


        # ___________________________ extract_unannotated_seqs ____________________________________#

    subparser_extract_unannotated_seqs = subparsers.add_parser(
        """extract_unannotated_seqs""",
        usage="""datafunk extract_unannotated_seqs --input-tree <file.tree> --input-metadata <file.csv> --output-tree <file.tree> --output-metadata <file.csv>""",
        description="""extract sequences with an empty cell in a specified cell in a metadata table""",
        help="""extract sequences with an empty cell in a specified cell in a metadata table""")

    subparser_extract_unannotated_seqs._action_groups.pop()
    required_extract_unannotated_seqs = subparser_extract_unannotated_seqs.add_argument_group('required arguments')

    required_extract_unannotated_seqs.add_argument('--input-fasta',
                        help='fasta file to extract sequences from',
                        required=True,
                        dest='fasta_in',
                        metavar='input.fasta')
    required_extract_unannotated_seqs.add_argument('--input-metadata',
                        help='metadata whose columns and rows will be checked',
                        required=True,
                        dest='metadata_in',
                        metavar='input.csv')
    required_extract_unannotated_seqs.add_argument('--null-column',
                        help='metadata column which will be checked as empty',
                        required=True,
                        dest='null_column')
    required_extract_unannotated_seqs.add_argument('--index-column',
                        help='metadata column to match to fasta file',
                        required=True,
                        dest='index_column')
    required_extract_unannotated_seqs.add_argument('--output-fasta',
                        help='fasta file to write',
                        required=True,
                        dest='fasta_out',
                        metavar='output.fasta')


    subparser_extract_unannotated_seqs.set_defaults(func=datafunk.subcommands.extract_unannotated_seqs.run)

    # _________________________________ AA_finder ____________________________#

    subparser_AA_finder = subparsers.add_parser(
        "AA_finder",
        description="Query a codon position for amino acids",
        help="Query a codon position for amino acids",
        usage="datafunk AA_finder -i <input.fasta> --codons-file <codons.csv> --genotypes-table <results.csv>",
    )
    subparser_AA_finder._action_groups.pop()
    required_AA_finder = subparser_AA_finder.add_argument_group('required arguments')

    required_AA_finder.add_argument('-i', '--input-fasta',
                        help='Alignment (to Wuhan-Hu-1) in Fasta format to type',
                        required=True,
                        dest='fasta_in',
                        metavar='input.fasta')
    required_AA_finder.add_argument('--codons-file',
                        help='Input CSV file with codon locations to parse. Format is: name,1-based start position of codon (dont include a header line), eg: d614g,23402',
                        required=True,
                        dest='AA_file',
                        metavar='codons.csv')
    required_AA_finder.add_argument('--genotypes-table',
                        help='CSV file with amino acid typing results to write. Returns the amino acid at each position in --codons-file for each sequence in --input-fasta, or "X" for missing data',
                        required=True,
                        dest='genotypes_file',
                        metavar='results.csv')

    subparser_AA_finder.set_defaults(func=datafunk.subcommands.AA_finder.run)

    # _________________________________ bootstrap ____________________________#

    subparser_bootstrap = subparsers.add_parser(
        "bootstrap",
        description="bootstrap an alignment",
        help="bootstrap an alignment",
        usage="datafunk bootstrap -i <input.fasta> --output-prefix boot -n 100",
    )
    subparser_bootstrap._action_groups.pop()
    required_bootstrap = subparser_bootstrap.add_argument_group('required arguments')
    optional_bootstrap = subparser_bootstrap.add_argument_group('optional arguments')

    required_bootstrap.add_argument('-i', '--input-fasta',
                        help='Alignment in fasta format to bootstrap',
                        required=True,
                        dest='fasta_in',
                        metavar='input.fasta')
    optional_bootstrap.add_argument('-p', '--output-prefix',
                        help='Prefix for output files (default is "bootstrap_")',
                        required=False,
                        default="bootstrap_",
                        dest='output_prefix',
                        metavar='boot')
    optional_bootstrap.add_argument('-n',
                        help='Number of bootstraps to generate (default is 1)',
                        required=False,
                        default=1,
                        metavar='1',
                        type=int)

    subparser_bootstrap.set_defaults(func=datafunk.subcommands.bootstrap.run)

    # ___________________________________________________________________________#

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
