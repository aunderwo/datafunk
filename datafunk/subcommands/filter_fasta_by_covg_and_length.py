from datafunk.filter_fasta_by_covg_and_length import *

def run(options):
    filter_sequences(options.input_fasta, options.output_fasta, options.failed_qc_output_fasta, options.min_covg, options.min_length)
