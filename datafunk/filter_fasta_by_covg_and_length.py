from Bio import SeqIO
import os 

def coverage_unaligned_non_N(sequence):
    unaligned_seq = sequence.replace("-", "")
    high_covg_seq = unaligned_seq.replace("N", "")
    return(len(unaligned_seq), float(len(high_covg_seq)) / len(unaligned_seq))


def sequence_has_low_coverage(sequence, coverage_threshold):
    unaligned_len, fraction_N = coverage_unaligned_non_N(sequence)
    if unaligned_len == 0 or  fraction_N < coverage_threshold / 100.0:
        return True
    return False


def length_unaligned_sequence(sequence):
    unaligned_seq = sequence.replace("-", "")
    if len(unaligned_seq) == 0:
        return 0
    i = 0
    while i < len(unaligned_seq) and unaligned_seq[i] == "N":
        i += 1
    j = len(unaligned_seq)
    while j > 0 and unaligned_seq[j-1] == "N":
        j -= 1
    unaligned_seq = unaligned_seq[i:j]
    return(len(unaligned_seq))

def sequence_too_short(sequence, length_threshold):
    unaligned_sequence_length = length_unaligned_sequence(sequence)
    if unaligned_sequence_length < length_threshold:
        return True
    return False

def filter_sequences(inpath, outpath, failed_outpath, min_covg=None, min_length=None):
    record_dict = SeqIO.index(inpath, "fasta")

    if outpath is None:
        outpath = inpath.replace(".fa",".filtered.fa")

    low_covg_seqs = []
    short_seqs = []

    if failed_outpath:
        failed_metrics_path = f'{os.path.splitext(failed_outpath)[0]}.metrics.tsv'
        with open(outpath, "w") as out_handle, open(failed_outpath, "w") as failed_out_handle, open(failed_metrics_path, "w") as failed_metrics_handle:
            # write header
            failed_metrics_handle.write(f'sample\tfraction_Ns\tlength\n')
            for seq_name in record_dict:
                record = record_dict[seq_name]
                if min_covg and sequence_has_low_coverage(str(record.seq), min_covg):
                    low_covg_seqs.append(seq_name)
                    failed_out_handle.write(f'>{record.id }\n{str(record.seq)}\n')
                    failed_metrics_handle.write(f'{record.id}\t{coverage_unaligned_non_N(str(record.seq))}\t{sequence_has_low_coverage(str(record.seq))}')
                    continue
                if min_length and sequence_too_short(str(record.seq), min_length):
                    short_seqs.append(seq_name)
                    failed_out_handle.write(f'>{record.id }\n{str(record.seq)}\n')
                    failed_metrics_handle.write(f'{record.id}\t{coverage_unaligned_non_N(str(record.seq))}\t{sequence_has_low_coverage(str(record.seq))}')
                    continue
                out_handle.write('>' + record.id + '\n')
                out_handle.write(str(record.seq) + '\n')
    else:
        with open(outpath, "w") as out_handle:
            for seq_name in record_dict:
                record_seq = str(record_dict[seq_name].seq)
                if min_covg and sequence_has_low_coverage(record_seq, min_covg):
                    low_covg_seqs.append(seq_name)
                    continue
                if min_length and sequence_too_short(record_seq, min_length):
                    short_seqs.append(seq_name)
                    continue
                record = record_dict[seq_name]
                out_handle.write('>' + record.id + '\n')
                out_handle.write(str(record.seq) + '\n')
                # SeqIO.write(record_dict[seq_name], out_handle, "fasta")

            if min_covg:
                print("#Low coverage sequences:")
                for seq_name in low_covg_seqs:
                    print(seq_name)
            if min_length:
                print("#Short/truncated sequences:")
                for seq_name in short_seqs:
                    print(seq_name)

    record_dict.close()
