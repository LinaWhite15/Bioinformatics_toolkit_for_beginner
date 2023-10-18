from modules.fastq_toolkit import *
from modules.protein_toolkit import *
import os

def read_fastq(input_path: str) -> Dict[str, Tuple[str, str,str]]:
    with open(input_path, 'r') as file:
        titles = []
        contigs = []
        comms = []
        quals = []
        for line in file:
            if line.startswith('@SRX'):
                name = line.strip('\n')
                titles.append(name)
                contig = file.readline().strip('\n')
                contigs.append(contig)
                comm = file.readline().strip('\n')
                comms.append(comm)
                qual = file.readline().strip('\n')
                quals.append(qual)
        val = list(zip(contigs, comms, quals))
        fastq = dict(zip(titles, val))
    return fastq


def write_fastq(filtered_fasq: dict, output_filename: str) -> None:
    if not os.path.exists('fastq_filtrator_resuls'):
        os.mkdir('fastq_filtrator_resuls')
    with open(f'fastq_filtrator_resuls.fastq', 'w') as output_file:
        for title, val in filtered_seqs.items():
            output_file.write(f'{title}\n')
            output_file.write(f'{val[0]}\n')
            output_file.write(f'{val[1]}\n')
            output_file.write(f'{val[2]}\n')


def print_result(result: list, corrupt_seqs: list):
    len_seq, len_corr_seq = len(result), len(corrupt_seqs)
    len_seqs = len_seq + len_corr_seq
    success = ["+" for _ in range(len_seqs)]
    if not len_corr_seq:
        print(f"All {len_seqs} sequence(s) processed successfully")
    elif len_corr_seq:
        for i in corrupt_seqs:
            success[i[0]] = "-"
        print(f'Processing result: [{"".join(success)}]\n')
        print(f"{len_seq} sequence(s) out of {len_seq + len_corr_seq} given have been processed successfully.")
        print(f"{len_corr_seq} has been recognized as corrupted, i.e. non-protein")


OPERATIONS = {"content_check": aa_content_check, "seq_length": seq_length, "protein_formula": protein_formula, "protein_mass": protein_mass, "charge": aa_chain_charge}


def protein_processing(*args, abbreviation: int = 1):
    """
    This function makes it possible to process protein sequences to identify them, determine their length, molecular weight, amino acid composition and charge.
    """
    *seqs, operation = args
    if operation not in OPERATIONS:
        raise ValueError(f'Unknown operation `{operation}`. Please, select from: "content_check", "seq_length", "protein_formula", "protein_mass", "charge"')

    result, corrupt_seqs = [], []
    for seq_index, seq in enumerate(seqs):
        if abbreviation == 3:
            seq = decomposition(seq)
        is_seq_valid = check_seq(seq, abbreviation)
        if is_seq_valid:
            if abbreviation == 3:
                seq = seq_transform(seq)
            result.append(OPERATIONS[operation](seq))
        elif not is_seq_valid:
            corrupt_seqs.append((seq_index, seq))

    print_result(result, corrupt_seqs)

    res_len, cor_seq_len = len(result), len(corrupt_seqs)
    result = result[0] if res_len >= 1 else result
    corrupt_seqs = corrupt_seqs[0] if cor_seq_len >= 1 else corrupt_seqs
    return result, corrupt_seqs


def filter_fastq(input_path: str, gc_bounds=(0,100), length_bounds=(0, 2**23), quality_threshold=0, output_filename = ''):
    for name, seq in seqs.items():
       gc_count = GC_content(seq)
       seq_length = length_seq(seq)
       average_quality_score = quality_seq(name)

    if output_filename is None:
        output_filename = os.path.split(input_path)[-1]

    seqs = read_fastq(input_path)
    filtered_seqs = {}

    lowergc, uppergc = gc_bounds
    lowerlength, upperlength = length_bounds

    for name, seq in seqs.items():
      if lowergc <= gc_count <= uppergc and lowerlength <= seq_length <= upperlength:
        if quality_threshold <= average_quality_score:
          filtered_seqs[name] = seq[name]

    write_fastq(filtered_seqs, output_filename)