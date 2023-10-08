from typing import Iterable

PROT_SET_1 = frozenset('ARNDCEQGHILKMFPSTWYV')
PROT_SET_3 = frozenset({'Ala','Arg', 'Asn', 'Asp', 'Cys',
                        'Gln', 'Glu', 'Gly', 'His', 'Ile',
                        'Leu', 'Lys', 'Met', 'Phe', 'Pro',
                        'Ser', 'Thr', 'Trp', 'Tyr', 'Val'})
AA_TR_DICT = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
              'Cys': 'C', 'Gln': 'E', 'Glu': 'Q', 'Gly': 'G',
              'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
              'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
              'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}
AA_UNIPROT_CONTENT = {
"A": 9.03, "R": 5.84, "N": 3.79, "D": 5.47, "C": 1.29,
"Q": 3.80, "E": 6.24, "G": 7.27, "H": 2.22, "I": 5.53,
"L": 9.85, "K": 4.93, "M": 2.33, "F": 3.88, "P": 4.99,
"S": 6.82, "T": 5.55, "W": 1.30, "Y": 2.88, "V": 6.86
}
AA_CHARGES = {
"A": 0, "R": 1, "N": 0, "D": -1, "C": 0,
"Q": 0, "E": -1, "G": 0, "H": 1, "I": 0,
"L": 0, "K": 1, "M": 0, "F": 0, "P": 0,
"S": 0, "T": 0, "W": 0, "Y": 0, "V": 0
}


def aa_content_check(seq: str) -> dict:
    seq_content = dict.fromkeys(AA_UNIPROT_CONTENT.keys(), 0)
    for AAcd in seq.upper():
        seq_content[AAcd] = seq_content[AAcd] + 1

    seq_length = len(seq)
    for AAcd, occurence in seq_content.items():
        seq_content[AAcd] = 100 * occurence / seq_length

    return seq_content


def Mann_Whitney_U(seq1: Iterable[int], seq2: Iterable[int]) -> bool:
    """
    Mann-Whitney U-test. Used to compare aminoacids composition in sequence with average composition provided by Uniprot.
    Used as a second step in `check_seq` function if sequence is 1-letter abbreviation.
    """
    len_seq1, len_seq2 = len(seq1), len(seq2)
    r1, r2 = dict.fromkeys(map(str, seq1), 0), dict.fromkeys(map(str, seq2), 0)

    r = sorted(list(seq1) + list(seq2))
    r_dict = dict.fromkeys(map(str, r), ())

    for index, value in enumerate(r):
        value = str(value)
        r_dict[value] = r_dict[value] + (index + 1,)
    for elem in r_dict:
        r_dict[elem] = sum(r_dict[elem]) / len(r_dict[elem])

    for value in seq1:
        value = str(value)
        r1[value] = r1[value] + r_dict[value]

    for value in seq2:
        value = str(value)
        r2[value] = r2[value] + r_dict[value]

    u1 = (len_seq1 * len_seq2) + len_seq1 * (len_seq1 + 1) / 2 - sum(r1.values())
    u2 = (len_seq1 * len_seq2) + len_seq2 * (len_seq2 + 1) / 2 - sum(r2.values())

    u_stat = min(u1, u2)

    if u_stat <= 127:
        return False
    return True


def decomposition(seq):
    len_seq, dec_seq = len(seq), []
    for i in range(0, len_seq, 3):
        dec_seq.append(seq[i:i+3].lower().capitalize())
    return dec_seq


def check_seq(seq: str, abbreviation: int = 1) -> bool:
    """
    Checks whether the string is protein.
    """
    if abbreviation == 3:
        exit_code = set(seq).issubset(PROT_SET_3)
    elif abbreviation == 1:
        seq_set = set(seq.upper())
        exit_code = seq_set.issubset(PROT_SET_1)
        if exit_code:
            seq_content, uniprot_content = aa_content_check(seq).values(), AA_UNIPROT_CONTENT.values()
            seq_Mann_Whitney_U = Mann_Whitney_U(seq_content, uniprot_content) if len(seq_set) == 20 else True
            exit_code = seq_Mann_Whitney_U

    return exit_code


def seq_transform(seq: list):
    seq_tr = ''
    for aa in seq:
        seq_tr += AA_TR_DICT[aa]

    return seq_tr


def seq_length(seq: str) -> int:
    return len(seq)


def protein_mass(seq: str):
    "Counts molecular weight of the protein"
    count = 0
    for amino in seq:
        if amino == "A":
            count += 89
        elif amino == "R":
            count += 174
        elif amino == "N":
            count += 132
        elif amino == "V":
            count += 117
        elif amino == "H":
            count += 155
        elif amino == "G":
            count += 75
        elif amino == "Q":
            count += 146
        elif amino == "E":
            count += 147
        elif amino == "I":
            count += 131
        elif amino == "L":
            count += 131
        elif amino == "K":
            count += 146
        elif amino == "M":
            count += 149
        elif amino == "P":
            count += 115
        elif amino == "S":
            count += 105
        elif amino == "Y":
            count += 181
        elif amino == "T":
            count += 119
        elif amino == "W":
            count += 204
        elif amino == "F":
            count += 165
        else:
            count += 133
    return count


def protein_formula(seq: str):
    "Returns molecular formula of the protein"
    fС = 0
    fH = 0
    fN = 0
    fO = 0
    fS = 0
    for amino in seq:
        if amino == "A":
            fС += 3
            fH += 7
            fN += 1
            fO += 2
            fS += 0
        elif amino == "R":
            fС += 6
            fH += 14
            fN += 4
            fO += 2
            fS += 0
        elif amino == "N":
            fС += 4
            fH += 8
            fN += 2
            fO += 3
            fS += 0
        elif amino == "V":
            fС += 5
            fH += 11
            fN += 1
            fO += 2
            fS += 0
        elif amino == "H":
            fС += 6
            fH += 9
            fN += 3
            fO += 2
            fS += 0
        elif amino == "G":
            fС += 2
            fH += 5
            fN += 1
            fO += 2
            fS += 0
        elif amino == "Q":
            fС += 5
            fH += 10
            fN += 2
            fO += 3
            fS += 0
        elif amino == "E":
            fС += 5
            fH += 9
            fN += 1
            fO += 4
            fS += 0
        elif amino == "I":
            fС += 6
            fH += 13
            fN += 1
            fO += 2
            fS += 0
        elif amino == "L":
            fС += 6
            fH += 13
            fN += 1
            fO += 2
            fS += 0
        elif amino == "K":
            fС += 6
            fH += 14
            fN += 2
            fO += 2
            fS += 0
        elif amino == "M":
            fС += 5
            fH += 11
            fN += 1
            fO += 2
            fS += 1
        elif amino == "P":
            fС += 5
            fH += 9
            fN += 1
            fO += 2
            fS += 0
        elif amino == "S":
            fС += 3
            fH += 7
            fN += 1
            fO += 3
            fS += 0
        elif amino == "Y":
            fС += 9
            fH += 11
            fN += 1
            fO += 3
            fS += 0
        elif amino == "T":
            fС += 4
            fH += 9
            fN += 1
            fO += 3
            fS += 0
        elif amino == "W":
            fС += 11
            fH += 12
            fN += 2
            fO += 2
            fS += 0
        elif amino == "F":
            fС += 9
            fH += 11
            fN += 1
            fO += 2
            fS += 0
        else:
            fС += 4
            fH += 7
            fN += 1
            fO += 4
            fS += 0
    if fS == 0:
        aa_formula = f'С: {fС}, H: {fH}, N: {fN}, O:{fO}'
    else:
        aa_formula = f'С: {fС}, H: {fH}, N: {fN}, O:{fO}, S: {fS}'
    return aa_formula


def aa_chain_charge(seq: str, aa_charges: dict = AA_CHARGES) -> int:
    aa_charge = 0
    for AAcd in seq.upper():
        aa_charge += aa_charges[AAcd]

    return aa_charge