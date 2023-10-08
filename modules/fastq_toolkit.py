def GC_content (seq:str) -> float:
    "Defines GC-content of sequence"

    gc = 0
    for nucl in seq:
        if nucl == "G" or nucl == "C":
            gc += 1
    gc_count = gc / len(seq) * 100

    return gc_count

def length_seq (seq: str) -> float:
    "Defines length of sequence"

    seq_length = len(seq)

    return seq_length

ascii_dict = {
    "!": 0, '"': 1, "#": 2, "$":3, "%":4,
    "&": 5, "'": 6, "(": 7, ")":8, "*": 9,
    "+": 10, ",": 11, "-": 12, ".": 13, "/": 14,
    "0": 15, "1": 16, "2": 17, "3": 18, "4": 19,
    "5": 20, "6": 21, "7": 22, "8": 23, "9": 24,
    ":": 25, ";": 26, "<": 27, "=": 28, ">": 29,
    "?": 30, "@": 31, "A": 32, "B": 33, "C": 34,
    "D": 35, "E": 36, "F": 37, "G": 38, "H": 39,
    "I": 40
}

def quality_seq (seqs: str) -> int:
    """
    Converts symbolic quality metrics to numeric ones
    and counts average quality score of FastQ sequence
    """
    quality_score = 0

    for char in seqs:
        quality_score += ascii_dict[char]

    average_quality_score = quality_score / len(seqs)

    return average_quality_score
