# Which DNA patterns play the role of molecular clocks
# A transcription factor regulates a gene by bind- ing to a specific short DNA interval called a regulatory motif,
# or transcription factor binding site, in the gene’s upstream region, a 600-1000 nucleotide-long region preceding
# the start of the gene

nucleotides = ["A", "C", "G", "T"]


# Chapter 1 function
def pattern_to_number(pattern):
    pattern = pattern.upper()
    number = 0
    for i in range(len(pattern)):
        number += nucleotides.index(pattern[i]) * 4 ** (len(pattern) - i - 1)
    return number


# Chapter 1 function
def number_to_pattern(index, k):
    text = []
    number = index
    while number // 4 != 0:
        text.insert(0, nucleotides[number % 4])
        number = number // 4

    if number != 0:
        text.insert(0, nucleotides[number])

    if len(text) != k:
        for i in range(k - len(text)):
            text.insert(0, "A")
    return "".join(text)


# Chapter 1 function
def hamming_distance(p, q):
    distance = 0
    if len(p) != len(q):
        return "Unequal size"
    for i in range(len(q)):
        if p[i] != q[i]:
            distance += 1
    return distance


# Chapter 1 function
def approximate_pattern_matching(pattern, text, d):
    pattern_positions = []
    for i in range(len(text) - len(pattern) + 1):
        current_check = text[i:(i + len(pattern))]
        if hamming_distance(current_check, pattern) <= d:
            pattern_positions.append(i)
    return pattern_positions


# Chapter 1 function
def neighbours(pattern, d):
    if d == 0:
        return pattern
    elif len(pattern) == 1:
        return {"A", "C", "G", "T"}

    first_symbol = pattern[0]
    suffix = pattern[1:]
    neighbourhood = []
    suffix_neighbours = neighbours(suffix, d)
    for text in suffix_neighbours:
        if hamming_distance(suffix, text) < d:
            for nucleotide in nucleotides:
                neighbourhood.append(nucleotide + text)
        else:
            neighbourhood.append(first_symbol + text)

    neighbourhood = set(neighbourhood)

    return neighbourhood


def hamming_diff_sites(pattern, string):
    if len(pattern) == len(string):
        return hamming_distance(pattern, string)
    else:
        final_distance = len(string)
        current_distance = 0
        for i in range(len(string) - len(pattern) + 1):
            current_distance = hamming_distance(string[i:(len(pattern) + i)], pattern)
            if final_distance > current_distance:
                final_distance = current_distance

    return final_distance


# Score function that doesn't work
def score_wrong(motifs):
    score_to_return = 0
    scores = 0

    for i in range(len(motifs[0])):
        current_column_dict = {"A": 0,
                               "C": 0,
                               "G": 0,
                               "T": 0}
        for j in range(len(motifs)):
            current_column_dict[motifs[j][i]] += 1
        most_abundant_nucleotide = max(current_column_dict.values())
        score_to_return = [int(value) for key, value in current_column_dict.items() if
                           value != most_abundant_nucleotide]
        scores += sum(score_to_return)

    return scores
