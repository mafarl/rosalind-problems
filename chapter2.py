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


# A brute force approach for solving the Implanted Motif Problem is based on the observation that any (k, d)-motif
# must be at most d mismatches apart from some k-mer appearing in one of the strings of DnaA brute force approach
# for solving the Implanted Motif Problem is based on the observation that any (k, d)-motif must be at most d mismatches
# apart from some k-mer appearing in one of the strings of Dna
# Rather slow for large k and d
def motif_enumeration(dna, k, d):
    # Split each DNA string into list
    dna_list = dna.splitlines()

    patterns = []
    for i in range(len(dna) - k + 1):
        neighbourhood = neighbours(dna[i:(i + k)].upper(), d)

        if type(neighbourhood) == set:
            for pattern in neighbourhood:
                count = 0

                for string in dna_list:
                    if approximate_pattern_matching(pattern, string, d):
                        count += 1

                if count == len(dna_list):
                    patterns.append(pattern)
        # When d = 0, neighbours returns a string (just one pattern since should be the same ie, 0 mismatches)
        else:
            count = 0

            for string in dna_list:
                if approximate_pattern_matching(neighbourhood, string, d):
                    count += 1

            if count == len(dna_list):
                patterns.append(neighbourhood)

    # Removing duplicated by converting to set
    patterns = str(set(patterns))

    return patterns


def hamming_diff_sites(pattern, string):
    if len(pattern) == len(string):
        return hamming_distance(pattern, string)
    else:
        final_distance = len(string)
        current_distance = 0
        for i in range(len(string) - len(pattern) + 1):
            current_distance = hamming_distance(string[i:(len(pattern)+i)], pattern)
            if final_distance > current_distance:
                final_distance = current_distance

    return final_distance


# ind a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern
# d(Pattern, Dna) - the sum of distances between Pattern and all strings in Dna
# call it median string
def median_string(dna, k):
    median = ""
    distance = 0

    # Generate all possible k-mers
    for i in range(4**k):
        pattern = number_to_pattern(i, k)

        # Used to count hamming distance
        count = 0
        for string in dna:
            count += hamming_diff_sites(pattern, string)
        if i == 0:
            distance = count
            median = pattern
        else:
            if distance > count:
                distance = count
                median = pattern

    # This code considers only k-mers that occur IN the dna
    """
        for i in range(len(dna)):
        for j in range(len(dna[0]) - k + 1):
            pattern = dna[i][j:(j+k)]

            count = 0
            for string in dna:
                count += hamming_diff_sites(pattern, string)
            if i == 0:
                distance = count
                median = pattern
            else:
                if distance > count:
                    distance = count
                    median = pattern
    """

    return median


# Given a profile matrix Profile, we can evaluate the probability of every k-mer in a string Text and
# find a Profile-most probable k-mer in Text, i.e., a k-mer that was most likely to have been generated
# by Profile among all k-mers in Text.
def profile_most_probable_k_mer(text, k, matrix):

    new_matrix = []
    for m in range(len(matrix)):
        a = "".join(matrix[m])
        a = a.split()
        row = map(float, a)
        new_matrix.append(list(row))
    print(new_matrix)

    final_probability = 0
    final_pattern = ""

    for i in range(len(text) - k + 1):
        pattern = text[i:(i + k)]
        current_probability = 1

        for j in range(len(pattern)):
            current_probability *= new_matrix[nucleotides.index(pattern[j])][j]

        if current_probability > final_probability:
            final_probability = current_probability
            final_pattern = pattern

    return final_pattern


def greedy_motif_search(dna, k, t):

    for i in range(len(dna) - k + 1):
        pattern = dna[i:(i+k)]

