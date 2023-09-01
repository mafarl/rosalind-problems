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
            current_distance = hamming_distance(string[i:(len(pattern) + i)], pattern)
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
    for i in range(4 ** k):
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
    
    # Converting matrix elements to floats if they are strings
    # new_matrix = []
    # for m in range(len(matrix)):
    #    a = "".join(matrix[m])
    #    a = a.split()
    #    row = map(float, a)
    #    new_matrix.append(list(row))

    new_matrix = matrix
    final_probability = -1
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


# Has each nucleotide of a motif in a list rather than list with probabilities for each nucleotide
def build_profile_matrix(motifs):
    """
        Finds the profile matrix for given list of motifs

        motifs: list of motif sequences (list)

        Returns: the profile matrix for motifs (list)
    """

    # profile_matrix = [[] for _ in range(len(motifs[0]))]
    profile_matrix = [[] for _ in range(4)]

    for i in range(len(motifs[0])):

        current_column = {"A": 0,
                          "C": 0,
                          "G": 0,
                          "T": 0}

        for motif in motifs:
            current_column[motif[i]] += 1
        counter = 0
        for key, value in current_column.items():
            profile_matrix[nucleotides.index(key)].append(value / 10)
            counter += 1
    return profile_matrix


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


def score(motifs):
    columns = [''.join(seq) for seq in zip(*motifs)]
    max_count = sum([max([c.count(nucleotide) for nucleotide in 'ACGT']) for c in columns])
    return len(motifs[0])*len(motifs) - max_count


# For the Subtle Motif problem returns score(motifs) = 70, while in book - score(motifs) = 58
# => need to improve
def greedy_motif_search(dna, k, t):
    # Initialise best score as high as possible
    import math
    best_score = math.inf
    # Collection of motifs from each string of dna
    best_motifs = []

    for string in dna:
        dna = [a.upper() for a in dna]
        best_motifs.append(string[:k])

    for i in range(len(dna[0]) - k + 1):
        # Array with current motifs for selected k-mer in dna[0]
        # print(dna[0][i:(i + k)])
        motifs = [dna[0][i:(i + k)]]

        for j in range(1, t):
            """
            profile_to_add = []
            for m in range(len(profile_matrix)):
                profile_to_add.append(profile_matrix[m][i:(i+k)])
            """
            profile_matrix = build_profile_matrix(motifs)
            current_motif = profile_most_probable_k_mer(dna[j], k, profile_matrix)
            motifs.append(current_motif)

        if score(motifs) < best_score:
            best_score = score(motifs)
            best_motifs = motifs

    return best_motifs


# Laplace’s Rule of Succession - substituting zero-probabilities with pseudocounts
# (adds 1 to each element of COUNT(Motifs)
def laplaces_rule(motifs, number):
    profile_matrix = [[] for _ in range(4)]

    for i in range(len(motifs[0])):

        current_column = {"A": 0,
                          "C": 0,
                          "G": 0,
                          "T": 0}

        for motif in motifs:
            current_column[motif[i]] += 1
        counter = 0
        for key, value in current_column.items():
            profile_matrix[nucleotides.index(key)].append((value + 1) / (4 + number))
            counter += 1
    return profile_matrix


# Gives score(motifs) = 39, while 41 in book
def greedy_motif_search_with_pseudocounts(dna, k, t):
    # Initialise best score as high as possible
    import math
    best_score = math.inf
    # Collection of motifs from each string of dna
    best_motifs = []

    for string in dna:
        dna = [a.upper() for a in dna]
        best_motifs.append(string[:k])

    number = 0

    for i in range(len(dna[0]) - k + 1):
        # Array with current motifs for selected k-mer in dna[0]
        # print(dna[0][i:(i + k)])
        motifs = [dna[0][i:(i + k)]]

        for j in range(1, t):

            # Changing this part - instead of building profile matrix, use Laplace's rule to build a profile matrix
            profile_matrix = laplaces_rule(motifs, number)
            number += 1
            current_motif = profile_most_probable_k_mer(dna[j], k, profile_matrix)
            motifs.append(current_motif)

        if score(motifs) < best_score:
            best_score = score(motifs)
            best_motifs = motifs

    return best_motifs


"""
    This part of code test how accurate greedy_motif_search is
    For the Subtle Motif problem returns score(motifs) = 70, while in book - score(motifs) = 58
    => need to improve
    
output = greedy_motif_search_with_pseudocounts(["atgaccgggatactgatAgAAgAAAGGttGGGggcgtacacattagataaacgtatgaagtacgttagactcggcgccgccg",
                           "acccctattttttgagcagatttagtgacctggaaaaaaaatttgagtacaaaacttttccgaatacAAtAAAAcGGcGGGa",
                           "tgagtatccctgggatgacttAAAAtAAtGGaGtGGtgctctcccgatttttgaatatgtaggatcattcgccagggtccga",
                           "gctgagaattggatgcAAAAAAAGGGattGtccacgcaatcgcgaaccaacgcggacccaaaggcaagaccgataaaggaga",
                           "tcccttttgcggtaatgtgccgggaggctggttacgtagggaagccctaacggacttaatAtAAtAAAGGaaGGGcttatag",
                           "gtcaatcatgttcttgtgaatggatttAAcAAtAAGGGctGGgaccgcttggcgcacccaaattcagtgtgggcgagcgcaa",
                           "cggttttggcccttgttagaggcccccgtAtAAAcAAGGaGGGccaattatgagagagctaatctatcgcgtgcgtgttcat",
                           "aacttgagttAAAAAAtAGGGaGccctggggcacatacaagaggagtcttccttatcagttaatgctgtatgacactatgta",
                           "ttggcccattggctaaaagcccaacttgacaaatggaagatagaatccttgcatActAAAAAGGaGcGGaccgaaagggaag",
                           "ctggtgagcaacgacagattcttacgtgcattagctcgcttccggggatctaatagcacgaagcttActAAAAAGGaGcGGa"], 15, 10)

print(score(output))
give = []
this = [''.join(seq) for seq in zip(*output)]
for i in range(len(this)):
    current_column_dict = {"A": 0,
                           "C": 0,
                           "G": 0,
                           "T": 0}
    for j in range(len(this[0])):
        current_column_dict[this[i][j]] += 1
    give.append([key for key, value in current_column_dict.items() if value == max(current_column_dict.values())])

result = ""
for letter in give:
    result += letter[0]
"""