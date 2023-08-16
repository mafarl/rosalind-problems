nucleotides = ["A", "C", "G", "T"]


# just finding repeating patterns in the text
# we say pattern is a most frequent k-mer in Text if it maximises COUNT(Text,Pattern) among all k-mers
# last k-mer of text begins at |text| - k
def pattern_count(text, pattern):
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if text[i:(i + len(pattern))].upper() == pattern.upper():
            count += 1
    return count


# nucleotide_array = [i for i in range(len(text))]
# complexity = O( |text|^2 * k )
def frequent_words(text, k):
    count = [0 for i in range(len(text))]
    frequent_patterns = []
    for i in range(len(text) - k + 1):
        pattern = text[i:(i + k)]
        count[i] = pattern_count(text, pattern)
    max_count = max(count)
    for i in range(len(text) - k + 1):
        if count[i] == max_count:
            frequent_patterns.append(text[i:(i + k)])
    frequent_patterns = set(frequent_patterns)
    return frequent_patterns


# we store all possible k-mers in an array (eg, AA, AC, AG, AT, CA...)
# so each k-mer has an index in its alphabetical order
def pattern_to_number(pattern):
    pattern = pattern.upper()
    number = 0
    for i in range(len(pattern)):
        number += nucleotides.index(pattern[i]) * 4 ** (len(pattern) - i - 1)
    return number


# vice versa - have pattern's index, need to convert it to string
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


# given a region, convert every k-mer into a number and store in the alphabetical frequency array
def computing_frequencies(text, k):
    frequency_array = [0 for i in range(4 ** k)]
    for i in range(len(text) - k + 1):
        pattern = text[i:(i + k)]
        j = pattern_to_number(pattern)
        frequency_array[j] += 1
    return frequency_array


def faster_frequent_words(text, k, t):
    frequent_patterns = []
    frequency_array = computing_frequencies(text, k)
    # max_count = max(frequency_array)
    for i in range(4 ** k):
        """
        # Code for a pattern with the maximum count
            if frequency_array[i] == max_count:
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)
        """
        if frequency_array[i] >= t:
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)
    return frequent_patterns


# DnaA boxes are not necessarily the same, they can be reverse patterns of each other
def reverse_complement(text):
    nucleotides_here = {"A": "T",
                        "C": "G",
                        "G": "C",
                        "T": "A"}
    reverse_text = []
    for letter in text:
        reverse_text.insert(0, letter)

    reverse_complement_text = []
    for letter in reverse_text:
        reverse_complement_text.append(nucleotides_here[letter])
    return "".join(reverse_complement_text)


# finding if there are any other (close) places in the genome where the pattern is repeated
def pattern_matching(pattern, genome):
    positions = []
    for i in range(len(genome) - len(pattern) + 1):
        if genome[i:(i + len(pattern))].upper() == pattern:
            positions.append(i)
    return positions


# The pseudocode below slides a window of length L down Genome
# After computing the frequency array for the current window,
# it identifies (L, t)-clumps simply by finding which k-mers occur at least t times within the window
def clump_finding(genome, k, L, t):
    frequent_patterns = []
    for i in range(len(genome) - L + 1):
        text = genome[i:(i + L)]
        current_patterns = faster_frequent_words(text, k, t)
        for pattern in current_patterns:
            frequent_patterns.append(pattern)

    frequent_patterns = list(set(frequent_patterns))
    return frequent_patterns


# SKEWi(Genome) as the difference between the total number of occurrences of G and the total number of occurrences of C
# in the first i nucleotides of Genome
# return indexes of the nucleotides at which skew difference is the lowest
# The position of the skew minimum is often only a rough indicator of oriC position
def minimum_skew(genome):
    g_c_difference = 0
    skew = [0 for i in range(len(genome))]
    for i in range(len(genome)):
        if genome[i].upper() == nucleotides[1]:
            g_c_difference -= 1
        elif genome[i].upper() == nucleotides[2]:
            g_c_difference += 1
        skew[i] = g_c_difference

    skew.insert(0, 0)
    minimum_indexes = [i for i in range(len(skew)) if skew[i] == min(skew)]
    return minimum_indexes


# The number of mismatches between strings p and q is called the Hamming distance between these strings
def hamming_distance(p, q):
    distance = 0
    if len(p) != len(q):
        return "Unequal size"
    for i in range(len(q)):
        if p[i] != q[i]:
            distance += 1
    return distance


# a k-mer Pattern appears as a substring of Text with at most d mismatches
# if there is some k-mer substring Pattern’ of Text having d or fewer mismatches with Pattern,
# i.e., hamming_distance(Pattern, Pattern') <= d
# (observation that a DnaA box may appear with slight variations leads to the following generalization
# of the Pattern Matching Problem)
def approximate_pattern_matching(pattern, text, d):
    pattern_positions = []
    for i in range(len(text) - len(pattern) + 1):
        current_check = text[i:(i + len(pattern))]
        if hamming_distance(current_check, pattern) <= d:
            pattern_positions.append(i)
    return pattern_positions


# A most frequent k-mer with up to d mismatches in Text is maximising length of array output in
# approximate_pattern_matching
# inefficient version since generate all possible k-mers and check them
def frequent_words_with_mismatches(text, k, d):
    frequent_k_mers = [0 for i in range(4 ** k)]

    for i in range(len(frequent_k_mers)):
        current_pattern = number_to_pattern(i, k)
        length = len(approximate_pattern_matching(current_pattern, text, d))
        frequent_k_mers[i] = length

    frequent_patterns = []
    for i in range(len(frequent_k_mers)):
        if frequent_k_mers[i] == max(frequent_k_mers):
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)

    return frequent_patterns


# generate the 1-neighbourhood of Pattern (hamming distance = 1)
def immediate_neighbours(pattern):
    neighbourhood = [pattern]
    for i in range(len(pattern)):
        symbol = pattern[i]
        for nucleotide in nucleotides:
            if nucleotide != symbol.upper():
                pattern_array = [*pattern]
                pattern_array[i] = nucleotide
                neighbour = "".join(pattern_array)
                neighbourhood.append(neighbour)

    neighbourhood = set(neighbourhood)

    return neighbourhood


# Generate the d-neighborhood NEIGHBORS(Pattern, d),
# the set of all k-mers whose Hamming distance from Pattern does not exceed d
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


# iterative version of neighbours
def iterative_neighbours(pattern, d):
    neighbourhood = [pattern]
    for j in range(d):
        text = pattern[j:]
        result = immediate_neighbours(text)
        for i in result:
            neighbourhood.append(i)
        neighbourhood = list(set(neighbourhood))

    return neighbourhood


# we set CLOSE(i) = 1 whenever Pattern = number_to_pattern(i,k) is close to some k-mer in Text
def faster_frequent_words_with_mismatches(text, k, d):
    frequent_patterns = []
    frequency_array = [0 for i in range(4 ** k)]
    close = [0 for i in range(4 ** k)]

    for i in range(len(text) - k + 1):
        neighbourhood = neighbours(text[i:(i + k)], d)

        if type(neighbourhood) == set:
            for pattern in neighbourhood:
                index = pattern_to_number(pattern)
                close[index] = 1
        else:
            index = pattern_to_number(neighbourhood)
            close[index] = 1

    for j in range(len(close)):
        if close[j] == 1:
            pattern = number_to_pattern(j, k)
            frequency_array[j] = len(approximate_pattern_matching(pattern, text, d))

    frequent_patterns = []
    for i in range(len(frequency_array)):
        if frequency_array[i] == max(frequency_array):
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)

    return frequent_patterns


# in addition to pattern's neighbours
def frequent_words_with_mismatches_reverse_complements(text, k, d):
    frequent_patterns = []
    frequency_array = [0 for i in range(4 ** k)]
    close = [0 for i in range(4 ** k)]

    for i in range(len(text) - k + 1):
        neighbourhood = neighbours(text[i:(i + k)], d)

        if type(neighbourhood) == set:
            for pattern in neighbourhood:
                index = pattern_to_number(pattern)
                close[index] = 1
        else:
            index = pattern_to_number(neighbourhood)
            close[index] = 1

    # what i think should be added because it should look at BOTH Text and its Reverse Complement
    # (i.e. not just looking at Text, and not just looking at the Reverse Complement of Text, but looking
    # at both)
    for m in range(len(close)):
        if close[m] == 1:
            pattern = number_to_pattern(m, k)
            reverse_pattern = reverse_complement(pattern)
            number = pattern_to_number(reverse_pattern)
            close[number] = 1

    for j in range(len(close)):
        if close[j] == 1:
            pattern = number_to_pattern(j, k)
            frequency_array[j] = len(approximate_pattern_matching(pattern, text, d)) + \
                                 len(approximate_pattern_matching(reverse_complement(pattern), text, d))

    frequent_patterns = []
    for i in range(len(frequency_array)):
        if frequency_array[i] == max(frequency_array):
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)

    return frequent_patterns


print(frequent_words_with_mismatches_reverse_complements("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1))
