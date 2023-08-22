# Which DNA patterns play the role of molecular clocks


nucleotides = ["A", "C", "G", "T"]


def hamming_distance(p, q):
    distance = 0
    if len(p) != len(q):
        return "Unequal size"
    for i in range(len(q)):
        if p[i] != q[i]:
            distance += 1
    return distance


def approximate_pattern_matching(pattern, text, d):
    pattern_positions = []
    for i in range(len(text) - len(pattern) + 1):
        current_check = text[i:(i + len(pattern))]
        if hamming_distance(current_check, pattern) <= d:
            pattern_positions.append(i)
    return pattern_positions


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


print(motif_enumeration("AACAA\n"
                        "AAAAA\n"
                        "AAAAA", 3, 0))
