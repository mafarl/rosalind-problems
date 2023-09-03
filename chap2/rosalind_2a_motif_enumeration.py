from chapter2 import neighbours, approximate_pattern_matching


# A brute force approach for solving the Implanted Motif Problem is based on the observation that any (k, d)-motif
# must be at most d mismatches apart from some k-mer appearing in one of the strings of DnaA brute force approach
# for solving the Implanted Motif Problem is based on the observation that any (k, d)-motif must be at most d mismatches
# apart from some k-mer appearing in one of the strings of Dna
# Rather slow for large k and d
def motif_enumeration(dna, k, d):
    patterns = []
    for i in range(len(dna) - k + 1):
        neighbourhood = neighbours(dna[i:(i + k)], d)

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
    patterns = set(patterns)
    patterns = list(patterns)

    return patterns


with open('../data/2a.txt') as input_data:
    k, d = map(int, input_data.readline().split())
    dna_list = [line.strip() for line in input_data.readlines()]
    dna = "".join(dna_list)

output = motif_enumeration(dna, k, d)

# Print and save the answer.
print('\n'. join(output))
with open('../output/Assignment_02A.txt', 'w') as output_data:
    output_data.write(' '.join(output))
