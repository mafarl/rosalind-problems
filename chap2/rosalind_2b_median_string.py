from chapter2 import number_to_pattern, hamming_diff_sites


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


with open('../data/2b.txt') as input_data:
    k = int(input_data.readline())
    dna_list = [line.strip() for line in input_data.readlines()]

output = median_string(dna_list, k)

# Print and save the answer.
print(output)
with open('../output/Assignment_02B.txt', 'w') as output_data:
    output_data.write(output)
