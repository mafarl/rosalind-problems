from chapter2 import nucleotides


# Given a profile matrix Profile, we can evaluate the probability of every k-mer in a string Text and
# find a Profile-most probable k-mer in Text, i.e., a k-mer that was most likely to have been generated
# by Profile among all k-mers in Text.
def profile_most_probable_k_mer(text, k, matrix):
    # Converting matrix elements to floats if they are strings
    # new_matrix = []
    # for m in range(len(matrix)):
    #    a = "".join(matrix[m])
    #    a = a.split())    #    row = map(float, a)
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


with open('../data/2c.txt') as input_data:
    dna = input_data.readline().strip()
    k = int(input_data.readline())
    profile = [list(map(float, line.strip().split())) for i, line in enumerate(input_data.readlines())]

output = profile_most_probable_k_mer(dna, k, profile)

# Print and save the answer.
print(output)
with open('../output/Assignment_02C.txt', 'w') as output_data:
    output_data.write(output)
