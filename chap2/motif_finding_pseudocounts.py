from chapter2 import nucleotides, profile_most_probable_k_mer, score


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