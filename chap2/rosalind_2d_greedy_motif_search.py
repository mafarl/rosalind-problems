from chapter2 import nucleotides
from rosalind_2c_profile_most_probable_k_mer import profile_most_probable_k_mer


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


if __name__ == '__main__':
    with open('../data/2d.txt') as input_data:
        k, t = map(int, input_data.readline().split())
        dna_list = [line.strip() for line in input_data.readlines()]

    best_motifs = greedy_motif_search(dna_list, k, t)

    # Print and save the answer.
    print('\n'.join(best_motifs))
    with open('../output/Assignment_02D.txt', 'w') as output_data:
        output_data.write('\n'.join(best_motifs))

