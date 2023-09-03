from rosalind_2c_profile_most_probable_k_mer import profile_most_probable_k_mer
from rosalind_2d_greedy_motif_search import score
from rosalind_2e_greedy_motif_search_pseudocounts import laplaces_rule
import random


def random_number_generator(N):
    return random.randrange(N+1)


def randomised_motif_search(dna, k, t):
    m = []

    for string in dna:
        m.append(random_number_generator(len(dna[0])-k))

    motifs = [dna[i][number:(number+k)] for i, number in enumerate(m)]
    best_motifs = motifs
    while True:
        profile = laplaces_rule(motifs, 0)
        motifs = [profile_most_probable_k_mer(dna[i], k, profile) for i in range(t)]
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs


if __name__ == '__main__':
    with open('../data/2f.txt') as input_data:
        k, t = map(int, input_data.readline().split())
        dna_list = [line.strip() for line in input_data.readlines()]

    # Initialize the best scoring motifs as a score higher than the highest possible score.
    best_motifs = [k * t, None]

    # Repeat the radomized motif search 1000 times.
    for repeat in range(1000):
        current_motifs = randomised_motif_search(dna_list, k, t)
        if score(current_motifs) < best_motifs[0]:
            best_motifs[0] = score(current_motifs)
            best_motifs[1] = current_motifs

    # Print and save the answer.
    print('\n'.join(best_motifs[1]))
    with open('../output/Assignment_02F.txt', 'w') as output_data:
        output_data.write('\n'.join(best_motifs[1]))
