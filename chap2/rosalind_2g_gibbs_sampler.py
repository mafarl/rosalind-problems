import random
from rosalind_2d_greedy_motif_search import score
from rosalind_2e_greedy_motif_search_pseudocounts import laplaces_rule
from rosalind_2c_profile_most_probable_k_mer import profile_most_probable_k_mer


# A more cautious iterative algorithm that discards a single k-mer from the current set of motifs
# at each iteration and decides to either keep it or replace it with a new one
def gibbs_sampler(dna, k, t, N):

    """
    Starts with randomly chosen k-mers in each of t DNA sequences,
    but it makes a random rather than a deterministic choice at each iteration.
    It uses randomly selected k-mers Motifs = (Motif1, . . . , Motift) to come up with another
    (hopefully better scoring) set of k-mers

    Randomly selects an integer i between 1 and t and then randomly changes a single k-mer Motifi
    """

    # Randomly generating initial k-mers
    numbers = [random.randint(0, len(dna[0]) - k) for i in range(t)]
    motifs = [dna[j][number:(number + k)] for j, number in enumerate(numbers)]

    best_motifs = [motifs, score(motifs)]

    for m in range(1, N):
        a = random.randint(0, t)
        profile = laplaces_rule([motif for j, motif in enumerate(motifs) if j != a], 0)
        motifs = [profile_most_probable_k_mer(dna[a], k, profile) if j == a else motif for j, motif in enumerate(motifs)]
        if score(motifs) < best_motifs[1]:
            best_motifs = [motifs, score(motifs)]

    return best_motifs


if __name__ == '__main__':
    with open('../data/2g.txt') as input_data:
        k, t, N = map(int, input_data.readline().split())
        dna_list = [line.strip().upper() for line in input_data.readlines()]

    # Initialize the best scoring motifs as a score higher than the highest possible score.
    best_motifs = [None, k * t]

    # Repeat the randomised motif search 1000 times.
    for repeat in range(20):
        current_motifs = gibbs_sampler(dna_list, k, t, N)
        if current_motifs[1] < best_motifs[1]:
            best_motifs[0] = current_motifs[0]
            best_motifs[1] = current_motifs[1]

    # Print and save the answer.
    print('\n'.join(best_motifs[0]))
    with open('../output/Assignment_02G.txt', 'w') as output_data:
        output_data.write('\n'.join(best_motifs[0]))

