things = []


# Works only if no two patterns fit to be next k-mers
# aka this program can't select the right pattern if there are several options
def string_reconstruction(k, patterns):

    # Find the first pattern (no other pattern ends with first k-1 nucleotides of this pattern)
    # Then find a pattern that starts with last k-1 nucleotides of the first pattern

    for j, pattern in enumerate(patterns):
        not_first = False
        if len(things) == 0:
            for i, patt in enumerate(patterns):
                if i != j and patt[1:] == pattern[:(k-1)]:
                    not_first = True
                    break
        else:
            # If new pattern ends with same nucleotides as previous starts - no match
            # mb create a list here where we'll store all the options and then try one by one seeing if we succeed
            if things[-1][:(k-1)] == pattern[1:]:
                not_first = True

        if not not_first:
            things.append(pattern)
            patterns.remove(pattern)
            if len(patterns) != 1:
                string_reconstruction(k, patterns)
            else:
                things.append(patterns[0])
                return 0


def try_two(k, patterns):
    repeats = []
    for j, pattern in enumerate(patterns):
        not_first = False
        if len(things) == 0:
            for i, patt in enumerate(patterns):
                if i != j and patt[1:] == pattern[:(k - 1)]:
                    not_first = True
                    break
        else:
            # If new pattern ends with same nucleotides as previous starts - no match
            # mb create a list here where we'll store all the options and then try one by one seeing if we succeed
            if things[-1][:(k - 1)] == pattern[1:]:
                not_first = True

        if not not_first:
            repeats.append(pattern)

    for repeat in repeats:
        things.append(repeat)
        patterns.remove(repeat)
        if len(patterns) != 1:
            string_reconstruction(k, patterns)
        else:
            things.append(patterns[0])
            return 0


collection = "AAT ATG ATG ATG CAT CCA GAT GCC GGA GGG GTT TAA TGC TGG TGT"
collection = collection.split(" ")
string_reconstruction(3, collection)
print(things)