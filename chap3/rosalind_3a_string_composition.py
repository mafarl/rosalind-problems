def string_composition(text, k):

    composition = sorted([text[i:(i+k)] for i in range(len(text) - k + 1)])

    return composition


with open('../data/3a.txt') as input_data:
    k = int(input_data.readline().strip())
    text = input_data.readline().strip()

output = string_composition(text, k)

# Print and save the answer.
print('\n'.join(output))
with open('../output/Assignment_03A.txt', 'w') as output_data:
    output_data.write('\n'.join(output))
