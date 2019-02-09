from collections import Counter
import math
def computePWM(dataset, pseudocount_value, alphabet):
    sequence_length = len(dataset[0])
    for i in range(1, len(dataset)):
        if len(dataset[i]) != sequence_length:
            raise UnalignedSequencesError
    num_sequences = len(dataset)
    
    pwm = dict()
    for amino_acid in alphabet:
        pwm[amino_acid] = list()
    num_amino_acids = len(alphabet)
    for position in range(0, sequence_length):
        column = list([sequence[position] for sequence in dataset])
        counts = Counter(column)
        for amino_acid in alphabet:
            #if that column doesn't contain the given amino acid, then counts[amino_acid] simply returns 0, as it should 
            count = counts[amino_acid]            
            pwm[amino_acid].append(math.log2((count + pseudocount_value)/(num_sequences + pseudocount_value*num_amino_acids)))
    return pwm

def likelihood(pwm, sequence):
    total = 0.0
    print('sequence')
    print(sequence)
    for i in range(0, len(sequence)):
        total += pwm[sequence[i]][i]
    return total
