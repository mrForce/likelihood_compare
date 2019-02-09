from pwm import *
import argparse
import statistics
import math
import matplotlib.pyplot as plt


def smooth(data, window_size):
    smoothed = []
    print('data')
    print(data)
    for i in range(0, len(data) - window_size):
        print(data[i:(i + window_size)])
        smoothed.append(statistics.mean(data[i:(i + window_size)]))
    return smoothed

def absolute_difference(data, smoothed):
    total = 0.0
    for i in range(0, len(smoothed)):
        total += abs(smoothed[i] - data[i])
    return total/len(smoothed)

parser = argparse.ArgumentParser(description='Give this a Q-value threshold and set of peptides, Q-values. Construct a PWM with peptides under the threshold, then plot the sum of absolute differences; select the window size, and plot the log likelihood')
parser.add_argument('file', help='a tab seperated file with peptide in first column, the q-value in second')
parser.add_argument('threshold', type=float, help='peptides below this threshold used to construct PWM')
pseudocount_value = 1


args = parser.parse_args()
data = []
alphabet = set()
with open(args.file, 'r') as f:
    for line in f:
        fields = line.split('\t')
        peptide = fields[0].strip()
        alphabet = alphabet.union(set(peptide))
        q_value = float(fields[1].strip())
        data.append((peptide, q_value))
data.sort(key= lambda x: x[1])

high_confidence_peptides = [x[0] for x in filter(lambda x: x[1] <= args.threshold, data)]

pwm = computePWM(high_confidence_peptides, pseudocount_value, list(alphabet))

likelihoods = [likelihood(pwm, x[0]) for x in data]

window_sizes = []
absolute_differences = []

for x in range(5, int(len(data)/10)):
    window_sizes.append(x)
    absolute_differences.append(absolute_difference(likelihoods, smooth(likelihoods, x)))

plt.plot(window_sizes, absolute_differences)
plt.show()

window_size = int(input('window size: '))
smoothed = smooth(likelihoods, window_size)
plt.plot(smoothed)
plt.show()
