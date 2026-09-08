"""Inspectable smoothed autoregressive sequence distribution; no efficacy oracle."""
import math
import random
from collections import Counter, defaultdict

ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'


def fit(sequences, order, tau=20.0):
    counts = defaultdict(lambda: [0] * 20)
    lengths = [1] * 43
    for sequence in sequences:
        lengths[len(sequence) - 8] += 1
        for i, aa in enumerate(sequence):
            index = ALPHABET.index(aa)
            for k in range(order + 1):
                context = ('^' * k + sequence[:i])[-k:] if k else ''
                counts[context][index] += 1
    return {'version': 1, 'alphabet': ALPHABET, 'order': order, 'tau': tau,
            'counts': dict(sorted(counts.items())), 'length_weights': lengths,
            'training_sequences': len(sequences)}


class Model:
    def __init__(self, weights):
        self.weights = weights
        self.order = weights['order']
        self.cache = {}

    def probabilities(self, context):
        if context in self.cache:
            return self.cache[context]
        counts = self.weights['counts'].get(context, [0] * 20)
        n = sum(counts)
        if not context:
            values = [(c + 1) / (n + 20) for c in counts]
        else:
            shorter = self.probabilities(context[1:])
            tau = self.weights['tau']
            values = [(c + tau * p) / (n + tau) for c, p in zip(counts, shorter)]
        self.cache[context] = values
        return values

    def nll(self, sequence):
        loss = 0.0
        for i, aa in enumerate(sequence):
            context = ('^' * self.order + sequence[:i])[-self.order:] if self.order else ''
            loss -= math.log2(self.probabilities(context)[ALPHABET.index(aa)])
        return loss / len(sequence)

    def sample(self, rng):
        length = rng.choices(range(8, 51), weights=self.weights['length_weights'], k=1)[0]
        sequence = ''
        for _ in range(length):
            context = ('^' * self.order + sequence)[-self.order:] if self.order else ''
            sequence += rng.choices(ALPHABET, weights=self.probabilities(context), k=1)[0]
        return sequence

    def features(self, sequence):
        length = len(sequence)
        return [self.nll(sequence),
                (sequence.count('K') + sequence.count('R') - sequence.count('D') - sequence.count('E')) / length,
                sum(aa in 'AILMFWVY' for aa in sequence) / length]
