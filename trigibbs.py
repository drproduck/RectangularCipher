import numpy as np
from ios import *
from collections import Counter


# pt = plaintext
# ct = ciphertext
# suppose ciphertext is generated from plaintext following a Dirichlet process with base dsitribution a sufficiently sparse cateforical, and hyperparameter \alpha
# p(ct | pt) = \alpha * p(ct | pt) + n(ct,pt) / \alpha

class Cache():
    def __init__(self, plain_seq, cipher_seq):
        counters = dict()
        for p, c in zip(plain_seq, cipher_seq):
            if p not in counters:
                counters[p] = Counter()
            counters[p][c] += 1
        self.counters = counters

    def get_counter(self, plain_char):
        if plain_char not in self.counters:
            self.counters[plain_char] = Counter()
        return self.counters[plain_char]

    def decrement(self, plain_char, cipher_char):
        cipher_counts = self.get_counter(plain_char)
        if cipher_char not in cipher_counts:
            return 0
        else:
            cipher_counts[cipher_char] -= 1
            return cipher_counts[cipher_char]

    def increment(self, plain_char, cipher_char):
        cipher_counts = self.get_counter(plain_char)
        cipher_counts[cipher_char] += 1
        return cipher_counts[cipher_char]

    def get_counts(self, plain_char, cipher_char):
        """
        return count of plain character and count of cipher character given the plain character
        """
        if plain_char not in self.counters:
            return 0, 0
        else:
            return sum(self.counters[plain_char].values()), self.counters[plain_char][cipher_char]


class GibbsSampler:
    def __init__(self, cipher_seq, bigram):
        self.T = len(cipher_seq)
        self.cipher_seq = cipher_seq
        plain_seq = self.random_plaintext(self.T)
        self.cache = Cache(plain_seq, cipher_seq)
        self.plain_seq = plain_seq
        self.bigram = bigram
        self.alpha = 0.01
        self.prior = 1.0 / 26
        self.plain_dict = np.arange(26, dtype=int)

    def random_plaintext(self, length):
        letter_pos = dict()
        counts = Counter()
        for pos, c in enumerate(self.cipher_seq):
            counts[c] += 1
            if c in letter_pos:
                letter_pos[c] += [pos]
            else:
                letter_pos[c] = [pos]
        common_letters = 'ETAOINSHRDLCUM'
        l = len(common_letters)
        plain_seq = []
        for i, (c, _) in enumerate(reversed(counts.most_common())):
            for pos in letter_pos[c]:
                plain_seq.append(common_letters[l - 1 - i if l - 1 - i >= 0 else 0])
        plain_seq = tointseq(plain_seq)
        return plain_seq

    def get_gap_probability(self, pos):
        """
        return probability of character at position in the senquence
        """

        if pos == 0:
            before = 1
        else:
            before = np.zeros(len(self.plain_dict))
            for i, plain_char in enumerate(self.plain_dict):
                before[i] = self.bigram[self.plain_seq[pos - 1], plain_char]

        if pos == self.T - 1:
            after = 1
        else:
            after = np.zeros(len(self.plain_dict))
            for i, plain_char in enumerate(self.plain_dict):
                after[i] = self.bigram[plain_char, self.plain_seq[pos + 1]]

        return before * after

    def get_emission_probability(self, pos):
        cipher_char = self.cipher_seq[pos]
        p = []
        for plain_char in self.plain_dict:
            pn_count, cr_count = self.cache.get_counts(plain_char, cipher_char)
            p.append((self.alpha * self.prior + cr_count) / (self.alpha + pn_count))

        return p

    def sample_gap(self, pos):
        old_plain_char = self.plain_seq[pos]
        cipher_char = self.cipher_seq[pos]

        # decrement
        self.cache.decrement(old_plain_char, cipher_char)

        # sample
        emit = self.get_emission_probability(pos)
        gap = self.get_gap_probability(pos)
        probs = []
        for i in range(len(emit)):
            probs += [emit[i] * gap[i]]
        probs = probs / sum(probs)
        sampled_char = np.random.choice(self.plain_dict, p=probs)

        # increment
        self.cache.increment(sampled_char, cipher_char)

        # update
        self.plain_seq[pos] = sampled_char

        return sampled_char

    def sample_seq(self):
        for i in range(self.T):
            self.sample_gap(i)

    def sample(self, n_iter):
        for i in range(n_iter):
            print('iteration', i)
            self.sample_seq()
            print(tocharseq(self.plain_seq))

    def get_likelihood(self):
        pass


def main():
    plain_seq = tointseq(get_plaintext('408plaincleaned'))
    cipher_seq = get_ciphertext('408ciphercleaned')
    bigram = get_bigram('data/processed/bigram.npy')
    sampler = GibbsSampler(cipher_seq, bigram)
    sampler.sample(5000)


if __name__ == '__main__':
    main()
