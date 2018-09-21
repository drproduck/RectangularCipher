import numpy as np
from ios import *
from collections import Counter
from suffstat import *
from transition import *
from time import time

# pt = plaintext
# ct = ciphertext
# suppose ciphertext is generated from plaintext following a Dirichlet process with base dsitribution a sufficiently sparse cateforical, and hyperparameter \alpha
# p(ct | pt) = \alpha * p(ct | pt) + n(ct,pt) / \alpha

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
        plain_seq = random_plaintext(self.T)
        self.cache = Cache(plain_seq, cipher_seq)
        self.plain_seq = plain_seq
        self.bigram = bigram
        self.alpha = 0.01
        self.prior = 1.0 / 26
        self.plain_dict = np.arange(26, dtype=int)

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


class BlockSampler:
    def __init__(self, n_hidden, n_obs, hidden_seq_list, obs_seq_list, transition_kernel):
        self.transition = transition_kernel
        self.emission = EmissionKernel(n_hidden, n_obs)
        assert (len(hidden_seq_list) == len(obs_seq_list))
        for hidden_seq, obs_seq in zip(hidden_seq_list, obs_seq_list):
            self.emission.increment(hidden_seq, obs_seq)
        self.n_hidden = n_hidden
        self.n_obs = n_obs
        self.hidden_seq_list = hidden_seq_list
        self.obs_seq_list = obs_seq_list

    def emit(self, h, o, dp_param=0.01):
        """
        rich-get-richer emission
        :param h: hidden
        :param o: observed
        :param dp_param: dirichlet process scale parameter
        :return:
        """
        h_emit = self.emission.get(h)
        return (dp_param * 1.0 / self.n_hidden + h_emit[o]) / (dp_param + h_emit.sum())

    def transit(self, h1, h2):
        return self.transition.get_transition_prob(h2, h1)

    def forward_backward_bigram_blocksample(self, hidden_seq, obs_seq):
        """
        first, calculate p(h_{t-1}=r,h_t = s| d_1^{t})
         = p(h_{t-1}=r,h_t=s,d_t| d_1^{t-1}) / C(h_{t-1},h_t)
         = p(h_{t-1}=r|d_1^{t-1}) p(h_t=s| h_{t-1}=r) p(d_t| h_t=s) / C(h_{t-1},h_t)

         Second, calculate p(h_t=r| d_1^t) (is P[i-1,:,r].sum())

         Third, sample p(\bold{h}| d_1^n) = p(h_n| d_1^n) \prod_t p(h_{n-t}| h_{n-t+1}^n, d_1^n)
         = p(h_n| d_1^n) \prod_t p(h_{n-t}| h_{n-t+1}, d_1^{n-t+1})

         where p(h_{n-t}=r| h_{n-t+1}, d_1^{n-t+1}) ~ p(h_{n-t}=r, h_{n-t+1}| d_1^{n-t+1}) see first step
        :param hidden_seq:
        :param obs_seq:
        :return:
        """
        l = len(hidden_seq)
        P = np.zeros([l - 1, self.n_hidden, self.n_hidden], dtype=float)
        self.emission.decrement(hidden_seq, obs_seq)
        for i in range(l - 1):
            if i == 0:
                for r in range(self.n_hidden):
                    for s in range(self.n_hidden):
                        P[i, r, s] = self.transit(r, s) * self.emit(s, obs_seq[i])
            else:
                for r in range(self.n_hidden):
                    for s in range(self.n_hidden):
                        P[i, r, s] = P[i - 1, :, r].sum() * self.transit(r, s) * self.emit(s, obs_seq[i])

            # normalize
            P[i, :, :] /= P[i, :, :].sum(axis=(0, 1))

        # sampling
        pn = P[l - 2, :, :].sum(axis=0)
        hidden_seq[l - 1] = np.random.choice(self.n_hidden, p=pn)
        for i in reversed(range(l - 1)):
            pi = P[i, :, hidden_seq[i + 1]]
            hidden_seq[i] = np.random.choice(self.n_hidden, p=pi / pi.sum())

        return hidden_seq

    def sample(self, n_iter):
        for it in range(n_iter):
            for i, (h_seq, o_seq) in enumerate(zip(self.hidden_seq_list, self.obs_seq_list)):
                t0 = time()
                seq = self.forward_backward_bigram_blocksample(h_seq, o_seq)
                print('iter',it,'sentence',i,'time',time() - t0)
                print(tocharseq(seq))

def main():
    plain_seq = tointseq(get_plaintext('408plaincleaned'))
    cipher_seq = get_ciphertext('408ciphercleaned')
    bigram = get_bigram('data/processed/bigram.npy')
    # sampler = GibbsSampler(cipher_seq, bigram)
    transition = BigramKernel(bigram)
    # sampler = BlockSampler(26, len(set(plain_seq)))
    # sampler.sample(5000)


if __name__ == '__main__':
    main()
