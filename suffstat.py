from collections import Counter


class Cache:
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

    def _increment_char(self, plain_char, cipher_char):
        cipher_counts = self.get_counter(plain_char)
        cipher_counts[cipher_char] += 1
        return cipher_counts[cipher_char]

    def increment(self, plain, cipher):
        try:
            for p, c in zip(plain, cipher):
                self._increment_char(p, c)
        except TypeError:
            self._increment_char(plain, cipher)

    def _decrement_char(self, plain_char, cipher_char):
        cipher_counts = self.get_counter(plain_char)
        if cipher_char not in cipher_counts:
            return 0
        else:
            cipher_counts[cipher_char] -= 1
            return cipher_counts[cipher_char]

    def decrement(self, plain, cipher):
        try:
            for p, c in zip(plain, cipher):
                self._decrement_char(p, c)
        except TypeError:
            self._decrement_char(plain, cipher)

    def get_counts(self, plain_char, cipher_char):
        """
        return count of plain character and count of cipher character given the plain character
        """
        if plain_char not in self.counters:
            return 0, 0
        else:
            return sum(self.counters[plain_char].values()), self.counters[plain_char][cipher_char]
