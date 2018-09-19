from numpy import loadtxt, load


def get_plaintext(file):
    plain_seq = []
    f = open(file, 'r')
    for line in f:
        line = line.strip().split()
        plain_seq.extend(line)
    return plain_seq


def get_ciphertext(file):
    cipher_seq = []
    f = open(file, 'r')
    for line in f:
        line = line.strip().split()
        cipher_seq.extend(line)
    return cipher_seq


def get_bigram(file):
    mat = load(file)
    assert (abs(mat[0, :].sum() - 1) < 0.0001)
    return mat


def get_trigram(file, fmt='npy'):
    if fmt == 'npy':
        mat = load(file)
        assert (abs(mat[0, 0, :].sum() - 1) < 0.0001)
    return mat


def toint(c):
    return ord(c) - 65


def tochar(i):
    return chr(i + 65)


def tointseq(charseq):
    seq = []
    for c in charseq:
        seq += [toint(c)]
    return seq


def tocharseq(intseq):
    seq = []
    for i in intseq:
        seq += [tochar(i)]
    return seq

