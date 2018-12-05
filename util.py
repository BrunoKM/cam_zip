from math import floor, log2

def convert_fraction_to_binary(x: float, num_bits: int):
    """Take a real number x in the interval (0., 1.) and convert the fractional part to binary"""

    assert x >= 0. and x < 1.

    bits = []
    for pos in range(num_bits):
        # multiply myf by 2 (shifting it "left" in binary)
        x *= 2
        if x >= 1.:
            bits.append(1)
            x -= 1.
        else:
            bits.append(0)
    return bits

def convert_int_to_binary(x: int):
    return list(map(lambda b: int(b), list("{0:b}".format(x))))

def decimal_int_to_binary_list(x: int):
    return list(map(lambda s: int(s), bin(x)[2:]))


class elias_gamma_coding(object):
    @staticmethod
    def encode(x: int):
        """Encode the integer using Elias Gamma Coding"""
        N = floor(log2(x))
        output_seq = [0 for i in range(N)]

        x_binary = convert_int_to_binary(x)
        assert len(x_binary) == N + 1
        output_seq.extend(x_binary)
        return output_seq

    @staticmethod
    def decode(input_seq):
        """The input sequence can be longer than necessary (contains both the code for integer and more)"""
        N = 0
        while input_seq[N] == 0:
            N += 1
        x_binary = list(input_seq[N:2 * N + 1])
        x = int(''.join(map(lambda x: str(x), x_binary)), 2)
        return x, N
