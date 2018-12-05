import string
import numpy as np

from util import decimal_int_to_binary_list, elias_gamma_coding
from math import floor, ceil
from collections import OrderedDict
from sys import stdout as so
from bisect import bisect


class ArithmeticAdaptive(object):
    def __init__(self, precision=32):
        self.update_precision(precision)

    def update_precision(self, precision):
        self.precision = precision
        self.one = int(2**precision - 1)
        self.quarter = int(ceil(self.one / 4.))
        self.half = 2 * self.quarter
        self.threequarters = 3 * self.quarter
        return

    def encode(self, input_seq: bytes):
        """Encode input bytes sequence including the number of symbols to ignore"""
        print("Encoding the sequence")
        compressed_seq = self.encode_seq(input_seq)
        decompressed_seq = self.decode_seq(compressed_seq)

        # Calculate the number of symbols to ignore during decoding
        symbols_to_ignore = len(decompressed_seq) - len(input_seq)
        # Find the Elias gamma coding representation of the number
        symbols_to_ignore_seq = elias_gamma_coding.encode(symbols_to_ignore)

        return symbols_to_ignore_seq + compressed_seq

    def decode(self, input_seq: bytes):
        """Encode input bytes sequence including the number of symbols to ignore"""
        print("Decoding the sequence")

        # Calculate the number of symbols to ignore during decoding
        symbols_to_ignore, N = elias_gamma_coding.decode(input_seq)

        compressed_seq = input_seq[2 * N + 1:]  # Ignore the elias gamma code at beginning
        decompressed_seq = self.decode_seq(compressed_seq)
        # Ignore the number of symbols determined
        return decompressed_seq[:-symbols_to_ignore]

    def encode_seq(self, input_seq: bytes):
        """Encode an input bytes sequence"""
        input_seq_len = len(input_seq)

        self.symbol_freq = OrderedDict({i: 1 for i in range(256)})
        self.alphabet = list(self.symbol_freq.keys())
        self.output_seq = []  # Reset to an empty output sequence

        self.lo, self.hi = 0, self.one # initialise lo and hi to be [0, 1.0)
        self.straddle = 0 # initialise the straddle counter to 0

        for k, symbol in enumerate(input_seq):
            self._update_prob_maps()
            if k % 100 == 0: self.display_progress(current_symbol=k, seq_len=input_seq_len)

            self._encode_single(symbol)

            # Update the frequency of symbol
            self.symbol_freq[symbol] += 1

        # after processing all input symbols, flush bits in the 'straddle' pipeline
        self.straddle += 1 # adding 1 to straddle for "good measure" (ensures prefix-freeness)
        if self.lo < self.quarter: # the position of lo determines the dyadic interval that fits
            self.output_seq.append(0)
            self.output_seq += [1 for s in range(self.straddle)] # output a zero followed by "straddle" ones
        else:
            self.output_seq.append(1)
            self.output_seq += [0 for s in range(self.straddle)] # output a 1 followed by "straddle" zeros
        return self.output_seq

    def _update_prob_maps(self):
        freqs = self.symbol_freq.values()
        total_freq = sum(freqs)
        prob_list = list(map(lambda freq: float(freq) / total_freq, freqs))

        cum_prob_list = [0.]
        for p in prob_list[:-1]: # Do not append the last cumulative probability with value 1.0
            cum_prob_list.append(cum_prob_list[-1] + p)

        self.cum_prob_list = cum_prob_list
        self.prob_list = prob_list
        self.prob_map = OrderedDict(zip(self.alphabet, prob_list))
        self.cum_prob_map = OrderedDict(zip(self.alphabet, cum_prob_list))
        return

    def _encode_single(self, symbol):
        lohi_range = self.hi - self.lo + 1  # The added 1 is necessary to avoid rounding issues

        # 2) narrow the interval end-points [lo,hi) to the new range
        self.lo += ceil(lohi_range * self.cum_prob_map[symbol])
        self.hi = self.lo + floor(lohi_range * self.prob_map[symbol])

        if self.lo == self.hi:
            raise NameError('Zero interval!')

        # Re-scale the interval
        while True:
            if self.hi < self.half: # if lo < hi < 1/2
                # output a 0 followed by 'straddle' ones (if any)
                self.output_seq.append(0)  # append a zero to the output list y
                self.output_seq += [1 for s in range(self.straddle)]  # extend by a sequence of 'straddle' ones
                self.straddle = 0  # zero the straddle counter
            elif self.lo >= self.half: # if hi > lo >= 1/2
                # stretch the interval by 2 and substract 1
                self.output_seq.append(1)  # append a 1 to the output list y
                self.output_seq += [0 for s in range(self.straddle)]  # extend 'straddle' zeros
                self.straddle = 0  # reset the straddle counter

                self.lo -= self.half
                self.hi -= self.half
            elif self.lo >= self.quarter and self.hi < self.threequarters: # if 1/4 < lo < hi < 3/4
                # we can increment the straddle counter and stretch the interval around the half way point.
                self.straddle += 1  # increment straddle
                self.lo -= self.quarter
                self.hi -= self.quarter
            else:
                break # we break the infinite loop if the interval has reached an un-stretchable state
            # now we can stretch the interval (for all 3 conditions above) by multiplying by 2
            self.lo *= 2
            self.hi *= 2
            self.hi += 1  # Add 1 for good measure
        return

    def decode_seq(self, input_seq):
        input_seq = list(input_seq).copy()
        input_seq_len = len(input_seq)

        self.symbol_freq = OrderedDict({i: 1 for i in range(256)})
        self.alphabet = list(self.symbol_freq.keys())
        self.output_seq = []  # Reset to an empty output sequence

        self.lo, self.hi = 0, self.one # initialise lo and hi to be [0, 1.0)

        input_seq.extend(self.precision*[0]) # dummy zeros to make sure we know value to max. precision possible
        self.value = int(''.join(str(sym) for sym in input_seq[:self.precision]), 2)
        self.input_position = self.precision # position where currently reading input seq.

        end_of_file = False
        while not end_of_file:
            if self.input_position % 100 == 0: self.display_progress(self.input_position, input_seq_len)
            self._update_prob_maps()

            try:
                self._decode_stretch_interval(input_seq)
            except IndexError:
                # End of input sequence reached
                end_of_file = True
                # Specify the largest dyadic interval within the final interval
                self.value_lo = self.value
                self.value_hi = self.value + 1
            next_symbol = self._decode_single()
            self.output_seq.append(next_symbol)

            # Update the probabilities and frequencies
            self.symbol_freq[next_symbol] += 1

        # Append the remaining symbols
        while True:
            self._decode_stretch_interval_end_of_seq()
            next_symbol = self._decode_single()
            if (self.lo >= self.value_lo or self.hi < self.value_hi):
                # If the dyadic interval not within [lo, hi] -> stop
                break
            else:
                # While the dyadic interval is still within [lo, hi] continue
                self.output_seq.append(next_symbol)
                # Update the probabilities and frequencies
                self.symbol_freq[next_symbol] += 1

        return self.output_seq

    def _decode_single(self):
        lohi_range = self.hi - self.lo + 1
        a = bisect(self.cum_prob_list, (self.value - self.lo) / lohi_range) - 1
        next_symbol = self.alphabet[a]

        self.lo = self.lo + int(ceil(self.cum_prob_list[a] * lohi_range))
        self.hi = self.lo + int(floor(self.prob_list[a]*lohi_range))
        if (self.lo == self.hi):
            raise NameError('Zero interval!')
        return next_symbol

    def _decode_stretch_interval(self, input_seq):
        while True:
            if self.hi < self.half:
                # do nothing
                pass
            elif self.lo >= self.half:
                self.lo = self.lo - self.half
                self.hi = self.hi - self.half
                self.value = self.value - self.half
            elif self.lo >= self.quarter and self.hi < self.threequarters:
                self.lo = self.lo - self.quarter
                self.hi = self.hi - self.quarter
                self.value = self.value - self.quarter
            else:
                break
            self.lo = 2 * self.lo
            self.hi = 2 * self.hi + 1
            self.value = 2 * self.value + input_seq[self.input_position]
            self.input_position += 1
            if self.input_position == len(input_seq):
                raise IndexError('Reached end of sequence')
        return

    def _decode_stretch_interval_end_of_seq(self):
        while True:
            if self.hi < self.half:
                pass
            elif self.lo >= self.half:
                self.lo = self.lo - self.half
                self.hi = self.hi - self.half
                self.value_lo = self.value_lo - self.half
                self.value_hi = self.value_hi - self.half
            elif self.lo >= self.quarter and self.hi < self.threequarters:
                self.lo = self.lo - self.quarter
                self.hi = self.hi - self.quarter
                self.value_lo = self.value_lo - self.quarter
                self.value_hi = self.value_hi - self.quarter
            else:
                break
            self.lo = 2 * self.lo
            self.hi = 2 * self.hi + 1
            self.value_lo = 2 * self.value_lo
            self.value_hi = 2 * self.value_hi
        return

    def display_progress(self, current_symbol: int, seq_len: int):
        so.write('Adaptive Arithmetic Progress: {}%  \r'.format(floor((current_symbol / seq_len) * 100)))
        so.flush()
        return


class ArithmeticAdaptiveWithStats(ArithmeticAdaptive):
    def __init__(self, precision=32):
        super(ArithmeticAdaptiveWithStats, self).__init__(precision=precision)
        self.prob_vs_time = []
        return

    def _update_prob_maps(self):
        freqs = self.symbol_freq.values()
        total_freq = sum(freqs)
        prob_list = list(map(lambda freq: float(freq) / total_freq, freqs))

        cum_prob_list = [0.]
        for p in prob_list[:-1]: # Do not append the last cumulative probability with value 1.0
            cum_prob_list.append(cum_prob_list[-1] + p)

        self.cum_prob_list = cum_prob_list
        self.prob_list = prob_list
        self.prob_map = OrderedDict(zip(self.alphabet, prob_list))
        self.cum_prob_map = OrderedDict(zip(self.alphabet, cum_prob_list))

        self.prob_vs_time.append(np.array(cum_prob_list))
        return
    
