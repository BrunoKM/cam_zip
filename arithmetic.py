from math import floor, ceil
from collections import OrderedDict
from sys import stdout as so
from bisect import bisect


def encode(input_seq, prob_map):
    """Encodes input_message using probabilities char_probs"""
    precision = 32
    one = int(2**precision - 1)
    quarter = int(ceil(one/4))
    half = 2*quarter
    threequarters = 3*quarter

    # Remove probabilities of zero
    prob_map = OrderedDict({key: p for key, p in prob_map.items() if p > 0.})

    # Compute cumulative probability
    cum_prob_list = [0.]
    for p in list(prob_map.values())[:-1]: # Do not append the last cumulative probability with value 1.0
        cum_prob_list.append(cum_prob_list[-1] + p)

    cum_prob_map = {symbol: cum_prob for symbol, cum_prob in zip(prob_map.keys(), cum_prob_list)}


    y = [] # initialise output list
    lo,hi = 0, one # initialise lo and hi to be [0, 1.0)
    straddle = 0 # initialise the straddle counter to 0


    for k, symbol in enumerate(input_seq):
        if k % 100 == 0:
            so.write('Arithmetic encoded %d%%    \r' % int(floor(k/len(input_seq) * 100)))
            so.flush()

        lohi_range = hi - lo + 1  # The added 1 is necessary to avoid rounding issues

        # 2) narrow the interval end-points [lo,hi) to the new range [f,f+p]
        # within the old interval [lo,hi]
        lo += ceil(lohi_range * cum_prob_map[symbol])
        hi = lo + floor(lohi_range * prob_map[symbol])

        if lo == hi:
            raise NameError('Zero interval!')

        # Now we need to re-scale the interval if its end-points have bits in common,
        # and output the corresponding bits where appropriate. We will do this with an
        # infinite loop, that will break when none of the conditions for output / straddle
        # are fulfilled
        while True:
            if hi < half: # if lo < hi < 1/2
                # stretch the interval by 2 and output a 0 followed by 'straddle' ones (if any)
                y.append(0)  # append a zero to the output list y
                y += [1 for s in range(straddle)]  # extend by a sequence of 'straddle' ones
                straddle = 0  # zero the straddle counter
            elif lo >= half: # if hi > lo >= 1/2
                # stretch the interval by 2 and substract 1, and output a 1 followed by 'straddle'
                # zeros (if any) and zero straddle after that. Again, HOLD OFF on doing the stretching
                # as this will be done after the if statement, but note that 2*interval - 1 is equivalent
                # to 2*(interval - 1/2), so for now just substract 1/2 from the interval upper and lower
                # bound (and don't forget that when we say "1/2" we mean the integer "half" we defined
                # above: this is an integer arithmetic implementation!
                y.append(1)  # append a 1 to the output list y
                y += [0 for s in range(straddle)]  # extend 'straddle' zeros
                straddle = 0  # reset the straddle counter
                # Scale by 2x and subtract on from lo and hi
                lo -= half
                hi -= half
            elif lo >= quarter and hi < threequarters: # if 1/4 < lo < hi < 3/4
                # we can increment the straddle counter and stretch the interval around
                # the half way point. This can be impemented again as 2*(interval - 1/4),
                # and as we will stretch by 2 after the if statement all that needs doing
                # for now is to subtract 1/4 from the upper and lower bound
                straddle += 1  # increment straddle
                lo -= quarter
                hi -= quarter
            else:
                break # we break the infinite loop if the interval has reached an un-stretchable state
            # now we can stretch the interval (for all 3 conditions above) by multiplying by 2
            lo *= 2
            hi *= 2
            hi += 1
            # ...  multiply hi by 2 and add 1 (I DON'T KNOW WHY +1 IS NECESSARY BUT IT IS. THIS IS MAGIC.
            #      A BOX OF CHOCOLATES FOR ANYONE WHO GIVES ME A WELL ARGUED REASON FOR THIS... It seems
            #      to solve a minor precision problem.)

    # termination bits
    # after processing all input symbols, flush any bits still in the 'straddle' pipeline
    straddle += 1 # adding 1 to straddle for "good measure" (ensures prefix-freeness)
    if lo < quarter: # the position of lo determines the dyadic interval that fits
        y.append(0)
        y += [1 for s in range(straddle)] # output a zero followed by "straddle" ones
    else:
        y.append(1)
        y += [0 for s in range(straddle)] # output a 1 followed by "straddle" zeros
    return(y)

def decode(compressed_seq, prob_map, num_symbols):
    precision = 32
    one = int(2**precision - 1)
    quarter = int(ceil(one/4))
    half = 2*quarter
    threequarters = 3*quarter

    prob_map = OrderedDict({key: p for key, p in prob_map.items() if p > 0.})

    alphabet = list(prob_map.keys())
    prob_list = list(prob_map.values())

    # Compute cumulative probability
    cum_prob_list = [0.]
    for p in list(prob_map.values())[:-1]: # Do not append the last cumulative probability with value 1.0
        cum_prob_list.append(cum_prob_list[-1] + p)

    compressed_seq.extend(precision*[0]) # dummy zeros to prevent index out of bound errors
    output_seq = num_symbols * [0] # initialise all zeros

    # initialise by taking first 'precision' bits from y and converting to a number
    value = int(''.join(str(sym) for sym in compressed_seq[:precision]), 2)
    input_seq_position = precision # position where currently reading y
    lo, hi = 0, one

    for output_seq_position in range(num_symbols):
        if output_seq_position % 100 == 0:
            so.write('Arithmetic decoded %d%%    \r' % int(floor(output_seq_position / num_symbols * 100)))
            so.flush()

        lohi_range = hi - lo + 1
        a = bisect(cum_prob_list, (value-lo) / lohi_range) - 1
        output_seq[output_seq_position] = alphabet[a]

        lo = lo + int(ceil(cum_prob_list[a] * lohi_range))
        hi = lo + int(floor(prob_list[a]*lohi_range))
        if (lo == hi):
            raise NameError('Zero interval!')

        while True:
            if hi < half:
                # do nothing
                pass
            elif lo >= half:
                lo = lo - half
                hi = hi - half
                value = value - half
            elif lo >= quarter and hi < threequarters:
                lo = lo - quarter
                hi = hi - quarter
                value = value - quarter
            else:
                break
            lo = 2*lo
            hi = 2*hi + 1
            value = 2 * value + compressed_seq[input_seq_position]
            input_seq_position += 1
            if input_seq_position == len(compressed_seq):
                raise NameError('Reached end of sequence prematurely')

    return(output_seq)
