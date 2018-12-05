from math import log2, ceil
from collections import OrderedDict
from util import convert_fraction_to_binary
from bisect import insort

def shannon_fano(prob_map):
    # Begin by sorting the probabilities in decreasing order, as required
    # in Shannon's paper.
    prob_map_sorted = OrderedDict(sorted([(key, prob_map[key]) for key in prob_map.keys() if prob_map[key] > 0.], key=lambda tup: tup[1], reverse=True))
    probs_sorted = list(prob_map_sorted.values())

    # Compute the cumulative probability distribution
    cum_prob_list = [0.]

    for p in probs_sorted[:-1]: # Do not append the last cumulative probability with value 1.0
        cum_prob_list.append(cum_prob_list[-1] + p)

    # Convert the cumulative prob. list into a dictionary
    cum_prob_map = dict([(key, cum_prob_list[i]) for i, key in enumerate(prob_map_sorted.keys())])

    # Assign the codewords to each symbol
    code_map = {}
    for symbol in prob_map_sorted.keys(): # For each symbol to encode
        prob = prob_map_sorted[symbol]
        cum_prob = cum_prob_map[symbol]

        # Compute the codeword length
        code_length = ceil(log2(1. / prob))

        # Compute the codeword which is the binary representation of cum_prob
        code = convert_fraction_to_binary(cum_prob, num_bits=code_length)
        code_map[symbol] = code
    return code_map


def huffman(prob_map):
    # Convert the prob_map to an ordered dictionary to make sure the order is preserved
    prob_map = OrderedDict(prob_map)

    # create an xtree with all the source symbols (to be the leaves) initially orphaned
    xtree = [[-1,[], symbol] for symbol in prob_map.keys()]
    # label the probabilities with a "pointer" to their corresponding nodes in the tree
    # in the process, we convert the probability vector from a Python dictionary to a
    # list of tuples (so we can modify it)
    node_prob_list = [(node_idx, p) for node_idx, p in enumerate(prob_map.values())]

    while len(node_prob_list) > 1:
        node_prob_list = sorted(node_prob_list, key=lambda tup: tup[1]) # sort in ascending order

        # Smallest probaility nodes:
        least_prob_nodes = list(map(lambda tup: tup[0], node_prob_list[:2]))

        # Append a new node to the tree with no parent (parent = -1) and children
        # being the two nodes with the smallest probabilities
        # the leaves are labeled according to the symbols in the probability vector.
        # label the remaining tree nodes with numbers starting at len(prob_map)
        new_node = [-1, least_prob_nodes, str(len(xtree))]
        new_node_prob = sum(map(lambda tup: tup[1], node_prob_list[:2]))

        xtree.append(new_node)

        # Assign parent to the nodes pointed to by the smallest probabilities
        # The parent idx is the current last node
        new_node_idx = len(xtree) - 1
        for idx in least_prob_nodes:
            xtree[idx][0] = new_node_idx

        # remove the two nodes with the smallest probability
        node_prob_list.pop(0)
        node_prob_list.pop(0)

        # Add the new combined node the the node prob. list
        node_prob_list.append((new_node_idx, new_node_prob))
    return xtree


def bits2bytes(x):
    n = len(x)+3
    r = (8 - n % 8) % 8
    prefix = format(r, '03b')
    x = ''.join(str(a) for a in x)
    suffix = '0'*r
    x = prefix + x + suffix
    x = [x[k:k+8] for k in range(0,len(x),8)]
    y = []
    for a in x:
        y.append(int(a,2))

    return y

def bytes2bits(y):
    x = [format(a, '08b') for a in y]
    r = int(x[0][0:3],2)
    x = ''.join(x)
    x = [int(a) for a in x]
    for k in range(3):
        x.pop(0)
    for k in range(r):
        x.pop()
    return x


def vl_encode(input_sequence, code_map):
    encoded_sequence = []
    for symbol in input_sequence:
        encoded_sequence.extend(code_map[symbol])
    return encoded_sequence


def vl_decode(y, xt):
    x = []
    root = [k for k in range(len(xt)) if xt[k][0]==-1]
    if len(root) != 1:
        raise NameError('Tree with no or multiple roots!')
    root = root[0]
    leaves = [k for k in range(len(xt)) if len(xt[k][1]) == 0]

    n = root
    for k in y:
        if len(xt[n][1]) < k:
            raise NameError('Symbol exceeds alphabet size in tree node')
        if xt[n][1][k] == -1:
            raise NameError('Symbol not assigned in tree node')
        n = xt[n][1][k]
        if len(xt[n][1]) == 0: # it's a leaf!
            x.append(xt[n][2])
            n = root
    return x
