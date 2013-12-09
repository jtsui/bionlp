import itertools
import progressbar
import pprint

pr = pprint.PrettyPrinter(indent=2)


def serialize_rxn(substrates, products):
    '''
    Serializes a list of substrates and products (either chemical ids or names) 
    into a readable and hashable string. 
    input: substrates = [28248, 2003], products = [10153, 10559]
    output: '2003 + 28248 => 10153 + 10559'
    '''
    return '%s => %s' % (' + '.join([str(x) for x in sorted(substrates)]),
                         ' + '.join([str(x) for x in sorted(products)]))


def deserialize_rxn(rxn):
    '''
    Deserializes reaction string into a tuple of substrate list, product list
    input: '2003 + 28248 => 10153 + 10559'
    output: substrates = [28248, 2003], products = [10153, 10559]
    '''
    subs = [y for y in rxn.split(' => ')[0].split(' + ')]
    prods = [y for y in rxn.split(' => ')[1].split(' + ')]
    return subs, prods


def merge_dols(dol1, dol2):
    '''
    Merges two dictionaries of lists by concatenating list values
    input: dol1 = {'a':[1], 'b':[1]}, dol2 = {'b':[2], 'c':[2]}
    output: {'a': [1], 'c': [2], 'b': [1, 2]}
    '''
    result = dict(dol1, **dol2)
    result.update((k, dol1[k] + dol2[k])
                  for k in set(dol1).intersection(dol2))
    return result


def invert_dol(dol):
    '''
    Inverts a dictionary of lists appending values that have the same key
    input: {'a':[1], 'b':[1,2,3]}
    output: {1: ['a', 'b'], 2: ['b'], 3: ['b']}
    '''
    inv_map = {}
    for k, alist in dol.iteritems():
        for v in alist:
            inv_map[v] = inv_map.get(v, [])
            inv_map[v].append(k)
    return inv_map


def pbar(size):
    '''Returns a new progress bar'''
    bar = progressbar.ProgressBar(maxval=size,
                                  widgets=[progressbar.Bar('=', '[', ']'),
                                           ' ', progressbar.Percentage(),
                                           ' ', progressbar.ETA(),
                                           ' ', progressbar.Counter(),
                                           '/%s' % size])
    return bar


def lines_in_file(file_name):
    '''Opens a file and returns the number of lines. Similar to wc -l'''
    num_lines = 0
    with open(file_name, 'rb') as f_in:
        for line in f_in:
            num_lines += 1
    return num_lines


def flatten_list(alist):
    '''Flattens a list of lists into a list of depth 1'''
    out_list = []
    for l in alist:
        if isinstance(l, list):
            out_list.extend(flatten_list(l))
        else:
            out_list.append(l)
    return out_list


def grouper(iterable, i=2):
    '''takes s and returns (s0,s1), (s1,s2), (s2, s3), etc'''
    # a list of i iterators
    iters = itertools.tee(iterable, i)
    # loop through the iterators
    for j in range(i):
        # advance next j times
        for k in range(j):
            next(iters[j], None)
    # * unpacks the list turning each element to an arg
    return itertools.izip(*iters)


def main():
    '''Put test code here'''
    pass

if __name__ == '__main__':
    main()
