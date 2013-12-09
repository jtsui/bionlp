import json
from utils import *
from parse_utils import *
import chemtagger
import chem_canonicalizer
from pymongo import *
import sys

COMMON = set(['a', 'purified'])


def clean_chemicals():
    '''
    input: train_chemicals.json
    output: train_clean_chemicals.json
    loads chemicals and removes common words and duplicates, and add plurals
    '''
    chemicals = json.load(open('../data/train_chemicals.json'))
    clean_chemicals = {}
    bar, i = pbar(len(chemicals)), 0
    print 'Cleaning chemical names'
    bar.start()
    for chem_id, names in chemicals.iteritems():
        clean_names = [x.lower() for x in names if x.lower() not in COMMON]
        # add plural versions of chemicals
        clean_names = flatten_list([[x, pluralize(x)] for x in clean_names])
        clean_names = set(clean_names)
        clean_chemicals[chem_id] = list(clean_names)
        i += 1
        bar.update(i)
    bar.finish()
    json.dump(clean_chemicals,
              open('../data/train_clean_chemicals.json', 'wb'), indent=2)
    print 'Result dumped to ../data/train_clean_chemicals.json'


def generate_sentences():
    '''
    input: train_abstracts.json
    output: train_sentences.json
    loads abstracts and splits into training set of sentences as a map with
    keys <pmid>-<sentence id>
    '''
    abstracts = json.load(open('../data/train_abstracts.json'))
    bar, i = pbar(len(abstracts)), 0
    print 'Generating list of sentences from abstracts'
    bar.start()
    sentences = {}
    for pmid, data in abstracts.iteritems():
        j = 0
        for sentence in splitSentences(data['abstract']):
            sentences['%s-%s' % (pmid, j)] = sentence
            j += 1
        i += 1
        bar.update(i)
    bar.finish()
    json.dump(sentences, open('../data/train_sentences.json', 'wb'),
              indent=2, sort_keys=True)
    print 'Result dumped to ../data/train_sentences.json'


def tag_sentences():
    '''
    input: train_sentences.json
    output: train_tag_sentences.json
    Tags the chemicals in each sentence using ChemicalTagger 
    (http://chemicaltagger.ch.cam.ac.uk/). This method communicates with 
    ChemicalTagger through a custom REST API running on pathway.berkeley.edu
    '''
    sentences = json.load(open('../data/train_sentences.json'))
    bar, i = pbar(len(sentences)), 0
    print 'Tagging chemicals in sentences'
    bar.start()
    chemicals = {}
    for sid, sentence in sentences.iteritems():
        chems = chemtagger.get_compounds(sentence)
        if chems:
            chemicals[sid] = chems
        i += 1
        bar.update(i)
    bar.finish()
    json.dump(chemicals, open('../data/train_tag_sentences.json', 'wb'),
              indent=2, sort_keys=True)
    print 'Result dumped to ../data/train_tag_sentences.json'


def match_criteria(substrates, products):
    '''
    returns True or False if substrates and products list matches our matching
    criteria. match only if there is at least one substrate and at least one 
    product and the lists cannot be the same
    '''
    return (len(substrates) > 0 and
            len(products) > 0 and
           (len(substrates) > 1 or
            len(products) > 1 or
            substrates != products))
    # return len(substrates) > 0 or len(products) > 0


def find_reactants(sentence, substrate_set, product_set):
    '''
    Finds the susbstrates and products in the sentence. Returns the serialized
    reaction or None if not found in sentence
    '''
    found_subs, found_prods = [], []
    sentence_lower = sentence.lower()
    sentence_words = sentence.split()
    sentence_lower_words = sentence_lower.split()
    for substrate in sorted(substrate_set, key=lambda x: len(x)):
        if ' ' in substrate:
            try:
                i = sentence_lower.index(substrate)
                found_subs.append(sentence[i:i + len(substrate)])
                break
            except ValueError:
                continue
        else:
            try:
                i = sentence_lower_words.index(substrate)
                found_subs.append(sentence_words[i])
                break
            except ValueError:
                continue
    for product in sorted(product_set, key=lambda x: len(x)):
        if ' ' in product:
            try:
                i = sentence_lower.index(product)
                found_prods.append(sentence[i:i + len(product)])
                break
            except ValueError:
                continue
        else:
            try:
                i = sentence_lower_words.index(product)
                found_prods.append(sentence_words[i])
                break
            except ValueError:
                continue
    if match_criteria(found_subs, found_prods):
        return serialize_rxn(found_subs, found_prods)
    return None


def match_name():
    '''
    input: train_sentences.json, train_clean_chemicals.json, 
           train_abstracts.json
    output: train_chemicals.json 
    '''
    sentences = json.load(open('../data/train_sentences.json'))
    chemicals = json.load(open('../data/train_clean_chemicals.json'))
    abstracts_rxns = json.load(open('../data/train_abstracts.json'))
    bar, i = pbar(len(sentences)), 0
    print 'Getting matches by name'
    bar.start()
    matches = {}
    for sid, sentence in sentences.iteritems():
        pmid = sid.split('-')[0]
        reactions = abstracts_rxns[pmid]['reactions']
        reactants = defaultdict(set)
        for rxn_id, reaction in reactions.iteritems():
            sub_ids = set(reaction['substrates'])
            prod_ids = set(reaction['products'])
            sub_set = set([y for x in sub_ids for y in chemicals[str(x)]])
            prod_set = set([y for x in prod_ids for y in chemicals[str(x)]])
            rxn_found = find_reactants(sentence, sub_set, prod_set)
            if rxn_found:
                rxn_readable = serialize_rxn(sub_ids, prod_ids)
                reactants[rxn_found].add(rxn_readable)
        sentence_reactants = dict([(x, list(y)) for x, y in reactants.items()])
        if len(sentence_reactants) > 0:
            matches[sid] = {'reactants': sentence_reactants,
                            'sentence': sentence,
                            }
        i += 1
        bar.update(i)
    bar.finish()
    json.dump(matches, open('../data/train_match_name.json', 'wb'),
              indent=2, sort_keys=True)
    print 'Results dumped to ../data/train_match_name.json'


def match_inchi():
    '''
    input: train_sentences.json, train_tag_sentences.json, 
           train_abstracts.json, chemid_inchi_map.json
    output: train_match_inchi.json 
    '''
    chemical_inchi_map = json.load(open('../data/chemid_inchi_map.json'))
    sentences = json.load(open('../data/train_sentences.json'))
    chemicals = json.load(open('../data/train_tag_sentences.json'))
    abstracts_rxns = json.load(open('../data/train_abstracts.json'))
    bar, i = pbar(len(sentences)), 0
    print 'Getting matches by inchi'
    bar.start()
    matches = {}
    for sid, sent in sentences.iteritems():
        i += 1
        bar.update(i)
        pmid = sid.split('-')[0]
        chems = chemicals.get(sid, [])
        inchis = chem_canonicalizer.names_to_inchi(chems)
        inchi_set = set(inchis.keys())
        rxns = abstracts_rxns[pmid]['reactions']
        sentence_reactants = defaultdict(set)
        for rxn_id, reaction in rxns.iteritems():
            substrate_set = set([chemical_inchi_map[str(x)]
                                for x in reaction['substrates']])
            product_set = set([chemical_inchi_map[str(x)]
                               for x in reaction['products']])
            s_intersect = inchi_set.intersection(substrate_set)
            p_intersect = inchi_set.intersection(product_set)
            if match_criteria(s_intersect, p_intersect):
                key = serialize_rxn([inchis[x] for x in s_intersect],
                                    [inchis[x] for x in p_intersect])
                val = serialize_rxn(reaction['substrates'],
                                    reaction['products'])
                sentence_reactants[key].add(val)
        if len(sentence_reactants) > 0:
            reactants = dict([(x, list(y))
                             for x, y in sentence_reactants.items()])
            matches[sid] = {'sentence': sent,
                            'reactants': reactants,
                            }
    bar.finish()
    json.dump(matches, open('../data/train_match_inchi.json', 'wb'),
              indent=2, sort_keys=True)
    print 'Results dumped to ../data/train_match_inchi.json'


def remove_substrings(chemicals):
    '''
    Compress chemical lists by removing names that are substrings of others
    Example: ['ent-kaurenoic acid', u'kaurenoic acid'] -> ['ent-kaurenoic acid']
    '''
    compressed = []
    for a in chemicals:
        is_substring = False
        for b in chemicals:
            if a != b and a in b:
                is_substring = True
        if not is_substring:
            compressed.append(a)
    return compressed


def combine():
    '''
    input: train_match_inchi.json, train_match_name.json
    output: train_match.json 
    Combines the matches by inchi and matches by name and merges chemical 
    names that map to the same reaction
    '''
    match_inchi = json.load(open('../data/train_match_inchi.json'))
    match_name = json.load(open('../data/train_match_name.json'))
    combined = match_name
    for sid, data in match_inchi.iteritems():
        if sid in combined:
            combined[sid]['reactants'] = merge_dols(
                combined[sid]['reactants'], data['reactants'])
            merged_reactants = defaultdict(set)
            for rxn, reactants in invert_dol(combined[sid]['reactants']).items():
                if len(reactants) == 1:
                    merged_reactants[reactants[0]].add(rxn)
                    continue
                # for same reaction, merge two sets of reactants
                subs = remove_substrings(set([y for x in reactants
                                              for y in x.split(' => ')[0].split(' + ')]))
                prods = remove_substrings(set([y for x in reactants
                                               for y in x.split(' => ')[1].split(' + ')]))
                merged_reactants[serialize_rxn(subs, prods)].add(rxn)
            combined[sid]['reactants'] = dict([(x, list(y))
                                               for x, y in merged_reactants.items()])
        else:
            combined[sid] = data
    json.dump(combined, open('../data/train_match.json', 'wb'),
              indent=2, sort_keys=True)


def stats():
    '''
    input: train_abstracts.json, train_sentences.json, train_match_inchi.json,
           train_match_name.json, train_match.json 
    generates a summary report after loading all intermediary match data
    '''
    abstracts = json.load(open('../data/train_abstracts.json'))
    sentences = json.load(open('../data/train_sentences.json'))
    match_inchi = json.load(open('../data/train_match_inchi.json'))
    match_name = json.load(open('../data/train_match_name.json'))
    match = json.load(open('../data/train_match.json'))
    print 'Overview:'
    print '  Abstracts: %s' % len(abstracts)
    print '  Sentences: %s' % len(sentences)
    print 'Match by name:'
    abstracts_name = set([x.split('-')[0] for x in match_name.keys()])
    print '  Sentences: %s' % len(match_name)
    print '  Abstracts: %s' % len(abstracts_name)
    print 'Match by inchi:'
    abstracts_inchi = set([x.split('-')[0] for x in match_inchi.keys()])
    print '  Sentences: %s' % len(match_inchi)
    print '  Abstracts: %s' % len(abstracts_inchi)
    print 'Intersection:'
    print '  Sentences: %s' % len(set(match_inchi.keys()).intersection(set(match_name.keys())))
    print '  Abstracts: %s' % len(abstracts_inchi.intersection(abstracts_name))
    print 'Combined:'
    print '  Sentences: %s' % len(match)
    print '  Abstracts: %s' % len(set([x.split('-')[0] for x in match.keys()]))


def abstract_stats():
    chemical_inchi_map = json.load(open('../data/chemid_inchi_map.json'))
    chemicals = json.load(open('../data/train_tag_sentences.json'))
    clean_chemicals = json.load(open('../data/train_clean_chemicals.json'))
    abstracts_rxns = json.load(open('../data/train_abstracts.json'))
    chemicals_abstract = defaultdict(list)
    for sid, chems in chemicals.iteritems():
        chemicals_abstract[sid.split('-')[0]].extend(chems)
    print 'Matching abstracts by inchi'
    match_inchi = set()
    bar, i = pbar(len(abstracts_rxns)), 0
    bar.start()
    for pmid, data in abstracts_rxns.iteritems():
        rxns = data['reactions']
        abstract = data['abstract']
        i += 1
        bar.update(i)
        chems = chemicals_abstract.get(pmid, [])
        inchi_set = set(chem_canonicalizer.names_to_inchi(chems))
        matched = False
        for rxn_id, reaction in rxns.iteritems():
            substrate_set = set([chemical_inchi_map[str(x)]
                                for x in reaction['substrates']])
            product_set = set([chemical_inchi_map[str(x)]
                               for x in reaction['products']])
            s_intersect = inchi_set.intersection(substrate_set)
            p_intersect = inchi_set.intersection(product_set)
            if match_criteria(s_intersect, p_intersect):
                matched = True
                break
        if matched:
            match_inchi.add(pmid)
    bar.finish()
    print 'Matching abstracts by name'
    match_name = set()
    bar, i = pbar(len(abstracts_rxns)), 0
    bar.start()
    for pmid, data in abstracts_rxns.iteritems():
        reactions = data['reactions']
        abstract = data['abstract']
        i += 1
        bar.update(i)
        matched = False
        for rxn_id, reaction in reactions.iteritems():
            sub_ids = set(reaction['substrates'])
            prod_ids = set(reaction['products'])
            sub_set = set(
                [y for x in sub_ids for y in clean_chemicals[str(x)]])
            prod_set = set(
                [y for x in prod_ids for y in clean_chemicals[str(x)]])
            rxn_found = find_reactants(abstract, sub_set, prod_set)
            if rxn_found:
                matched = True
                break
        if matched:
            match_name.add(pmid)
    bar.finish()
    print 'Abstracts: %s' % len(abstracts_rxns)
    print 'Match by name: %s' % len(match_name)
    print 'Match by inchi: %s' % len(match_inchi)
    print 'Intersection: %s' % len(match_name.intersection(match_inchi))
    print 'Union: %s' % len(match_name.union(match_inchi))


def main():
    '''
    Commands are in order of dependencies.
    '''
    if len(sys.argv) == 2:
        this_file, command = sys.argv
        if command == 'sentences':
            generate_sentences()  # generates train_sentences.json
            return
        elif command == 'clean_chemicals':
            clean_chemicals()  # generates train_clean_chemicals.json
            return
        elif command == 'match_name':
            match_name()  # generates train_match_name.json
            return
        elif command == 'tag_sentences':
            tag_sentences()  # generates train_tag_sentences.json
            return
        elif command == 'match_inchi':
            match_inchi()  # train_match_inchi.json
            return
        elif command == 'combine':
            combine()  # generates train_match.json
            return
        elif command == 'stats':
            stats()
            return
        elif command == 'abstract_stats':
            abstract_stats()
            return
    print 'Wrong number of arguments. Usage: python rapier.py [sentences, chemicals, matchname, tag, matchinchi]'


if __name__ == '__main__':
    main()
