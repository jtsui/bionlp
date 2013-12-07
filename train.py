import json
from utils import *
from parse_utils import *
import chemtagger
import chem_canonicalizer
from pymongo import *
import sys

COMMON = set(['a', 'purified'])


def clean_chemicals():
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


def find_reactants(sentence, substrate_set, product_set):
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
    if (len(found_subs) > 0 and len(found_prods) > 0 and
        (len(found_subs) > 1 or len(found_prods) > 1 or
         found_subs != found_prods)):  # chemical can be substrate and product
        return '%s => %s' % (' + '.join(found_subs), ' + '.join(found_prods))
    return None


def match_name():
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
            rxn_readable = '%s => %s' % (' + '.join([str(x) for x in sorted(list(sub_ids))]),
                                         ' + '.join([str(x) for x in sorted(list(prod_ids))]))
            sub_set = set([y for x in sub_ids for y in chemicals[str(x)]])
            prod_set = set([y for x in prod_ids for y in chemicals[str(x)]])
            rxn_found = find_reactants(sentence, sub_set, prod_set)
            if rxn_found:
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
    chemical_inchi_map = json.load(open('../data/chemid_inchi_map.json'))
    sentences = json.load(open('../data/train_sentences.json'))
    chemicals = json.load(open('../data/train_tag_sentences.json'))
    abstracts_rxns = json.load(open('../data/train_abstracts.json'))
    bar, i = pbar(len(sentences)), 0
    print 'Getting matches by inchi'
    bar.start()
    matches = {}
    for sid, sent in sentences.iteritems():
        pmid = sid.split('-')[0]
        chems = chemicals.get(sid, [])
        inchis = chem_canonicalizer.names_to_inchi(chems)
        inchi_set = set(inchis.keys())
        rxns = abstracts_rxns[pmid]['reactions']
        if len(inchis) < 2:
            continue
        sentence_reactants = defaultdict(set)
        for rxn_id, reaction in rxns.iteritems():
            substrate_set = set([chemical_inchi_map[str(x)]
                                for x in reaction['substrates']])
            product_set = set([chemical_inchi_map[str(x)]
                               for x in reaction['products']])
            s_intersect = inchi_set.intersection(substrate_set)
            p_intersect = inchi_set.intersection(product_set)
            # chemical can be substrate and product
            if (len(s_intersect) > 0 and len(p_intersect) > 0 and
                (len(s_intersect) > 1 or len(p_intersect) > 1 or
                 s_intersect != p_intersect)):
                key = '%s => %s' % (' + '.join([inchis[x] for x in s_intersect]),
                                    ' + '.join([inchis[x] for x in p_intersect]))
                val = '%s => %s' % (' + '.join([str(x) for x in sorted(list(reaction['substrates']))]),
                                    ' + '.join([str(x) for x in sorted(list(reaction['products']))]))
                sentence_reactants[key].add(val)
        if len(sentence_reactants) > 0:
            reactants = dict([(x, list(y))
                             for x, y in sentence_reactants.items()])
            matches[sid] = {'sentence': sent,
                            'reactants': reactants,
                            }
        i += 1
        bar.update(i)
    bar.finish()
    json.dump(matches, open('../data/train_match_inchi.json', 'wb'),
              indent=2, sort_keys=True)
    print 'Results dumped to ../data/train_match_inchi.json'


def combine():
    match_inchi = json.load(open('../data/train_match_inchi.json'))
    match_name = json.load(open('../data/train_match_name.json'))
    combined = match_name
    for sid, data in match_inchi.iteritems():
        if sid in combined:
            combined[sid]['reactants'].update(data['reactants'])
            pr.pprint(combined[sid])
            import pdb
            pdb.set_trace()
        else:
            combined[sid] = data
    json.dump(combined, open('../data/train_match.json', 'wb'),
              indent=2, sort_keys=True)


def stats():
    abstracts = json.load(open('../data/train_abstracts.json'))
    sentences = json.load(open('../data/train_sentences.json'))
    match_inchi = json.load(open('../data/train_match_inchi.json'))
    match_name = json.load(open('../data/train_match_name.json'))
    match = json.load(open('../data/train_match.json'))
    print 'Overview:'
    print '  Abstracts: %s' % len(abstracts)
    print '  Sentences: %s' % len(sentences)
    print 'Sentence with at least 1 substrate and at least 1 product:'
    print '  Match by name'
    abstracts_name = set([x.split('-')[0] for x in match_name.keys()])
    print '    Sentences: %s' % len(match_name)
    print '    Abstracts: %s' % len(abstracts_name)
    print '  Match by inchi'
    abstracts_inchi = set([x.split('-')[0] for x in match_inchi.keys()])
    print '    Sentences: %s' % len(match_inchi)
    print '    Abstracts: %s' % len(abstracts_inchi)
    print '  Intersection'
    print '    Sentences: %s' % len(set(match_inchi.keys()).intersection(set(match_name.keys())))
    print '    Abstracts: %s' % len(abstracts_inchi.intersection(abstracts_name))
    print '  Combined'
    print '    Sentences: %s' % len(match)
    print '    Abstracts: %s' % len(set([x.split('-')[0] for x in match.keys()]))


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
    print 'Wrong number of arguments. Usage: python rapier.py [sentences, chemicals, matchname, tag, matchinchi]'


if __name__ == '__main__':
    main()
