import json
from utils import *
from parse_utils import *
import chemtagger
import chem_canonicalizer
from pymongo import *
import sys

COMMON_WORDS = set(['a', 'purified'])


def clean_chemicals():
    chemicals = json.load(open('../data/train_chemicals.json'))
    clean_chemicals = {}
    bar, i = pbar(len(chemicals)), 0
    print 'Cleaning chemical names'
    bar.start()
    for chem_id, names in chemicals.iteritems():
        clean_names = [x.lower() for x in names]
        # add plural versions of chemicals
        clean_names = flatten_list([[x, pluralize(x)] for x in clean_names])
        clean_names = set(clean_names)
        clean_names.difference_update(COMMON_WORDS)
        clean_chemicals[chem_id] = list(clean_names)
        i += 1
        bar.update(i)
    bar.finish()
    json.dump(clean_chemicals,
              open('../data/train_chemicals_cleaned.json', 'wb'), indent=2)


def readable(rxn_rep, chemicals):
    sub_ids = rxn_rep.split(' => ')[0].split(' + ')
    prod_ids = rxn_rep.split(' => ')[1].split(' + ')
    prod_readable = '\n\t+\n'.join([', '.join(chemicals[x]) for x in prod_ids
                                    if len(chemicals[x]) > 0])
    sub_readable = '\n\t+\n'.join([', '.join(chemicals[x]) for x in sub_ids
                                   if len(chemicals[x]) > 0])
    print '%s\n\t => \n%s' % (sub_readable, prod_readable)


def generate_training():
    abstracts = json.load(open('../data/train_abstracts.json'))
    bar, i = pbar(len(abstracts)), 0
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


def tag_chemicals():
    sentences = json.load(open('../data/train_sentences.json'))
    bar, i = pbar(len(sentences)), 0
    bar.start()
    chemicals = {}
    for sid, sentence in sentences.iteritems():
        chems = chemtagger.get_compounds(sentence)
        if chems:
            chemicals[sid] = chems
        i += 1
        bar.update(i)
    bar.finish()
    json.dump(chemicals, open('../data/train_sentences_tagged.json', 'wb'),
              indent=2, sort_keys=True)


def find_reaction(sentence, substrate_set, product_set):
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


def find_reactants(sentence, reactions, chemicals):
    reactants = defaultdict(set)
    for rxn_id, reaction in reactions.iteritems():
        sub_ids = set(reaction['substrates'])
        prod_ids = set(reaction['products'])
        if sub_ids == prod_ids:
            continue
        sub_names = [y for x in sub_ids for y in chemicals[str(x)]]
        prod_names = [y for x in prod_ids for y in chemicals[str(x)]]
        if len(sub_names) == 0 or len(prod_names) == 0:
            continue
        sub_set = set(sub_names)
        prod_set = set(prod_names)
        rxn_readable = '%s => %s' % (' + '.join([str(x) for x in sub_ids]),
                                     ' + '.join([str(x) for x in prod_ids])),
        rxn_found = find_reaction(sentence, sub_set, prod_set)
        if rxn_found:
            reactants[rxn_found].add(rxn_readable)
    return dict([(x, list(y)) for x, y in reactants.items()])


def match_names():
    sentences = json.load(open('../data/train_sentences.json'))
    chems = json.load(open('../data/train_chemicals_cleaned.json'))
    abstracts_rxns = json.load(open('../data/train_abstracts.json'))
    bar, i = pbar(len(sentences)), 0
    print 'Getting matches by name'
    bar.start()
    matches = {}
    for sid, sentence in sentences.iteritems():
        pmid = sid.split('-')[0]
        rxns = abstracts_rxns[pmid]['reactions']
        sentence_reactants = find_reactants(sentence, rxns, chems)
        if len(sentence_reactants) > 0:
            matches[sid] = {'reactants': sentence_reactants,
                            'sentence': sentence,
                            }
        i += 1
        bar.update(i)
    bar.finish()
    json.dump(matches, open('../data/train_matches_name.json', 'wb'),
              indent=2, sort_keys=True)
    print 'Results dumped to ../data/train_matches_name.json'


def match_inchis():
    chemical_inchi_map = json.load(open('../data/chemid_inchi_map.json'))
    sentences = json.load(open('../data/train_sentences.json'))
    chemicals = json.load(open('../data/train_sentences_tagged.json'))
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
            substrate_set = set([chemical_inchi_map(x)
                                for x in reaction['substrates']])
            product_set = set([chemical_inchi_map(x)
                               for x in reaction['products']])
            s_intersect = inchi_set.intersection(substrate_set)
            p_intersect = inchi_set.intersection(product_set)
            # chemical can be substrate and product
            if (len(s_intersect) > 0 and len(p_intersect) > 0 and
                (len(s_intersect) > 1 or len(p_intersect) > 1 or
                 s_intersect != p_intersect)):
                key = '%s => %s' % (' + '.join([inchis[x] for x in s_intersect]),
                                    ' + '.join([inchis[x] for x in p_intersect]))
                val = '%s => %s' % (' + '.join([str(x) for x in reaction['substrates']]),
                                    ' + '.join([str(x) for x in reaction['products']]))
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
    json.dump(matches, open('../data/train_matches_inchi.json', 'wb'),
              indent=2, sort_keys=True)
    print 'Results dumped to ../data/train_matches_inchi.json'


def stats():
    abstracts = json.load(open('../data/train_abstracts.json'))
    sentences = json.load(open('../data/train_sentences.json'))
    match_inchi = json.load(open('../data/train_matches_inchi.json'))
    match_name = json.load(open('../data/train_matches_name.json'))
    print 'Abstracts: %s' % len(abstracts)
    print 'Sentences: %s' % len(sentences)
    print 'Match by name'
    abstracts_name = set([x.split('-')[0] for x in match_name.keys()])
    print '  Sentences: %s' % len(match_name)
    print '  Abstracts: %s' % len(abstracts_name)
    print 'Match by inchi'
    abstracts_inchi = set([x.split('-')[0] for x in match_inchi.keys()])
    print '  Sentences: %s' % len(match_inchi)
    print '  Abstracts: %s' % len(abstracts_inchi)
    print 'Intersection'
    print '  Sentences: %s' % len(set(match_inchi.keys()).intersection(set(match_name.keys())))
    print '  Abstracts: %s' % len(abstracts_inchi.intersection(abstracts_name))
    combined = {}
    for sid, data in match_inchi.iteritems():
        if sid in combined:
            combined[sid]['reactants'].update(data['reactants'])
        else:
            combined[sid] = data
    print 'Combined'
    print '  Sentences: %s' % len(combined)
    print '  Abstracts: %s' % len(set([x.split('-')[0] for x in combined.keys()]))


def main():
    if len(sys.argv) == 2:
        this_file, command = sys.argv
        if command == 'sentences':
            generate_training()  # generates train_sentences.json
            return
        elif command == 'chemicals':
            clean_chemicals()  # generates train_chemicals_cleaned.json
            return
        elif command == 'matchname':
            match_names()  # generates train_matches_name.json
            return
        elif command == 'tag':
            tag_chemicals()  # generates train_sentences_tagged.json
            return
        elif command == 'matchinchi':
            match_inchis()  # train_matches_inchi.json
            return
        elif command == 'stats':
            stats()  # train_matches_inchi.json
            return
    print 'Wrong number of arguments. Usage: python rapier.py [sentences, chemicals, matchname, tag, matchinchi]'


if __name__ == '__main__':
    main()
