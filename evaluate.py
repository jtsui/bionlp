import re
import sys
import json
import patterns
from utils import *
from nltk.stem.lancaster import LancasterStemmer


def sanitize_chemicals(chemicals):
    chemicals = set(chemicals)
    return [x for x in chemicals if len(x) > 1]


def extract(sentence, chemicals, stemmer, pattern_set):
    reactants = []
    if chemicals is None:
        return reactants
    chemicals = sanitize_chemicals(chemicals)
    chems = [y for x in chemicals for y in x.split()]
    stemmed_sentence = ' '.join([x if x in chems else stemmer.stem(x)
                                 for x in sentence.split()])
    tagged_sentence = stemmed_sentence
    for chemical in chemicals:
        tagged_sentence = tagged_sentence.replace(
            chemical, '$%s$cmp' % chemical)
    for pattern_id in pattern_set:
        for pattern in pattern_set[pattern_id]:
            matches = re.findall(pattern, tagged_sentence)
            if matches:
                reactants.append((pattern_id, matches))
    return reactants


def match_reaction(tagged_reactions, reactants):
    for tagged_reaction in tagged_reactions:
        substrates, products = deserialize_rxn(tagged_reaction)
        matched_substrate = False
        for substrate in substrates:
            if substrate in reactants:
                matched_substrate = True
                break
        matched_product = False
        for product in products:
            if product in reactants:
                matched_product = True
                break
        if matched_product and matched_substrate:
            return True
    return False


def evaluate():
    stemmer = LancasterStemmer()
    pattern_set = patterns.expand_patterns()
    chemicals = json.load(open('../data/train_tag_sentences.json'))
    sentences = json.load(open('../data/train_sentences.json'))
    training_set = json.load(open('../data/train_match.json'))
    tp = []  # Pattern match and subset matches with BRENDA tagged reaction
    tn = []  # Pattern match but not BRENDA tagged or tagged but not a subset
    fp = []  # Not a pattern match and not tagged by BRENDA
    fn = []  # Not a pattern match but tagged by BRENDA
    bar, i = pbar(len(sentences)), 0
    print 'Evaluating patterns'
    bar.start()
    for sid, sentence in sentences.iteritems():
        reactants = extract(sentence, chemicals.get(sid), stemmer, pattern_set)
        if reactants:  # pattern match
            if sid in training_set:
                tagged_reactions = training_set[sid]['reactants'].keys()
                found_match = False
                for pattern_id, reactions in reactants:
                    if found_match:
                        break
                    for reaction in reactions:
                        if match_reaction(tagged_reactions, reaction):
                            tp.append((sid, [(pattern_id, [reaction])]))
                            found_match = True
                            break
                if not found_match:
                    fp.append((sid, reactants))
            else:
                fp.append((sid, reactants))
        else:  # not a pattern match
            if sid in training_set:  # tagged by BRENDA
                fn.append(sid)
            else:  # not tagged by BRENDA
                tn.append(sid)
        i += 1
        bar.update(i)
    bar.finish()
    json.dump(tp, open('../data/evaluate_tp.json', 'wb'), indent=2)
    json.dump(fp, open('../data/evaluate_fp.json', 'wb'), indent=2)
    json.dump(tn, open('../data/evaluate_tn.json', 'wb'), indent=2)
    json.dump(fn, open('../data/evaluate_fn.json', 'wb'), indent=2)


def stats():
    pass


def main():
    '''
    Commands are in order of dependencies.
    '''
    if len(sys.argv) == 2:
        this_file, command = sys.argv
        if command == 'run':
            evaluate()
            return
        elif command == 'stats':
            stats()
            return
    print 'Wrong number of arguments. Usage: python evaluate.py [run, stats]'


if __name__ == '__main__':
    main()
