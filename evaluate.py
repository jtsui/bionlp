import sys
import json
import patterns
from utils import *
from collections import defaultdict
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
    for chem in chemicals:
        tagged_sentence = tagged_sentence.replace(chem, '$%s$chem' % chem)
    grouped_sentence = patterns.group_list(tagged_sentence)
    for pattern_id in pattern_set:
        for pattern in pattern_set[pattern_id]:
            groups = pattern.findall(grouped_sentence)
            matches = patterns.expand_chems(groups)
            for match in matches:
                if match and len(match) > 1:
                    reactants.append((pattern_id, match))
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
    # Pattern match and subset matches with BRENDA tagged reaction
    tp = defaultdict(list)
    # Pattern match but no subset match with BRENDA tagged reaction
    tm = defaultdict(list)
    tn = []  # Not a pattern match and not tagged by BRENDA
    fp = defaultdict(list)  # Pattern match but not BRENDA tagged
    fn = []  # Not a pattern match but tagged by BRENDA
    bar, i = pbar(len(sentences)), 0
    print 'Evaluating patterns'
    bar.start()
    for sid, sentence in sentences.iteritems():
        reactants = extract(sentence, chemicals.get(sid), stemmer, pattern_set)
        if reactants and sid in training_set:
            tagged_reactions = training_set[sid]['reactants'].keys()
            found_match = False
            for pattern_id, reactant_list in reactants:
                # pattern match and subset match with BRENDA tag
                if match_reaction(tagged_reactions, reactant_list):
                    tp[sid].append([pattern_id, reactant_list])
                    found_match = True
                    break
            # pattern match and tagged by BRENDA but not a subset match
            if not found_match:
                tm[sid].extend(reactants)
        elif not reactants and sid in training_set:
            fn.append(sid)  # not a pattern match and tagged by BRENDA
        elif reactants and sid not in training_set:
            # pattern match but not tagged by BRENDA
            fp[sid].extend(reactants)
        else:
            tn.append(sid)  # not a pattern match and not tagged by BRENDA
        i += 1
        bar.update(i)
    bar.finish()
    json.dump(tp, open('../data/evaluate_tp.json', 'wb'),
              indent=2, sort_keys=True)
    json.dump(tm, open('../data/evaluate_tm.json', 'wb'),
              indent=2, sort_keys=True)
    json.dump(fp, open('../data/evaluate_fp.json', 'wb'),
              indent=2, sort_keys=True)
    json.dump(sorted(tn), open('../data/evaluate_tn.json', 'wb'), indent=2)
    json.dump(sorted(fn), open('../data/evaluate_fn.json', 'wb'), indent=2)


def stats():
    tp = json.load(open('../data/evaluate_tp.json'))
    fp = json.load(open('../data/evaluate_fp.json'))
    tn = json.load(open('../data/evaluate_tn.json'))
    fn = json.load(open('../data/evaluate_fn.json'))
    tm = json.load(open('../data/evaluate_tm.json'))
    sentences = json.load(open('../data/train_sentences.json'))
    pos_labels = json.load(open('../data/train_match.json'))
    tp_num = len(tp) * 1.0
    fp_num = len(fp) * 1.0
    tn_num = len(tn) * 1.0
    fn_num = len(fn) * 1.0
    tm_num = len(tm) * 1.0
    print 'STATS'
    print '-----------------------------'
    width = 20
    print '%s: %i' % ('True posistive'.rjust(width), tp_num)
    print '%s: %i\t(Not a subset match)' % ('True positive miss'.rjust(width), tm_num)
    print '%s: %i' % ('False positive'.rjust(width), fp_num)
    print '%s: %i' % ('True negative'.rjust(width), tn_num)
    print '%s: %i' % ('False negative'.rjust(width), fn_num)
    print '%s: %.4f' % ('Precision'.rjust(width), (tp_num / (tp_num + fp_num + tm_num)))
    print '%s: %.4f' % ('Recall'.rjust(width), (tp_num / (tp_num + tm_num + fn_num)))
    print '-----------------------------'
    width = 25
    print '%s: %i' % ('Number of sentences'.rjust(width), len(sentences))
    print '%s: %s' % ('Check sentences'.rjust(width), ((tp_num + fp_num + tn_num + tm_num + fn_num) == len(sentences)))
    print '%s: %i' % ('Number of positive labels'.rjust(width), len(pos_labels))
    print '%s: %s' % ('Check positive matches'.rjust(width), ((tp_num + tm_num + fn_num) == len(pos_labels)))
    print '-----------------------------'
    print 'DONE'


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
