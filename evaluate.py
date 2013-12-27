import sys
import json
import patterns
import chemtagger
from utils import *
from tabulate import tabulate
from collections import defaultdict


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
    sentences = json.load(open('../data/train_sentences.json'))
    training_set = json.load(open('../data/train_match.json'))
    # Pattern match and subset matches with BRENDA tagged reaction
    tp = defaultdict(list)
    # Pattern match but no subset match with BRENDA tagged reaction
    tm = defaultdict(list)
    # Not a pattern match and not tagged by BRENDA
    tn = []
    # Pattern match but not BRENDA tagged
    fp = defaultdict(list)
    # Not a pattern match but tagged by BRENDA
    fn = []
    bar, i = pbar(len(sentences)), 0
    print 'Evaluating patterns'
    bar.start()
    for sid, sentence in sentences.iteritems():
        reactants = patterns.extract(sid, sentence)
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
    table = [['Precision', (tp_num / (tp_num + tm_num + fp_num))],
             ['Recall', (tp_num / (tp_num + tm_num + fn_num))]]
    print tabulate(table, [], tablefmt="simple")
    table = [['', 'BRENDA Tagged', 'Not BRENDA Tagged', ''],
             ['Pattern Match', '%i tp | %s tm' %
                 (tp_num, tm_num), '%i fp' % fp_num, '%i' % (tp_num + tm_num + fp_num)],
             ['No Pattern Match', '%i fn' %
                 fn_num, '%s tn' % tn_num, '%i' % (fn_num + tn_num)],
             ['', '%i' % (fn_num + tp_num + tm_num), '%i' % (fp_num + tn_num), '%i' %
                 (tp_num + tm_num + fn_num + fp_num + tn_num)],
             ]
    print tabulate(table, [], tablefmt="grid")
    table = [['True Positive', tp_num],
             ['True Positive Miss', tm_num],
             ['False Positive', fp_num],
             ['False Negative', fn_num],
             ['True Negative', tn_num],
             ]
    print tabulate(table, [], tablefmt="plain")
    table = [['Number of Sentences', len(sentences)],
             ['Check Sentences',
                 (tp_num + fp_num + tn_num + tm_num + fn_num) == len(sentences)],
             ['Number of Positive Labels', len(pos_labels)],
             ['Check Positive Matches',
                 (tp_num + tm_num + fn_num) == len(pos_labels)]
             ]
    print tabulate(table, [], tablefmt="simple")


def main():
    '''
    Commands are in order of dependencies.
    '''
    if len(sys.argv) == 2:
        this_file, command = sys.argv
        if command == 'run':
            evaluate()
            chemtagger.save_map()
            return
        elif command == 'stats':
            stats()
            return
    print 'Wrong number of arguments. Usage: python evaluate.py [run, stats]'


if __name__ == '__main__':
    main()
