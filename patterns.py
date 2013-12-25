import re
from utils import *
from collections import defaultdict
from nltk.stem.lancaster import LancasterStemmer

stemmer = LancasterStemmer()
chem = r'\$([^\$.]*?)\$chem'
chem_list = r'{([^{}.]*?)}chem_list'
list_pattern = re.compile(r'(%s)(,?\s*(?:and)?\s*%s)*' % (chem, chem))

trans = 'from|to|into|by|are|yield'
# pattern 1 : [trigger1] <0,1> chem_list <0,1> [transition] chem_list
trig1 = 'reduc|convert|produc|form|oxid|transform|bioconvert|synthes|react|interconvert'
patterns1 = [r'(?:%s) %s (?:%s) %s',
             r'(?:%s) \w* %s (?:%s) %s',
             r'(?:%s) \w* %s \w* (?:%s) %s',
             r'(?:%s) %s \w* (?:%s) %s',
             ]

# pattern2 : chem_list <0,1> [trigger2] <0,1> [transition] <0,1> chem_list
trig2 = 'convert|oxid|produc|interconvert'
patterns2 = [r'%s (?:%s) (?:%s) %s',
             r'%s \w* (?:%s) (?:%s) %s',
             r'%s (?:%s) \w* (?:%s) %s',
             r'%s (?:%s) (?:%s) \w* %s',
             r'%s \w* (?:%s) \w* (?:%s) %s',
             r'%s (?:%s) \w* (?:%s) \w* %s',
             r'%s \w* (?:%s) (?:%s) \w* %s',
             r'%s \w* (?:%s) \w* (?:%s) \w* %s',
             ]

# pattern3: chem_list [trigger3] <0,1> chem_list
trig3 = 'yield'
patterns3 = [r'%s (?:%s) %s',
             r'%s (?:%s) \w* %s',
             ]

# pattern4: [trigger4] <0,1> chem_lists
trig4 = 'convert|interconvert'
patterns4 = [r'(?:%s) %s',
             r'(?:%s) \w* %s',
             ]

# pattern5: chem_list [transition5] <0,1> [trigger5] <0,1> chem_list
trig5 = 'produc|metabolit'
trans5 = 'is|ar'
patterns5 = [r'%s (?:%s) (?:%s) %s',
             r'%s (?:%s) \w* (?:%s) %s',
             r'%s (?:%s) \w* (?:%s) \w* %s',
             r'%s (?:%s) (?:%s) \w* %s',
             ]


def expand_patterns():
    patterns = defaultdict(list)
    for pattern in patterns1:
        patterns[1].append(
            re.compile(pattern % (trig1, chem_list, trans, chem_list)))
    for pattern in patterns2:
        patterns[2].append(
            re.compile(pattern % (chem_list, trig2, trans, chem_list)))
    for pattern in patterns3:
        patterns[3].append(re.compile(pattern % (chem_list, trig3, chem_list)))
    for pattern in patterns4:
        patterns[4].append(re.compile(pattern % (trig4, chem_list)))
    for pattern in patterns5:
        patterns[5].append(
            re.compile(pattern % (chem_list, trans5, trig5, chem_list)))
    return patterns


def group_list(sentence):
    return re.sub(list_pattern, '{\g<0>}chem_list', sentence)


def main():
    patterns = expand_patterns()
    tests = [
        (0, 1, ['cholesterol', 'chlorohydrins'],
         'These results indicate that myeloperoxidase can convert $cholesterol$chem to $chlorohydrins$chem and epoxides by a reaction involving HOCl'),
        (1, 2, ['prostaglandins', 'COX-2'],
            'We investigated the possible role of $prostaglandins$chem produced by $COX-2$chem in the immunosuppression observed during Trypanosoma cruzi infection'),
        (2, 3, ['5-formyl-6-methoxy-2,3', '6- and 5-methyl derivatives'],
         'Reduction of 5-methoxy-6-formyl(Ia)- and $5-formyl-6-methoxy-2,3$chem yielded $6- and 5-methyl derivatives$chem Ib and IVb, respectively'),
        (3, 4, ['ester', 'amide'],
         'It results in the mutual conversion of $ester$chem and $amide$chem bonds'),
        (4, 1, ['cholesterol', 'another chemical', 'chlorohydrins'],
         'These results indicate that myeloperoxidase can convert $cholesterol$chem and $another chemical$chem to $chlorohydrins$chem and epoxides by a reaction involving HOCl'),
        (5, 1, ['a', 'b', 'c', 'd', 'e'],
         'The reaction converts $a$chem, $b$chem, and $c$chem to $d$chem and $e$chem.'),
        (6, 1, ['a', 'b', 'c', 'd', 'e'],
         'The reaction converts $a$chem, $b$chem, and $c$chem to $d$chem, $e$chem.'),
        (7, 1, ['a', 'b', 'c', 'd', 'e'],
         'The reaction converts $a$chem, $b$chem, and $c$chem to $d$chem, $e$chem but not $x$chem.'),
        (8, 1, [],
         'The reaction converts $a$chem notatransition $c$chem.'),
        (9, 1, [],
         'The reaction xxxxxxx $a$chem to $c$chem.'),
    ]
    for test_index, pattern_num, output, sent in tests:
        stemmed_sent = ' '.join(map(lambda word: word if word.startswith(
            '$') else stemmer.stem(word), sent.split()))
        fail = True
        grouped_sent = group_list(stemmed_sent)
        for pattern in patterns[pattern_num]:
            result = flatten_list_tuple(pattern.findall(grouped_sent))
            expanded_result = [x for y in result for x in re.findall(chem, y)]
            if expanded_result == output:
                print 'Test %s passed.' % test_index
                fail = False
                break
        if fail:
            print grouped_sent
            print 'TEST %s FAILED!' % test_index


if __name__ == '__main__':
    main()
