import re
import itertools
from collections import defaultdict
from nltk.stem.lancaster import LancasterStemmer


transitions = ['from', 'to', 'into', 'by', 'are', 'yield']

# pattern 1 : [trigger1] <0,1> cmp <0,1> [transition] cmp
patterns1 = [r'%s \$([^\$.]*?)\$cmp %s \$([^\$.]*?)\$cmp',
             r'%s \w* \$([^\$.]*?)\$cmp %s \$([^\$.]*?)\$cmp',
             r'%s \w* \$([^\$.]*?)\$cmp \w* %s \$([^\$.]*?)\$cmp',
             r'%s \$([^\$.]*?)\$cmp \w* %s \$([^\$.]*?)\$cmp',
             ]
triggers1 = ['reduc', 'convert', 'produc', 'form', 'oxid',
             'transform', 'bioconvert', 'synthes', 'react',
             'interconvert']

# pattern2 : cmp <0,1> [trigger2] <0,1> [transition] <0,1> cmp
triggers2 = ['convert', 'oxid', 'produc', 'interconvert']
patterns2 = [r'\$([^\$.]*?)\$cmp %s %s \$([^\$.]*?)\$cmp',
             r'\$([^\$.]*?)\$cmp \w* %s %s \$([^\$.]*?)\$cmp',
             r'\$([^\$.]*?)\$cmp %s \w* %s \$([^\$.]*?)\$cmp',
             r'\$([^\$.]*?)\$cmp %s %s \w* \$([^\$.]*?)\$cmp',
             r'\$([^\$.]*?)\$cmp \w* %s \w* %s \$([^\$.]*?)\$cmp',
             r'\$([^\$.]*?)\$cmp %s \w* %s \w* \$([^\$.]*?)\$cmp',
             r'\$([^\$.]*?)\$cmp \w* %s %s \w* \$([^\$.]*?)\$cmp',
             r'\$([^\$.]*?)\$cmp \w* %s \w* %s \w* \$([^\$.]*?)\$cmp',
             ]

# pattern3: cmp [trigger3] <0,1> cmp
triggers3 = ['yield']
patterns3 = [r'\$([^\$.]*?)\$cmp %s \$([^\$.]*?)\$cmp',
             r'\$([^\$.]*?)\$cmp %s \w* \$([^\$.]*?)\$cmp',
             ]

# pattern4: [trigger4] <0,1> cmp and cmp
triggers4 = ['convert', 'interconvert']
patterns4 = [r'%s \$([^\$.]*?)\$cmp and \$([^\$.]*?)\$cmp',
             r'%s \w* \$([^\$.]*?)\$cmp and \$([^\$.]*?)\$cmp',
             ]

# pattern5: cmp [transition5] <0,1> [trigger5] <0,1> cmp
triggers5 = ['produc', 'metabolit']
transitions5 = ['is', 'ar']
patterns5 = [r'\$([^\$.]*?)\$cmp %s %s \$([^\$.]*?)\$cmp',
             r'\$([^\$.]*?)\$cmp %s \w* %s \$([^\$.]*?)\$cmp',
             r'\$([^\$.]*?)\$cmp %s \w* %s \w* \$([^\$.]*?)\$cmp',
             r'\$([^\$.]*?)\$cmp %s %s \w* \$([^\$.]*?)\$cmp',
             ]


def grow_patterns():
    patterns = defaultdict(list)
    for pattern in patterns1:
        for vals in itertools.product(triggers1, transitions):
            patterns[1].append(pattern % vals)
    for pattern in patterns2:
        for vals in itertools.product(triggers2, transitions):
            patterns[2].append(pattern % vals)
    for pattern in patterns3:
        for vals in triggers3:
            patterns[3].append(pattern % vals)
    for pattern in patterns4:
        for vals in triggers4:
            patterns[4].append(pattern % vals)
        for pattern in patterns5:
            for vals in itertools.product(transitions5, patterns5):
                patterns[5].append(pattern % vals)
    return patterns


def main():
    patterns = grow_patterns()
    stemmer = LancasterStemmer()
    tests = [
        (0, 1, [('cholesterol', 'epoxides')],
         'These results indicate that myeloperoxidase can convert $cholesterol$cmp to $chlorohydrins$cmp and epoxides by a reaction involving HOCl'),
        (1, 2, [('prostaglandins', 'COX-2')],
            'We investigated the possible role of $prostaglandins$cmp produced by $COX-2$cmp in the immunosuppression observed during Trypanosoma cruzi infection'),
        (2, 3, [
            ('5-formyl-6-methoxy-2,3-diphenylbenzofuran (IVa)', '6- and 5-methyl derivatives')],
         'Reduction of 5-methoxy-6-formyl(Ia)- and $5-formyl-6-methoxy-2,3-diphenylbenzofuran (IVa)$cmp yielded $6- and 5-methyl derivatives$cmp Ib and IVb, respectively'),
        (3, 4, [('ester', 'amide')],
         'It results in the mutual conversion of $ester$cmp and $amide$cmp bonds'),
    ]
    for test_index, pattern_num, output, sent in tests:
        stemmed_sent = ' '.join(map(lambda word: word if word.startswith(
            '$') else stemmer.stem(word), sent.split()))
        for pattern in patterns[pattern_num]:
            result = re.findall(pattern, stemmed_sent)
            if result == output:
                print 'Test %s passed.' % test_index

if __name__ == '__main__':
    main()