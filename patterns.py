import chemtagger
from utils import *
import regex as re
from collections import defaultdict
from nltk.stem.lancaster import LancasterStemmer


STEMMER = LancasterStemmer()
CHEM = r'\$([^\$.]*?)\$chem'
CHEM_LIST = r'{([^{}.]*?)}chem_list'
LIST_PATTERN = re.compile(
    r'(%s)(\(?,?\s*(?:and)?(?:with)?\s*%s\)?)*' % (CHEM, CHEM))
PATTERNS = defaultdict(list)

# pattern0: chem_list <--> chem_list
PATTERN0 = [r'%s[\<\-\>\s]+%s']

TRANS = 'from|to|into|by|are|yield'
# pattern 1 : [trigger1] <0,1> chem_list <0,1> [transition] chem_list
TRIG1 = 'phosphoryl|condens|hydrolys|metabol|reduc|convert|produc|form|oxid|transform|bioconvert|synthes|react|interconvert'
PATTERN1 = [r'(?:%s) %s (?:%s) %s',
            r'(?:%s) \w* %s (?:%s) %s',
            r'(?:%s) \w* %s \w* (?:%s) %s',
            r'(?:%s) %s \w* (?:%s) %s',
            ]
# pattern2 : chem_list <0,1> [trigger2] <0,1> [transition] <0,1> chem_list
TRIG2 = 'convert|oxid|produc|interconvert'
PATTERN2 = [r'%s (?:%s) (?:%s) %s',
            r'%s \w* (?:%s) (?:%s) %s',
            r'%s (?:%s) \w* (?:%s) %s',
            r'%s (?:%s) (?:%s) \w* %s',
            r'%s \w* (?:%s) \w* (?:%s) %s',
            r'%s (?:%s) \w* (?:%s) \w* %s',
            r'%s \w* (?:%s) (?:%s) \w* %s',
            r'%s \w* (?:%s) \w* (?:%s) \w* %s',
            ]
# pattern3: chem_list [trigger3] <0,1> chem_list
TRIG3 = 'yield'
PATTERN3 = [r'%s (?:%s) %s',
            r'%s (?:%s) \w* %s',
            ]
# pattern4: [trigger4] <0,1> chem_lists
TRIG4 = 'convert|interconvert'
PATTERN4 = [r'(?:%s) of %s',
            ]
# pattern5: chem_list [transition5] <0,1> [trigger5] <0,1> chem_list
TRIG5 = 'produc|metabolit'
TRANS5 = 'is|ar'
PATTERN5 = [r'%s (?:%s) (?:%s) %s',
            r'%s (?:%s) \w* (?:%s) %s',
            r'%s (?:%s) \w* (?:%s) \w* %s',
            r'%s (?:%s) (?:%s) \w* %s',
            ]


def expand_patterns():
    global PATTERNS
    if PATTERNS:
        return PATTERNS
    for pattern in PATTERN0:
        PATTERNS[0].append(
            re.compile(pattern % (CHEM_LIST, CHEM_LIST), re.IGNORECASE))
    for pattern in PATTERN1:
        PATTERNS[1].append(
            re.compile(pattern % (TRIG1, CHEM_LIST, TRANS, CHEM_LIST), re.IGNORECASE))
    for pattern in PATTERN2:
        PATTERNS[2].append(
            re.compile(pattern % (CHEM_LIST, TRIG2, TRANS, CHEM_LIST), re.IGNORECASE))
    for pattern in PATTERN3:
        PATTERNS[3].append(
            re.compile(pattern % (CHEM_LIST, TRIG3, CHEM_LIST), re.IGNORECASE))
    for pattern in PATTERN4:
        PATTERNS[4].append(
            re.compile(pattern % (TRIG4, CHEM_LIST), re.IGNORECASE))
    for pattern in PATTERN5:
        PATTERNS[5].append(
            re.compile(pattern % (CHEM_LIST, TRANS5, TRIG5, CHEM_LIST), re.IGNORECASE))
    return PATTERNS


def group_list(sentence):
    return re.sub(LIST_PATTERN, '{\g<0>}chem_list', sentence)


def expand_chems(matches):
    if not matches:
        return []
    if not isinstance(matches[0], tuple):
        matches = [tuple(matches)]
    result = []
    for group in matches:
        expanded = [c for rxn in group for c in re.findall(CHEM, rxn)]
        result.append(expanded)
    return result


def sanitize_chemicals(chemicals):
    chemicals = set(chemicals)
    return [x for x in chemicals if len(x) > 1]


def extract(sid, sentence):
    reactants = []
    chemicals = chemtagger.get_compounds(sid, sentence)
    if chemicals is None:
        return reactants
    chemicals = sanitize_chemicals(chemicals)
    chems = [y for x in chemicals for y in x.split()]
    stemmed_sentence = ' '.join([x if x in chems else STEMMER.stem(x)
                                 for x in sentence.split()])
    tagged_sentence = ' %s ' % stemmed_sentence
    for chem in sorted(chemicals, key=len, reverse=True):
        tagged_sentence = tagged_sentence.replace(
            ' %s ' % chem, ' $%s$chem ' % chem)
    tagged_sentence = tagged_sentence.strip()
    grouped_sentence = group_list(tagged_sentence)
    for pattern_id in expand_patterns():
        for pattern in expand_patterns()[pattern_id]:
            groups = pattern.findall(grouped_sentence, overlapped=True)
            matches = expand_chems(groups)
            for match in matches:
                if match and len(match) > 1:
                    reactants.append((pattern_id, match))
    return reactants


def main():
    tests = [
        [
            0,
            [(1, ['serine', 'glycine'])],
            '10347152-2',
            'As for flux through serine hydroxymethyltransferase and GCS, the conversion of serine to glycine occurred fairly rapidly, followed by GCS-mediated slow decarboxylation of the accumulated glycine',
        ],
        [
            1,
            [(1, ['3-phosphohydroxypyruvate', 'l-phosphoserine'])],
            '10024454-0',
            'Phosphoserine aminotransferase (PSAT; EC 2.6.1.52), a member of subgroup IV of the aminotransferases, catalyses the conversion of 3-phosphohydroxypyruvate to l-phosphoserine'
        ],
        [
            2,
            [(2, ['methionine', 'methanethiol'])],
            '10482527-7',
            'We therefore propose that in P. putida methionine is converted to methanethiol and then oxidized to methanesulfonate',
        ],
        [
            3,
            [(2, ['ornithine', 'N-alpha-acetylornithine']),
             (5, ['ornithine', 'N-alpha-acetylornithine'])],
            '10692366-1',
            'Only two exceptions had been reported-the Enterobacteriaceae and Myxococcus xanthus (members of the gamma and delta groups of the class Proteobacteria, respectively)-in which ornithine is produced from N-alpha-acetylornithine by a deacylase, acetylornithinase (AOase) (argE encoded)',
        ],
        [
            4,
            [(4, ['l-lysine', 'l-beta-lysine'])],
            '17944492-3',
            'The energetics is here reported for the action of lysine 2,3-aminomutase (LAM), which catalyzes the interconversion of l-lysine and l-beta-lysine',
        ],
        [
            5,
            [(4, ['estrone', 'estradiol'])],
            '8013376-0',
            'Estradiol 17 beta-hydroxysteroid dehydrogenase (17 beta HSD) mediates the interconversion of estrone and estradiol in endocrine-responsive tissues such as the breast',
        ],
        [
            6,
            [(5, ['5-Hydroxymethyltryptophan', '5-hydroxy-4-methyltryptophan', '5-methyltryptophan'])],
            '10587452-4',
            '5-Hydroxymethyltryptophan and 5-hydroxy-4-methyltryptophan are the products from 5-methyltryptophan',
        ],
        [
            7,
            [(6, ['maltose', 'trehalose'])],
            '18505459-2',
            'We show that TreS from Mycobacterium smegmatis, as well as recombinant TreS produced in Escherichia coli, has amylase activity in addition to the maltose <--> trehalose interconverting activity (referred to as MTase)',
        ],
    ]
    for index, output, sid, sent in tests:
        res = extract(sid, sent)
        if res == output:
            print 'Test %s passed.' % index
        else:
            print 'TEST %s FAILED!\nResult: %s\n' % (index, res)
    chemtagger.save_map()


if __name__ == '__main__':
    main()
