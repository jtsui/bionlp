import requests
import json
from utils import *
import xmltodict

URL = 'http://pathway.berkeley.edu:27329/tag'


def get_tree(text):
    data = {'paper': text}
    headers = {'Content-type': 'application/json'}
    result = requests.get(URL, data=json.dumps(data), headers=headers)
    try:
        return xmltodict.parse(result.text)
    except Exception:
        return []


def parse_tree(tree, compounds):
    if isinstance(tree, list):
        for element in tree:
            parse_tree(element, compounds)
    if isinstance(tree, dict):
        for key in tree:
            if key == 'OSCAR-CM':
                if isinstance(tree[key], list):
                    compounds.append(' '.join(tree[key]))
                else:
                    compounds.append(tree[key])
            else:
                parse_tree(tree[key], compounds)


def get_compounds(text):
    tree = get_tree(text)
    cmps = []
    parse_tree(tree, cmps)
    return cmps


def main():
    tests = [
        "Aminoimidazole ribonucleotide (AIR) synthetase (PurM) catalyzes the conversion of formylglycinamide ribonucleotide (FGAM) and ATP to AIR, ADP, and P(i), the fifth step in de novo purine biosynthesis",
        "Induction tests with cells grown on limonene revealed that the oxygen consumption rates with limonene-1,2-epoxide, limonene-1,2-diol, 1-hydroxy-2-oxolimonene, and carveol were high",
        "Acetyl-coenzyme A (acetyl-CoA) synthetase (ADP forming) represents a novel enzyme in archaea of acetate formation and energy conservation (acetyl-CoA + ADP + P(i) --> acetate + ATP + CoA)",
        "It was also found that this enzyme catalyzed hydrolysis of such fructooligosaccharides as 1-kestose, nystose, and 1-fructosylnystose by liberating fructose",
        "In the presence of NADH, it could catalyze the stereospecific reduction of racemic acetoin ((3R/3S)- acetoin; K(m) = 4.5 mm, k(cat) = 98,000 min(-)(1)) to (2R,3R)-2,3-butanediol and meso-butanediol, respectively",
        "5-exo-Bromocamphor is readily metabolized by the P-450cam mixed function oxidase to 5-ketocamphor at rates and yields similar to that of the normal substrate, camphor, suggesting abstraction of the endo-hydrogen of 5-exo-bromocamphor and oxygen addition to produce a transient 5-bromo-5-hydroxycamphor intermediate",
        "Lactate monooxygenase catalyzes the oxidation of L-lactate with molecular oxygen to acetate, CO2, and water",
        "Therefore, the data provide firm evidence for the concept that delta\"mu Na+ is the primary driving force for the synthesis of ATP in P. modestum",
        "S-Ribosylhomocysteinase (LuxS) catalyzes the cleavage of the thioether linkage in S-ribosylhomocysteine (SRH) to produce homocysteine (Hcys) and 4,5-dihydroxy-2,3-pentanedione (DPD), the precursor of type II bacterial autoinducer (AI-2)",
    ]
    for test in tests:
        print get_compounds(test)

if __name__ == '__main__':
    main()
