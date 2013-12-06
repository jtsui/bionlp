import nltk
import string
import itertools
from collections import defaultdict
import unicodedata
from utils import *
import nltk.data


PUNCT = set(string.punctuation.replace('-', '').replace('+', '')
            .replace('[', '').replace(']', '').replace('(', '').replace(')', ''))

# Pubmed helper methods

sent_detector = None


def cleanAbstract(paper):
    """returns cleaned and normalized pubmed abstract text"""
    article = paper.get('MedlineCitation', {}).get('Article', {})
    abstract = article.get('Abstract', {}).get('AbstractText', '')
    if isinstance(abstract, list):
        abstract = ''.join(abstract)
    if isinstance(abstract, unicode):
        abstract = unicodedata.normalize(
            'NFKD', abstract).encode('ascii', 'ignore')
    return abstract


def cleanJournal(paper):
    """returns the cleaned journal of a pubmed paper"""
    journal = paper.get('MedlineCitation', {}).get(
        'Article', {}).get('Journal', {}).get('Title', '')
    if isinstance(journal, list):
        journal = ''.join(journal)
    if isinstance(journal, unicode):
        journal = unicodedata.normalize(
            'NFKD', journal).encode('ascii', 'ignore')
    return journal.lower()


def cleanTitle(paper):
    """returns the cleaned title of a pubmed paper"""
    article = paper.get('MedlineCitation', {}).get('Article', {})
    title = article.get('ArticleTitle', {})
    if isinstance(title, list):
        title = ''.join(title)
    if isinstance(title, unicode):
        title = unicodedata.normalize('NFKD', title).encode('ascii', 'ignore')
    return title.lower()

# Chemical db helper methods


def grab_names(names_entry):
    """Returns all names and synonyms of a chemical"""
    names = []
    names += [x.strip() for x in names_entry.get('brenda', [])]
    names += [x.strip() for x in names_entry.get('synonyms', [])]
    names += [x.strip() for x in flatten_list([x['values'] for x in
                                               names_entry.get('pubchem', [])
                                               if 'values' in x])]
    names += [x.strip() for x in names_entry.get('genbank', [])]
    names += [x.strip() for x in names_entry.get('pubmed', [])]
    return names

# String helper methods


def isNumber(s):
    """returns true if string s is a number"""
    try:
        float(s)
        return True
    except ValueError:
        return False

# Sentence helper methods


def getPos(sents, strip=True):
    """returns tuple of tagged tokens, dictionary of POS tags to words"""
    if strip:
        sents = [stripPunct(s) for s in sents]
    toks = [nltk.word_tokenize(sent) for sent in sents]
    tagged = [nltk.pos_tag(tok) for tok in toks]
    flat = sum(tagged, [])
    posDict = defaultdict(list)
    for wd, pos in flat:
        posDict[pos].append(wd)
    return tagged, posDict

# Paragraph helper methods


def stripPunct(text):
    """Returns string w/o punctuations"""
    newtext = ''
    for ch in text:
        if ch not in PUNCT:
            newtext += ch
        elif ch == '/' or ch == '-' or ch == ':':
            newtext += ' '
    return newtext


def getMispelled(dic, text):
    """returns a list of mispelled words from a paragraph"""
    return [wd for wd in stripPunct(text).split(' ')
            if wd != '' and not dic.check(wd)]


def nextLetter(text, index):
    """return next letter that comes after index. helper for splitSentences"""
    index += 1
    next = ''
    while index < len(text) and text[index] == ' ':
        index += 1
        next = text[index]
    return next


def splitSentences(text):
    """
    Splits a paragraph into sentences.

    Deal with cases like:
    The isopenicillin N synthase (cyclase) of Streptomyces lactamdurans (syn.
    Nocardia lactamdurans) has been purified to near homogeneity as judged by
    SDS-PAGE and isoelectric focusing.
    Some ambiguity when sentence ends in things like "L." referring to a
    chemical name
    """
    split = []
    parens = False
    brackets = False
    part = ''
    for index in range(len(text)):
        i = text[index]
        if i == '[':
            brackets = True
        elif i == ']':
            brackets = False
        elif i == '(':
            parens = True
        elif i == ')':
            parens = False
        elif i == '.' and not parens and not brackets:
            # split the sentence if period does not come after a single letter
            # or if the next word after the period is capitalized.
            if text[index - 2] != ' ' or nextLetter(text, index).isupper():
                split.append(part)
                part = ''
                continue
        part = part + i
    # remove leading space
    for i in range(len(split)):
        if split[i].startswith(' '):
            split[i] = split[i][1:]
    return split


def getTree(chunker, tagged):
    """returns a tree given list of tagged words using the chunker"""
    if tagged == []:
        return None
    (words, tags) = zip(*tagged)
    chunks = chunker.tag(tags)
    wtc = itertools.izip(words, chunks)
    lines = [' '.join([w, t, c]) for (w, (t, c)) in wtc if c]
    tree = nltk.chunk.conllstr2tree('\n'.join(lines))
    return tree


def shake(tree, fruit):
    """returns a list of all branches of the tree with the POS tag type"""
    leaves = []
    for subtree in tree.subtrees(filter=lambda t: t.node == fruit):
        leaves.append(subtree.leaves())
    return leaves


ABERRANT_PLURAL_MAP = {
    'appendix': 'appendices',
    'barracks': 'barracks',
    'cactus': 'cacti',
    'child': 'children',
    'criterion': 'criteria',
    'deer': 'deer',
    'echo': 'echoes',
    'elf': 'elves',
    'embargo': 'embargoes',
    'focus': 'foci',
    'fungus': 'fungi',
    'goose': 'geese',
    'hero': 'heroes',
    'hoof': 'hooves',
    'index': 'indices',
    'knife': 'knives',
    'leaf': 'leaves',
    'life': 'lives',
    'man': 'men',
    'mouse': 'mice',
    'nucleus': 'nuclei',
    'person': 'people',
    'phenomenon': 'phenomena',
    'potato': 'potatoes',
    'self': 'selves',
    'syllabus': 'syllabi',
    'tomato': 'tomatoes',
    'torpedo': 'torpedoes',
    'veto': 'vetoes',
    'woman': 'women',
}

VOWELS = set('aeiou')


def pluralize(singular):
    """
    from http://code.activestate.com/recipes/577781-pluralize-word-convert-singular-word-to-its-plural/
    Return plural form of given lowercase singular word (English only). Based on
    ActiveState recipe http://code.activestate.com/recipes/413172/
    >>> pluralize('')
    ''
    >>> pluralize('goose')
    'geese'
    >>> pluralize('dolly')
    'dollies'
    >>> pluralize('genius')
    'genii'
    >>> pluralize('jones')
    'joneses'
    >>> pluralize('pass')
    'passes'
    >>> pluralize('zero')
    'zeros'
    >>> pluralize('casino')
    'casinos'
    >>> pluralize('hero')
    'heroes'
    >>> pluralize('church')
    'churches'
    >>> pluralize('x')
    'xs'
    >>> pluralize('car')
    'cars'
    """
    if not singular:
        return ''
    plural = ABERRANT_PLURAL_MAP.get(singular)
    if plural:
        return plural
    root = singular
    try:
        if singular[-1] == 'y' and singular[-2] not in VOWELS:
            root = singular[:-1]
            suffix = 'ies'
        elif singular[-1] == 's':
            if singular[-2] in VOWELS:
                if singular[-3:] == 'ius':
                    root = singular[:-2]
                    suffix = 'i'
                else:
                    root = singular[:-1]
                    suffix = 'ses'
            else:
                suffix = 'es'
        elif singular[-2:] in ('ch', 'sh'):
            suffix = 'es'
        else:
            suffix = 's'
    except IndexError:
        suffix = 's'
    plural = root + suffix
    return plural


def main():
    """Put test code here"""
    pass

if __name__ == '__main__':
    main()
