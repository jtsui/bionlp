import os
import parse_utils
import json
import cirpy
import urllib2
import time
import sys
from random import choice
from utils import *

MAP_PATH = 'smiles_map.json'
# words that have a smiles but are not chemicals
BLACKLIST_COMMON = set(['heat', 'I', 'II', 'region', 'for', 'of', 'sp',
                       'cds', 'str', 'complete', 'genome', 'draft', 'and',
                       'the', 'str', 'factor', 'RNA', 'leader', 'spectrum',
                       'compound', 'acid', 'or'])


if os.path.exists(MAP_PATH):
    CHEM_SMILES_MAP = json.load(open(MAP_PATH))
    print 'Loaded existing smiles map.'
else:
    CHEM_SMILES_MAP = {}
    print 'No smiles map found. Starting from scratch.'


def get_smiles(chem):
    if chem in BLACKLIST_COMMON or parse_utils.isNumber(chem):
        return None
    if chem in CHEM_SMILES_MAP:
        return CHEM_SMILES_MAP[chem]
    else:
        smiles = query_smiles(chem)
        CHEM_SMILES_MAP[chem] = smiles
        return smiles


def query_smiles(chem):
    print 'Query for smiles'
    try:
        return cirpy.resolve(chem, 'smiles')
    except urllib2.URLError:
        print 'Sleeping for smiles'
        time.sleep(1)
        query_smiles(chem)
    except Exception as e:
        print e
        return None


def save_map():
    json.dump(CHEM_SMILES_MAP, open(MAP_PATH, 'wb'), indent=2)
    print 'Smiles map saved successfully.'


def verify_map():
    """
    Randomly selects entries from smiles map and checks if it is valid.

    This is a tool to verify the correctness of the map and estimate
    whether it needs to be updated via update_map().
    This method does not modify the map in any way.
    """
    smiles_map = CHEM_SMILES_MAP.items()
    i, j = 0, 0
    while True:
        name, smiles = choice(smiles_map)
        actual_smiles = query_smiles(name)
        if actual_smiles != smiles:
            j += 1
            data = (name, smiles, actual_smiles)
            print 'INVALID\n\tname:%s\n\tsmiles:\t%s\n\tactual:\t%s' % data
        i += 1
        if i % 50 == 0:
            print '%s/%s invalid smiles' % (j, i)


def update_map():
    """
    Updates every entry in CHEM_SMILES_MAP by performing a new query.

    A new cirpy query is performed for every chemical and old smiles
    are overwritten with the new ones.
    After processing, the map is now fully synchronized with cirpy.
    """
    bar, i = pbar(len(CHEM_SMILES_MAP)), 0
    bar.start()
    for name, smiles in CHEM_SMILES_MAP.items():
        actual_smiles = query_smiles(name)
        CHEM_SMILES_MAP[name] = actual_smiles
        if actual_smiles != smiles:
            print 'updated'
        i += 1
        bar.update(i)
    bar.finish()
    save_map()
    print 'Smiles map update complete. It is now in sync with cirpy.'


def main():
    if len(sys.argv) != 2:
        print 'Wrong number of arguments. Usage: python smiles_map.py [verify, update]'
        return
    the_file, command = sys.argv
    if command == 'verify':
        verify_map()
    elif command == 'update':
        update_map()
    else:
        print 'Invalid command. Usage: python smiles_map.py [verify, update]'
        return


if __name__ == '__main__':
    main()
