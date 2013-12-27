import os
import sys
import json
import time
import shlex
import cirpy
import urllib2
import subprocess
import parse_utils
from utils import *
from random import choice
from indigo.indigo import *
from indigo.indigo_inchi import *

INDIGO = Indigo()
INDIGO_INCHI = IndigoInchi(INDIGO)
MAP_PATH = '../data/inchi_map.json'
if os.path.exists(MAP_PATH):
    CHEM_INCHI_MAP = json.load(open(MAP_PATH))
    print 'Loaded existing inchi map.'
else:
    CHEM_INCHI_MAP = {}
    print 'No inchi map found. Starting from scratch.'


def get_canonical_inchi(chem):
    return canonicalize_inchi(get_inchi(chem))


def canonicalize_inchi(inchi):
    """Returns canonicalized version of inchi used in our db based on indigoinchi"""
    if inchi is None:
        return inchi
    try:
        o = INDIGO_INCHI.loadMolecule(str(inchi))
        o.clearAlleneCenters()
        o.clearCisTrans()
        o.clearStereocenters()
        return INDIGO_INCHI.getInchi(o)
    except UnicodeEncodeError:
        return inchi
    except IndigoException:
        return inchi


def verify_canonicalize():
    """Verify that canonicalize_inchi() has same output as java version"""
    inchis = [x for x in CHEM_INCHI.values() if x is not None]
    os.chdir('/Users/jeff/Documents/synthesis/Act/Installer')
    i, j = 0, 0
    while i < 100:
        inchi = choice(inchis)
        command = str(
            'java -jar mongoactsynth.jar -quiet -config ../MongoActSynth/war/config.xml -exec CONSISTENT_INCHI -target "%s"' % inchi)
        mine = canonicalize_inchi(inchi)
        theirs = subprocess.check_output(shlex.split(command)).strip()
        if mine != theirs:
            j += 1
        i += 1
    print '%s/%s did not match' % (j, i)


def get_inchi(chem):
    """
    Returns the stdinchi of the chem via CHEM_INCHI_MAP.

    Use cirpy if lookup fails and update the map.
    """
    chem = chem.lower()
    if chem in CHEM_INCHI_MAP:
        return CHEM_INCHI_MAP[chem]
    else:
        inchi = query_inchi(chem)
        CHEM_INCHI_MAP[chem] = inchi
        return inchi


def save_map():
    json.dump(CHEM_INCHI_MAP, open(MAP_PATH, 'wb'), indent=2)
    print 'Inchi map saved successfully.'


def query_inchi(chem):
    """Returns the stdinchi of the chem via cirpy"""
    print 'Query for inchi'
    try:
        return cirpy.resolve(chem, 'stdinchi')
    except urllib2.URLError:
        print 'Sleeping for inchikey'
        time.sleep(1)
        query_inchi(chem)


def generate_map(port=30000):
    """
    Generates or adds to an existing map by walking through a chemical db.

    If a chemical already exists in the map, cirpy.resolve is not called.
    """
    try:
        chemicals = Connection(
            'pathway.berkeley.edu', port).actv01['chemicals']
        chem_iter = chemicals.find(timeout=False)
        bar, i = pbar(chemicals.find().count()), 0
        bar.start()
        for chem in chem_iter:
            for name in parse_utils.grab_names(chem['names']):
                get_inchi(name)
            i += 1
            bar.update(i)
        bar.finish()
        save_map()
    except Exception as e:
        print e
        save_map()


def verify_map():
    """
    Randomly selects entries from inchi map and checks if it is valid.

    This is a tool to verify the correctness of the map and estimate
    whether it needs to be updated via update_map().
    This method does not modify the map in any way.
    """
    inchi_map = CHEM_INCHI_MAP.items()
    i, j = 0, 0
    while True:
        name, inchi = choice(inchi_map)
        actual_inchi = query_inchi(name)
        if actual_inchi != inchi:
            j += 1
            data = (name, inchi, actual_inchi)
            print 'INVALID\n\tname:%s\n\tinchi:\t%s\n\tactual:\t%s' % data
        i += 1
        if i % 50 == 0:
            print '%s/%s invalid inchis' % (j, i)


def update_map():
    """
    Updates every entry in CHEM_INCHI_MAP by performing a new query.

    This differs from generate_map in that a query is performed for every
    chemical and old inchis are overwritten with the new ones.
    After processing, the map is now fully synchronized with cirpy.
    """
    bar, i = pbar(len(CHEM_INCHI_MAP)), 0
    bar.start()
    for name, inchi in CHEM_INCHI_MAP.items():
        actual_inchi = query_inchi(name)
        CHEM_INCHI_MAP[name] = actual_inchi
        if actual_inchi != inchi:
            print 'updated'
        i += 1
        bar.update(i)
    bar.finish()
    save_map()
    print 'Inchi map update complete. It is now in sync with cirpy.'


def add_map(chemicals):
    bar, i = pbar(len(chemicals)), 0
    bar.start()
    for chem in chemicals:
        get_inchi(chem)
        i += 1
        bar.update(i)
    bar.finish()
    save_map()


def main():
    if len(sys.argv) == 2:
        this_file, command = sys.argv
        if command == 'generate':
            generate_map()
            return
        elif command == 'verify':
            verify_map()
            return
        elif command == 'update':
            update_map()
            return
    elif len(sys.argv) == 3:
        this_file, command, names_file = sys.argv
        if command == 'add':
            chemicals = json.load(open(names_file))
            if type(chemicals) == type([]):
                add_map(chemicals)
                return
    else:
        print 'Invalid arguments. Usage: python inchi_map.py [generate, verify, update] or python inchi_map.py [add] [filename of chemicals]'
        return


if __name__ == '__main__':
    main()
