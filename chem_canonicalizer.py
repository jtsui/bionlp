import inchi_map
import smiles_map
from indigo.indigo import *
from indigo.indigo_inchi import *

INDIGO = Indigo()
INDIGO_INCHI = IndigoInchi(INDIGO)


def smiles_to_inchi(smiles):
    o = INDIGO.loadMolecule(smiles)
    o.clearAlleneCenters()
    o.clearCisTrans()
    o.clearStereocenters()
    return INDIGO_INCHI.getInchi(o)


def inchi_to_smiles(inchi):
    return INDIGO_INCHI.loadMolecule(inchi).smiles()


def names_to_structures(names):
    '''
    Converts a list of chemical names to (smiles, canonical inchi, name).
    Names that can not be canonicalized will be filtered out of the list.
    '''
    structs = []
    for name in names:
        smiles = smiles_map.get_smiles(name)
        inchi = inchi_map.get_canonical_inchi(name)
        if not smiles and not inchi:
            continue
        elif not inchi:
            inchi = smiles_to_inchi(smiles)
        elif not smiles:
            smiles = inchi_to_smiles(inchi)
        structs.append((smiles, inchi, name))
    return structs


def names_to_inchi(names):
    inchis = {}
    for name in names:
        inchi = inchi_map.get_canonical_inchi(name)
        inchis[inchi] = name
    return inchis
