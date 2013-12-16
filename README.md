# bionlp

This project is a regular expression based information extraction system that extracts chemical reactions from biological publications. The code here depends on several datasets and REST APIs that are not publicly available. Examples of the structure and contents of these datasets are included in the report.

Review the report for project overview and evaluation results: https://docs.google.com/document/d/110WzRsLv0m4Aof3l3V1BoJllMLr_t6i6ML1tuyv46qU/edit


#### Where to start
* train.py - methods to generate the training set
* patterns.py - the set of regular expression patterns along with basic test cases
* evaluate.py - applies the patterns on the training set and computes statistics

#### Utilities
* chem_canonicalizer.py - uses Indigo to convert between InChI and SMILES notation
* smiles_map.py - a cache of chemical names to SMILES, generated via CIR
* smiles_inchi.py - a cache of chemical names to InChI, generated via CIR
* chemtagger.py - simple library to interface with ChemicalTagger API (custom API wrapper around ChemicalTagger library http://chemicaltagger.ch.cam.ac.uk/)
