"""
Resolve additional pubchem information.
"""
from settings import RESULTS_PATH, PARAMETER_JSON
import json
import os
import logging
from kinan.inchi import inchi_to_inchikey
from pprint import pprint
from pubchem import get_pubchem_collections, get_pubchem_synonyms

# ---------------------------------------------------
# Dictionaries for information lookup
# ---------------------------------------------------

from pubchem import COLLECTIONS, MIRIAM_COLLECTIONS

# chebi to inchi dictionary
chebi2inchi = {}
with open(os.path.join("..", "resources", "chebiId_inchi.tsv"), "r") as f:
    for line in f.readlines():
        tokens = [t.strip() for t in line.split("\t")]
        chebi2inchi[tokens[0]] = tokens[1]
    # delete header
    del chebi2inchi["CHEBI_ID"]

# inchi to inchikey dictionary
inchi2inchikey = {}
for inchi in chebi2inchi.values():
    inchikey = inchi_to_inchikey(inchi)
    if inchikey:
        inchi2inchikey[inchi] = inchikey
    else:
        logging.error(f"inchikey could not be generated for: {inchi}")


def get_chebis_from_parameters():
    """ Get all chebis from parameter entries.
    :return: set of chebis
    """
    chebis = set()
    with open(PARAMETER_JSON, "r") as fpars:
        pdata = json.load(fpars)

        for key, pentry in pdata.items():
            for kin_key in ['KM', 'KI', 'TN', 'KKM', 'SA']:

                if kin_key in pentry["data"]:
                    for item in pentry["data"][kin_key]:

                        if "chebi" in item:
                            chebi_id = item["chebi"].split("_")[1]
                            chebis.add(chebi_id)

    return chebis


def get_pubchem_info(chebis, force=False):
    """ Get all pubchem information for given chebi

    :param chebis:
    :return:
    """

    chebi_map = {}
    N_chebis = len(chebis)
    for k, chebi_id in enumerate(sorted(chebis)):
        print(f"{chebi_id}\t{k}/{N_chebis}")

        chebi_path = os.path.join(RESULTS_PATH, "pubchem", f"{chebi_id}.json")

        if not force:
            # only do query if file does not yet exist
            if os.path.exists(chebi_path):
                with open(chebi_path, "r") as f:
                    chebi_map[chebi_id] = json.load(f)
                continue

        map = {'chebi': {chebi_id}}

        # inchi
        if chebi_id not in chebi2inchi:
            logging.warning(f"No inchi for chebi: {chebi_id}")
        else:
            inchi = chebi2inchi[chebi_id]
            map['inchi'] = f"{inchi}"

            # inchikey
            if inchi not in inchi2inchikey:
                logging.error(f"No inchikey for chebi/inchi: {chebi_id}/{inchi}")
            else:
                inchikey = inchi2inchikey[inchi]
                map['inchikey'] = f"{inchikey}"

                # get additional information from pubchem
                synonyms = get_pubchem_synonyms(inchikey)
                if "error" in synonyms:
                    # catch errors of form: {'error': " The query InChIkey
                    # 'YAJCHEVQCOHZDC-QMMNLEPNSA-N' is not present in 'UniChem.'}
                    logging.error(f"{chebi_id}: {synonyms}")
                    synonyms = {}

                for d in synonyms:
                    src_id = d["src_id"]
                    src_compound_id = d["src_compound_id"]

                    # only add collections on identifiers.org
                    if src_id in MIRIAM_COLLECTIONS:
                        collection = MIRIAM_COLLECTIONS[src_id]

                        # fix kegg collection
                        if collection == "kegg":
                            if collection.startswith("C"):
                                collection = "kegg.compound"
                            elif collection.startswith("G"):
                                collection = "kegg.glycan"

                        # store collections as sets (brenda, chebi and others have multiple entries)
                        if collection in map:
                            map[collection].add(src_compound_id)
                        else:
                            map[collection] = {src_compound_id}
        # term fixes
        map["chebi"] = {f"CHEBI:{term}" for term in map["chebi"]}
        # storing results
        _serialize_json(map, chebi_path)
        chebi_map[chebi_id] = map

    return chebi_map


def _serialize_json(data, path):
    """Serialize dict to json."""

    def set_default(obj):
        if isinstance(obj, set):
            return list(obj)
        raise TypeError

    with open(path, 'w') as fp:
        json.dump(data, fp, indent=2, default=set_default)

# get all chebis
chebis = get_chebis_from_parameters()
# chebis = ["5931",]

# get information for chebis
chebi_map = get_pubchem_info(chebis, force=False)
path_pubchem = os.path.join(RESULTS_PATH, "pubchem_mapping.json")
_serialize_json(data=chebi_map, path=path_pubchem)
