"""
UniChem

Resolving
https://www.ebi.ac.uk/unichem/info/webservices


Name:  Get src_compound_ids from InChIKey
Description:  Obtain a list of all src_compound_ids (from all sources) which are CURRENTLY assigned to a query InChIKey
Extension of base url:  inchikey
Number of required input parameters:  1
Input:  /InChIKey
Output:  list of two element arrays, containing 'src_compound_id' and 'src_id'.
Example:  https://www.ebi.ac.uk/unichem/rest/inchikey/AAOVKJBEBIDNHE-UHFFFAOYSA-N


"""

import logging
import requests
import os
import json
from pprint import pprint

COLLECTIONS = None
MIRIAM_COLLECTIONS = {
    '2': "drugbank",
    '3': "pdb",
    '5': "pubchem.substance",
    '6': "kegg",  #KEGG LIGAND is a composite DB consisting of COMPOUND, 'GLYCAN, REACTION, RPAIR, RCLASS, and ENZYME DBs, whose 'entries are identified by C, G, R, RP, RC, and EC 'numbers, respectively.
    '7': "chebi",
    '9': "zinc",
    '17': "pharmgkb.drug",
    '18': "hmdb",
    '21': "pubchem.substance",
    '22': "pubchem.compound",
    '25': "lincs.smallmolecule",
    '27': "vmhmetabolite",
    '31': "bindingdb",
    '32': "comptox",
    '33': "lipidmaps",
    '36': "metabolights",
    '37': "brenda.ligand",
    '38': "rhea",
    '41': "swisslipid",
}


def get_pubchem_collections():
    """ Available collections in PubChem.

    Name:  Get all src_ids
    Description:  Obtain all src_ids currently in UniChem
    Extension of base url:  src_ids
    Number of required input parameters:  0
    Input:   - none -
    Output:  list of 'src_id's.
    Example:  https://www.ebi.ac.uk/unichem/rest/src_ids/

    Name:  Get source infomation
    Description:  Obtain all information on a source by querying with a source id (src_id).
    Extension of base url:  sources
    Number of required input parameters:  1
    Input:  /src_id
    Output:  list containing:
    src_id (the src_id for this source),
    src_url (the main home page of the source),
    name (the unique name for the source in UniChem, always lower case),
    name_long (the full name of the source, as defined by the source),
    name_label (A name for the source suitable for use as a 'label' for the source within a web-page. Correct case setting for source, and always less than 30 characters),
    description (a description of the content of the source),
    base_id_url_available (an flag indicating whether this source provides a valid base_id_url for creating cpd-specific links [1=yes, 0=no]).
    base_id_url (the base url for constructing hyperlinks to this source [append an identifier from this source to the end of this url to create a valid url to a specific page for this cpd], unless aux_for_url=1),
    aux_for_url (A flag to indicate whether the aux_src field should be used to create hyperlinks instead of the src_compound_id [1=yes, 0=no]
    Example:  https://www.ebi.ac.uk/unichem/rest/sources/1

    :return:
    """
    logging.info("get_unichem_collections")
    url_source_ids = "https://www.ebi.ac.uk/unichem/rest/src_ids/"
    response = requests.get(url_source_ids)
    src_ids_data = response.json()
    data = {}
    for item in src_ids_data:
        src_id = item['src_id']
        url_src_id = f"https://www.ebi.ac.uk/unichem/rest/sources/{src_id}"
        response = requests.get(url_src_id)
        src_json = response.json()
        data[src_id] = src_json

    return data


COLLECTIONS = get_pubchem_collections()


def get_pubchem_synonyms(inchikey):
    """
    :return:
    """
    url = f"https://www.ebi.ac.uk/unichem/rest/inchikey/{inchikey}"

    response = requests.get(url)
    return response.json()




if __name__ == "__main__":
    from pprint import pprint
    print("Loading chebi information")

    # get_substance_info(brenda_ligand_id=100)
    data = get_pubchem_collections()

