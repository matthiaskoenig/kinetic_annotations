"""
Using brendapy to parse all kinetic parameters.
Stores resulting information as json.
"""
import os
from collections import Counter
from copy import deepcopy
import json
import logging
from pprint import pprint

from brendapy import BrendaParser

BRENDA_PARSER = BrendaParser()


def parse_all_parameters():
    """ Parses all KM, KI, SA, TN and KKM parameters from BRENDA.

    Data is mapped on ontologies in the process. Unmapped information
    is tracked.

    :return:
    """
    # Keep track of information which cannot be mapped to ontologies
    unmapped_organisms = Counter()
    unmapped_substances = Counter()
    unmapped_tissues = Counter()

    def find_mapped_substances(items):
        """ Find entries with substance information mapped to chebi."""
        mapped = []
        unmapped = []
        for item in items:
            if "substrate" in item:
                if "chebi" in item:
                    mapped.append(item)
                else:
                    unmapped.append(item['substrate'])

        return mapped, unmapped

    results = {}
    for ec in BRENDA_PARSER.keys():
        proteins = BRENDA_PARSER.get_proteins(ec)
        for pkey, p in proteins.items():

            # unmapped organisms
            if not p.taxonomy:
                unmapped_organisms.update([p.organism])

            # unmapped tissues
            if p.ST:
                unmapped_tissues.update(
                    [item['data'] for item in p.ST if "bto" not in item]
                )

            # resolve parameters (track unmapped substances
            data = {}
            if p.KM:
                mapped, unmapped = find_mapped_substances(p.KM)
                unmapped_substances.update(unmapped)
                data['KM'] = mapped
            if p.KI:
                mapped, unmapped = find_mapped_substances(p.KI)
                unmapped_substances.update(unmapped)
                data['KI'] = mapped
            if p.TN:
                mapped, unmapped = find_mapped_substances(p.TN)
                unmapped_substances.update(unmapped)
                data['TN'] = mapped
            if p.KKM:
                data['KKM'] = p.KKM
            if p.SA:
                data['SA'] = p.SA

            for key in ['KM', 'KI', 'TN', 'KKM', 'SA']:
                if key in data:
                    # deepcopy all data which have values
                    items = []
                    for item in data[key]:
                        if "value" in item:
                            # resolve pubmeds
                            pubmeds = []
                            info = deepcopy(item)
                            for ref_id in info["refs"]:
                                if ref_id in p.references:
                                    pmid = p.references[ref_id].get("pubmed")
                                    if pmid:
                                        pubmeds.append(pmid)
                            info['pubmeds'] = pubmeds
                            items.append(info)

                    if items:
                        data[key] = items
                    else:
                        del data[key]

            # parameter for export (full annotation)
            if p.taxonomy and data:
                ec_str = p.ec.replace('.', '_')
                eid = f"EC{ec_str}__PR{p.protein_id}"
                entry = {
                    'ec': p.ec,
                    'protein_id': p.protein_id,
                    'organism': p.organism,
                    'taxonomy': p.taxonomy,
                    'tissues': p.tissues,
                }
                if p.uniprot:
                    entry['uniprot'] = p.uniprot
                entry['data'] = data
                results[eid] = entry

    return results, unmapped_organisms, unmapped_substances, unmapped_tissues


def _serialize_json(data, path):
    """Serialize dict to json."""

    def set_default(obj):
        if isinstance(obj, set):
            return list(obj)
        raise TypeError

    with open(path, 'w') as fp:
        json.dump(data, fp, indent=2, default=set_default)


def _serialize_counter(counter, path):
    """Serialize Counter to text file."""
    with open(path, "w") as fout:
        for key, count in counter.most_common():
            fout.write(f"{count} : {key}\n")


def parse_kinetic_parameters(results_path, parameters_path):
    print("-" * 80)
    print("Parsing kinetic parameters")
    print("-" * 80)

    results, unmapped_organisms, unmapped_substances, unmapped_tissues = parse_all_parameters()

    _serialize_json(
        results,
        path=parameters_path
    )
    _serialize_counter(
        unmapped_organisms,
        path=os.path.join(results_path, "unmapped_organisms.txt")
    )
    _serialize_counter(
        unmapped_tissues,
        path=os.path.join(results_path, "unmapped_tissues.txt")
    )
    _serialize_counter(
        unmapped_substances,
        path=os.path.join(results_path, "unmapped_substances.txt")
    )


if __name__ == "__main__":
    from settings import RESULTS_PATH, PARAMETER_JSON
    parse_kinetic_parameters(results_path=RESULTS_PATH, parameters_path=PARAMETER_JSON)
