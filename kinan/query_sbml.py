"""
Loads BRENDA parameters.json (created with `parameters.py`) and writes
all parameters as SBML.

- units are annotated as SBML units
- tissues are annotated via BRENDA Tissue Ontology (bto)
- organisms are annotated via NCBI Taxononmy (taxonomy)
- compounds (for KM, KI and TN) are annotated to ChEBI (chebi)
- parameters are annotated to SBO terms (sbo)
- literature is annotated via pubmed ids

All annotations use BQB.IS qualifiers with exception of pubmed ids which use BQB.IS_DESCRIBED_BY.
"""
import os
import logging
import json
from copy import deepcopy
from sbmlutils.units import *
from sbmlutils.factory import *
from sbmlutils.modelcreator import templates
from sbmlutils.annotation.miriam import *
from sbmlutils.annotation.sbo import *

from settings import RESULTS_PATH, RESOURCE_PATH

# ---------------------------------------------------------------------------------------------------------------------
mid = 'brenda_query'
version = 1
notes = Notes([
    """
    <h1>Model for querying kinetic parameters</h1>
    <h2>Description</h2>
    
    <p>This file contains a minimal glycolysis model for querying kinetic parameters.</p>
    
    <p>
    The SBML was created using
    <a href="https://github.com/matthiaskoenig/sbmlutils">sbmlutils</a>
    (<a href="https://doi.org/10.5281/zenodo.597149">https://doi.org/10.5281/zenodo.597149</a>).
    </p>
    """,
    templates.terms_of_use
])
creators = [
    Creator(familyName='Koenig',
            givenName='Matthias',
            email='koenigmx@hu-berlin.de',
            organization='Humboldt-University Berlin, Institute for Theoretical Biology',
            site="https://livermetabolism.com")
]
# FIXME: support model annotations
# annotations

# ---------------------------------------------------------------------------------------------------------------------
# Units
# ---------------------------------------------------------------------------------------------------------------------
model_units = ModelUnits(time=UNIT_s, extent=UNIT_mmole, substance=UNIT_mmole,
                         length=UNIT_m, area=UNIT_m2, volume=UNIT_KIND_LITRE)
units = [
    UNIT_s,
    UNIT_mmole,
    UNIT_m,
    UNIT_m2,
    Unit('mM', [(UNIT_KIND_MOLE, 1, -3, 1.0), (UNIT_KIND_LITRE, -1.0)]),
    Unit('per_s', [(UNIT_KIND_SECOND, -1.0, 0, 1)]),
    Unit('per_mM_s', [(UNIT_KIND_MOLE, 1, -3, -1.0), (UNIT_KIND_LITRE, 1.0), (UNIT_KIND_SECOND, -1.0, 0, 1)]),
]

# -----------------------------------------------------------------------------
# Compartments
# -----------------------------------------------------------------------------
compartments = [
    Compartment('ext', name='external compartment (blood)', value=1.0, unit=UNIT_KIND_LITRE,
                constant=True),
    Compartment('pancreas', name='pancreas', value=1.0, unit=UNIT_KIND_LITRE,
                constant=True),
    Compartment('cyto', name='pancreas cytosol', value=1.0, unit=UNIT_KIND_LITRE,
                constant=True),
    Compartment('pm', name='plasma membrane', spatialDimensions=2, value=1.0, unit="m2",
                constant=True),
]

# -----------------------------------------------------------------------------
# Species
# -----------------------------------------------------------------------------

species = [
    Species('glc_ext', name="D-glucose", compartment="ext",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False,
            constant=True, boundaryCondition=True),
    Species('lac_ext', name="lactate", compartment="ext",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False,
            constant=True, boundaryCondition=True),
    Species('adp', name="ADP", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('atp', name="ATP", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('bpg13', name="1,3-biphosphoglycerate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('dhap', name="dihydroxyacetone phosphate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('f26bp', name="d-fructose 2,6-bisphosphate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('f6p', name="fructose 6-phosphate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('fbp', name="fructose 1,6-bisphosphate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('g6p', name="glucose 6-phosphate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('glc', name="d-glucose", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('grap', name="glyceraldehyde 3-phosphate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('nad', name="NAD", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('nadh', name="NADH", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('lac', name="lactate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('phos', name="phosphate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('pep', name="phosphoenolpyruvate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('pg2', name="2-phosphoglycerate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('pg3', name="3-phosphoglycerate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('pyr', name="pyruvate", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False),
    Species('h', name="hydrogen ion", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False,
            constant=True, boundaryCondition=True),
    Species('h2o', name="water", compartment="cyto",
            initialConcentration=1.0, substanceUnit='mmole', hasOnlySubstanceUnits=False,
            constant=True, boundaryCondition=True),
]

# -----------------------------------------------------------------------------
# Reactions
# -----------------------------------------------------------------------------
reactions = [

    # ------------------------------------------
    # Glycolysis
    # ------------------------------------------
    Reaction(
        sid="GLCt1r",
        name="Glucose Transport (Uniport)",
        equation="glc_ext <-> glc",
        compartment='pm'
    ),

    Reaction(
        sid="HEX1",
        name="Hexokinase (D-Glucose:ATP) (HEX1)",
        equation="glc + atp <-> g6p + adp + h [g6p]",
        compartment='cyto'
    ),

    Reaction(
        sid="PGI",
        name="Glucose-6-Phosphate Isomerase",
        equation="g6p <-> f6p",
        compartment='cyto'
    ),

    Reaction(
        sid="PFK",
        name="Phosphofructokinase",
        equation="f6p + atp <-> adp + fbp + h [pyr, f26bp]",
        compartment='cyto'
    ),

    Reaction(
        sid="FBA",
        name="Fructose-Bisphosphate Aldolase",
        equation="fbp <-> dhap + grap",
        compartment='cyto'
    ),

    Reaction(
        sid="TPI",
        name="Triose-Phosphate Isomerase",
        equation="dhap <-> grap",
        compartment='cyto'
    ),

    Reaction(
        sid="GAPD",
        name="Glyceraldehyde-3-Phosphate Dehydrogenase",
        equation="grap + phos + nad <-> bpg13 + nadh + h [atp]",
        compartment='cyto'
    ),

    Reaction(
        sid="PGK",
        name="Phosphoglycerate Kinase",
        equation="adp + bpg13 <-> atp + pg3",
        compartment='cyto'
    ),

    Reaction(
        sid="PGM",
        name="Phosphoglycerate Mutase",
        equation="pg3 <-> pg2",
        compartment='cyto'
    ),

    Reaction(
        sid="ENO",
        name="Enolase",
        equation="pg2 <-> h2o + pep",
        compartment='cyto'
    ),

    Reaction(
        sid="PK",
        name="Pyruvate Kinase",
        equation="pep + adp + h <-> pyr + atp",
        compartment='cyto'
    ),

    Reaction(
        sid="LDH_L",
        name="L-Lactate Dehydrogenase",
        equation="pyr + nadh + h <-> lac + nad",
        compartment='cyto'
    ),

    Reaction(
        sid="L_LACt2r",
        name="L-Lactate Reversible Transport",
        equation="lac_ext <-> lac",
        compartment='pm'
    ),

    Reaction(
        sid="ATPuse",
        name="ADP consumption",
        equation="atp + h2o -> adp + phos + h",
        compartment='cyto'
    ),
]

# ------------------------------------------------
# Annotations
# ------------------------------------------------
import os
with open(os.path.join(RESOURCE_PATH, 'mass_charge.json'), 'r') as f:
    mass_charge_data = json.load(f)
with open(os.path.join(RESOURCE_PATH, 'species_annotations.json'), 'r') as f:
    species_annotations_data = json.load(f)
with open(os.path.join(RESOURCE_PATH, 'reactions_annotations.json'), 'r') as f:
    reactions_annotations_data = json.load(f)


for s in species:
    sid = s.sid

    # mass charge balance
    item = mass_charge_data.get(sid, None)
    if not item:
        logging.error(f"No mass_charge data for metabolite: {sid}")
        continue
    charge = item['charge']
    formula = item['chemicalFormula']
    if charge is not None:
        s.charge = charge
    if formula is not None:
        s.chemicalFormula = formula

    # annotations
    if s.annotations is None:
        s.annotations = []
    for collection, annotations in species_annotations_data[sid].items():
        for a in annotations:
            s.annotations.append(
                (BQB(a[0]), a[1])
            )

for r in reactions:
    sid = r.sid

    # annotations
    if r.annotations is None:
        r.annotations = []

    reaction_annotations = reactions_annotations_data.get(sid, None)
    if not reaction_annotations:
        logging.error(f"No reaction_annotations for reaction: {sid}")
        continue
    for collection, annotations in reaction_annotations.items():
        for a in annotations:
            r.annotations.append(
                (BQB(a[0]), a[1])
            )


if __name__ == "__main__":
    from sbmlutils.modelcreator.creator import create_model
    create_model(modules=['query_sbml'],
                 target_dir=RESULTS_PATH,
                 filename="brenda_query.xml",
                 create_report=True,
                 validate=True)

    import os
    import cobra
    from cobra.io import read_sbml_model
    model_path = os.path.join(RESULTS_PATH, "brenda_query.xml")
    print(f"*** MASS/CHARGE Balance: {model_path}")
    print('-' * 80)

    model = read_sbml_model(model_path)

    for r in model.reactions:
        mass_balance = r.check_mass_balance()
        if len(mass_balance) > 0:
            print(r)
            print(mass_balance)
            print('-' * 80)

