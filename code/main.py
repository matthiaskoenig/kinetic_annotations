from settings import RESULTS_PATH, PARAMETER_JSON
from parameters import parse_kinetic_parameters
from sbmlutils.modelcreator.creator import create_model

# parameters json
parse_kinetic_parameters(results_path=RESULTS_PATH, parameters_path=PARAMETER_JSON)

# parameters SBML
create_model(modules=['parameters_sbml'],
             target_dir=RESULTS_PATH,
             create_report=False,
             validate=False)


