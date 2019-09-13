import os

RESULTS_PATH = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "results")
)
RESOURCE_PATH = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "resources")
)
if not os.path.exists(RESULTS_PATH):
    raise IOError(f"Path does not exist: {RESULTS_PATH}")


PARAMETER_JSON = os.path.join(RESULTS_PATH, "brenda_parameters.json")  # created via parameters.py
PARAMETER_SBML = os.path.join(RESULTS_PATH, "brenda_parameters.xml")   # created via parameters_sbml.py
