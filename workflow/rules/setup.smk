import os
import re
import tempfile
import functools as ft

from bids.layout import parse_file_entities
from snakebids import bids, generate_inputs, filter_list
from snakemake import utils as sutils
from snakemake.exceptions import WorkflowError, IncompleteCheckpointException
from snakemake.io import checkpoint_target
import pandas as pd

from pathlib import Path
from snakeboost import Tar, Pyscript, XvfbRun, PipEnv, Boost, Datalad, Env
import snakeboost.bash as sh
from templateflow import api as tflow
from lib.participants import filter_participants
from lib.path_store import PathStore

def get_labels(label):
    cli = config.get(label, None)
    vals = cli if isinstance(cli, int) else cli.split(",") if cli is not None else None
    return vals

def get_participants(participant_file):
    subs = pd.read_csv(participant_file, sep='\t')
    return subs['participant_id'].map(lambda s: s[4:])

participant_label = get_labels("participant_label")
exclude_participant_label = (
    get_labels("exclude_participant_label") if participant_label is None else None
)
if participant_label is None and exclude_participant_label is None:
    try:
        participant_label = list(get_participants(
            Path(config['bids_dir'], 'derivatives', 'snakedwi-0.1.0', 'participants.tsv')
        ))
    except FileNotFoundError:
        pass



tmpdir = eval(workflow.default_resources._args.get("tmpdir"), {"system_tmpdir": tempfile.gettempdir()})

if workflow.run_local:
    workflow.shadow_prefix = os.environ.get("SLURM_TMPDIR")

###
# Input Globals
###
inputs = generate_inputs(
    bids_dir=config['bids_dir'],
    pybids_inputs=config['pybids_inputs'],
    derivatives=True,
    participant_label=participant_label,
    exclude_participant_label=exclude_participant_label,
    pybids_database_dir=config.get("pybids_database_dir"),
    pybids_reset_database=config.get("pybids_reset_database"),
)
del inputs['surf']

wildcards = inputs.input_wildcards['preproc_dwi']

###
# Output Globals
###

work = Path(tmpdir) / "sn_prepdwi_recon"
shared_work = Path(config['output_dir'])/'work'/'prepdwi_recon'
output = Path(config['output_dir'])/"prepdwi_recon"
qc = Path(output)/"qc"

# Unique ID for easy naming in temporary files
uid = '.'.join(wildcards.values())

def shell_uid(sample):
    return '.'.join(
        re.sub(r'^\{', '{wildcards.', val)
        for val in inputs.input_wildcards[sample].values()
    )

###
# bids Partials
###
bids_output_dwi = ft.partial(bids, root=output, space="T1w", datatype="dwi", **wildcards)
bids_output_anat = ft.partial(bids, root=output, space="orig", datatype="anat", **wildcards)

###
# Utility functions
###
tar = Tar(work)
xvfb_run = XvfbRun()
boost = Boost(work, logger, debug=True)
datalad = Datalad(config['bids_dir'])
env = Env()
pyscript = Pyscript(workflow.basedir)

# Patch shells in workflow, allowing implicit use of boost in every shell: ...
workflow.shellcmd = lambda *cmd: Workflow.shellcmd(workflow, boost(*cmd))
base_shell = shell
class BoostShell(base_shell):
    def __new__(cls, *cmd, iterable=False, read=False, bench_record=None,  **kwargs):
        import inspect
        locals().update(inspect.currentframe().f_back.f_locals)
        formatted = (
            sutils.format(boost(*cmd), stepout=2, **kwargs)
            .replace('{', '{{')
            .replace('}', '}}')
        )
        return base_shell(
            formatted,
            iterable=iterable,
            read=read,
            bench_record=bench_record,
            **kwargs,
        )
shell = BoostShell

###
# Pipenvs
###
wma_env = PipEnv(
    packages = [
        'whitematteranalysis',
        'vtk==8.1.2',
        '/scratch/knavynde/snakeboost'
    ],
    flags = config.get("pip-flags", ""),
    root = work
)

dipy_env = PipEnv(
    packages = [
        'dipy',
        'pandas',
        'numpy',
        'scipy',
        'snakeboost',
        'nibabel',
    ],
    flags = config.get("pip-flags", ""),
    root = work
)

r1_env = PipEnv(
    packages = [
        'numpy',
        '/scratch/knavynde/snakeboost',
        'nibabel',
    ],
    flags = config.get("pip-flags", ""),
    root = work
)

convert_env = PipEnv(
    packages = [
        'dipy',
        'fury',
        'vtk==8.1.2',
        'more-itertools',
        '/scratch/knavynde/snakeboost',
        'networkx',
    ],
    flags = config.get("pip-flags", ""),
    root = work
)

parcellation_env = PipEnv(
    packages = [
        'intersection',
    ],
    flags = config.get("pip-flags", ""),
    root = work,
)

rich_club_env = PipEnv(
    packages = [
        "more-itertools",
        "networkx",
        "numpy",
        "pandas",
        "/scratch/knavynde/snakeboost",
        "plotly",
        "colour",
        "ipywidgets",
        "xarray",
    ],
    flags = config.get("pip-flags", ""),
    root = work,
)

nodal_props_venv = PipEnv(
    packages = [
        "numpy",
        "pandas",
        "networkx",
        "colour",
        "/scratch/knavynde/snakeboost",
        "plotly",
        "attrs",
    ],
    flags = config.get("pip-flags", ""),
    root = work,
)

design_matrix_env = PipEnv(
    packages = [
        "more-itertools",
        "numpy",
        "pandas",
        "scikit-learn",
        "/scratch/knavynde/snakeboost",
    ],
    flags = config.get("pip-flags", ""),
    root = work,
)

test_env = PipEnv(
    packages = [
        'colorama',
        'numpy',
        'black',
        '../snakeboost/'
    ],
    root = work
)
