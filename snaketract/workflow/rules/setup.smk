import os
import re
import tempfile
import functools as ft
import itertools as it
import subprocess as sp

from bids.layout import parse_file_entities
from snakebids import bids, generate_inputs, filter_list
import pandas as pd
import numpy as np

from pathlib import Path
from snakeboost import Tar, Pyscript, XvfbRun, PipEnv, Boost, Datalad, Env
import snakeboost.bash as sh
from templateflow import api as tflow
from lib.path_store import PathStore

# def get_labels(label):
#     cli = config.get(label, None)
#     vals = cli if isinstance(cli, int) else cli.split(",") if cli is not None else None
#     return vals

# def get_participants(participant_file):
#     subs = pd.read_csv(participant_file, sep='\t')
#     return subs['participant_id'].map(lambda s: s[4:])

# participant_label = get_labels("participant_label")
# exclude_participant_label = (
#     get_labels("exclude_participant_label") if participant_label is None else None
# )
# if participant_label is None and exclude_participant_label is None:
#     try:
#         participant_label = list(get_participants(
#             Path(config['bids_dir'], 'derivatives', 'snakedwi-0.1.0', 'participants.tsv')
#         ))
#     except FileNotFoundError:
#         pass



tmpdir = eval(workflow.default_resources._args.get("tmpdir"), {"system_tmpdir": tempfile.gettempdir()})

if workflow.run_local:
    workflow.shadow_prefix = os.environ.get("SLURM_TMPDIR")

###
# Input Globals
###
inputs = generate_inputs(
    bids_dir=config['bids_dir'],
    pybids_inputs=config['pybids_inputs'],
    derivatives=config.get("derivatives", False),
    participant_label=config.get('participant_label'),
    exclude_participant_label=config.get('exclude_participant_label'),
    pybidsdb_dir=config.get("pybidsdb_dir"),
)

wildcards = inputs.input_wildcards['preproc_dwi']

###
# Output Globals
###

work = Path(tmpdir) / "snaketract"
shared_work = Path(config['output_dir'])/'work'/'prepdwi_recon'
output = Path(config['output_dir'])
qc = Path(output)/"qc"

DIFFUSION_PARAMS = ["FA", "MD", "RD", "L1"]
VALID_WEIGHTS = DIFFUSION_PARAMS + ["R1"]

# Unique ID for easy naming in temporary files

def _uid(comp = None, entities = None):
    parts = comp.wildcards.values() if comp is not None else []
    wrapped_entities = (f"{{{entity}}}" for entity in entities or [])
    return ".".join(it.chain(parts, wrapped_entities))

def tempout(rulename: str, comp = None, extension = "", *entities):
    uid = _uid(comp, entities)
    return temp(work.joinpath(uid, rulename).with_suffix(extension))

def log(rulename: str, comp = None, *entities):
    uid = _uid(comp, entities)
    return Path("code", "logs", uid, rulename).with_suffix(".log")

def benchmark(rulename: str, comp = None, *entities):
    uid = _uid(comp, entities)
    return Path("code", "benchmarks", uid,  rulename).with_suffix(".tsv")

def resource(path):
    return os.path.join(workflow.basedir, "..", "resources", path)

def script(path):
    return os.path.join(workflow.basedir, "scripts", path)

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


###
# Pipenvs
###

r1_env = PipEnv(
    packages = [
        'numpy',
        '/scratch/knavynde/snakeboost',
        'nibabel',
    ],
    flags = config.get("pip-flags", ""),
    root = work
)
