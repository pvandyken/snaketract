import os
import re
import tempfile
import functools as ft

from snakebids import bids, generate_inputs
from snakemake import utils as sutils

from pathlib import Path
from snakeboost import Tar, Pyscript, ScriptDict, XvfbRun, PipEnv, Boost, Datalad, Env
import snakeboost.bash as sh

def get_labels(label):
    cli = config.get(label, None)
    vals = cli.split(",") if cli is not None else None
    if vals is None:
        try:
            with open(f"{label}.txt") as f:
                return f.read().splitlines()
        except FileNotFoundError:
            pass
    return vals

participant_label = get_labels("participant_label")
exclude_participant_label = (
    get_labels("exclude_participant_label") if participant_label is None else None
)


###
# Evaluate env vars passed through CLI
###
def eval_environ(s):
    s = s.strip("'\"") if s else s
    if s and len(s) > 1 and s[0] == "$":
        if s[1] == "$":
            return s[1:]
        return os.environ.get(s[1:])
    return s

tmpdir = eval_environ(workflow.default_resources._args.get("tmpdir"))

if tmpdir is not None:
    workflow.default_resources.set_resource("tmpdir", tmpdir)

workflow.shadow_prefix = eval_environ(workflow.shadow_prefix)

###
# Input Globals
###
inputs = generate_inputs(
    bids_dir=config['bids_dir'],
    pybids_inputs=config['pybids_inputs'],
    derivatives=True,
    participant_label=participant_label,
    exclude_participant_label=exclude_participant_label,
    use_bids_inputs=True,
    pybids_database_dir=config.get("pybids_database_dir"),
    pybids_reset_database=config.get("pybids_reset_database"),
)

wildcards = inputs.input_wildcards['preproc_dwi']

###
# Output Globals
###

work = Path(tmpdir or tempfile.gettempdir()) / "sn_prepdwi_recon"
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

test_env = PipEnv(
    packages = [
        'colorama',
        'numpy',
        'black',
        '../snakeboost/'
    ],
    root = work
)
