from typing import List
from pathlib import Path
from lib.utils import hash_name


PYTHON_VENV_CREATE_ERR = "[ERROR] (jobid={jobid}): Error creating python environment"
TIMED_OUT_ERR = "[ERROR] (jobid={jobid}): Script timed out when waiting for Python"


class PipEnv:
    def __init__(
        self,
        packages: List[str],
        root: Path,
        flags: str = "",
        requirements: List[str] = [],

    ):
        name = hash_name(str(sorted(packages)) + str(sorted(requirements)))
        self._dir = root/'__snakemake_venvs__'/name
        self._venv_lock = self._dir/".venv_lock"
        self.venv = self._dir/"venv"
        self.bin = self.venv/"root"
        self.python_path = self.venv/"bin"/"python"
        if not flags:
            self._flags = ""
        else:
            self._flags = flags
            if not isinstance(flags, str):
                raise TypeError("flags attribute in PipEnv must be a string")

        self._packages = ' '.join(packages)
        self._requirements = '-r ' + ' -r '.join(requirements) if requirements else ""


    @property
    def get_venv(self):
        install_prefix = f"{self.python_path} -m pip install {self._flags}"
        install_cmd = " && ".join(
            filter(None, [
                f"{install_prefix} --upgrade pip",
                f"{install_prefix} {self._packages}" if self._packages else "",
                f"{install_prefix} {self._requirements}" if self._requirements else ""
            ])
        )
        return (
            f"mkdir -p {self._dir} && ("
                f" echo '[[ -x {self.python_path} ]] ||"
                "("
                    f"virtualenv --no-download {self.venv} && "
                    f"{install_cmd}"
                ") || ("
                    f"echo '\"'\"'{PYTHON_VENV_CREATE_ERR}'\"'\"' 1>&2 && exit 1"
                f")' | flock -w 900 {self._dir} bash "
            ")"
        )

    def python(self, cmd: str):
        return f"{self.get_venv} && {self.python_path} {cmd}"

    def script(self, cmd: str):
        stripped = cmd.strip()
        return self.python(f"{self.venv}/bin/{stripped}")

    def make_venv(self, cmd: str):
        return f"{self.get_venv} && {cmd}"
