from typing import List
from pathlib import Path


PYTHON_VENV_CREATE_ERR = "[ERROR] (jobid={jobid}): Error creating python environment"
TIMED_OUT_ERR = "[ERROR] (jobid={jobid}): Script timed out when waiting for Python"


class PipEnv:
    def __init__(
        self,
        name: str,
        packages: List[str],
        root: Path,
        flags: str,
        requirements: List[str] = [],

    ):
        self._dir = root/name
        self._venv_lock = root/name/".venv_lock"
        self._venv = root/name/"venv"
        self._python = self._venv/"bin"/"python"
        self._flags = flags
        self._packages = ' '.join(packages)
        self._requirements = '-r ' + ' -r '.join(requirements) if requirements else ""


    @property
    def get_venv(self):
        install_prefix = f"{self._python} -m pip install {self._flags}"
        install_cmd = " && ".join(
            filter(None, [
                f"{install_prefix} --upgrade pip",
                f"{install_prefix} {self._packages}" if self._packages else "",
                f"{install_prefix} {self._requirements}" if self._requirements else ""
            ])
        )
        return (
            "("
                f"[[ -d {self._dir} ]] || mkdir {self._dir}"
            ") && ("
                f"[[ -x {self._python} ]] || ("
                    f"[[ ! -e {self._venv_lock} ]] && ( "
                        "("
                            f"touch {self._venv_lock} &&"
                            f"virtualenv --no-download {self._venv} && "
                            f"{install_cmd} && "
                            f"rm {self._venv_lock}"
                        ") || ("
                            f"echo '{PYTHON_VENV_CREATE_ERR}' 1>&2 && exit 1"
                        ") "
                    ") || ("
                        f"echo 'while [[ -e {self._venv_lock} ]]; do sleep 5; done' | "
                            "timeout 15m bash || "
                        "("
                            f"echo '{TIMED_OUT_ERR}' 1>&2 && "
                            "exit 1"
                        ")"
                    ")"
                ")"
            ") && "
        )

    def python(self, cmd: str):
        return f"{self.get_venv} {self._python} {cmd}"

    def script(self, cmd: str):
        stripped = cmd.strip()
        return self.python(f"{self._venv}/bin/{stripped}")
