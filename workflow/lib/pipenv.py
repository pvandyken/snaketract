from typing import List
from pathlib import Path

class PipEnv:
    def __init__(
        self,
        name: str,
        packages: List[str],
        requirements: List[str],
        root: Path,
        flags: str,
        
    ):
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
                f"{install_prefix} --upgrade pip"
                f"{install_prefix} {self._packages}" if self._packages else "",
                f"{install_prefix} {self._requirements}" if self._requirements else ""
            ])
        )
        return (
            "("
                f"[[ -x {self._python} ]] || ( "
                    f"virtualenv --no-download {self._venv} && "
                    f"{install_cmd}"
                ")"
            ")"
        )

    def python(self, cmd: str):
        return f"{self.get_venv} && {self._python}"

    def script(self, cmd: str):
        stripped = cmd.strip()
        return f"{self.python(cmd)} {self._venv}/bin/{stripped}"

        


