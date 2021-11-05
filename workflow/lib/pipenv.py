from typing import List
from pathlib import Path

class PipEnv:
    def __init__(
        self,
        name: str,
        packages: List[str],
        root: Path,
        flags: str,
        
    ):
        self.venv = root/name/"venv"
        self.python = self.venv/"bin"/"python"
        self.flags = flags
        self.packages = ' '.join(packages)

    @property
    def script(self):
        return (
            "("
                f"[[ -x {self.python} ]] || ( "
                    f"virtualenv --no-download {self.venv} && "
                    f"{self.python} -m pip install {self.flags} --upgrade pip && "
                    f"{self.python} -m pip install {self.flags} {self.packages} "
                ")"
            f") && {self.python} {self.venv}/bin/"
        )


