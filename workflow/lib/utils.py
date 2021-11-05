from typing import Dict


def xvfb_run(config: Dict, cmd: str):
    if config.get('x11_srv', False):
        escaped = cmd.replace("'", "'\"'\"'")
        return f"echo '{escaped}' | xvfb-run -a bash"
    return cmd
