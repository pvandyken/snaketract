from typing import Dict, List
from uuid import uuid4


class XvfbRun:
    def __init__(self, active: bool):
        self.active = active

    def __call__(self, cmd: str):
        if self.active:
            return f"echo '{_quote_escape(cmd)}' | xvfb-run -a bash"
        return cmd

def _quote_escape(text: str):
    return text.replace("'", "'\"'\"'")

def _silent_mv(src: str, dest: str):
    """This was written out of concern for mv affecting the file timestamp, but it
    doesn't seem to. Leaving this here for now, but should eventually be removed if
    we never encounter problems."""
    return (
        #f"timestamp=$(stat -c %y {src}) && "
        f"mv {src} {dest}"
        #f"touch -hd \"$timestamp\" {dest}"
    )

def _cp_timestamp(src: str, dest: str):
    return (
        f"timestamp=$(stat -c %y {src}) && "
        f"touch -hd \"$timestamp\" {dest}"
    )


def hash_name(name: str):
        return f"$(realpath '{_quote_escape(name)}' | md5sum | awk '{{{{print $1}}}}')"


def _rm_if_exists(path: str, recursive=False):
    if recursive:
        flag = "-rf"
    else:
        flag = ""
    return f"( [[ ! -e {path} ]] || rm {flag} {path} )"


class Tar:
    def __init__(self, root: str):
        self.root = root

    def __call__(
        self,
        cmd: str,
        inputs: List[str]=[],
        outputs: List[str]=[],
        modify: List[str]=[]
    ):
        input_pre, input_posts = [*zip(
            *(
                (
                    self._open_tar(src, f"{self.root}/{hash_name(src)}"),
                    self._close_tar(src)

                ) for src in inputs
            )
        )] if inputs else ( [], [] )

        output_pre, output_posts = [*zip(
            *(
                (
                    f"( [[ ! -e {dest_stowed} ]] || ("
                        "echo "
                            "\"Found stashed tar file: "
                            f"'{dest_stowed}' "
                            "while atempting to generate the output: "
                            f"'{dest}' "
                            "Please rename this file, remove it, or manually change "
                            "its extension back to .tar.gz. If this file should not "
                            "have been processed, you may with to run `snakemake "
                            "--touch` to enforce correct timestamps for files\" && "
                        "false"
                    ")) && "
                    f"{_rm_if_exists(dest)} && "
                    f"{_rm_if_exists(tmpdir, True)} && "
                    f"mkdir -p {tmpdir} && "
                    f"ln -s {tmpdir} {dest}",

                    f"{self._save_tar(dest, tmpdir)}"

                ) for dest in outputs if (
                    (dest_stowed := self._stowed(dest)) and
                    (tmpdir := f"{self.root}/{hash_name(dest)}")
                )
            )
        )] if outputs else ([], [])

        modify_pre, modify_success, modify_fail = [*zip(
            *(
                (
                    self._open_tar(tar, tmpdir),
                    f"{self._save_tar(tar, tmpdir)}",
                    self._close_tar(tar)

                ) for tar in modify if (
                    (tmpdir := f"{self.root}/{hash_name(tar)}")
                )
            )
        )] if modify else ([], [], [])

        pre_script = " && ".join((
            *input_pre,
            *output_pre,
            *modify_pre
        ))

        post_success = (
            f"&& {s}" if (s:=" && ".join((
                *output_posts,
                *input_posts,
                *modify_success
            ))) else ""
        )

        post_fail = (
            f"|| ({s} && exit 1)" if (s:=" && ".join((
                *input_posts,
                *modify_fail
            ))) else ""
        )

        return f"{pre_script} && {cmd} {post_success} {post_fail}"



    def _open_tar(self, tar: str, mount: str):
        stowed = self._stowed(tar)
        return (
            f"([[ -d {mount} ]] && ( "
                f"[[ -e {stowed} ]] || {_silent_mv(tar, stowed)}"
            ") || ("
                f"mkdir -p {mount} && "
                f"([[ -e {stowed} ]] && ("
                    f"echo \"Found stowed tarfile: '{stowed}''. Extracting...\" && "
                    f"tar -xzf {stowed} -C {mount} && "
                    f"{_rm_if_exists(tar)}"
                ") || ("
                    f"echo \"Extracting and stowing tarfile: '{tar}'\" && "
                    f"tar -xzf {tar} -C {mount} && "
                    f"{_silent_mv(tar, stowed)} "
                "))"
            ")) && "
            f"ln -s {mount} {tar} && {_cp_timestamp(stowed, tar)}"
        )

    def _close_tar(self, tar: str):
        return f"rm {tar} && {_silent_mv(self._stowed(tar), tar)}"

    def _save_tar(self, tar: str, mount: str):
        stowed = self._stowed(tar)
        return (
            f"rm {tar} && "
            f"echo \"Packing tar file: {tar}\" && "
            f"tar -czf {tar} -C {mount} . && "
            f"{_rm_if_exists(stowed)}"
        )

    def _stowed(self, tar: str):
        return tar + ".__unpacked"

def display(cmd: str):
    return (
        f"$(: \033[1A\033[K\n\033[K\n\033[0;32m Command:\033[K\n\033[K\033[1;36m){cmd}$(: \n\033[1;30m)"
    )