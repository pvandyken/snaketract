#!/usr/bin/env python3
import argparse
import yaml
import itertools as it
import os
from pathlib import Path
import re
import sys

def reverse_split(s: str, reverse: bool = False):
    if reverse:
        l = s.split()
        l.reverse()
        return l
    return s.split()

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("pattern_file", type=Path)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--dry", "-n", action="store_true")
    parser.add_argument("--reverse", "-r", action="store_true")

    args = parser.parse_args(sys.argv[1:])
    with args.pattern_file.open() as f:
        scheme = yaml.load(f)

    templates = [reverse_split(template, args.reverse) for template in scheme["files"]]
    wildcards = {
        name: f"({match})" for name, match in scheme.get("wildcards", {}).items() 
    }

    sub_tags = {
        name: f"\\1{i+1}" for name, i in enumerate()
    }
    print([
        ( regex.format(**wildcards), sub.format(**wildcards) )
        for regex, sub in templates
    ])
    exit()

    patterns = [
        ( re.compile(regex.format(**wildcards)), sub.format(**wildcards) )
        for regex, sub in templates
    ]
    os.chdir(args.directory)

    for path, (regex, replace) in it.product(Path().rglob("*"), patterns):
        if regex.match(str(path)):
            new = regex.sub(replace, str(path))
            if args.dry:
                print(str(path) + " -> " + new)
            else:
                path.rename(new)

if __name__ == "__main__":
    main()
    exit(0)