#!/usr/bin/env python3
import argparse
import itertools as it
import os
from pathlib import Path
import re
import sys

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("pattern_file", type=Path)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--dry", "-n", action="store_true")
    parser.add_argument("--reverse", "-r", action="store_true")

    args = parser.parse_args(sys.argv[1:])
    with args.pattern_file.open() as f:
        patterns = f.read()
        lines = patterns.split("\n")

    os.chdir(args.directory)

    def sp(s: str):
        if args.reverse:
            l = s.split()
            l.reverse()
            return l
        return s.split()

    patterns = [(re.compile(sp(line)[0]), sp(line)[1]) for line in filter(None, lines)]

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