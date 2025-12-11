#!/usr/bin/env python3
import os
import sys
import subprocess
import argparse

def cmp_meshes(f1, f2):
    """Default operation: call external tool cmp_meshes."""
    try:
        result = subprocess.run(
            ["cmp_meshes", f1, f2],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        # print(f"[cmp_meshes] {os.path.basename(f1)}:")
        print(result.stdout, end='')
        if result.stderr:
           print(result.stderr, file=sys.stderr, end='')
    except FileNotFoundError:
        print("Error: cmp_meshes executable not found.", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Apply operations to matching OFF files in two directories."
    )
    parser.add_argument("d1", help="Directory 1")
    parser.add_argument("d2", help="Directory 2")
    parser.add_argument(
        "--op",
        help="Operation to apply on matching files. Default: cmp_meshes",
        default="cmp_meshes"
    )
    args = parser.parse_args()

    d1 = args.d1
    d2 = args.d2

    if not os.path.isdir(d1) or not os.path.isdir(d2):
        print("Error: both arguments must be directories", file=sys.stderr)
        sys.exit(1)

    # Operation dispatch
    if args.op == "cmp_meshes":
        operation = cmp_meshes
    else:
        print(f"Error: unknown operation {args.op}", file=sys.stderr)
        sys.exit(1)

    # Collect off files (case-insensitive)
    def off_files(directory):
        return {
            f for f in os.listdir(directory)
            if f.lower().endswith(".off") and os.path.isfile(os.path.join(directory, f))
        }

    off1 = off_files(d1)
    off2 = off_files(d2)

    # Files present in both
    common = sorted(off1 & off2)
    # Files present only in one
    only1 = sorted(off1 - off2)
    only2 = sorted(off2 - off1)

    # Report missing files
    if only1:
        print("Files present in d1 but not in d2:")
        for f in only1:
            print("  ", f)

    if only2:
        print("Files present in d2 but not in d1:")
        for f in only2:
            print("  ", f)

    # Apply operation on matching files
    for f in common:
        f1 = os.path.join(d1, f)
        f2 = os.path.join(d2, f)
        print(f"=== Processing {f} ===")
        operation(f1, f2)


if __name__ == "__main__":
    main()
