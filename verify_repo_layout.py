#!/usr/bin/env python3
"""
verify_repo_layout.py - covalent5 reproducibility archive self-check.

Walks the deposited tree, compares it against the directory inventory documented in
README.md, and exits non-zero on any mismatch. Written after the first revision, in
which the archive was populated through the GitHub web upload interface and was
silently truncated at 100 files.

Usage:
    python verify_repo_layout.py
    python verify_repo_layout.py -o repo_inventory.txt
"""
from __future__ import annotations
import argparse, os, sys

REQUIRED = {
    "data":           {"min_files": 10, "must_contain": ["anchor_ki.csv", "descriptors.csv",
                                                         "features.csv", "molecule_manifest.csv"]},
    "dft":            {"min_files": 40, "must_contain": ["job_manifest.csv"]},
    "docking":        {"min_files": 20, "must_contain": ["docking_results.csv", "4EY7.pdb"]},
    "qm_cluster":     {"min_files":  5, "must_contain": ["scan_alkoxide.out",
                                                         "cluster_reactant.xyz"]},
    "conf_stability": {"min_files": 20, "must_contain": ["conf_stability_summary.csv",
                                                         "conf_manifest.csv"]},
    "experimental":   {"min_files":  1, "must_contain": []},
}

JOB_FILES = {"opt.inp", "opt.out", "opt.xyz", "sp.inp", "sp.out"}
JOB_ROOTS = ("dft", "conf_stability")

# Enumerated stereoisomers that carry no measured kinetics and were deliberately
# not computed. Declared here so the check reports them as documented, not as loss.
DECLARED_NOT_COMPUTED = {
    "BSAR_01_P5R", "DETAB_01_P3S", "FEN_01_P3S",
    "GF_01_P1R", "METH_01_P2S", "VR_01_P8R",
}


def walk(root):
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d != ".git"]
        out += [os.path.join(dirpath, f) for f in filenames]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("-o", "--out")
    a = ap.parse_args()
    os.chdir(a.root)

    lines = []
    def say(s=""):
        print(s); lines.append(s)

    failures = 0
    say("covalent5 reproducibility archive - layout verification")
    say("=" * 62)
    say(f"{'directory':<18}{'files':>8}   status")
    say("-" * 62)

    for name, spec in REQUIRED.items():
        if not os.path.isdir(name):
            say(f"{name:<18}{'--':>8}   FAIL  directory absent")
            failures += 1
            continue
        files = walk(name)
        base = {os.path.basename(f) for f in files}
        problems = []
        if len(files) < spec["min_files"]:
            problems.append(f"only {len(files)} files, expected >= {spec['min_files']}")
        missing = [m for m in spec["must_contain"] if m not in base]
        if missing:
            problems.append("missing " + ", ".join(missing))
        if problems:
            say(f"{name:<18}{len(files):>8}   FAIL  " + "; ".join(problems))
            failures += 1
        else:
            say(f"{name:<18}{len(files):>8}   ok")

    say("-" * 62)
    incomplete = []; declared = []
    for root in JOB_ROOTS:
        if not os.path.isdir(root):
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            if not any(f.endswith(".inp") for f in filenames):
                continue
            if dirnames:
                continue
            missing = JOB_FILES - set(filenames)
            if missing:
                if os.path.basename(dirpath) in DECLARED_NOT_COMPUTED:
                    declared.append(dirpath)
                else:
                    incomplete.append((dirpath, sorted(missing)))
    if incomplete:
        say(f"FAIL  {len(incomplete)} job directory(ies) incomplete:")
        for d, m in incomplete[:20]:
            say(f"      {d}  missing {', '.join(m)}")
        if len(incomplete) > 20:
            say(f"      ... and {len(incomplete) - 20} more")
        failures += 1
    else:
        say("ok    all job directories carry a complete file set")
    if declared:
        say(f"info  {len(declared)} enumerated stereoisomer(s) declared not computed "
            "(no measured kinetics; excluded from all analyses):")
        for d in declared:
            say(f"      {d}")

    total = len(walk("."))
    say("-" * 62)
    say(f"total tracked-tree files: {total}")
    if total == 100:
        say("WARN  exactly 100 files - this is the GitHub web-upload cap. "
            "Confirm the archive was pushed with git, not uploaded through the browser.")
    say("=" * 62)
    say("RESULT: PASS" if failures == 0 else f"RESULT: FAIL ({failures} check(s) failed)")

    if a.out:
        with open(a.out, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"\nreport written to {a.out}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
