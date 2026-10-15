"""Publishes the executed notebooks, outputs included, on a separate branch so that GitHub shows them rendered while main keeps stripped notebooks.

    python tools/publishNotebooks.py            # commit the local executed notebooks to the 'notebooks' branch and force-push it
    python tools/publishNotebooks.py --no-push  # only update the local branch

main never holds outputs: nbstripout strips them on every commit (.gitattributes). This script writes the executed working copies as git objects
without filters (git hash-object --no-filters), in a commit with no parent, on the branch 'notebooks', replacing it each time: the branch is always one
commit, so outputs never accumulate in history. Run the notebooks first (jupyter nbconvert --execute --inplace, or with BICLUSTERING_USE_CACHE=1);
the commit message records the main commit they were run from.
"""
import argparse
import pathlib
import subprocess
import sys

root = pathlib.Path(__file__).resolve().parent.parent
branch = "notebooks"
notebooks = ["examples/channelBudgetShowcase.ipynb", "examples/vetoPatternShowcase.ipynb"]
readme = """# Executed notebooks

Generated branch: the notebooks of VetoBiclustering with their outputs, so that GitHub renders them. It holds a single commit, replaced each time by
`tools/publishNotebooks.py`; the sources live on `main`, where notebook outputs are stripped. Run from main commit {head}.

- [channelBudgetShowcase.ipynb](examples/channelBudgetShowcase.ipynb): channel budgets on the regenerated RoIs and the Geant4 photons
- [vetoPatternShowcase.ipynb](examples/vetoPatternShowcase.ipynb): a first look on the February 2023 dump
"""


def git(*args, input=None):
    return subprocess.run(["git", *args], cwd=root, input=input, capture_output=True, text=True, check=True).stdout.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--no-push", action="store_true")
    args = parser.parse_args()
    head = git("rev-parse", "--short", "HEAD")
    entries = {}
    for path in notebooks:
        if not (root / path).exists():
            sys.exit(f"missing {path}")
        directory, name = path.rsplit("/", 1)
        entries.setdefault(directory, []).append((name, git("hash-object", "-w", "--no-filters", path)))
    subtrees = [f"040000 tree {git('mktree', input=''.join(f'100644 blob {sha}\t{name}\n' for name, sha in files))}\t{directory}\n"
                for directory, files in entries.items()]
    readmeBlob = git("hash-object", "-w", "--stdin", input=readme.format(head=head))
    tree = git("mktree", input="".join(subtrees) + f"100644 blob {readmeBlob}\tREADME.md\n")
    commit = git("commit-tree", tree, "-m", f"Executed notebooks from main {head}")
    git("update-ref", f"refs/heads/{branch}", commit)
    print(f"branch {branch} -> {commit[:7]} (notebooks executed from {head})")
    if not args.no_push:
        git("push", "--force", "origin", f"{branch}:{branch}")
        print(f"pushed {branch}")


if __name__ == "__main__":
    main()
