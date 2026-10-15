# Contributing

## Workflow
- Large changes go on a branch (or a fork) and are merged when reviewed; small, atomic commits, one per topic.
- Every number in the code carries its source (a constant, a paper, or the word "Placeholder"); camelCase names, no unit suffixes. See [docs/MANUAL.md](docs/MANUAL.md), section "Conventions".
- Run `python tests/test_sanityChecks.py` before committing. Notebooks are run with `jupyter nbconvert --to notebook --execute --inplace` and must exit 0.

## Crediting AI assistance
This repository follows the standard of the mplhep package and of the Linux kernel ([`Documentation/process/coding-assistants.rst`](https://docs.kernel.org/process/coding-assistants.html)):

- **An AI is never a co-author.** Do not add a `Co-authored-by:` trailer for an AI: it asserts authorship and a copyright interest, and the work stays the human's.
- **An AI never signs off.** Do not let a tool add `Signed-off-by:` on your behalf: it certifies the [Developer Certificate of Origin](https://developercertificate.org/), a legal statement only a human can make.
- **Credit the assistance with a trailer in the commit**, in the form `Assisted-by: <harness>:<model>`, for example `Assisted-by: claude-code:claude-opus-5-5`. If a squash message is proposed in a pull request body, put the trailer there too.

The rules are checked by [`tools/checkCommitTrailers.py`](tools/checkCommitTrailers.py). Enable the commit-msg hook once per clone, and audit a range of history when needed:

```
git config core.hooksPath .githooks
python tools/checkCommitTrailers.py --range origin/main..HEAD
```
