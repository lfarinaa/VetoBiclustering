"""Checks commit messages against the rules for crediting AI assistance (the standard of the mplhep package and of the Linux kernel, Documentation/process/coding-assistants.rst):

- an AI is never a co-author: no "Co-authored-by:" trailer naming an AI (it asserts authorship and a copyright interest, and the work stays the human's);
- an AI never signs off: no "Signed-off-by:" naming an AI or a tool (it certifies the Developer Certificate of Origin, which only a human can);
- AI assistance is credited with "Assisted-by: <harness>:<model>", for example "Assisted-by: claude-code:claude-sonnet-5-5".

Use as a commit-msg hook (see CONTRIBUTING.md) or to audit history:

    python tools/checkCommitTrailers.py --message .git/COMMIT_EDITMSG      # one message, exit code 1 if it breaks a rule
    python tools/checkCommitTrailers.py --range origin/main..HEAD          # every commit in a range, listing the offenders
"""
import argparse
import re
import subprocess
import sys

aiNames = re.compile(r"claude|anthropic|openai|chatgpt|\bgpt|copilot|gemini|codex|cursor|\bllm\b|\bai\b|assistant|bot\b", re.IGNORECASE)   # Names that mark a trailer as belonging to an AI or a tool.
assistedByFormat = re.compile(r"^Assisted-by:\s*[A-Za-z0-9._-]+:[A-Za-z0-9._-]+(\s*\(.*\))?\s*$")   # Assisted-by: <harness>:<model>.


def findProblems(message):
    """The rules a commit message breaks, as a list of sentences (empty if it follows them)."""
    problems = []
    for line in message.splitlines():
        if re.match(r"(?i)^co-authored-by:", line) and aiNames.search(line):
            problems.append(f"an AI cannot be a co-author: '{line.strip()}' (use 'Assisted-by: <harness>:<model>')")
        if re.match(r"(?i)^signed-off-by:", line) and aiNames.search(line):
            problems.append(f"an AI cannot sign off: '{line.strip()}' (the sign-off certifies the Developer Certificate of Origin, which only a human can)")
        if re.match(r"(?i)^assisted-by:", line) and not assistedByFormat.match(line):
            problems.append(f"malformed credit: '{line.strip()}' (the form is 'Assisted-by: <harness>:<model>', for example 'Assisted-by: claude-code:claude-sonnet-5-5')")
    return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--message", help="path of a commit message file (a commit-msg hook gets it as its argument)")
    parser.add_argument("--range", help="a git revision range to audit, for example origin/main..HEAD")
    arguments = parser.parse_args()
    if arguments.message:
        problems = findProblems(open(arguments.message, encoding="utf-8").read())
        for problem in problems:
            print(f"commit message rejected: {problem}", file=sys.stderr)
        return 1 if problems else 0
    if arguments.range:
        hashes = subprocess.run(["git", "rev-list", arguments.range], capture_output=True, text=True, check=True).stdout.split()
        offenders = []
        for commit in hashes:
            body = subprocess.run(["git", "log", "-1", "--format=%B", commit], capture_output=True, text=True, check=True).stdout
            problems = findProblems(body)
            if problems:
                offenders.append((commit, body.splitlines()[0], problems))
        for commit, subject, problems in offenders:
            print(f"{commit[:9]} {subject[:80]}\n    " + "\n    ".join(problems))
        print(f"{len(offenders)} of {len(hashes)} commits in {arguments.range} break the rules")
        return 1 if offenders else 0
    parser.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
