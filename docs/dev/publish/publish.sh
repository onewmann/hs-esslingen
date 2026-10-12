#!/usr/bin/env bash
# Build the clean history for onewmann/2D-FFT-Code from the working branch.
#
#   bash docs/dev/publish/publish.sh          build and check only
#   bash docs/dev/publish/publish.sh --push   also push main and the tag
#
# Needs git and git-filter-repo (pip install git-filter-repo); on Windows run it
# in Git Bash. The new repository must exist and be empty before --push.
# The rewrite drops docs/dev and the 2025 GitLab template README from every
# commit, replaces the remaining mentions of the university, sets every author
# and committer to Oliver Neumann and removes the AI trailers from the messages.
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
SRC=${SRC:-https://github.com/onewmann/hs-esslingen}
BRANCH=${BRANCH:-claude/pensive-gates-pzcxla}
DEST=${DEST:-https://github.com/onewmann/2D-FFT-Code}
WORK=$(mktemp -d)

cp "$HERE/replacements.txt" "$HERE/messages.txt" "$HERE/mailmap.txt" "$WORK/"
git clone -q "$SRC" "$WORK/repo"
cd "$WORK/repo"
git checkout -q -B main "origin/$BRANCH"
git tag -a student-project-2025 8f0aec4 \
  -m "Student project scripts as uploaded. The AI-assisted rework starts after this commit." 2>/dev/null || true
for ref in $(git for-each-ref --format='%(refname)' refs/remotes refs/heads refs/tags); do
  case "$ref" in refs/heads/main|refs/tags/student-project-2025) ;; *) git update-ref -d "$ref" ;; esac
done
git remote remove origin

git filter-repo --force \
  --invert-paths --path docs/dev --path 2D-FFT/README.md \
  --replace-text ../replacements.txt \
  --replace-message ../messages.txt \
  --mailmap ../mailmap.txt

if git grep -i -q esslingen $(git rev-list --all); then
  echo "error: a file in the history still names the university" >&2; exit 1
fi
if git log --all --format='%B%n%an %ae%n%cn %ce' \
    | grep -i -q -e esslingen -e 'co-authored-by: claude' -e 'claude-session' -e 'noreply@anthropic'; then
  echo "error: a commit message or identity still needs cleaning" >&2; exit 1
fi
echo "clean history: $WORK/repo ($(git rev-list --count main) commits)"

if [ "${1:-}" = "--push" ]; then
  git remote add origin "$DEST"
  git push -u origin main
  git push origin student-project-2025
fi
