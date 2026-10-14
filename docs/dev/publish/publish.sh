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
# It also removes the location data (EXIF, XMP) of the 2025 photo, a Windows
# user folder in one 2025 script and the early README sentences about the
# time before the rework (replacements.txt, strip_metadata.py), and checks
# every blob for them afterwards (check_history.py, needs Python).
# It stops if the local repository that holds this script has uncommitted
# changes (outside docs/dev) or another HEAD than the commit it publishes;
# ALLOW_DIFF=1 turns that into a warning.
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
SRC=${SRC:-https://github.com/onewmann/hs-esslingen}
BRANCH=${BRANCH:-claude/pensive-gates-pzcxla}
DEST=${DEST:-https://github.com/onewmann/2D-FFT-Code}
TAG=student-project-2025
OWNER_NAME='Oliver Neumann'
OWNER_EMAIL='149961666+onewmann@users.noreply.github.com'
WORK=$(mktemp -d)

fail() { echo "error: $*" >&2; exit 1; }
differs() {
  if [ "${ALLOW_DIFF:-0}" = 1 ]; then echo "warning: $*" >&2; else fail "$* (ALLOW_DIFF=1 to go on)"; fi
}

PY=
for c in python3 python py; do                  # on Windows python3 can be a store stub that does not run
  if "$c" -c 'import sys' >/dev/null 2>&1; then PY=$c; break; fi
done
[ -n "$PY" ] || fail "no Python found (python3, python or py)"
cp "$HERE/replacements.txt" "$HERE/messages.txt" "$HERE/mailmap.txt" "$HERE/strip_metadata.py" \
  "$HERE/check_history.py" "$WORK/"
git clone -q "$SRC" "$WORK/repo"
cd "$WORK/repo"
git checkout -q -B main "origin/$BRANCH"
echo "source: $SRC $BRANCH at $(git rev-parse --short HEAD): $(git log -1 --format=%s)"
if ROOT=$(git -C "$HERE" rev-parse --show-toplevel 2>/dev/null); then
  if [ -n "$(git -C "$ROOT" status --porcelain -- . ':(exclude)docs/dev')" ]; then
    differs "$ROOT has uncommitted changes outside docs/dev; they would not be published"
  fi
  if [ "$(git -C "$ROOT" rev-parse HEAD)" != "$(git rev-parse HEAD)" ]; then
    differs "the local HEAD $(git -C "$ROOT" rev-parse --short HEAD) is not the published commit;" \
      "push it first, or run with SRC=$ROOT"
  fi
fi
# the tagger takes the committer identity: fixed here, not whatever this machine has configured
GIT_COMMITTER_NAME=$OWNER_NAME GIT_COMMITTER_EMAIL=$OWNER_EMAIL git tag -f -a "$TAG" 8f0aec4 \
  -m "Student project scripts as uploaded. The AI-assisted rework starts after this commit." >/dev/null
for ref in $(git for-each-ref --format='%(refname)' refs/remotes refs/heads refs/tags); do
  case "$ref" in refs/heads/main|refs/tags/$TAG) ;; *) git update-ref -d "$ref" ;; esac
done
git remote remove origin

git filter-repo --force \
  --invert-paths --path docs/dev --path 2D-FFT/README.md \
  --replace-text ../replacements.txt \
  --replace-message ../messages.txt \
  --mailmap ../mailmap.txt \
  --blob-callback "$(cat ../strip_metadata.py)"

if git grep -i -q esslingen $(git rev-list --all); then
  fail "a file in the history still names the university"
fi
"$PY" ../check_history.py || fail "a blob in the history still needs cleaning (see above)"
if { git log --all --format='%B%n%an %ae%n%cn %ce'; git for-each-ref refs/tags --format='%(contents)'; } \
    | grep -i -q -e esslingen -e 'docs/dev' -e 'co-authored-by: claude' -e 'claude-session' -e 'noreply@anthropic' \
      -e 'without ai' -e 'rework only'; then
  fail "a commit or tag message or identity still needs cleaning"
fi
[ "$(git for-each-ref --format='%(refname)' refs/heads refs/tags | tr '\n' ' ')" = "refs/heads/main refs/tags/$TAG " ] \
  || fail "expected exactly the branch main and the tag $TAG"
[ "$(git for-each-ref refs/tags --format='%(objecttype) %(taggername) %(taggeremail)')" = "tag $OWNER_NAME <$OWNER_EMAIL>" ] \
  || fail "the tag $TAG is not annotated by $OWNER_NAME <$OWNER_EMAIL>"
git merge-base --is-ancestor "$TAG" main || fail "the tag $TAG is not in the history of main"
git cat-file -e "$TAG:2D-FFT" 2>/dev/null || fail "the README links $TAG/2D-FFT, which the tag does not hold"
echo "clean history: $WORK/repo ($(git rev-list --count main) commits, tag $TAG)"

if [ "${1:-}" = "--push" ]; then
  git remote add origin "$DEST"
  git push --atomic -u origin main "$TAG"
fi
