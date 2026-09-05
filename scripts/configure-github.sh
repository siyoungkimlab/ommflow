#!/usr/bin/env bash
# Apply the repository settings that make main a protected, tested branch.
#
# Run once after the repository exists:
#     bash scripts/configure-github.sh [OWNER/REPO]
#
# Defaults to siyoungkimlab/ommflow. Safe to re-run; every call is a PUT or
# PATCH of the full desired state.
#
# Note: on GitHub Free, branch protection is available on public repositories.
# Private repositories need a paid plan, and the protection call below will
# fail with 403 on a private repo without one.
set -euo pipefail

REPO="${1:-siyoungkimlab/ommflow}"
BRANCH=main

command -v gh >/dev/null 2>&1 || { echo "gh CLI is required." >&2; exit 1; }
gh auth status >/dev/null 2>&1 || { echo "Run 'gh auth login' first." >&2; exit 1; }

echo "==> Merge settings for $REPO"
# Squash-only, so main gets one commit per pull request, and that commit's
# message is the pull request title and body.
gh api -X PATCH "repos/$REPO" \
    -F allow_squash_merge=true \
    -F allow_merge_commit=false \
    -F allow_rebase_merge=false \
    -f squash_merge_commit_title=PR_TITLE \
    -f squash_merge_commit_message=PR_BODY \
    -F delete_branch_on_merge=true \
    -F has_wiki=false \
    >/dev/null
echo "    squash-only; commit message taken from the pull request; branches deleted on merge"

echo "==> Protecting $BRANCH"
# "all checks" is the gate job in .github/workflows/tests.yml, which depends on
# every other job. Requiring it by name survives changes to the test matrix.
#
# enforce_admins is false so an owner can still repair a broken main; set it to
# true to hold yourself to the same rule as everyone else.
gh api -X PUT "repos/$REPO/branches/$BRANCH/protection" --input - >/dev/null <<'JSON'
{
  "required_status_checks": {
    "strict": false,
    "contexts": ["all checks"]
  },
  "enforce_admins": false,
  "required_pull_request_reviews": {
    "required_approving_review_count": 0,
    "dismiss_stale_reviews": true
  },
  "restrictions": null,
  "required_linear_history": true,
  "required_conversation_resolution": true,
  "allow_force_pushes": false,
  "allow_deletions": false
}
JSON
echo "    pull request required; 'all checks' must pass; no force pushes or deletions"

echo
echo "==> Verifying"
gh api "repos/$REPO/branches/$BRANCH/protection" \
    --jq '"    required checks : " + (.required_status_checks.contexts | join(", "))
        + "\n    pull request    : " + (if .required_pull_request_reviews then "required" else "NOT required" end)
        + "\n    force pushes    : " + (if .allow_force_pushes.enabled then "ALLOWED" else "blocked" end)
        + "\n    linear history  : " + (if .required_linear_history.enabled then "required" else "not required" end)'
gh api "repos/$REPO" \
    --jq '"    squash title    : " + .squash_merge_commit_title
        + "\n    squash body     : " + .squash_merge_commit_message
        + "\n    merge commits   : " + (if .allow_merge_commit then "ALLOWED" else "blocked" end)'
echo
echo "Done."
