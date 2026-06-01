#!/usr/bin/env bash
set -eu

repo_root="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

fail() {
  printf 'FAIL: %s\n' "$1" >&2
  exit 1
}

assert_installed() {
  install_dir="$1"
  test -f "$install_dir/ttest/SKILL.md" || fail "basic statistics skills are not flattened into the target directory"
  test -f "$install_dir/logistic-reg/SKILL.md" || fail "advanced statistics skills are not flattened into the target directory"
  test -f "$install_dir/ps-matching/SKILL.md" || fail "literature statistics skills are not flattened into the target directory"
  test -f "$install_dir/r-script/SKILL.md" || fail "r-script skill is not copied"
  test -f "$install_dir/quarto-report/SKILL.md" || fail "quarto-report skill is not copied"
  test -f "$install_dir/jupyter-notebook/SKILL.md" || fail "jupyter-notebook skill is not copied"
  test ! -d "$install_dir/basic-stats" || fail "category directories should not be copied as nested skills"
}

AGENT_SKILLS_DIR="$tmp_dir/local-skills" bash "$repo_root/install.sh" >/dev/null
assert_installed "$tmp_dir/local-skills"

cp "$repo_root/install.sh" "$tmp_dir/install.sh"
AGENT_SKILLS_DIR="$tmp_dir/external-skills" R_MED_STATS_SOURCE="$repo_root" bash "$tmp_dir/install.sh" >/dev/null
assert_installed "$tmp_dir/external-skills"

tar -czf "$tmp_dir/source.tar.gz" --exclude .git -C "$(dirname "$repo_root")" "$(basename "$repo_root")"
AGENT_SKILLS_DIR="$tmp_dir/archive-skills" R_MED_STATS_ARCHIVE_URL="file://$tmp_dir/source.tar.gz" bash "$tmp_dir/install.sh" >/dev/null
assert_installed "$tmp_dir/archive-skills"
