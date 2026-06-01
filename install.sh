#!/usr/bin/env bash
set -eu

repo_root="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
target_dir="${AGENT_SKILLS_DIR:-${CODEX_HOME:-$HOME/.codex}/skills}"
repo_ref="${R_MED_STATS_REF:-main}"
archive_url="${R_MED_STATS_ARCHIVE_URL:-https://github.com/LeiGao0203/R-Medical-Statistics-Skills/archive/refs/heads/$repo_ref.tar.gz}"
tmp_dir=""

cleanup() {
  if [ -n "$tmp_dir" ]; then
    rm -rf "$tmp_dir"
  fi
}
trap cleanup EXIT

if [ -n "${R_MED_STATS_SOURCE:-}" ]; then
  source_dir="$R_MED_STATS_SOURCE"
elif [ -d "$repo_root/basic-stats" ]; then
  source_dir="$repo_root"
else
  tmp_dir="$(mktemp -d)"
  archive_file="$tmp_dir/source.tar.gz"

  if command -v curl >/dev/null 2>&1; then
    curl -fsSL "$archive_url" -o "$archive_file"
  elif command -v wget >/dev/null 2>&1; then
    wget -qO "$archive_file" "$archive_url"
  else
    printf 'Install failed: curl or wget is required to download %s\n' "$archive_url" >&2
    exit 1
  fi

  tar -xzf "$archive_file" -C "$tmp_dir"
  source_dir="$(find "$tmp_dir" -mindepth 1 -maxdepth 1 -type d | head -n 1)"
fi

for required in advanced-stats basic-stats literature-stats r-script quarto-report jupyter-notebook; do
  if [ ! -d "$source_dir/$required" ]; then
    printf 'Install failed: missing %s under %s\n' "$required" "$source_dir" >&2
    exit 1
  fi
done

mkdir -p "$target_dir"

for category in advanced-stats basic-stats literature-stats; do
  find "$source_dir/$category" -mindepth 1 -maxdepth 1 -type d -exec cp -R {} "$target_dir/" \;
done

for skill in r-script quarto-report jupyter-notebook; do
  cp -R "$source_dir/$skill" "$target_dir/"
done

printf 'Installed R Medical Statistics Skills to %s\n' "$target_dir"
