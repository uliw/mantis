#!/usr/bin/env bash
# Local script to export org-mode tutorials to GitHub-Flavored Markdown for Sphinx.
# NOTE: This script is intended for LOCAL USE ONLY and should NOT be run in CI/RTD.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if ! command -v emacsclient >/dev/null 2>&1; then
    echo "Error: emacsclient is not found in PATH. Please run an Emacs server locally."
    exit 1
fi

echo "Exporting Tutorial org files to Markdown..."

files_to_export=("$@")
if [ ${#files_to_export[@]} -eq 0 ]; then
    files_to_export=(Tutorial_1.org Tutorial_2.org Tutorial_3.org Tutorial_4.org getting_help.org)
fi

for org_file in "${files_to_export[@]}"; do
    # Strip any directory prefix if passed like docs/Tutorial_1.org
    base_file="$(basename "$org_file")"
    if [ -f "$base_file" ]; then
        echo "Exporting $base_file ..."
        emacsclient -e "(progn (switch-to-buffer (find-file-noselect \"$base_file\")) (end-of-buffer) (org-gfm-export-to-markdown) (kill-buffer))"
    else
        echo "Warning: $base_file not found in $SCRIPT_DIR"
    fi
done

echo "Export completed successfully. Generated markdown files are ready in docs/"
