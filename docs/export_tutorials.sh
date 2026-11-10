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

for org_file in Tutorial_1.org Tutorial_2.org; do
    if [ -f "$org_file" ]; then
        echo "Exporting $org_file ..."
        emacsclient -e "(progn (switch-to-buffer (find-file-noselect \"$org_file\")) (end-of-buffer) (org-gfm-export-to-markdown) (kill-buffer))"
    else
        echo "Warning: $org_file not found in $SCRIPT_DIR"
    fi
done

echo "Export completed successfully. Generated markdown files are ready in docs/"
