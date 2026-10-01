#!/usr/bin/env bash
# ==============================================================================
# OASpepDB: Reference Proteome Download Wrapper
# Downloads UniProt Swiss-Prot and NCBI RefSeq Human Proteomes into Data/
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/../scripts/download_reference_proteomes.py"

if command -v python3 &>/dev/null; then
    python3 "${PYTHON_SCRIPT}" --out-dir "${SCRIPT_DIR}" "$@"
elif command -v python &>/dev/null; then
    python "${PYTHON_SCRIPT}" --out-dir "${SCRIPT_DIR}" "$@"
else
    echo "Error: Python 3 is required to run the automated reference downloader." >&2
    exit 1
fi
