#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
echo "=== trilattice : minimal 3-panel viewer ==="
python3 -m pip install -e ".[fast]" -q
python3 -m trilattice.cli animate examples/config.toml -T 1800 --minimal
