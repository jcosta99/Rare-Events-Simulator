#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_path="${SIMULATION_CONFIG:-$project_dir/simulations.toml}"
exec "${PYTHON_BIN:-python3}" "$project_dir/Simulator/launch_simulations.py" \
    "$config_path" "$@"
