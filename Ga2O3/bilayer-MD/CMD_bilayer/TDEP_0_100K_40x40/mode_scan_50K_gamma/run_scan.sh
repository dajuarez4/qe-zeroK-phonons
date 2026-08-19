#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
: > energies.raw
for datafile in structures/mode_*.data; do
    output=$(../bin/lmp-tabgap -var datafile "$datafile" -in in.single_point -log none)
    energy=$(printf '%s\n' "$output" | awk '/^ENERGY / {print $2}')
    printf '%s %s\n' "$datafile" "$energy" >> energies.raw
done
