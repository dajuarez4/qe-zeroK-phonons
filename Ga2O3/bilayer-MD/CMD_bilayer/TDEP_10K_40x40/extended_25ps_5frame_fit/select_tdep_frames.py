#!/usr/bin/env python3
"""Select five evenly spaced configurations from the full 26-frame TDEP data."""

from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
NAT = 16000
SELECT = [0, 6, 12, 19, 25]


def subset_blocks(name):
    source = HERE / f"{name}.all26"
    current = HERE / name
    if not source.exists():
        shutil.move(current, source)
    wanted = set(SELECT)
    with source.open() as inp, current.open("w") as out:
        for iframe in range(26):
            block = [inp.readline() for _ in range(NAT)]
            if any(line == "" for line in block):
                raise RuntimeError(f"Unexpected end of {source} at frame {iframe}")
            if iframe in wanted:
                out.writelines(block)


for filename in ("infile.positions", "infile.forces"):
    subset_blocks(filename)

stat_all = HERE / "infile.stat.all26"
stat_now = HERE / "infile.stat"
if not stat_all.exists():
    shutil.move(stat_now, stat_all)
rows = stat_all.read_text().splitlines()
with stat_now.open("w") as out:
    for new_index, old_index in enumerate(SELECT, start=1):
        fields = rows[old_index].split()
        fields[0] = str(new_index)
        out.write(" ".join(fields) + "\n")

(HERE / "infile.meta").write_text(
    f"{NAT} # N atoms\n{len(SELECT)} # N configurations\n"
    "1000.0 # original trajectory sampling interval, fs\n10.0 # temperature, K\n"
)
print("Selected frames:", SELECT)
