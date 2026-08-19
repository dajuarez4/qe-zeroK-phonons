#!/usr/bin/env python3
"""Create 60--100 K inputs from the validated 50 K calculation templates."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "tdep_50K_40x40"

for temperature in range(60, 101, 10):
    folder = ROOT / f"tdep_{temperature}K_40x40"
    folder.mkdir(exist_ok=True)

    md = (SOURCE / "in.md_50K").read_text()
    md = md.replace("50K", f"{temperature}K")
    md = md.replace("50 K", f"{temperature} K")
    md = md.replace("50.0", f"{temperature}.0")
    md = md.replace("505050", str(temperature * 10000 + temperature * 100 + temperature))
    (folder / f"in.md_{temperature}K").write_text(md)

    converter = (SOURCE / "prepare_tdep_50K.py").read_text()
    converter = converter.replace("50K", f"{temperature}K")
    converter = converter.replace("50 K", f"{temperature} K")
    converter = converter.replace("50.0", f"{temperature}.0")
    (folder / f"prepare_tdep_{temperature}K.py").write_text(converter)
    print(f"Prepared inputs for {temperature} K")
