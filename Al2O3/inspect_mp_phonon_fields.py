import os
from mp_api.client import MPRester

with MPRester(os.environ["MP_API_KEY"]) as mpr:
    print("AVAILABLE PHONON FIELDS:")
    fields = sorted(mpr.materials.phonon.available_fields)
    for f in fields:
        print(f)

    print("\nPOSSIBLY RELEVANT FIELDS:")
    for f in fields:
        low = f.lower()
        if any(k in low for k in ["band", "dos", "phon", "freq", "q", "branch", "label"]):
            print(f)
