import json
import os
from inspect import signature

from mp_api.client import MPRester


MATERIAL_ID = "mp-1143"


def save_object(obj, filename):
    if obj is None:
        print(f"{filename}: no data returned")
        return

    # Prefer the MSONable path when available. This avoids the noisy
    # Pydantic serializer warnings seen for the band-structure object.
    if hasattr(obj, "as_dict"):
        text = json.dumps(obj.as_dict(), indent=2, default=str)
    elif hasattr(obj, "model_dump"):
        text = json.dumps(obj.model_dump(mode="json"), indent=2, default=str)
    elif hasattr(obj, "model_dump_json"):
        text = obj.model_dump_json(indent=2)
    elif hasattr(obj, "json"):
        try:
            text = obj.json(indent=2)
        except TypeError:
            text = json.dumps(obj, indent=2, default=str)
    else:
        text = json.dumps(obj, indent=2, default=str)

    with open(filename, "w") as f:
        f.write(text)
    print(f"saved {filename}")


def relevant_names(obj):
    keys = ("phonon", "band", "dos", "freq")
    return sorted(name for name in dir(obj) if any(k in name.lower() for k in keys))


api_key = os.environ["MP_API_KEY"]

with MPRester(api_key) as mpr:
    print("MPRester phonon-related names:")
    for name in relevant_names(mpr):
        print(" ", name)

    print("\nmaterials.phonon names containing phonon/band/dos/freq:")
    for name in relevant_names(mpr.materials.phonon):
        print(" ", name)

    candidates = [
        "get_phonon_bandstructure_by_material_id",
        "get_phonon_dos_by_material_id",
        "get_phonon_bandstructure_from_material_id",
        "get_phonon_dos_from_material_id",
    ]

    print("\nDirect retrieval attempts:")
    found = False
    for name in candidates:
        if not hasattr(mpr, name):
            print(f"  {name}: not available")
            continue

        found = True
        method = getattr(mpr, name)
        print(f"  {name}{signature(method)}")

        try:
            result = method(MATERIAL_ID)
        except Exception as exc:
            print(f"    call failed: {exc}")
            continue

        print(f"    returned type: {type(result).__name__}")
        suffix = "band" if "band" in name else "dos"
        save_object(result, f"{MATERIAL_ID}-{suffix}.json")

    if not found:
        print("  No direct phonon bandstructure/DOS helper was found on this MPRester.")
