import json
import os
from mp_api.client import MPRester

api_key = os.environ["MP_API_KEY"]

with MPRester(api_key) as mpr:
    docs = mpr.materials.phonon.search(material_ids=["mp-1143"])

print("number of docs:", len(docs))

for i, doc in enumerate(docs):
    data = doc.model_dump()
    with open(f"mp-1143-phonon-{i}.json", "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"saved mp-1143-phonon-{i}.json")
    print("keys:", sorted(data.keys()))