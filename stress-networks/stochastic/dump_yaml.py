import pandas as pd
import yaml
import os

# Data from the table
data = {
    "Case ID": [
        "Case_0",
        "Case_1a",
        "Case_1b",
        "Case_1c",
        "Case_2a",
        "Case_2b",
        "Case_2c",
        "Case_3a",
        "Case_3b",
    ],
    "sigma_xx": [0e6, 1e6, 1e6, 3e6, 5e6, 5e6, 15e6, 5e6, 20e6],
    "sigma_yy": [0e6, 1e6, 3e6, 1e6, 5e6, 15e6, 5e6, 20e6, 5e6],
    "sigma_zz": [0e6, 1e6, 1e6, 1e6, 5e6, 5e6, 5e6, 5e6, 5e6],
}

df = pd.DataFrame(data)

# Where to write the YAML files
out_dir = "cases_yaml"
os.makedirs(out_dir, exist_ok=True)

for _, row in df.iterrows():
    case_dict = {
        "case_id": row["Case ID"],
        "sigma_xx": float(row["sigma_xx"]),
        "sigma_yy": float(row["sigma_yy"]),
        "sigma_zz": float(row["sigma_zz"]),
    }

    out_path = os.path.join(out_dir, f"{row['Case ID']}.yaml")
    with open(out_path, "w") as f:
        yaml.dump(case_dict, f, default_flow_style=False, sort_keys=False)

    print(f"Wrote {out_path}")
    