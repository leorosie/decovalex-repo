import pandas as pd
import subprocess 

# Data from the table
data = {
    "Case ID": [
        "Case_0",
        "Case_1a",
        "Case_1b",
        "Case_1c",
        "Case_2a",
        "Case_2b",
        "Case_2c"
    ],
    "sigma_xx": [
        0 *1e6,
        1 *1e6,
        1 *1e6,
        3 *1e6,
        5 *1e6,
        5 *1e6,
        15 *1e6
    ],
    "sigma_yy": [
        0 *1e6,
        1 *1e6,
        3 *1e6,
        1 *1e6,
        5 *1e6,
        15 *1e6,
        5 *1e6
    ],
    "sigma_zz": [
        0 *1e6,
        1 *1e6,
        1 *1e6,
        1 *1e6,
        5 *1e6,
        5 *1e6,
        5 *1e6
    ]
}

# Convert the data into a pandas DataFrame
df = pd.DataFrame(data)

# Display the DataFrame
print(df)

# Print each row individually
print("\nRows of the DataFrame:")
for index, row in df.iterrows():
    #print(f"Row {index}: {row.to_dict()}")
    #print(f"{row['Case ID']}")
    cmd = f"python3.11 apply_stress_4frac.py {row['Case ID']} {row['sigma_xx']} {row['sigma_yy']} {row['sigma_zz']}"
    print(cmd) 
    subprocess.call(cmd, shell = True)
