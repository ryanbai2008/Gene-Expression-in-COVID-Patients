import pandas as pd

counts_df = pd.read_excel("data/data.xlsx")
truncated_counts_df = pd.read_csv("data/truncated normal counts.csv")
patient_IDs = truncated_counts_df["Sample"].tolist()

# Filter and select columns properly
counts_df = counts_df.loc[counts_df["Patient_ID"].isin(patient_IDs), ["Patient_ID", "Status"]]

# Convert Status to binary values
counts_df["Status"] = counts_df["Status"].apply(lambda x: 1 if x in ["COVID-19", "Coronavirus"] else 0)

# Save to CSV
counts_df.to_csv("data/traits.csv", index=False)
