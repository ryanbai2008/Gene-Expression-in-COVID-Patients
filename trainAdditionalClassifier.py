import pandas as pd
# import numpy as np
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
from sklearn.metrics import classification_report
from sklearn.feature_selection import VarianceThreshold
import re


#loads gene expression data (placeholder)

df = pd.read_csv("./data2/normalized_counts_DESeq2(2).csv", index_col=0)

# Transpose so rows = samples, columns = genes
X = df.T.copy()

# Normalize sample names (e.g., nCWES001d1 → nCWES001)
def normalize_sample_name(name):
    return re.sub(r'd\d+$', '', name)

print(X[150:220])
X['Normalized'] = X.index.to_series().apply(normalize_sample_name)
print(X[150:220])

# Load metadata
metadata = pd.read_excel("data/data.xlsx")
metadata = metadata.set_index('Patient_ID')  # Ensure 'Patient_ID' is the index
print(metadata['Status'].unique())

# Merge metadata based on normalized names
X = X.merge(metadata, left_on='Normalized', right_index=True, how='inner')
group_labels = X['Status'].copy()

# Drop extra columns (Normalized + Status)
X.drop(columns=['Normalized', 'Status'], inplace=True)

# Ensure numeric, fill NaNs
X = X.apply(pd.to_numeric, errors='coerce').fillna(0)

# Filter for high variance genes
selector = VarianceThreshold(threshold=0.37)
X_highvar = selector.fit_transform(X)
selected_genes = X.columns[selector.get_support()]
X_highvar_df = pd.DataFrame(X_highvar, index=X.index, columns=selected_genes)


# #test labels
# np.random.seed(42)
# y = np.random.choice([0, 1, 2], size=X.shape[0])

#loads labels
labels_df = pd.read_csv("endotype_labels.csv")  # must have 'Sample', 'Endotype'
labels_df = labels_df.set_index("Sample")

#match order
y = labels_df.loc[X_highvar_df.index]["Endotype"].values

#train/test split
X_train, X_test, y_train, y_test = train_test_split(X_highvar_df, y, stratify=y, random_state=42)

#trains XGBoost
model = XGBClassifier(eval_metric="mlogloss", random_state=42)
model.fit(X_train, y_train)

#evaluates
y_pred = model.predict(X_test)
print(classification_report(y_test, y_pred))