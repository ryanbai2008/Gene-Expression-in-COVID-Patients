import pandas as pd
# import numpy as np
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
from sklearn.metrics import classification_report

#loads gene expression data (placeholder)
X = pd.read_csv("placeholder.csv", index_col=0).T


# #test labels
# np.random.seed(42)
# y = np.random.choice([0, 1, 2], size=X.shape[0])

#loads labels
labels_df = pd.read_csv("endotype_labels.csv")  # must have 'Sample', 'Endotype'
labels_df = labels_df.set_index("Sample")

#match order
y = labels_df.loc[X.index]["Endotype"].values

#train/test split
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=42)

#trains XGBoost
model = XGBClassifier(eval_metric="mlogloss", random_state=42)
model.fit(X_train, y_train)

#evaluates
y_pred = model.predict(X_test)
print(classification_report(y_test, y_pred))