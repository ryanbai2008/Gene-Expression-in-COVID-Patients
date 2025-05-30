import pandas as pd
import numpy as np

# pip install pandas openpyxl
df = pd.read_excel("data/data.xlsx")
print(df.head()) #head shows the first 5 rows of the dataframe, but if argument than it can do more



#sheet name -- tells what sheet to read
df = pd.read_excel("data/data.xlsx", sheet_name="Sheet1")
df

# ------------------------------
# Label Data
# ------------------------------

# Label data is the data that we want to predict (by disease)
mappingdisease_map = {
    'Control': 0,
    'COVID-19': 1,
    'Coronavirus': 1,
    'Influenza_A': 2,
    'Influenza_B': 3,
    'Sepsis': 4,
    'Septic_shock': 5,
    'Co-infection (Bac/viral)': 6,
    'Co-infection (Bac/ViralFungal)': 7,
    'Co-infection (Bac/Viral/Fungal)': 7,
    'Co-infection (Viral/Fungal)': 8,
    'Co-infection (Viral/viral)': 9,
    'Co-infection (viral/viral)': 9,
    'Co-infection (viral/viral/fungal)': 10,
}

df['label'] = df['Status'].map(mappingdisease_map) #add labels to each coulumns
print(df[['Status']].head())
print(df)
df.sort_values(by='label', inplace=True) #sort by the label 
print(df[['Status', 'label']].head())
print(df)


# ------------------------------
# Transform Data
# ------------------------------
data = df.drop(["RecruitmentSite"], axis=1) #drops the columns that are not needed
data = data.drop(["Gender"], axis=1)
data = data.drop(["Age"], axis=1)
data = data.drop(["vaccination"], axis=1)
print(data) 

# data
#data["patient_ID"].plot(kind="bar", figsize=(20, 10), title="Patient ID Distribution")

