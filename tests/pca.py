#Script for principal component analysis (PCA) on a dataset.
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import LabelEncoder


#Read the dataset
dataset = pd.read_excel('data/data.xlsx')
#dataset = dataset.drop(["Gender", "Age", "Patient_ID", "Status"], axis=1)


le = LabelEncoder()
dataset['Status'] = le.fit_transform(dataset['Status'])

#Selects features and labels, X contains the first 13 columns and y the 14th.
X = dataset.iloc[:, 4:].values
y = dataset['Status'].values

#Splits data into training and testing subjects randomly with test size 20% and training 80%
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

sc = StandardScaler()
#Fit the data
X_train = sc.fit_transform(X_train)  # calculates mean and standard deviation, then scales it
X_test = sc.transform(X_test) # Uses same mean and std for test data

# Apply PCA into training and testing set 
pca = PCA(n_components=2) #reduces components to 2 principal components
X_train = pca.fit_transform(X_train) #finds components of training data
X_test = pca.transform(X_test) #projects test data onto components

# how much variance in each principal component
explained_varience = pca.explained_variance_ratio_

classifier = LogisticRegression(random_state=0) # use logistic regression for best fit line
classifier.fit(X_train, y_train)

y_pred = classifier.predict(X_test)

cm = confusion_matrix(y_test, y_pred)

# Predicting the training set
# result through scatter plot
X_set, y_set = X_train, y_train
X1, X2 = np.meshgrid(np.arange(start=X_set[:, 0].min() - 1,
                               stop=X_set[:, 0].max() + 1, step=0.01),
                     np.arange(start=X_set[:, 1].min() - 1,
                               stop=X_set[:, 1].max() + 1, step=0.01))

plt.contourf(X1, X2, classifier.predict(np.array([X1.ravel(),
                                                  X2.ravel()]).T).reshape(X1.shape), alpha=0.75,
             cmap=ListedColormap(('yellow', 'white', 'aquamarine')))

plt.xlim(X1.min(), X1.max())
plt.ylim(X2.min(), X2.max())

for i, j in enumerate(np.unique(y_set)):
    plt.scatter(X_set[y_set == j, 0], X_set[y_set == j, 1],
                color=ListedColormap(('red', 'green', 'blue'))(i), label=j)

plt.title('Logistic Regression (Training set)')
plt.xlabel('PC1')  # for Xlabel
plt.ylabel('PC2')  # for Ylabel
plt.legend()  # to show legend

# show scatter plot
plt.show()