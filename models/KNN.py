import numpy as np
import sys
import os

# this is not generalized version of KNN btw, designed only to work on the dataset that we have for this project
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from data.getData import *

# might take a bit to run
X, y = getData()

N, D = X.shape
X_train = X[0 : int(0.9 * N)]
y_train = y[0 : int(0.9 * N)]
X_test = X[int(0.9 * N):]
y_test = y[int(0.9 * N):]

def classify(X_train, y_train, X_test, y_test):
    classification_accuracy_list = np.zeros(len(X_test))
    for index, XyList_test in enumerate(zip(X_test, y_test)):
        data_test, test_class = XyList_test
        out = None
        least_distance = np.sum(np.abs(data_test - X_train[0]))
        for XyList_train in zip(X_train, y_train):
            data_train, train_class = XyList_train
            distance = np.sum(np.abs(np.abs(data_train - data_test)))
            if distance < least_distance:
                least_distance = distance
                out = train_class
        classification_accuracy_list[index] = (out == test_class)
    
    # oops this part repetitive so might be a bit inefficient
    classification_accuracy = np.sum(classification_accuracy_list) / len(classification_accuracy_list)
    return f"classification accuracy is {classification_accuracy}"

# test
print(classify(X_train, y_train, X_test, y_test))


