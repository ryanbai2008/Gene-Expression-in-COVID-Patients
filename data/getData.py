import pandas as pd
import numpy as np
import os

def getData(gene_expression_folder_name=r'./GSE217948_RAW', severity_folder_name=r'./[trimmed] data.xlsx'):
    # setup folder names list
    GENE_EXPRESSION_FOLDER_NAME = gene_expression_folder_name
    SEVERITY_FOLDER_NAME = severity_folder_name

    # check folder names are valid
    assert not GENE_EXPRESSION_FOLDER_NAME == None
    assert not SEVERITY_FOLDER_NAME == None

    # get dimensions & setup numpy array for data
    file_list = os.listdir(GENE_EXPRESSION_FOLDER_NAME)
    n = len(file_list)
    with open(os.path.join(GENE_EXPRESSION_FOLDER_NAME,'GSM6730791_HC22.counts.txt'), 'r') as file_content:
        d = len(file_content.readlines())
    data = np.zeros((n,d))

    # set actual values for the data array (np array)
    for file_number in range(n):
        file_name = file_list[file_number]
        with open(os.path.join(GENE_EXPRESSION_FOLDER_NAME, file_name), 'r') as file_content:
            added_data = []
            for line in file_content:
                magnitude_expression = line.split()[1]
                added_data.append(magnitude_expression)
            data[file_number] = added_data

    # set up severity values and identifiers
    identifiers = []
    for identifier in file_list:
        identifiers.append(identifier[11:].split('.')[0])
    severity = np.full(n, fill_value=-1, dtype = float)

    # get the gene patient severity values
    severity_file = pd.read_excel(SEVERITY_FOLDER_NAME)
    for i in range(len(severity_file['Sample identifier'])):
        identifier = severity_file['Sample identifier'][i]
        indexAppear = identifiers.index(identifier)
        severity[indexAppear] = severity_file['WHO severity'][i]

    data_array = []
    severity_array = []
    for i in range(len(severity)):
        severity_magnitude = severity[i]
        if not severity_magnitude == -1:
            data_array.append(data[i])
            severity_array.append(severity_magnitude)

    return data_array, severity_array
