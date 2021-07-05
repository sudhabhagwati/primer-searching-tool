"""
Usage:
    python generate_database.py <path to *.json file> 

Author:
    Sudha Bhagwati - 05.07.2021
"""

## imports

import sys
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import matplotlib.pyplot as plt
import numpy as np
import json
from progress.spinner import MoonSpinner


## utility functions

'''
input: row of the dataframe
output: ID as string
'''
def extractID(row):
    idx = row['protein'].find("ID=")
    return row['protein'][idx+3:].split(";")[0]

'''
input: row of the dataframe
output: product name as string
'''
def extractProductName(row):
    idx = row['protein'].find("Product=")
    return row['protein'][idx+8:].split(";")[0]

# extract exons from the excel data
# returns dictionary with key being id and value being a list of 'start:end' string of exons
def extractExons(data):
    exons = {}
    for index, row in data.iterrows():
        if row['type'] == 'mRNA':
            ID = extractID(row)
            exons[ID] = {}
            productname = extractProductName(row)
            exons[ID][productname] = []
        if row['type'] == 'CDS':
            exons[ID][productname].append(str(row['start'])+':'+str(row['end']))
    return exons

def extracttempExons(data):
    exons = {}
    for index, row in data.iterrows():
        if row['type'] == 'mRNA':
            ID = extractID(row)
            exons[ID] = []
        if row['type'] == 'CDS':
            exons[ID].append(str(row['start'])+':'+str(row['end']))
    return exons

def findlinenumber(ID, lines):
    for idx, line in enumerate(lines):
        if '>' in line:
            if ID in line:
                return idx
            
def findsequence(ID, lines):
    idx = findlinenumber(ID, lines)
    sequence = ""
    for line in lines[idx+1:]:
        if '>' in line:
            break
        sequence += line
    return sequence

def getsequence(lines, molid, startend, indexes):
    st = 0
    seq = findsequence(molid, lines)
    for idx in indexes:
        if startend in idx:
            return seq[st:st+(int(idx.split(':')[1])-int(idx.split(':')[0]))+1]
        else:
            st += (int(idx.split(':')[1])-int(idx.split(':')[0]))+1

## Database generation

xlsfile = sys.argv[1]
fastafile = sys.argv[2]

print("\n Using " + sys.argv[1] + " and " + sys.argv[2] + " to generate output.json file \n")
print("This will take about 5 mins....please wait while output.json is being created!\n")

# read the genomics data from the excel file
df = pd.read_excel(sys.argv[1], sheet_name='genomic')
extractedExons = extractExons(df)
extractedtempExons = extracttempExons(df)

# read FASTA file contents
file1 = open(sys.argv[2], 'r')
Lines = file1.readlines()
file1.close()
Lines = [line.strip() for line in Lines]

output = {}
with MoonSpinner('Processing…') as bar:
	for molid in extractedExons.keys():
	    output[molid] = {}
	    for protein in extractedExons[molid]:
	        output[molid][protein] = {}
	        for startend in extractedExons[molid][protein]:
	            start = int(startend.split(':')[0])
	            end = int(startend.split(':')[1])
	            output[molid][protein][startend] = getsequence(Lines, molid, startend, extractedtempExons[molid])
	    bar.next()

print("Please wait while the file is being saved!\n")
# save the output to json file
with open('output.json', 'w') as fp:
    json.dump(output, fp, sort_keys=True, indent=4)

print("Succesfully created and saved output.json!\n")
