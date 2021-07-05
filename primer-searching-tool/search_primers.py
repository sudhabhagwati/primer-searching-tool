"""
Usage:
    python search_primers.py <path to *.xls file> <path to *.fa file>

Author:
    Sudha Bhagwati - 05.07.2021
"""


## imports

import sys
from pydna.dseq import Dseq
from pydna.amplify import pcr
from pydna.primer import Primer
from pydna.design import primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.tm import *
import arlpy.plot as aplt
import numpy as np
import json
from progress.spinner import MoonSpinner
import csv

## utility functions

# input : list of keys to be sorted (string)
# output : list of sorted keys (string)
def sortKeys(keys):
    tempkeys = []
    sortedkeys =[]
    for key in keys:
        tempkeys.append(int(key.split(":")[0]))
        tempkeys.append(int(key.split(":")[1]))
    tempkeys = sorted(tempkeys)
    for i in range(len(tempkeys)):
        if i % 2 == 0:
            sortedkeys.append(str(tempkeys[i])+":"+str(tempkeys[i+1]))
        else:
            continue
    return sortedkeys

"""
# Input: template sequence
# Output: primers
"""
def getPrimer(dna, temperature=60.0):
    ampl = primer_design(Dseqrecord(dna), target_tm=temperature)
    return ampl

"""
# Input: template sequence, product length, temperature
# Output: list containing forward and reverse primers
"""
def get_primers_list(dna, product_length, temperature=60.0):
    list_primers_info = []
    for index in range((len(dna)-product_length)+1):
        start_index = index
        end_index = index + product_length
        ampl = getPrimer(dna[start_index:end_index], temperature)
        primers_info = (start_index, end_index, ampl.forward_primer, ampl.reverse_primer)
        list_primers_info.append(primers_info)
    return list_primers_info


"""
# Input: primer sequence
# Output: True if repeats present, False otherwise
"""
def checkRepeats(seq):
    if ("aaaaa" in seq) or ("ttttt" in seq) or ("ggggg" in seq) or ("ccccc" in seq) :
        return True
    else:
        return False

"""
# Input: primer sequence
# Output: True if 3' end with gc, False otherwise
"""
def isgcAtend(seq):
    if (seq[-2:] == "gc"):
        return True

## Read from database

jsonfile = sys.argv[1]
geneid = input("Enter gene ID: ")
genename = input("Enter gene name: ")

# read database (json file)
with open(jsonfile, 'r') as fp:
    data = json.load(fp)

# example search
exon_dict = data[geneid][genename]

# get the sequence
sortedkeys = sortKeys(exon_dict.keys())
sequence = ''
for item in sortedkeys:
    sequence +=  exon_dict[item]

print("This will take a while (about 10 mins)..\n")

## STEP1: Select primers based on product length, temperature and gc
# select all such primers with temperature setting 60 <= t <= 63
# select all such primers which have product length 75 <= pl <= 150
# select all such primers which have 50 <= gc <= 60
firstscreen = []
with MoonSpinner('Peforming the first screening') as bar:
	for pl in range(75, 151):
	    for temp in range(60, 64):
	        primers_list = get_primers_list(sequence, pl, temp)
	        for j in range(len(primers_list)):
	            if (primers_list[j][2].gc()) >= 50 and (primers_list[j][3].gc() >= 50
	            and primers_list[j][2].gc()) < 60 and (primers_list[j][3].gc() < 60):
	                firstscreen.append(primers_list[j])
	bar.next()
print("Number of primers after first screening: " + str(len(firstscreen)) + "\n")

## STEP2: Select primers based on primer length
# select all such primers which have 18 >= primer length <= 25
secondscreen = []
with MoonSpinner('Peforming the second screening') as bar:
	for i in range(len(firstscreen)):
	    if (len(firstscreen[i][2]) >= 18 and len(firstscreen[i][2]) <= 25
	       and len(firstscreen[i][3]) >= 18 and len(firstscreen[i][3]) <= 25):
	        secondscreen.append(firstscreen[i])
	bar.next()
print("Number of primers after second screening: " + str(len(secondscreen)) + "\n")

## STEP3: Discard primers with >4 repeats of single base
thirdscreen = []
with MoonSpinner('Peforming the third screening') as bar:
	for i in secondscreen:
	    # check repeats in fp and rp
	    if (not checkRepeats(i[2].seq)) and (not checkRepeats(i[3].seq)):
	        thirdscreen.append(i)
	bar.next()
print("Number of primers after third screening: " + str(len(thirdscreen)) + "\n")

## STEP4: Prefer primers with GC clamp at 3' end
### both fp and rp have gc
allgc = []
for i in thirdscreen:
    if (isgcAtend(i[2].seq)) and (isgcAtend(i[3].seq)):
        allgc.append(i)

### one of fp or rp have gc
onegc = []
for i in thirdscreen:
    if (isgcAtend(i[2].seq) and not isgcAtend(i[3].seq)) or (not isgcAtend(i[2].seq) and isgcAtend(i[3].seq)):
        onegc.append(i)

### none of them have gc
nogc = []
for i in thirdscreen:
    if (not isgcAtend(i[2].seq)) and (not isgcAtend(i[3].seq)):
        nogc.append(i)

print("Hold on a moment while we save the primers found to a file! \n")

## save to files
with open('allgc.csv','w') as result_file:
    wr = csv.writer(result_file, dialect='excel')
    wr.writerows(allgc)

with open('onegc.csv','w') as result_file:
    wr = csv.writer(result_file, dialect='excel')
    wr.writerows(onegc)

with open('nogc.csv','w') as result_file:
    wr = csv.writer(result_file, dialect='excel')
    wr.writerows(nogc)

# with open("allgc.txt", "w") as output:
#     output.write(str(allgc))

# with open("onegc.txt", "w") as output:
#     output.write(str(onegc))

# with open("nogc.txt", "w") as output:
#     output.write(str(nogc))
