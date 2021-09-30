#!/usr/bin/env python

"""Download AMR data from NCBI

"""
import pandas as pd
from Bio import Entrez
import xmltodict
from collections import Counter

Entrez.email = "adam.rivers@usda.gov"
handle = Entrez.esearch(db="biosample", retmax=10000, term='"antibiogram"[filter] AND "salmonella enterica"[organism]') # or esearch, efetch, ...
record = Entrez.read(handle)
handle.close()

handle2 = Entrez.efetch(db="biosample", retmax=10000, id=record['IdList'], rettype="xml", retmode="xml")

#handle2 = Entrez.efetch(db="biosample", retmax=10, id=record['IdList'], rettype="xml", retmode="xml")
rec2 = handle2.read()

doc = xmltodict.parse(rec2)



def _id_parser(bslist_item):
    d1 = {}
    for rec in bslist_item["Ids"]["Id"]:
        if  "@db" in rec:
            d1[rec["@db"]] = rec["#text"]
    return d1


def _amr_parser(bslist_item):
    d2 = {}
    testtype = None
    try:
        for item in bslist_item["Description"]["Comment"]["Table"]["Body"]["Row"]:
            cell = item["Cell"]
            antibiotic = cell[0]
            if cell[1] == 'susceptible':
                bin = int(0)
            elif cell[1] == 'resistant':
                bin = int(1)
            else:
                bin = None
            ineq = cell[2]
            mic = cell[3].split("/")[0] # report only first value in combo antibiotic
            binname = "bin_" + antibiotic
            micname = "mic_" + antibiotic
            d2[binname] = bin
            d2[micname] = float(mic)
            testtype = cell[5]
            if testtype =="disk diffusion":
                d2 = {}
    except:
        raise
    return testtype, d2


def data_parser(doc):
    rclist = []
    bslist = doc["BioSampleSet"]["BioSample"]
    testtypelist = []
    for count, bslist_item in enumerate(bslist):
        try:
            recdict = {}
            recdict.update(_id_parser(bslist_item))
            testtype, amrdict = _amr_parser(bslist_item)
            testtypelist.append(testtype)
            if bool(amrdict):
                recdict.update(amrdict)
                rclist.append(recdict)
        except:
            print("exception in {}".format(count))
            pass
    print("Test Types:")
    print(set(testtypelist))
    return rclist

data = data_parser(doc)

pddataframe = pd.DataFrame.from_records(data)

pddataframe.to_csv("AMR_Data_NCBI.csv")
