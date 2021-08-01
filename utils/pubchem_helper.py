import pubchempy as pcp
import pandas as pd
import numpy as np


def get_compounds_by_name(name):
    cpd_list = pcp.get_compounds(name, 'name')
    cid_list = []
    for cpd in cpd_list:
        cid_list.append(cpd.cid)
    if len(cid_list) != 0:
        return cid_list[0]
    return -1


# Read a list of CID and output a list of compounds with CID, Name, SMILES, and AIDs
def get_compounds_from_cids(cid_list):
    num = len(cid_list)
    df_list = []
    for cid in cid_list:
        print(num, "compounds left")
        if cid < 0:
            df_list.append({"CID": cid, "Name": np.nan, "SMILES": np.nan, "AIDs": np.nan})
            continue
        compound_list = pcp.get_compounds(cid, 'cid')
        if len(compound_list) == 0:
            df_list.append({"CID": cid, "Name": np.nan,  "SMILES": np.nan, "AIDs": np.nan})
            continue
        if len(compound_list) > 1:
            print("Warning! CID {} has multiple corresponding compounds as listed below:".format(cid))
            for cmd in compound_list:
                print(cmd.name+"\n")
        compound = compound_list[0]
        smiles = compound.canonical_smiles
        name = compound.iupac_name
        aids = compound.aids
        df_list.append({"CID": cid, "Name": name,  "SMILES": smiles, "AIDs": aids})
        print("CID " + str(cid) + " is done")
        num -= 1
    return pd.DataFrame(df_list)


def similar_compounds_from_cid_simple(cid, similarity_percent=90, time_limit=60, number_returned=100):
    """
    input a compound and return a dataframe of similar compounds
    :param cid: the cid of the compound to do similarity search on
    :param similarity_percent: the threshold of similarity (e.g. 90, 95, 98)
    :param time_limit: how long the search can last
    :param number_returned: how many compounds to return
    :return: a pandas df containing three columns: cid, name, smiles
    """
    compound_list = pcp.get_compounds(cid, 'cid', searchtype='similarity', threshold=similarity_percent,
                                      MaxSeconds=time_limit,
                                      listkey_count=number_returned)
    cid_list = [compound.cid for compound in compound_list]
    return cid_list


# Read a list of AID and output the AID, Name, and Description
def get_assay_info(aid_list):
    df_list = []
    for aid in aid_list:
        if aid < 0:
            df_list.append({"AID": aid, "Name": np.nan, "Description": np.nan})
            continue
        assay_list = pcp.get_assays(aid, u'aid')
        if len(assay_list) == 0:
            df_list.append({"AID": aid, "Name": np.nan, "Description": np.nan})
            continue
        if len(assay_list) > 1:
            print("Warning! AID {} has multiple corresponding assays as listed below:".format(aid))
            for asy in assay_list:
                print(asy.name+"\n")
        assay = assay_list[0]
        name = assay.name
        # results = assay.results
        description = assay.description
        df_list.append({"AID": aid, "Name": name, "Description": description})
    return pd.DataFrame(df_list)