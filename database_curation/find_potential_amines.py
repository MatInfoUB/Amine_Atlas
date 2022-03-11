import pandas as pd
from utils import similar_compounds_from_cid_simple, get_compounds_from_cids, find_all_amine_types_simple
pd.options.mode.chained_assignment = None


# in order to do everything in one file
df_amine = pd.read_csv("data/existing_amines.csv")

df_amine.fillna(-9999, inplace=True)
df_amine['CID'] = df_amine['CID'].astype(int)
print(df_amine.head())
df_amine_output = get_compounds_from_cids(list(df_amine["CID"]))
df_to_concat = df_amine_output[['SMILES', 'AIDs']]
print(df_to_concat.head())
df_origin = pd.concat([df_amine, df_to_concat], axis=1)
print(df_origin.head(30))
df_origin.to_csv("data/existing_amines_SMILES_AIDs.csv", index=False)

# Part 1 fetch the cids
# df_origin = pd.read_csv("data/original_set_with_SMILES_and_AIDs.csv", usecols=["CID"])
cid_list = list(df_origin["CID"])
new_cid_set = set()
for cid in cid_list:
    new_list = similar_compounds_from_cid_simple(cid, similarity_percent=95, number_returned=10)
    for new_cid in new_list:
        new_cid_set.add(new_cid)
    print(str(cid)+" is done")
new_cid_list = list(new_cid_set)
df_new_cid = pd.DataFrame(new_cid_list, columns=['CID'])
# df_new_cid.to_csv("data/expanded_CID_list.csv", index=False)

# Part 2 get the SMILES
# df_new_cid = pd.read_csv("data/expanded_CID_list.csv")
new_cid_list = list(df_new_cid['CID'])
df = get_compounds_from_cids(new_cid_list)
# df.to_csv("data/expanded_set_unfiltered.csv", index=False)

# # Part 3 remove ions and non-amines
print(df.shape)
df = df.dropna(axis=0)
print(df.shape)
df_processed = df[~df['SMILES'].str.contains('\.')]
print(df_processed.shape)
df_processed = find_all_amine_types_simple(df_processed)
df_processed = df_processed[df_processed['fr_all_N'] > 0]

# # part 4 remove amines with no assays
print(df_processed.shape)
df_activity = df_processed[df_processed['AIDs'].astype(str) != '[]']
print(df_activity.shape)

# part 5 remove existing amines
df_true = df_activity[~df_activity["CID"].isin(cid_list)]
print(df_true.shape)
df_true.to_csv("data/potential_amines_SMILES_AIDs_unprocessed.csv", index=False)