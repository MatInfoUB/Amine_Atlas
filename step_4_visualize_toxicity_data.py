import requests
import pandas as pd
import numpy as np
import plotly.express as px
from plotly.offline import plot


def hit_ratio(row):
    active_num = row['Number_of_active']
    inactive_num = row['Number_of_inactive']
    assay_num = active_num + inactive_num
    if assay_num == 0:
        return -100
    if assay_num < 2:
        return -1
    return active_num/assay_num


if __name__ == "__main__":

    df_amines = pd.read_csv("processed_data/final_db_with_two_level_classifications.csv")
    df_assays = pd.read_csv("input_data/selected_bioassays_info_preprocessed.csv")
    assay_list = df_assays['AID'].to_list()
    for assay in assay_list:
        df = pd.read_csv('input_data/processed_AID_datatables/processed_AID_' + str(assay) + '.csv',
                         usecols=['PUBCHEM_CID', 'PUBCHEM_ACTIVITY_SCORE'])
        df = df.rename(columns={'PUBCHEM_CID': 'CID', 'PUBCHEM_ACTIVITY_SCORE': 'PUBCHEM_ACTIVITY_SCORE_AID_'+str(assay)})
        df_amines = df_amines.merge(df,how='left', on='CID')
    df_amines.set_index(["CID","RDKIT_SMILES","Type","UMAP-1","UMAP-2","Class","Subclass"], inplace=True)
    df_nan_to_zero = df_amines.copy().fillna(0)   # convert all nan to zero
    df_nan_to_minus1 = df_amines.copy().fillna(-1)  # convert all nan to -1
    assay_num = len(df_amines.copy().columns)
    df_amines['Number_of_inactive'] = assay_num - np.count_nonzero(df_nan_to_minus1, axis=1)
    df_amines['Number_of_active'] = np.count_nonzero(df_nan_to_zero, axis=1)
    df_amines['Hit_ratio'] = df_amines.apply(hit_ratio, axis=1)
    df_amines.reset_index(inplace=True)
    df_activity = df_amines[df_amines['Number_of_inactive'] + df_amines['Number_of_active'] > 0]
    df_hit_ratio = df_activity[df_activity['Hit_ratio'] > -1]
    df_hit_ratio_for_plot = df_hit_ratio.copy()
    df_hit_ratio_for_plot['Size'] = df_hit_ratio_for_plot['Hit_ratio'] + 0.1

    fig = px.scatter(df_hit_ratio_for_plot, x='UMAP-1', y='UMAP-2', color='Hit_ratio', size='Size', template='plotly_white',
                     hover_data=["CID", 'RDKIT_SMILES', 'Type', 'Class', 'Subclass'],
                     color_continuous_scale=px.colors.sequential.Bluered, opacity=0.7, width=1000, height=600)
    plot(fig, filename="output_html/fig6_amine_atlas_with_toxicity_data.html")
    df_hit_ratio.to_csv("output_csv/amine_atlas_with_toxicity_data.csv", index=False)

df_vap = pd.read_csv('input_data/vapor_office.csv')

df_vap = df_vap[~df_vap['volatile'].isin([-1])]
print(df_vap)
df_vap['volatile_lg'] = np.log10(df_vap['volatile'])

fig = px.scatter(df_vap, x='UMAP-1', y='UMAP-2', color='volatile_lg', template='plotly_white',
                 hover_data=["CID", 'RDKIT_SMILES', 'Type', 'Class', 'Subclass'],
                 color_continuous_scale=px.colors.sequential.Bluered, opacity=0.7, width=1000, height=600)
plot(fig, filename="output_html/fig21_amine_atlas_with_toxicity_data.html")

