import pandas as pd
from utils import find_rdkit_fragments, find_benzyl_amines
import plotly.express as px
from plotly.offline import plot


def first_class(row):
    if row['fr_Ar_N'] > 0 or row['fr_thiophene'] > 0:
        return 'Heterocyclic aromatic amines'
    if row['fr_benzene'] > 0:
        if row['fr_aniline'] > 0:
            return 'Anilines'
        return 'Phenylalkyl amines'
    if row['RDKIT_SMILES'].find("C1") != -1:
            return 'Cyclic aliphatic amines'
    return "Noncyclic aliphatic amines"


def second_class_heterocyclic(row):
    if row['fr_thiophene'] > 0:
        return 'Thiophenes'
    if row['fr_Nhpyrrole'] > 0:
        return 'Pyrroles'
    if row['fr_pyridine'] > 0:
        return 'Pyridines'
    if row['fr_imidazole'] > 0:
        return 'Imidazoles'
    if row['fr_thiazole'] > 0:
        return 'Thiazoles'
    return 'Other heterocyclic amines'


def second_class_phenylalkyl(row):
    if row['aromatic_piperazine'] > 0:
        return 'Aromatic piperazines'
    if row['2phenylpropyl'] > 0:
        return '2-Phenylpropylamines'
    if row['1phenylethyl'] > 0:
        return '1-Phenylethylamines'
    if row['amphetamine'] > 0:
        return 'Amphetamines'
    if row['n_alkylbenzyl'] > 0:
        return 'N-Alkyl benzylamines'
    if row['n_methylphenylethyl'] > 0:
        return 'N-Methyl phenylethylamines'
    if row['phenylbutyl'] > 0:
        return 'Phenylbutylamines'
    if row['phenylpropyl'] > 0:
        return 'Phenylpropylamines'
    if row['phenylethyl']>0:
        return 'Phenylethylamines'
    if row['phenylmethyl']>0:
        return 'Benzylamines'
    return 'Other'


if __name__ == "__main__":
    # Read UMAP result
    df_scores = pd.read_csv("processed_data/final_db_with_UMAP.csv")
    # Find RDKIT fragments
    df_scores_fingers = find_rdkit_fragments(df_scores.copy())
    df_scores_fingers = df_scores_fingers.loc[:, (df_scores_fingers != 0).any(axis=0)]
    # See above for "classify" function
    df_scores_fingers['Class'] = df_scores_fingers.apply(first_class, axis=1)

    # main class classification
    fig = px.scatter(df_scores_fingers, x='UMAP-1', y='UMAP-2', color='Class', template='plotly_white',
                     hover_data=["CID",'RDKIT_SMILES'], width=1000, height=600)
    plot(fig, filename="output_html/fig2_amine_atlas_in_classes.html")
    # fig.write_image("output_html/fig2.eps")    # static img output only workable in Mac or Linux
    df_filtered = df_scores_fingers[['CID','RDKIT_SMILES','Type','UMAP-1','UMAP-2','Class']]
    df_filtered.to_csv("output_csv/amine_atlas_in_classes.csv", index=False)

    # Second class: heterocyclic aromatics
    df_hetero = df_scores_fingers[df_scores_fingers['Class'] == 'Heterocyclic aromatic amines']
    df_hetero['Subclass'] = df_hetero.apply(second_class_heterocyclic, axis=1)
    fig = px.scatter(df_hetero, x='UMAP-1', y='UMAP-2', color='Subclass', template='plotly_white',
                     hover_data=["CID",'RDKIT_SMILES'], width=1000, height=600)
    plot(fig, filename="output_html/fig3_amine_atlas_in_subclass_heterocyclic.html")
    df_filtered_hetero = df_hetero[['CID','RDKIT_SMILES','Type','UMAP-1','UMAP-2','Class', 'Subclass']]
    df_filtered_hetero.to_csv("output_csv/amine_atlas_subclass_hetero.csv", index=False)

    # Second class: phenylalkyl amines
    df_phenylalkyl = df_scores_fingers[df_scores_fingers['Class'] == 'Phenylalkyl amines']
    df_phenylalkyl_groups = find_benzyl_amines(df_phenylalkyl)
    df_phenylalkyl_groups['Subclass'] = df_phenylalkyl_groups.apply(second_class_phenylalkyl, axis=1)
    fig = px.scatter(df_phenylalkyl_groups, x='UMAP-1', y='UMAP-2', color='Subclass', template='plotly_white',
                     hover_data=["CID",'RDKIT_SMILES'], width=1000, height=600)
    plot(fig, filename="output_html/fig4_amine_atlas_in_subclass_phenylalkyl.html")
    df_filtered_phenylalkyl_groups = df_phenylalkyl_groups[['CID', 'RDKIT_SMILES', 'Type', 'UMAP-1', 'UMAP-2', 'Class', 'Subclass']]
    df_filtered_phenylalkyl_groups.to_csv("output_csv/amine_atlas_subclass_phenylalkyl.csv", index=False)

    # Second class: aniline
    df_aniline = df_scores_fingers[df_scores_fingers['Class'] == 'Anilines']
    df_aniline['Subclass'] = 'Anilines'
    df_filtered_aniline = df_aniline[['CID', 'RDKIT_SMILES', 'Type', 'UMAP-1', 'UMAP-2', 'Class', 'Subclass']]

    # Second class: noncyclic aliphatics (manual classification)
    df_aliphatic_second = pd.read_csv("input_data/subclass_aliphatics_manual_classification.csv")
    df_aliphatic_noncyclic = df_aliphatic_second[df_aliphatic_second['Class'] == 'Noncyclic aliphatic amines'].copy()
    fig = px.scatter(df_aliphatic_noncyclic, x='UMAP-1', y='UMAP-2', color='Subclass', template='plotly_white',
                     hover_data=["CID", 'RDKIT_SMILES'], width=1000, height=600)
    plot(fig, filename="output_html/fig5_amine_atlas_in_subclass_noncyclic_aliphatic.html")

    df_new_summary = pd.concat([df_filtered_hetero, df_filtered_phenylalkyl_groups, df_filtered_aniline, df_aliphatic_second])
    df_new_summary.to_csv("processed_data/final_db_with_two_level_classifications.csv", index=False)










