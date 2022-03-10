import pandas as pd
from utils import find_rdkit_fragments, find_benzyl_amines, find_heterocyclic
import plotly.express as px
from plotly.offline import plot
pd.options.mode.chained_assignment = None


def first_class(row):
    if row['fr_Ar_N'] > 0 or row['fr_thiophene'] > 0:
        return 'Heteroaromatic amines'
    if row['fr_benzene'] > 0:
        if row['fr_aniline'] > 0:
            return 'Anilines'
        return 'Phenylalkyl amines'
    return "Aliphatic amines"


# Changed the rules for heteroaromatics subclasses 03/09/2022
def second_class_heterocyclic(row):
    if row['imidazole'] > 0:
        return 'Imidazoles'
    if row['indole'] > 0:
        return 'Indoles'
    if row['triazole'] > 0:
        return 'Triazoles'
    if row['thiophene'] > 0:
        return 'Thiophene amines'
    if row['fr_thiazole'] > 0:
        return 'Thiazoles'
    if row['pyridine'] > 0:
        return 'Pyridines'
    return 'Other heteroaromatic amines'


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

    # for UMAP, please comment for t-SNE
    version = "UMAP_"+ "15_025_0310"
    df_scores = pd.read_csv("processed_data/final_db_with_"+version+".csv")
    x_axis = 'UMAP-1'
    y_axis = 'UMAP-2'
    #
    # for t-SNE, please comment for UMAP
    # df_scores = pd.read_csv("processed_data/final_db_with_TSNE_70_random_0309.csv")
    # x_axis = 'TSNE-1'
    # y_axis = 'TSNE-2'
    # version = "0309_TSNE_70_random"

    # Find RDKIT fragments
    df_scores_fingers = find_rdkit_fragments(df_scores.copy())
    df_scores_fingers = df_scores_fingers.loc[:, (df_scores_fingers != 0).any(axis=0)]
    # See above for "classify" function
    df_scores_fingers['Class'] = df_scores_fingers.apply(first_class, axis=1)

    # main class classification
    fig = px.scatter(df_scores_fingers, x=x_axis, y=y_axis, color='Class', template='plotly_white',
                     hover_data=["CID",'RDKIT_SMILES'], width=1000, height=600)
    plot(fig, filename="output_html/fig2_amine_atlas_main_classes_"+version+".html")
    # fig.write_image("output_html/fig2.eps")    # static img output only workable in Mac or Linux

    # Second class: heterocyclic aromatics
    df_hetero = df_scores_fingers[df_scores_fingers['Class'] == 'Heteroaromatic amines']
    df_hetero_groups = find_heterocyclic(df_hetero)
    df_hetero_groups['Subclass'] = df_hetero_groups.apply(second_class_heterocyclic, axis=1)
    fig = px.scatter(df_hetero_groups, x=x_axis, y=y_axis, color='Subclass', template='plotly_white',
                     hover_data=["CID",'RDKIT_SMILES'], width=1000, height=600)
    plot(fig, filename="output_html/fig3_amine_atlas_in_subclass_hetero_"+version+".html")

    # Second class: phenylalkyl amines
    df_phenylalkyl = df_scores_fingers[df_scores_fingers['Class'] == 'Phenylalkyl amines']
    df_phenylalkyl_groups = find_benzyl_amines(df_phenylalkyl)
    df_phenylalkyl_groups['Subclass'] = df_phenylalkyl_groups.apply(second_class_phenylalkyl, axis=1)
    fig = px.scatter(df_phenylalkyl_groups, x=x_axis, y=y_axis, color='Subclass', template='plotly_white',
                     hover_data=["CID",'RDKIT_SMILES'], width=1000, height=600)
    plot(fig, filename="output_html/fig4_amine_atlas_in_subclass_phenylalkyl_"+version+".html")

    # second class: aliphatic
    df_aliphatic = df_scores_fingers[df_scores_fingers['Class'] == 'Aliphatic amines']
    print(df_aliphatic.head())
    df_aliphatic_classification = pd.read_csv("input_data/subclass_aliphatics_manual_classification_new.csv", usecols=['CID', 'Subclass'])
    print(df_aliphatic_classification.head())
    df_aliphatic_groups = pd.merge(df_aliphatic, df_aliphatic_classification, how="left", on="CID")
    fig = px.scatter(df_aliphatic_groups, x=x_axis, y=y_axis, color='Subclass', template='plotly_white',
                     hover_data=["CID", 'RDKIT_SMILES'], width=1000, height=600)
    plot(fig, filename="output_html/fig5_amine_atlas_in_subclass_aliphatic_"+version+".html")


















