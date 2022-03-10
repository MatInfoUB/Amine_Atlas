import requests
import pandas as pd
import numpy as np
import plotly.express as px
from plotly.offline import plot
import math


df_vap = pd.read_csv('input_data/selected_vapor_pressure.csv')

df_vap = df_vap[~df_vap['volatile'].isin([-1])]
print(df_vap)
df_vap['volatile_lg'] = np.log10(df_vap['volatile'])

fig = px.scatter(df_vap, x='UMAP-1', y='UMAP-2', color='volatile_lg', template='plotly_white',
                 hover_data=["CID", 'RDKIT_SMILES', 'Type', 'Class', 'Subclass'],
                 color_continuous_scale=px.colors.sequential.Bluered, opacity=0.7, width=1000, height=600)
plot(fig, filename="output_html/fig21_amine_atlas_with_toxicity_data.html")
