# Amine_Atlas
 This repository is for the manuscript "AI Informed Toxicity Screening of Amine Chemistries used in the Synthesis of Hybrid Organic-Inorganic Perovskites" which is currently under revision. The preprint of an earlier version is available at Authorea: https://doi.org/10.22541/au.163251254.47236038/v1 

The revised manuscript is comming soon! 

 
# Amine_Atlas

   - **Installation**
   - **How to Use Amine_Atlas**
    
**Installation**
    
    1.Do not worry about operating systems as Amine_Atlas has been tested on macOS Catalina, Linux AS-Map has been tested on macOS Catalina, 
      Linux(Ubuntu LTS 20.04), and Windows 10.
    2.It is essential to install Pycharm first if Pycharm is not installed on your computer.
      https://www.jetbrains.com/pycharm/download/#section=windows
      Click on the download page of the community edition, it is free and for pure Python development.
    3.Install Anaconda for Python 3.7 if your computer do not have Anaconda for Python 3.7.
      https://www.anaconda.com/products/individual
    4.In windows search box, find “Anaconda Prompt”, select “run with administrator”,input:
      "conda create -c rdkit -n amine rdkit python=3.7"
      This will setup a conda virtual environment "amine" that you will use to run Amine_Atlas.
    5.When the amine environment is properly setup, input:
      "conda activate amine"
      Now you are in the amine virtual environment.Let's install the basic machine learning and visualization libraries.
    6.To install rdkit, input:
      "conda install -c conda-forge rdkit"
    7.To install pubchempy, input:
      “conda install -c mcs07 pubchempy"
    8.Next, install the required libraries,input:
      "conda install pandas numpy scikit-learn matplotlib plotly seaborn requests nb_conda_kernels"
    9.Start exploring Amine_Atlas! After entering PyCharm, click File→Open and select the folder to import the project 
      in the popup window.After the Python project is started, you need to configure the Python corresponding to the project
      to run properly.Follow these steps:
      "File → settings→ → Project → Python Interpreter → Add → conda environments → Existing environment → ok
    10.Start exploring Amine_Atlas!If there is any problem in the installation, please contact ansuzjut@outlook.com. 
     
    Note: This project can be used not only with Pycharm but also with another IDE or Command Line, but it is important to 
          switch to the Conda Environment first!

  
**How to Use Amine_Atlas**
    
    1."step_1_calculate_mhfp.py": In order to compute molecular fingerprints.
    2."step_2_calculate_umap.py": Use UMAP to reduce the dimension of the data.
    3."step_3_classification_atlas_visualization.py": Classify and visualize perovskite amines.
    4."step_4_visualize_toxicity_data.py": Visualization of toxicity data for perovskite amines.

**Authors**
    
    An Su, Associate Professor of Research, Zhejiang University of Technology; Former postdoctoral associate, University at Buffalo
    Yingying Cheng & Haotian Xue, Graduate students, Zhejiang University of Technology
    Yuanbin She, Professor, Zhejiang University of Technology
    Krishna Rajan, Erich Bloch Chair & Empire Innovation Professor, University at Buffalo
