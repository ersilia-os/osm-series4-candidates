# Open Source Malaria Series 4 Candidates by Ersilia Open Source Initiative

## In a nutshell

1. We trained predictive models of activity based on the Series 4 experimental screening data from OSM.
2. We used reinforcement learning to populate the chemical space around Series 4 compounds.
3. We assembled and filtered results from multiple generative runs.


### Results

* [Explore](https://share.streamlit.io/ersilia-os/osm-series4-candidates/main/app.py) a set of selected 1000 candidates.
* [Download](https://github.com/ersilia-os/osm-series4-candidates/blob/main/postprocess/210530_EOSI_OSM_Series4_All.csv) a CSV file containing >100k candidates for a deeper exploration. This set contains clusters of highly similar molecules.
* Navigate the set of [>100,000 candidates](https://ersilia-os.github.io/osm-series4-candidates/tmap/osm_eosi_s4.html) using a tree map.
* See the [Colaboratory Notebook](https://colab.research.google.com/drive/1MK4UJP6Vw1FjaVaTu9CwpBb2XrjPw9xz?usp=sharing) of the generative model.
* Learn more in the OSM GitHub Issue [#33](https://github.com/OpenSourceMalaria/Series4_PredictiveModel/issues/33).

### Disclaimer

This work is exploratory and we did our best to achieve molecules with desirable properties. We will be happy to see run further generative batches based on feedback by OSM. We may optimise for synthetic accessibility, solubility, drug-likeness, etc., in order to achieve better Series 4 derivates.

## Introduction

This project aims to contribute new candidates to the [Open Source Malaria Series 4](https://github.com/OpenSourceMalaria/Series4) hit-to-lead optimisation. Series 4 compounds correspond to an already well-defined chemical space, with most compounds maintaining a core triazolopyrazine with NE and NW substituents.
Our approach to the issue of generating more derivatives with improved activity and solubility has been to couple different activity predictors to the [Reinvent 2.0](https://github.com/MolecularAI/Reinvent) generative model using reinforcement learning.
Below you can find an overview of the project as well as a summary of the different molecule batches we have generated in subsequent runs of the generative model.

## Project Overview
1. We have downloaded the experimental data for Series 4 compounds from the [Master List](https://github.com/OpenSourceMalaria/Series4/wiki/Sources-of-Data) in OSM GitHub (20/05/21).
2. The data has been used to train different Activity (IC50) predictors, which have subsequently been incorporated as reinforcement learning agents in generative models (23-27/05/21)
3. Six different batches of new Series 4 candidates have been generated as explained below (28-29/05/21). Not all batches were optimised for activity (see below) - some were aimed at exploring the chemical space.

## Candidate lists
### Batch00 and Batch01
We have tried to narrow the chemical space to series 4 compounds, maintaining the triazolopyrazine core with NW and NE substituents. This has been achieved by forcing matching to the core substructure and similarity to series 4 molecules. In addition, core structures with substituents at undesired atoms have been restricted. QED and solubility (SLogP), as well as synthetic (SA) and retrosynthetic accessibility (RA) have been included as scoring components. Finally, activity scoring has been achieved with a simple Random Forest Classifier and Regressor with optimised hyperparameters (AutoML).
Batch00 rendered 1867 compounds, and Batch01 rendered 1138 compounds. Due to high chemical constraints to maintain similarity to series 4 compounds, 500 of the generated molecules were common between batches, and, unfortunately, molecular weight (MW) of new molecules was lower (≈350g/mol) than those of active series 4 compounds (≈450g/mol).

### Batch02
With Batch02 we tried to explore the chemical space starting from Batch00 and Batch01 and relaxing the restrictions of similarity to series 4 compounds. The goal was to increase molecular weight (adding a molecular weight component) in order to diversify the substituents.
Batch02 includes 46349 compounds with molecular weight around 400-500 and more diversity (Tanimoto Similarity in pairwise comparison is ≈ 0.65). No activity constraints were added here.

### Batch03
Batch03 was an enlargement of Batch02 including more restrictions for desired physicochemical characteristics (QED, SLogP, SA and RA). It generated 56700 compounds.

### Batch04 and Batch05
Finally, taking as priors the diversified chemical space from series 4, we focused our efforts on obtaining compounds with improved predicted activity. To this end, we used a [Chemprop](https://github.com/chemprop/chemprop) model (Batch04) and the original AutoML predictor (Batch05).
Batch04 includes 8833 molecules and Batch05 includes 5352 molecules. We put strong scaffold constraints to keep diversity, which resulted in a small number of molecules predicted to be highly active.

## Results
We provide the following data for each candidate in this [file](https://github.com/ersilia-os/osm-series4-candidates/blob/main/postprocess/210530_EOSI_OSM_Series4_All.csv):
* EosID: new series 4 candidates have been annotated as EOS-{BATCH}-{NUMBER}
* Smiles: canonical smiles
* InChIKey
* MolWeight: molecular weight (g/mol)
* NumRings: number of rings
* Batch: batch of generation (0 to 5)
* Calculated scores for:
    * QED (0:not drug like - 1: max similarity to known drugs)
    * SLogP: octanol-water partition coefficient as a measure of solubility. Lower values indicate higher hydrophilicity
    * SAScore (0-1): synthetic accessibility (RDKit)
    * RAScore (0-1): retrosynthetic accessibility from the [Reymond Group](https://github.com/reymond-group/RAscore)
    * Tanimoto: similarity with existing series 4 molecules (0: new - 1: repeated molecule)
    * Triazolopyrazine core maintenance (0: no core or 1: core)
    * Activity against P. Falciparum. **Please note that scores are a relative ranking, not a probability**
        * ActivityClfAutoML (0: inactive - 1: active): Random Forest classifier
        * ActivityRegAutoML (0: inactive - 1: active): Random Forest regressor
        * ActivityClfGraph (0: inactive - 1: active): Chemprop classifier
        * ActivityRegGraph (0: inactive - 1: active): Chemprop regressor
        * ActivityAvg (0: inactive - 1: active): average of the four activity prediction scores
* Clustering: the candidates have been clustered with K-Means according to chemical similarity in:
    * Cluster100: 100 clusters
    * Cluster1000: 1000 clusters

In addition, we provide a [web-based interface](https://share.streamlit.io/ersilia-os/osm-series4-candidates/main/app.py) for the 1000 molecules (one from each cluster) with highest predicted activity (using average activity). If there is a particular cluster of interest please refer to the full excel file to select other molecules from the same cluster.

Finally, the 116728 molecules are represented in this [TreeMap](https://ersilia-os.github.io/osm-series4-candidates/tmap/osm_eosi_s4.html).


## About us

[Ersilia Open Source Initiative](https://ersilia.io) is small non-profit aimed at disseminating AI/ML for drug discovery. We have a special interest in infectious and neglected diseases of low- and middle-income countries.
