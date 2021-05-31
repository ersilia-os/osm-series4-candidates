# Open Source Malaria Series 4 Candidates by Ersilia Open Source Initiative

## In a nutshell

1. We trained predictive models of activity based on Series 4 experimental screening.
2. We used reinforcement learning to populate the chemical space around Series 4 compounds.
3. We assembled and filtered results from multiple generative runs.

### Results

* [Explore](https://share.streamlit.io/ersilia-os/osm-series4-candidates/main/app.py) a set of 1000 candidates selected
* [Download](https://github.com/ersilia-os/osm-series4-candidates/blob/main/postprocess/210530_EOSI_OSM_Series4_All.csv) a CSV file containing >100k candidates for a deeper exploration. This set has redundancy.
* Navigate a set the set of [1000 candidates](https://ersilia-os.github.io/osm-series4-candidates/tmap/map_1k/index.html) and the set of [>100,000 candidates](https://ersilia-os.github.io/osm-series4-candidates/tmap/map_100k/index.html) using a tree map.

### Disclaimer

This work is exploratory and we did our best to achieve molecules with desirable properties. We will be happy to see run further generative batches based on feedback by OSM. We may optimise for synthetic accessibility, solubility, drug-likeness, etc, in order to achieve better series 4 derivates.

## Introduction

This project aims to contribute new candidates to the [Open Source Malaria Series 4](https://github.com/OpenSourceMalaria/Series4). Series 4 compounds correspond to an already well-defined chemical space, with most compounds maintaining a core triazolopyrazine with NE and NW substituents.
Our approach to the issue of generating more derivatives with improved activity and solubility has been to couple different activity predictors to a generative model (Reinvent 2.0) using reinforcement learning.
Below you can find an overview of the project as well as a summary of the different molecule batches we have generated in subsequent runs of the Reinvent2.0

## Project Overview
1. We have downloaded the experimental data for Series 4 compounds from the [Master List](https://github.com/OpenSourceMalaria/Series4/wiki/Sources-of-Data) in OSM GitHub (20/05/21).
2. The data has been used to train different Activity (IC50) predictors, which have subsequently been incorporated as scoring components of the [REINVENT 2.0](https://github.com/MolecularAI/Reinvent) Reinforcement Learning generative model (23-27/05/21)
3. Six different batches of new series 4 molecules have been generated as explained below (28-29/05/21).

## New candidate lists
### Batch00 and Batch01
We have tried to narrow the chemical space to series 4 compounds, maintaining the triazolopyrazine core with NW and NE substituents. This has been achieved by forcing matching to the core substructure and similarity to series 4 molecules. In addition, core structures with substituents at undesired atoms have been restricted. QED Beauty and solubility (SLogP), as well as synthetic (SA) and retrosynthetic accessibility (RA) have been included as scoring components. Finally, activity scoring has been achieved with an autoML random forest classifier and regressor.
Batch00 rendered 1867 compounds (non-duplicates from original series 4), and Batch01 rendered 1138 compounds (non-duplicates from original series 4).
Due to high chemical constraints to maintain similarity to series 4 compounds, 500 of the generated molecules are common amongst both batches, and, interestingly, molecular weight (MW) of new molecules is lower (≈350g/mol) than those of active series 4 compounds (≈450g/mol)

### Batch02
With Batch02 we have tried to explore the chemical space starting from Batch00 and Batch01, eliminatig the restrictions of similarity to series 4 compounds. The goal is to increase molecular weight (adding a molecular weight component) in order to diversify the substituents.
Batch02 includes 46349 compounds with molecular weight around 400-500 and more diversity (Tanimoto Similarity in pairwise comparison is ≈ 0.65).

### Batch03
Batch03 is an enlargement of Batch02 including more restrictions for desired physicochemical characteristics (QED_Score, SLogP, SA and RA). It includes 56700 compounds.

### Batch04 and Batch05
Finally, with a diversified chemical space from series 4, we have focused our efforts on obtaining compounds with improved predicted activity. To this end, we have used a Chemprop model trained on the 100K datapoints from previous experiments (Batch04) and the original AutoML predictor (Batch05).
Batch04 includes 8833 molecules and Batch05 includes 5352 molecules.

## Results
An excel containing all candidates with the following information can be found [here]():
* EosID: new series 4 candidates have been annotated as EOS - #batch - #xxxxx
* Smiles: canonical smiles
* InChIKey
* MolWeight: molecular weight (g/mol)
* NumRings: number of rings
* Batch: batch of generation (0 to 5)
* Predicted scores for:
    * QED (0:not drug like - 1: max similarity to known drugs)
    * SLogP: octanol-water partition coefficient as a measure of solubility.Lower values indicate higher hidrophilicity
    * SAScore (0-1): synthetic accessibility (rdKit)
    * RAScore (0-1): retrosynthetic accessibility from the [Reymond Group](https://github.com/reymond-group/RAscore)
    * Tanimoto similarity with existing series 4 molecules (0: no similarity - 1: repeated molecule)
    * Triazolopyrazine core maintenance (0: no core -  1: core)
    * Activity against P. Falciparum:
        * ActivityClfAutoML (0: inactive - 1: active): autoML random forest classifier from sklearn
        * ActivityRegAutoML (0: inactive - 1: active): autoML random forest regressor from sklearn
        * ActivityClfGraph (0: inactive - 1: active): Chemprop-based random forest classifier
        * ActivityRegGraph (0: inactive - 1: active): Chemprop-based random forest regressor
        * ActivityAvg (0: inactive - 1: active): average of the four activity prediction scores
* Clustering: the candidates have been clustered according to chemical similarity in:
    * Cluster100: 100 clusters
    * Cluster1000: 1000 clusters

In addition, we provide a [web-based interface](https://share.streamlit.io/ersilia-os/osm-series4-candidates/main/app.py) for the 1000 molecules (one from each cluster) with highest predicted activity (using average activity). If there is a particular cluster of interest please refer to the full excel file to select other molecules from the same cluster.

Finally, the 116728 molecules are represented in this [TreeMap](https://ersilia-os.github.io/osm-series4-candidates/tmap/map_100k/index.html).

## About us

[Ersilia Open Source Initiative](https://ersilia.io) is small non-profit aimed at disseminating AI/ML for drug discovery. We have a special interest in infectious and neglected diseases of low- and middle-income countries.
