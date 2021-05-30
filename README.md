# osm-series4-candidates
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
Batch04 includes 8833 molecules and Batch05 includes

## Results
* Batch data in CSV:
  * [Batch 00](https://github.com/ersilia-os/osm-series4-candidates/blob/main/data/ReinventResults/ProcessedResults/processed0.csv)
  * [Batch 01](https://github.com/ersilia-os/osm-series4-candidates/blob/main/data/ReinventResults/ProcessedResults/processed1.csv)
  * [Batch 02](https://github.com/ersilia-os/osm-series4-candidates/blob/main/data/ReinventResults/ProcessedResults/processed2.csv)
  * [Batch 03](https://github.com/ersilia-os/osm-series4-candidates/blob/main/data/ReinventResults/ProcessedResults/processed3.csv)
  * [Batch 04](https://github.com/ersilia-os/osm-series4-candidates/blob/main/data/ReinventResults/ProcessedResults/processed4.csv)
* Streamlit app visualization:  https://share.streamlit.io/ersilia-os/osm-series4-candidates/main/app.py

* Tree map:
