# osm-series4-candidates
New candidates for the Open Source Malaria series 4 compounds


## Project Overview
1. We have downloaded the experimental data for Series 4 compounds from the [Master List](https://github.com/OpenSourceMalaria/Series4/wiki/Sources-of-Data) in OSM GitHub (20/05/21)
2. The data has been used to train different Activity (IC50) predictors, which have subsequently been incorporated as scoring components of the [REINVENT 2.0](https://github.com/MolecularAI/Reinvent) Reinforcement Learning generative model (23-27/05/21)
3. Five different batches of new series 4 molecules have been generated as explained below (28-29/05/21)


## New candidate lists
### Batch00
With this batch, we have tried to narrow the chemical space to series 4 compounds, maintaining the triazolopyrazine core with NW and NE substituents. This has been achieved by forcing matching to the core substructure and similarity to series 4 molecules. In addition, core structures with substituents at undesired atoms have been restricted. QED Beauty and solubility (SLogP), as well as synthetic (SA) and retrosynthetic accessibility (RA) have been included as scoring components. Finally, activity scoring has been achieved with an autoML random forest classifier and regressor.
Batch00 rendered 1867 unique compounds (non-duplicates from original series 4)

### Batch01
Batch01 has been generated upon batch00 prior and including the same scoring components, on an effor to generate compounds with improved activity.
Batch01 rendered 1138 unique compounds (non duplicates from original series 4), of whom 532 are shared with batch00. Both Batch00 and Batch01 have produced similar structures, with very few modifications in th substituents. Interestingly, molecular weight (MW) of new molecules is lower (≈350g/mol) than those of active series 4 compounds (≈450g/mol)

### Batch02
With Batch02 we have tried to explore the chemical space, giving eliminating the similarity to series 4 score component as well as activity prediction. The goal is to increase molecular weight (adding a molecular weight component) in order to diversify the substituents. It has been generated upon batch01 prior.
Batch02 includes 46349 compounds with molecular weight around 400-500 and more diversity (Tanimoto Similarity in pairwise comparison is ≈ 0.65). After molecule generation, we have predicted the activity using the autoML RF classifier previously trained, with highest scores of 0.7

### Batch03
Batch03 is an enlargement of batch02 including more restrictions for desired physicochemical characteristics (QED_Score, SLogP, SA and RA).

### Batch04
Finally, with a diversified chemical space from series 4, we have focused our efforts on obtaining compounds with improved predicted activity. To this end, we have first trained a [Chemprop model](https://github.com/chemprop/chemprop) using

## Results

* Batch data in CSV:
  * [Batch 00-01](https://github.com/ersilia-os/osm-series4-candidates/blob/main/data/batch00_01/processed_results/batch00_01.csv)
  * [Batch 02](https://github.com/ersilia-os/osm-series4-candidates/blob/main/data/batch02/processed_results/batch02.csv)
  * [Batch 03](https://github.com/ersilia-os/osm-series4-candidates/blob/main/data/batch03/processed_results/batch03.csv)
  * [Batch 04]()
* Streamlit app visualization:  https://share.streamlit.io/ersilia-os/osm-series4-candidates/main/app.py

* Tree map:
## Code and linked repositories
