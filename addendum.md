# Addendum

## Handling of Low-Hit-Count Assays

In our paper, we chose to include assays with very few hits (positive examples) in our training data. 
This decision was made to capture a realistic setting where some assays of interest may have limited positive examples due to the expense and difficulty of wet-lab experiments. 
Including these assays allowed us to evaluate the performance of multi-task learning approaches in this challenging scenario.

One of the main concerns raised about our analysis was how we handled these low-hit-count assays in our cross-validation procedure. 
Specifically, the concern was that for assays with fewer than 5 hits, at least one of the 5 cross-validation folds would necessarily contain no positive examples, leading to an undefined AUC for that fold. 
When aggregating AUCs across folds (by taking the median or mean), these undefined values (NaNs) would be ignored.

However this does not inflate performance estimates. Here's a detailed explanation of why our approach is valid:

Consider an example assay (175_561) which has 1 active compound out of 31 total compounds, and for which we report an AUC of 1.0 using the Cell Painting data modality. 
In a 5-fold cross-validation setup, the distribution of active and inactive compounds might look like this:

| Fold   | Inactives | Actives |
|--------|-----------|---------|
| Fold 1 | 8         | 1       |
| Fold 2 | 3         | 0       |
| Fold 3 | 11        | 0       |
| Fold 4 | 3         | 0       |
| Fold 5 | 5         | 0       |
| Total  | 30        | 1       |

In this case, AUC can only be evaluated for Fold 1, as it is the only fold containing an active compound. 
The other folds, having no actives, will produce undefined AUC values (NaNs) which are ignored when aggregating.

So, for this assay, the reported "cross-validation" AUC of 1.0 is in fact based on a single train-test split, where:

- The model is trained on Folds 2-5, which contain no active compounds. This is valid in our multi-task learning setup, as the model can still learn useful information from the other assays in these folds.
- The model is then evaluated on Fold 1, which contains 8 inactives and 1 active. With this setup, there are only 9 possible AUC values: 1.0, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0.0. Our reported AUC of 1.0 indicates that the single active compound was ranked above all 8 inactives.

This means that:

1. For the 36 assays in our dataset that have only 1 active compound, the reported AUC is based on a single train-test split rather than an average over 5 folds. More generally, this applies to any assay where only one of the cross-validation folds contains active compounds.
2. For assays where at least one fold (but not all folds) has no active compounds, the reported AUC is a cross-validation AUC in the sense that it averages over multiple train-test splits, but it is not a "5-fold" cross-validation AUC because some of the folds are ignored in the aggregation.

This approach does not inflate the reported performance estimates. 
The model still has to rank the single active compound above the inactives in the test set based on what it has learned from the other assays, which is a non-trivial task. 
The fact that we observe high AUCs in these cases suggests that the model is indeed learning useful cross-assay information. 
While this is a valid evaluation given the limited data, it's important to note that the uncertainty around these performance estimates is higher than it would be for a typical 5-fold cross-validation.

## Morphology Performance and Cell Health

While morphology-based models performed best on some assays with viability-related readouts, predicting hits for these assays likely involves more than identifying generally toxic compounds. For example, some assays sought compounds selectively affecting specific cell types or conditions, not just broadly cytotoxic ones. The X reactivation assay, despite using cell count, probably required specific effects beyond toxicity to enable cell survival (perhaps contingent on X reactivation). In non-mammalian assays (yeast, bacteria), extrapolating from mammalian toxicity is not trivial.

Thus, while cell health indicators likely contribute to the predictive morphological changes, they're probably not the whole story for these assays. A follow-up analysis comparing morphology-based models to models based on cell count alone could help tease apart the contributions of general toxicity versus more specific morphological effects, but with the caveat that "toxicity" itself may have very different implications across the diverse assay types in this dataset. 

<details>
<summary>Code</summary>
  
```r
# https://raw.githubusercontent.com/CaicedoLab/2023_Moshkov_NatComm/20fd97178796148c2b464ddc2c025a91654a6433/predictions/scaffold_median_AUC.csv
scaffold_median_AUC <- 
  read_csv("~/Documents/GitHub/2023_Moshkov_NatComm/predictions/scaffold_median_AUC.csv")

assay_metadata_expanded <-
  read_csv("~/Documents/GitHub/2023_Moshkov_NatComm/assay_data/assay_metadata_expanded.csv")

scaffold_median_AUC <- 
  scaffold_median_AUC %>% 
  filter(descriptor %in% c("mobc_es_op", "cp_es_op", "ge_es_op")) %>% 
  select(assay_id, auc_90, descriptor) %>% 
  pivot_wider(id_cols = assay_id, names_from = "descriptor", values_from = "auc_90") %>%
  filter(mobc_es_op & !cp_es_op & !ge_es_op) %>% 
  inner_join(assay_metadata_expanded, by = c("assay_id" = "PUMA_ASSAY_ID")) %>%
  arrange(ASSAY_DB_ID)

scaffold_median_AUC %>% 
  select(assay_id, ASSAY_NAME, ASSAY_DESC, OBS_NAME, ASSAY_DB_ID) %>% 
  write_csv("~/Desktop/MO_best.csv")
```
</details>

<details>
<summary>Data</summary>

`MO_best.csv`

```csv
assay_id,ASSAY_NAME,ASSAY_DESC,OBS_NAME,ASSAY_DB_ID
178_591,MITF,MITF.suppression,observation_1,CBIP:2007-01
184_606,STK33,NOMO.CTG.viability,observation_0,CBIP:2021-01
186_669,MitochondrialAbundance,Mitochondrial.Abundance,MitoIntensity/CellArea,CBIP:2026-01
187_618,STK33invitro,ADP-GLo,adp-glo,CBIP:2029-01
275_753,MLPCN X-reactivation,High content Imaging Assay,totalcells_sum,CBIP:7015-01
250_725,GPR85,Primary Screen,Fold Change over DMSO,CBIP:7048-01
280_759,MLPCN C. diff tox,C. difficile toxin inhibitors,amplexred,CBIP:7074-01
42_104,CellMetabolicProfiling.Absorb.MTTStain,MetabolismCellProfiling,TubesC.CCCP.48h,ScreenID:1020
8_13,AspulvinoneRegulation.Absorb.600,AspulvinoneUpregulation,ATCC10020.CsA-,ScreenID:1022
175_559,YeastGrowthInhibition.Absorb.600,YeastGrowthModifierScreen2005,BY4716._Basal.48h,ScreenID:1024
175_561,YeastGrowthInhibition.Absorb.600,YeastGrowthModifierScreen2005,BY4716.Sert14.48h,ScreenID:1024
44_106,CellMetabolicProfiling.FluorDye.AmplexRed,GlucoseMetabolismAssay,L6,ScreenID:1030
11_22,BacterialViability.Lumin.BTG,BacterialViabilityProfiling,SaureusR,ScreenID:1064
9_14,BacterialViability.Absorb.600,MycobacteriumTuberculosisViabilityInhibition,EFAE,ScreenID:1106
85_247,GSISyntheticLethality.FluorDye.Resazurin,GSISyntheticLethal,KoptK1.GSI,ScreenID:1111
5_5,AnthraxPathogenesis.Lumin.CTG,AnthraxLethalFactorInhibtion,TIB-67.LFPA,ScreenID:1135
17_28,BetaCateninSignaling.FluorDye.PosW2,BetaCatenin,HTB-96.bCAT-EGFP,ScreenID:1152
18_29,BetaCateninSignaling.FluorDye.TotalNuclei,BetaCatenin,HTB-96.bCAT-EGFP,ScreenID:1152
37_77,CardiovascularProfiling.FluorDye.JC1Mito,EndothelialCellProfiling2,UV0482._Basal.24h,ScreenID:914
```
  
</details>


<details>
  <summary>Table</summary>

In the table below `Cell Viability Focus` and `Description` are generated by GPT4o; they may contain errors and falsehoods
  |ASSAY_NAME                               |ASSAY_DESC                                     |OBS_NAME              |Our guess about viability                                                                |GPT4o Description, based on ASSAY_NAME, ASSAY_DESC, OBS_NAME                     |
|-----------------------------------------|-----------------------------------------------|----------------------|-----------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
|STK33                                    |NOMO.CTG.viability                             |observation_0         |No                                                                                       |Evaluates cell viability in the context of STK33 activity.                       |
|STK33invitro                             |ADP-GLo                                        |adp-glo               |No                                                                                       |Measures ADP levels to assess kinase activity of STK33.                          |
|MLPCN C. diff tox                        |C. difficile toxin inhibitors                  |amplexred             |No                                                                                       |Screens for inhibitors of C. difficile toxin affecting cell survival.            |
|CardiovascularProfiling.FluorDye.JC1Mito |Endothelial Cell Profiling                     |UV0482._Basal.24h     |No                                                                                       |Profiles mitochondrial function in endothelial cells using JC-1 dye.             |
|AnthraxPathogenesis.Lumin.CTG            |Anthrax Lethal Factor Inhibition               |TIB-67.LFPA           |No                                                                                       |Measures inhibition of anthrax lethal factor using a luminescent assay.          |
|MITF                                     |MITF suppression                               |observation_1         |No                                                                                       |Measures the suppression of MITF activity relevant to pigmentation or melanoma.  |
|MitochondrialAbundance                   |Mitochondrial Abundance                        |MitoIntensity/CellArea|No                                                                                       |Quantifies mitochondrial abundance within cells.                                 |
|MLPCN X-reactivation                     |High content Imaging Assay                     |totalcells_sum        |No                                                                                       |Uses high-content imaging to measure cell number in pathway reactivation studies.|
|AspulvinoneRegulation.Absorb.600         |Aspulvinone Upregulation                       |ATCC10020.CsA-        |No                                                                                       |Measures regulation of Aspulvinone                                               |
|BetaCateninSignaling.FluorDye.PosW2      |Beta Catenin                                   |HTB-96.bCAT-EGFP      |No                                                                                       |Measures beta-catenin signaling pathway activity using fluorescence.             |
|BetaCateninSignaling.FluorDye.TotalNuclei|Beta Catenin                                   |HTB-96.bCAT-EGFP      |No                                                                                       |Measures total nuclei count in beta-catenin signaling studies.                   |
|GSISyntheticLethality.FluorDye.Resazurin |GSISynthetic Lethal                            |KoptK1.GSI            |No, seems to be measuring differential viability in the presence of a particular compound|Measures synthetic lethality in GSI treatment using resazurin fluorescence       |
|GPR85                                    |Primary Screen                                 |Fold Change over DMSO |Not enough info                                                                          |Screens GPR85 activity comparing changes to a DMSO control.                      |
|CellMetabolicProfiling.Absorb.MTTStain   |Metabolism Cell Profiling                      |TubesC.CCCP.48h       |Probably                                                                                 |Profiles cell metabolism using MTT staining.                                     |
|CellMetabolicProfiling.FluorDye.AmplexRed|Glucose Metabolism Assay                       |L6                    |Probably                                                                                 |Profiles cell metabolism related to glucose using Amplex Red fluorescence.       |
|BacterialViability.Lumin.BTG             |Bacterial Viability Profiling                  |SaureusR              |Yes, but bacteria only                                                                   |Measures bacterial viability using a luminescent assay for Staphylococcus aureus.|
|BacterialViability.Absorb.600            |Mycobacterium Tuberculosis Viability Inhibition|EFAE                  |Yes, but bacteria only                                                                   |Measures viability inhibition of Mycobacterium tuberculosis.                     |
|YeastGrowthInhibition.Absorb.600         |Yeast Growth Modifier Screen                   |BY4716._Basal.48h     |Yes, but yeast only                                                                      |Screens for yeast growth inhibition using absorbance at 600 nm.                  |
|YeastGrowthInhibition.Absorb.600         |Yeast Growth Modifier Screen                   |BY4716.Sert14.48h     |Yes, but yeast only                                                                      |Screens for yeast growth inhibition under different conditions.                  |


</details>

## Leveraging Cross-Assay Associations

One key insight from our work is that multi-task learning allows models to make reasonable predictions even for assays where no positive examples were seen during training. 
This is possible because the model learns associations between related assays, and can draw on this learned structure to infer the likely activity of compounds in the held-out assay.

Intuitively, when an assay has no positive examples in the training set, the model can still learn to push the predicted probabilities of the (known) negatives to be very low. 
Then, when an unmarked positive example is encountered in the test set, its predicted probability may not be high in an absolute sense, but will likely be higher than the negative examples, allowing it to rank above them. 
This explains how the model can achieve perfect AUC even when no positive examples were available during training.

## Cross-Validation Procedure

To make maximal use of our limited data while avoiding overfitting, we used a nested cross-validation scheme for hyperparameter tuning and model evaluation:

```py
# Main Experiment
all_data = [fold1, fold2, fold3, fold4, fold5]

for fold in range(5):
  train_data = all_data.copy().drop(fold)
  hold_out = all_data[fold]

  # Hyper-param optimization:
  train, val, test = split(train_data, 80, 10, 10)
  candidate_params = BayesianOptimization(train, val)

  for p in candidate_params:
    model = train_model(train+val, p)
    acc = eval_model(test)
    results.add(acc, p)

  best_params = best(results)

  # Train fold model
  model = train_model(train_data, best_params)
  fold_results = eval_model(hold_out)
  main_results.add(fold_results, best_params)

# Final results
median(main_results)
```

In each outer fold, the held-out 20% is used only for final evaluation. 
The remaining 80% is further split to allow selection of hyperparameters without overfitting to the test set. 
While computationally expensive, this procedure makes efficient use of the limited data while avoiding bias in model selection.

## Evaluation Choices and Alternatives

In our paper, we made specific choices about how to evaluate the performance of our models. 
The starting point of our evaluation was the output of our cross-validation procedure: a 270x16170 matrix M' of predicted probabilities, which is the completion of the sparse 270x16170 binary matrix M of assay readouts. 
We also have a 16170-dimensional vector v indexing the cross-validation splits.

Our chosen evaluation strategy was to:

1. Compute AUCs per-assay, rather than globally: We felt that per-assay metrics better captured the key question of which assays could be reliably predicted by each data modality.
2. Calculate AUROC and AUPRC per split and per assay
3. Aggregate AUCs across folds by taking the median: We chose the median as it is more robust to outliers than the mean. For assays with positive examples in all folds, this aggregates across 5 values; for assays with fewer folds containing hits, fewer values are aggregated.
4. Summarize these median AUROCs and AUPRCs in various ways, such as the Venn diagram in Figure 2, using an AUC threshold of 0.9 to define "hits" for purposes of model comparison. This was an arbitrary threshold chosen to represent strong predictive performance. We believe most key conclusions would hold for other reasonable thresholds, but specific "hit counts" and Venn diagrams would change.

We believe these choices provide a fair and informative assessment of the models' performance, and that our main claims around the relative strengths of the different data modalities are robust to reasonable variations in these choices. 
However, we acknowledge that different choices could be made, and that these choices could lead to different conclusions. 
Some alternative evaluation strategies include:

1. Dropping assays with no positives in the training set, either at the start of evaluation or during evaluation
2. Grouping assays by category and reporting metrics per category
3. Drawing Venn diagrams based on different AUROC thresholds
4. Comparing the performance of our models to a baseline model based on cell count alone
5. Computing a single AUROC/AUPRC across all splits for each assay (ignoring v)
6. Including confidence intervals for AUROC/AUPRC and factoring these into the summary
7. Evaluating performance on the entire matrix M', not just the ~13% of points with readouts in M

All of these alternatives are valid and could provide interesting additional perspectives on the models' performance. 
We encourage interested readers to explore these and other variations, as our full code and data are available to facilitate such analyses.

While the specific numerical results and comparisons may change under different evaluation schemes, we believe our main qualitative conclusions are likely to hold: namely, that the different data modalities offer complementary information useful for predicting a variety of cellular assays. 
The robustness of these high-level conclusions is supported by the consistency of the patterns we observed across a range of evaluation choices. 
However, we acknowledge that the inherent sparsity of the data means that conclusions about specific assays or modalities may be more sensitive to the details of the evaluation procedure.

## Summary

We hope these clarifications are helpful to readers looking to better understand the details behind our paper's key conclusions. 
We believe that the PUMA dataset provides a valuable resource for further methodological explorations in multi-task learning and multimodal predictive modeling, and look forward to seeing future work that builds on and expands the analysis we have presented.
