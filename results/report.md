# Brain aging vs Alzheimer's disease: summary

Dataset: GSE48350  |  Region: hippocampus

## Group sizes
- young: 10
- old: 25
- AD: 19

## Differential expression and overlap
- Genes tested in both: 21655
- Aging DEGs: 137
- AD DEGs: 28
- Shared: 0 (concordant 0, discordant 0)
- Aging-only: 137
- AD-only: 28

- Overlap hypergeometric p: 1.00e+00
- Jaccard index: 0.000
- Spearman r of logFC (all common genes): -0.389

## Top pathways
**shared**
- no enrichment results

**aging_only**
- Immune System R-HSA-168256 (Reactome_2022, up, adj p 3.1e-14)
- Staphylococcus aureus infection (KEGG_2021_Human, up, adj p 1.2e-09)
- Innate Immune System R-HSA-168249 (Reactome_2022, up, adj p 1.2e-09)
- Leishmaniasis (KEGG_2021_Human, up, adj p 1.7e-09)
- Tuberculosis (KEGG_2021_Human, up, adj p 2.0e-09)

**ad_only**
- Actin Filament-Based Transport (GO:0099515) (GO_Biological_Process_2023, all, adj p 4.2e-02)
- Vesicle Transport Along Actin Filament (GO:0030050) (GO_Biological_Process_2023, all, adj p 1.6e-01)
- Retrograde Protein Transport, ER To Cytosol (GO:0030970) (GO_Biological_Process_2023, all, adj p 1.6e-01)
- Nuclear Migration (GO:0007097) (GO_Biological_Process_2023, all, adj p 1.6e-01)
- Endoplasmic Reticulum To Cytosol Transport (GO:1903513) (GO_Biological_Process_2023, all, adj p 1.6e-01)

## Top hubs (co-expression degree)
- **shared**: none
- **aging_only**: IL13RA1 (78), TMEM176B (78), HLA-DMA (76), FXYD5 (75), MSN (75)
- **ad_only**: CTD-2587H24.10 (16), HOXA7 (11), CYP2A7P1 (10), GIPC3 (10), LA16c-380H5.4 (10)

## Reading the result
AD shows a distinct molecular signature beyond aging.

## Limitations
- Single brain region per run.
- Single dataset; results need replication in an independent cohort.
- "Early AD" depends on how the dataset labels its AD cases; many are late-stage post-mortem diagnoses.
- Post-mortem tissue: agonal state, RNA quality and post-mortem interval are not modeled.
- Cell-type composition changes (neuron loss, gliosis) are not modeled and can drive expression differences.
- Co-expression hubs are correlational, not causal.
- Batch is not modeled; in GSE48350 all AD samples were submitted separately from the controls, so AD vs old is confounded with batch.
