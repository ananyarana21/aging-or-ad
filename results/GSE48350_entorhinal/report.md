# Brain aging vs Alzheimer's disease: summary

Dataset: GSE48350  |  Region: entorhinal

## Group sizes
- young: 11
- old: 18
- AD: 15

## Differential expression and overlap
- Genes tested in both: 21655
- Aging DEGs: 0
- AD DEGs: 15
- Shared: 0 (concordant 0, discordant 0)
- Aging-only: 0
- AD-only: 15

- Overlap hypergeometric p: 1.00e+00
- Jaccard index: 0.000
- Spearman r of logFC (all common genes): 0.063

## Top pathways
**shared**
- no enrichment results

**aging_only**
- no enrichment results

**ad_only**
- Ammonium Homeostasis (GO:0097272) (GO_Biological_Process_2023, all, adj p 4.5e-02)
- Cell Projection Assembly (GO:0030031) (GO_Biological_Process_2023, all, adj p 4.5e-02)
- Phospholipid Homeostasis (GO:0055091) (GO_Biological_Process_2023, all, adj p 4.5e-02)
- Ammonium Transmembrane Transport (GO:0072488) (GO_Biological_Process_2023, all, adj p 4.5e-02)
- Sequestering Of Extracellular Ligand From Receptor (GO:0035581) (GO_Biological_Process_2023, all, adj p 4.5e-02)

## Top hubs (co-expression degree)
- **shared**: none
- **aging_only**: none
- **ad_only**: CYP2A7P1 (11), LOC646588 (11), SLC25A46 (11), ANKIB1 (10), ZNF621 (10)

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
