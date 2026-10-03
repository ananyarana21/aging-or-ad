# Brain aging vs Alzheimer's disease: summary

Dataset: GSE5281  |  Region: hippocampus

## Group sizes
- young: 0
- old: 13
- AD: 10

## Differential expression and overlap
- AD DEGs: 3242 of 21655 genes tested

## Top pathways
**ad**
- Metabolism Of RNA R-HSA-8953854 (Reactome_2022, down, adj p 3.6e-30)
- Amyotrophic lateral sclerosis (KEGG_2021_Human, down, adj p 3.6e-24)
- Parkinson disease (KEGG_2021_Human, down, adj p 2.1e-23)
- HIV Infection R-HSA-162906 (Reactome_2022, down, adj p 2.7e-23)
- Cellular Responses To Stimuli R-HSA-8953897 (Reactome_2022, down, adj p 5.0e-23)

## Top hubs (co-expression degree)
- **ad**: PARK7 (1856), RTN3 (1830), ZBTB7A (1808), APLP2 (1784), ELAVL3 (1784)

## Reading the result
Aging comparison not available (no young controls in this dataset); only AD vs old controls is reported.

## Limitations
- Single brain region per run.
- Single dataset; results need replication in an independent cohort.
- "Early AD" depends on how the dataset labels its AD cases; many are late-stage post-mortem diagnoses.
- Post-mortem tissue: agonal state, RNA quality and post-mortem interval are not modeled.
- Cell-type composition changes (neuron loss, gliosis) are not modeled and can drive expression differences.
- Co-expression hubs are correlational, not causal.
- Batch is not modeled; in GSE48350 all AD samples were submitted separately from the controls, so AD vs old is confounded with batch.
