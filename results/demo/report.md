# Brain aging vs Alzheimer's disease: summary

Dataset: synthetic demo  |  Region: hippocampus

## Group sizes
- young: 20
- old: 20
- AD: 20

## Differential expression and overlap
- Genes tested in both: 2000
- Aging DEGs: 303
- AD DEGs: 351
- Shared: 150 (concordant 120, discordant 30)
- Aging-only: 153
- AD-only: 201

- Overlap hypergeometric p: 5.08e-46
- Jaccard index: 0.298
- Spearman r of logFC (all common genes): -0.082

## Top pathways
**shared**
- no enrichment results

**aging_only**
- no enrichment results

**ad_only**
- no enrichment results

## Top hubs (co-expression degree)
- **shared**: GENE0000 (119), GENE0001 (119), GENE0002 (119), GENE0003 (119), GENE0004 (119)
- **aging_only**: GENE0165 (124), GENE0202 (124), GENE0211 (122), GENE0271 (115), GENE0227 (113)
- **ad_only**: GENE0341 (198), GENE0322 (197), GENE0438 (197), GENE0478 (197), GENE0403 (196)

## Reading the result
Mixed: partial overlap with substantial AD-specific changes.

## Limitations
- Single brain region per run.
- Single dataset; results need replication in an independent cohort.
- "Early AD" depends on how the dataset labels its AD cases; many are late-stage post-mortem diagnoses.
- Post-mortem tissue: agonal state, RNA quality and post-mortem interval are not modeled.
- Cell-type composition changes (neuron loss, gliosis) are not modeled and can drive expression differences.
- Co-expression hubs are correlational, not causal.
- Batch is not modeled; in GSE48350 all AD samples were submitted separately from the controls, so AD vs old is confounded with batch.
