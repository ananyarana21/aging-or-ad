# aging-or-ad

Is Alzheimer's disease (AD) accelerated normal brain aging, or a distinct process? This project compares gene expression changes in brain tissue across aging and AD using public microarray datasets.

## Approach

1. Aging: young healthy vs old healthy samples.
2. AD: old healthy vs AD samples.
3. Compare the two sets of differentially expressed genes (overlap and direction of change).
4. Run pathway enrichment on shared, aging-only and AD-only gene sets.
5. Identify hub genes using co-expression networks (optionally STRING).
6. Write a summary report.

## Data

- GSE48350: main dataset (young, old and AD samples). Default region is the hippocampus.
- GSE5281: validation dataset for the AD comparison only (no young group).

Datasets are downloaded automatically with GEOparse and cached in `data/`.

## Setup

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Usage

```bash
python main.py --demo                               # run on synthetic data
python main.py                                      # default dataset and region
python main.py --gse GSE48350 --region entorhinal   # choose dataset and region
python main.py --string                             # add STRING network degree to hub genes
```

Defaults are set in `config.py`.

## Output

Results are written to `results/<run>/` and figures to `figures/<run>/`. The main summary is `report.md` in the results folder. Intermediate results are cached, so reruns skip finished steps.

## Project layout

- `main.py`: command line entry point
- `config.py`: dataset, thresholds and gene set settings
- `src/data.py`: data loading and sample grouping
- `src/de.py`: differential expression and volcano plots
- `src/compare.py`: aging vs AD comparison
- `src/enrich.py`: pathway enrichment
- `src/network.py`: co-expression and hub genes
- `src/report.py`: report generation

See `OVERVIEW.md` for a plain-language description of the question and findings.

## Caveats

AD and control samples may differ in processing batch, tissue is post-mortem, and neuron loss in AD affects expression. Findings should be validated on other datasets.
