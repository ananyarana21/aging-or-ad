# What this project is about (simple version)

## The question

As people get older, their brains change. Alzheimer's disease (AD) also changes the brain.

So we ask: **is Alzheimer's just "normal aging, but faster"? Or is it something different?**

## Why it matters

- If AD is **just fast aging**, then anything that slows normal brain aging might also help with AD.
- If AD is **something different**, it needs its own treatments and its own tests.

## How we check it

Every cell turns genes "up" or "down". We measure how active about 21,000 genes are in brain tissue from people who have died. Then we make two comparisons:

1. **Aging:** old healthy brains vs young healthy brains. Which genes change as we get old?
2. **Alzheimer's:** brains with AD vs old healthy brains. Which genes change because of AD?

Then we put the two lists side by side:

- **Same genes change in both, in the same direction** → AD looks like fast aging.
- **Different genes change** → AD is its own thing.
- **Some of each** → it's mixed.

## The data

- **Main data:** GSE48350, a public dataset. We use the **hippocampus**, the memory part of the brain and one of the first areas hit by AD.
  - 10 young healthy people, 25 old healthy people, 19 people with AD
- **Check data:** GSE5281, another public dataset. It has no young people, so it can only test the AD side.

## What the program does, step by step

1. **Download the data** and sort people into young / old / AD.
2. **Find changed genes** for aging, and for AD.
3. **Compare the two lists**: how many genes are shared, and do they move the same way?
4. **Look up what those genes do**, for example immunity or energy.
5. **Find "hub" genes**: the genes most connected to others, the likely key players.
6. **Write a summary** (`results/report.md`).

## What we found (hippocampus)

- **Aging** changed 137 genes, mostly **immune / inflammation** genes going up.
- **AD** changed 28 genes.
- **Zero genes overlapped.**
- So, in this data, **AD looks different from normal aging**.

## Be careful

- All the AD samples were processed at a different time from the healthy ones. Some of the AD "differences" could come from that, not from the disease.
- The brain tissue came from people after death, which adds noise.
- Brains with AD lose nerve cells, and that alone changes gene activity.
- These results need to be checked in other datasets.

## How to run it

```bash
source .venv/bin/activate
python main.py --demo   # quick test with fake data
python main.py          # the real analysis
```

Then open `results/report.md`.
