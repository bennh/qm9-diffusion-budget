# An Empirical Study of CDGS-Based Molecular Generation on QM9

This repository contains the code and experiment assets for a course project on **molecular graph generation with CDGS** on the **QM9** dataset.  
Using a fixed pretrained CDGS checkpoint, we study how different **diffusion sampling step budgets** affect generation quality, diversity, reliability, and runtime, and we further evaluate a lightweight **chemistry-aware post-processing** strategy.

---

## Overview

Diffusion-based molecular graph generators often face a trade-off between **sampling efficiency** and **sample quality**.  
In this project, we keep the **model, dataset, and checkpoint fixed**, and vary only the **sampling budget** to conduct a controlled empirical study.

We compare:

- different sampling budgets: `T ∈ {25, 50, 100, 200}`
- two evaluation settings:
  - `none`: raw generated outputs
  - `filterplus`: raw outputs followed by lightweight chemistry-aware filtering

The goal is to understand how fewer or more sampling steps influence:

- molecular validity and reliability
- distributional quality
- diversity and novelty
- sampling runtime

---

## Project Setting

### Model and Dataset

- **Model:** CDGS
- **Dataset:** QM9
- **Checkpoint:** pretrained `checkpoint_200`

### Sampling Budgets

- `T = 25`
- `T = 50`
- `T = 100`
- `T = 200`

### Main Evaluation Configuration

- `num_samples = 5000`
- `batch_size = 16`
- `save_graph = True`

---

## Evaluation Modes

### 1. `none`

Generated molecules are evaluated directly without additional filtering.

### 2. `filterplus`

Generated molecules are post-processed with lightweight chemistry-aware filtering before evaluation.

The `filterplus` pipeline includes:

- RDKit sanitization
- atom-type check
- connectivity check
- size constraint
- simple valence check

---

## Evaluation Metrics

Our analysis focuses on the following metrics:

- **Validity**
- **FCD/Test**
- **Novelty**
- **Unique ratio**
- **Sampling runtime**

Depending on the experiment setting, we also analyze retention statistics and failure cases after filtering.

---

## Repository Structure

```text
├── configs/                 # experiment configuration files
├── evaluation/              # evaluation utilities and metric code
├── models/                  # model components
├── scripts/                 # experiment and post-processing scripts
├── figures/                 # generated figures for the report
├── main.py                  # main entry point
├── run_lib.py               # training / evaluation pipeline
├── sampling.py              # sampling procedures
├── dpm_solvers.py           # diffusion solver implementations
├── datasets.py              # dataset utilities
├── losses.py                # loss definitions
├── sde_lib.py               # SDE utilities
├── utils.py                 # helper functions
├── visualize.py             # visualization utilities
├── results_summary.ipynb    # result aggregation and analysis
└── README.md
