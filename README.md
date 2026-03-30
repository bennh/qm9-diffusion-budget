# An Empirical Study of CDGS-Based Molecular Generation on QM9

A course project repository for studying **CDGS-based molecular graph generation** on **QM9** under different diffusion sampling step budgets, with additional **lightweight chemistry-aware post-processing**.

## Highlights

- Fixed pretrained **CDGS** checkpoint on **QM9**
- Controlled comparison across `T ∈ {25, 50, 100, 200}`
- Two evaluation settings:
  - **none**: raw generated outputs
  - **filterplus**: raw outputs + lightweight chemistry-aware filtering
- Main metrics:
  - `FCD/Test`
  - `Novelty`
  - `Unique ratio`
  - `Sampling runtime`
- Final report built from saved `.pkl`, `.log`, and summary files

---

## Project Goal

This project investigates how the number of diffusion sampling steps affects molecular generation behavior in a fixed CDGS setup on QM9.

We focus on:

- how generation quality changes with `T`
- how runtime changes with `T`
- how novelty and uniqueness behave across different step budgets
- whether simple post-processing can improve final output quality

---

## Experimental Settings

### Model and Dataset

- **Model**: CDGS
- **Dataset**: QM9
- **Checkpoint**: `ckpt = 200`

### Main Sampling Budgets

- `T = 25`
- `T = 50`
- `T = 100`
- `T = 200`

### Main Evaluation Configuration

- `num_samples = 5000`
- `batch_size = 16`
- `save_graph = True`

---

## Experiment Variants

### 1. none

Raw generated molecules are evaluated directly.

### 2. filterplus

Raw generated molecules are passed through lightweight chemistry-aware filtering before evaluation.

The final filter includes:

- sanitization
- atom-set check
- connectivity check
- size limit check
- valence check

---

## Repository Layout

```text
.
├── configs/
│   └── vp_qm9_cdgs.py
├── scripts/
│   └── filter_only_eval.py
├── exp/
│   └── vpsde_qm9_cdgs/
│       ├── checkpoints/
│       ├── eval_t25_none_5k/
│       ├── eval_t50_none_5k/
│       ├── eval_t100_none_5k/
│       ├── eval_t200_none_5k/
│       ├── eval_t25_filter_5k/
│       ├── eval_t50_filter_5k/
│       ├── eval_t100_filter_5k/
│       └── eval_t200_filter_5k/
├── eval_t25_none_5k.log
├── eval_t50_none_5k.log
├── eval_t100_none_5k.log
├── eval_t200_none_5k.log
├── main.py
├── run_lib.py
└── README.md
