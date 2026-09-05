# Auditing repetition steering and the confidence metric that validates it

Code and data for an anonymous NeurIPS 2026 AI4DD workshop submission.

This repository contains the notebooks that produced every number in the paper, the
saved artifacts they wrote, and the scripts that render the figures.

---

## What this is

A generate-and-screen pipeline for protein language models has two moving parts: an
intervention that shapes generation, and a cheap readout that decides what is worth
folding. We audit one published intervention (repetition steering) and one standard
readout (mean ESMFold pLDDT) with the same diagnostic — remove the thing the method is
supposed to depend on, and see whether anything changes.

We introduce no new method. The contribution is the audit.

---

## Layout

    notebooks/        experiments, run on Kaggle (T4). Executed, with output cells intact.
    analysis/         zero-GPU analysis scripts, saved CSV/JSON artifacts, figure renderers
    analysis/final/   the rendered figure PDFs used in the paper

Every notebook is numbered by the order it was run. Numbers in the paper trace to the
output cells of the notebooks listed below.

---

## Which notebook produced which result

### Steering arm

| Notebook | Result |
|---|---|
| `03-ai4dd-uccs-baseline-test` | the N=50 unsteered control arm (62.0%) |
| `34-ai4dd-residual-norm-audit` | per-model relative push; the 41,033x spread across conditions all labelled "1x" |
| `35-ai4dd-random-direction-control` | direction controls at matched norm; this run's own control is 52.0% |
| `36-ai4dd-foldability-direction-steering` | steering toward foldability, attempt 1 (co-linear, cosine +0.9819) |
| `37-ai4dd-dose-response-curve` | the eight-dose sweep and the fitted cliff; this run's own control is 64.0% |
| `39-ai4dd-orthogonalized-foldability-and-refold-fix` | attempt 2 (orthogonalized, cosine -0.0000), and the orthogonal-2x refold |
| `43-ai4dd-pep-alpha-mapping` | norm-ratio mapping of the published method onto our dose axis |
| `58-ai4dd-piggen-convention-and-vector-stability` | split-half reproducibility of the steering direction |

### Prediction arm

| Notebook | Result |
|---|---|
| `08-ai4dd-scaled-final-validation` | N=400 layer sweep and all matched baselines, both thresholds. Unsteered generation throughout |
| `59-ai4dd-riskcoverage-n1200` | headline AUC at N=1200 and the risk-coverage curve. Unsteered |
| `42-ai4dd-subclinical-detection-and-early-abort` | the sub-clinical null; early abort on the primary model |
| `44-ai4dd-early-abort-multimodel` | early abort across three generators |
| `61-ai4dd-sae-vs-raw-trained-properly` | sparse-vs-raw comparison at matched budget |

### Deployment

| Notebook | Result |
|---|---|
| `45-ai4dd-compute-cost-and-screening` | unbatched per-component timings and triage policies |
| `46-ai4dd-batched-throughput` | batched throughput, and the check that batching does not change pLDDT |

### External check against experimental measurements

| Notebook | Result |
|---|---|
| `50-rocklin-scramble-controlled-foldability-EXECUTED` | design-vs-scramble discrimination against measured protease stability |

---

## Reproducing the figures

    python analysis/fig2_steering.py     # steering: relative push, dose cliff, direction
    python analysis/fig3_gate.py         # repetition vs confidence; confidence vs stability
    python analysis/fig4_predictor.py    # risk-coverage; early abort across models
    python analysis/fig_appendix.py      # appendix figures

These are self-contained matplotlib transcriptions of values already in the saved
artifacts. They need no GPU.

The structure panel is rendered separately from `analysis/fig1/`.

---

## Notes

**One large artifact is omitted.** The per-sequence activation tensor from the
design-vs-scramble experiment (~82 MB) is not included. Re-running notebook 50
regenerates it. All derived results — predictions, per-pair outcomes, summary
statistics — are in `analysis/`.

**Controls are not interchangeable.** Three unsteered N=50 control arms appear across
these experiments, at 52.0%, 62.0% and 64.0%. They come from different runs and the
spread is consistent with binomial noise at that sample size. Every steering condition
is compared only with the control from its own run. The larger N=400 natural-collapse
estimate is 53.2%.

**Several AUC values are close in magnitude and are not comparable.** They come from
different pools and evaluation schemes: repeated stratified splits at N=1200 and N=400,
out-of-fold predictions for risk-coverage, a separate early-abort pool, and a smaller
compute pool. The paper states which is which at each use.

**Two claims were withdrawn during this work.** A sparse-autoencoder comparison was
retracted once we found the dictionary had been trained on too few sequences to support
it; an early-versus-late layer sensitivity claim was retracted once we found every model
had been over-pushed at its early layer, a dose confound rather than a depth effect.
Notebooks from the superseded analyses are retained for provenance and are not cited in
the paper.

**Hardware.** All GPU work ran on a single Tesla T4. Absolute timings are
hardware-specific; the ratios are the transferable quantities.