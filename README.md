# brana-tov

Modified Tolman–Oppenheimer–Volkoff solver and Bayesian inference code for compact stars
in an effective braneworld model, as used in

> S. I. dos Santos Júnior, R. C. R. de Lima and P. H. R. S. Moraes,
> *Equation-of-state dependence of compact stars in a large-tension braneworld:
> a multi-messenger Bayesian analysis*.

## What this does

Integrates the modified TOV system with a local quadratic brane correction and a
phenomenological closure for the non-local Weyl sector, parametrised by the brane tension
`λ`, the Weyl coupling `α_U` and the Weyl equation-of-state parameter `w_U`. Stellar
sequences are compared against NICER and gravitational-wave mass–radius posteriors through
an empirical kernel likelihood, sampled with an affine-invariant ensemble MCMC.

## Build

Requires a C++17 compiler with OpenMP.

```
make -j
```

This produces `build/brana_tov`. A CMake build is also provided
(`cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j`).

A quick check: `build/brana_tov sequence-custom inf 0.0 0.0 out.dat data/SLy.txt 20 0.03 0`
integrates the general-relativistic SLy sequence and should give
`M_max = 2.048 M_sun` (the exact value depends mildly on the
central-density grid resolution).

## Use

```
# a single stellar sequence at fixed parameters
build/brana_tov sequence-custom <lambda> <alpha_U> <w_U> out.dat <eos.txt> <n_rho> <dr> <closure>

# a production chain
build/brana_tov mcmc <steps> <walkers> <ensembles> <threads> <profile> <eos.txt> <outdir>
```

The profile `perturbative-gr` reproduces the configuration reported in the paper: kernel
bandwidths 1.2 km (GW) and 0.85 km (NICER), dataset weights 0.30 / 1.50 / 1.25, coverage
penalty `κ = 8`, and the priors given in Sec. 4.4 of the paper. Every run writes an
`effective_config.dat` listing the fully resolved configuration and a
`run_provenance.dat` recording the commit, compiler and build flags.

## Reproducing the paper

`scripts/rev2_chains.sh` and `scripts/rev2_chains_b.sh` run the seven production chains.
`scripts/rev2_summary.py` reduces them to the summary table, and
`scripts/build_eos_paper_assets.py` and `scripts/rev2_alpha_vs_R14.py` regenerate the
tables and figures. Python dependencies are in `requirements.txt`.

## Data

`data/SLy.txt` is included as a worked example. The other six equation-of-state tables come
from CompOSE (QHC19-A from Zenodo record 14808639) and are redistributed, together with the
posterior chains, in the Zenodo deposit accompanying the paper.

The observational posterior samples are **not** redistributed here. Obtain them from the
original analyses: GW170817 from the LIGO/Virgo data release, PSR J0740+6620 from
Dittmann et al. (2024), and PSR J1231−1411 from Salmi et al. (2024).

## Licence

MIT. See `LICENSE`.
