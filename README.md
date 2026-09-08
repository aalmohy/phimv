# PHIMV

MATLAB implementation and reproducibility files for

**Awad H. Al-Mohy, _Computing Linear Combinations of $\varphi$-Function Actions for Exponential Integrators_.**

Preprint: [arXiv:2509.26475](https://arxiv.org/abs/2509.26475) To appear in IMA Journal of Numerical Analysis.

## Overview

`phimv` computes one or several linear combinations of matrix
$\varphi$-function actions,

$w_i = \varphi_0(t_iA)v_0 + \sum_{j=1}^{p}\alpha_i^j\,\varphi_j(t_iA)v_j,\qquad i=1,\ldots,r,$
where
$\varphi_0(z)=e^z,\qquad \varphi_j(z)=\frac{e^z-\sum_{k=0}^{j-1}z^k/k!}{z^j},\quad j\ge 1$.

The routine is designed for exponential-integrator computations and supports

- explicit dense or sparse matrices;
- matrix-free, block-capable linear operators;
- a single linear combination or several combinations evaluated simultaneously;
- independent stage parameters $t_i$ and polynomial weights $\alpha_i$;
- user-controlled accuracy;
- reuse of the scaling and shift parameters when the same operator is used repeatedly.

The method uses scaling and recovering with truncated Taylor series and a
spectral shift.  The shift and scaling parameters are selected by
`find_optimal_param`.

## Main files

```text
phimv/
|
|-- phimv.m
|-- find_optimal_param.m
|-- setup.m
|-- README.md
|
|-- experiments/
|   |-- experiment1.m
|   |-- experiment2.m
|   |-- experiment3.m
|   |-- experiment4.m
|   |
|   |-- common/
|   |   `-- phi_func_ex.m
|   |
|   `-- results/
|
`-- external/                 # created locally by setup.m
```

The directory `external/` contains third-party comparison codes downloaded by
`setup.m`; these packages are not part of PHIMV itself and should normally not
be committed to this repository.

## Requirements

The core routines

```text
phimv.m
find_optimal_param.m
```

are designed for **MATLAB R2022b** and require no additional MATLAB toolbox.

The numerical experiments reported in the paper were carried out in
**MATLAB R2022b** on a PC with an Intel Core i7-7700T CPU @ 2.90 GHz and
16 GB of RAM.

Later MATLAB releases may also run the code, but timings can differ between
MATLAB versions and hardware.

### Optional requirements for the experiments

The comparison experiments use external software installed by `setup.m`:

- [`phi_funm`](https://github.com/xiaobo-liu/phi_funm), by Al-Mohy and Liu;
- [`bamphi`](https://github.com/francozivcovich/bamphi), by Caliari, Cassini, and Zivcovich;
- [`kiops`](https://gitlab.com/stephane.gaudreault/kiops), by Gaudreault, Rainwater, and Tokman;
- [`phipm`](https://calgo.acm.org/919.zip), by Niesen and Wright, ACM Algorithm 919;
- [`Anymatrix`](https://github.com/north-numerical-computing/anymatrix), by Higham and Mikaitis;
- [`matrices-expm`](https://github.com/xiaobo-liu/matrices-expm), by Xiaobo Liu, installed by `setup.m` as the Anymatrix collection `expm`.

Some high-precision reference computations use the
[Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/).
Advanpix is commercial software and is **not** downloaded by `setup.m`.

Experiment 3 can use the Parallel Computing Toolbox to impose hard timeouts.
If it is unavailable, the script falls back to serial execution with a soft
timeout.

## Installation

Clone or download the repository and start MATLAB in the repository root.

For PHIMV itself, add the repository root to the MATLAB path:

```matlab
addpath(pwd)
```

To install the external packages required by the numerical experiments, run

```matlab
setup
```

The packages are downloaded into the local `external/` directory.  The setup
script adds only the intended package directories to the MATLAB path so that,
in particular, bundled copies of solvers do not shadow the standalone
versions used in the experiments.

To remove and redownload the external packages, use

```matlab
setup(true)
```

The setup script also installs `matrices-expm` inside Anymatrix as the
collection `expm` and runs

```matlab
anymatrix('scan')
```

to register the collection.

## Basic usage

The calling sequence is

```matlab
[W,st,mv] = phimv(t,alpha,A,v0,V,tol,s,shift);
```

The last four inputs are progressively optional:

```matlab
W = phimv(t,alpha,A,v0);
W = phimv(t,alpha,A,v0,V);
[W,st,mv] = phimv(t,alpha,A,v0,V,tol);
[W,st,mv] = phimv(t,alpha,A,v0,V,tol,s,shift);
```

Here

- `t` is a scalar or a vector of stage parameters;
- `alpha` is a scalar or has `numel(t)` entries;
- `A` is an `n`-by-`n` matrix or a block-capable function handle;
- `v0` is the vector multiplying $\varphi_0(tA)=e^{tA}$;
- `V = [v1,...,vp]` contains the vectors multiplying
  $\varphi_1,\ldots,\varphi_p$;
- `tol` is the requested tolerance;
- `s` and `shift` are optional reusable parameters returned by
  `find_optimal_param`;
- `st` is the effective scaling count;
- `mv` counts the number of columns passed to the operator, including
  parameter-selection work when the parameters are computed internally.

Either `v0` or `V` may be empty, but not both.

### Example 1: explicit matrix

```matlab
n = 100;
A = gallery('tridiag',n,-1,2,-1);

v0 = randn(n,1);
V  = randn(n,3);

w = phimv(1,1,A,v0,V);
```

This computes

$w=e^A v_0+\varphi_1(A)v_1+\varphi_2(A)v_2+\varphi_3(A)v_3.$


### Example 2: matrix-free operator

For a large sparse matrix, the operator can be supplied as a function handle:

```matlab
Afun = @(X) A*X;

w = phimv(1,1,Afun,v0,V);
```

The function handle must be **block capable**: when `X` has several columns,
`Afun(X)` must return the corresponding block \(AX\).

### Example 3: several combinations in one call

```matlab
t     = [1/3, 1/2, 5/6];
alpha = t;

W = phimv(t,alpha,Afun,v0,V);
```

Column `i` of `W` is

$W(:,i)=e^{t_iA}v_0+\sum_{j=1}^{p}t_i^j\varphi_j(t_iA)v_j.$


The parameters `t` and `alpha` are intentionally separate.  For example,

```matlab
W = phimv([1/3,5/6],[1,1],Afun,v0,V);
```

uses different arguments of the \(\varphi\)-functions while keeping the
external polynomial weights equal to one.

## Reusing the scaling and shift

When `phimv` is called repeatedly for the same operator, tolerance, and
working precision, the scaling and shift can be computed once:

```matlab
tol = 1e-7;
m   = 61;

[s,shift] = find_optimal_param(Afun,n,m,tol);
```

and then reused:

```matlab
W1 = phimv(t1,alpha1,Afun,v0,V,tol,s,shift);
W2 = phimv(t2,alpha2,Afun,v0,V,tol,s,shift);
```

This is useful in exponential integrators, where the same spatial operator is
used at many time steps.

Do not reuse `s` and `shift` for a different operator, tolerance, or working
precision.

## Reproducing the numerical experiments

First install the external dependencies:

```matlab
setup
```

Then run the experiments from the `experiments` directory.

### Experiment 1

```matlab
run experiments/experiment1.m
```

This tests PHIMV against `phi_funm` on the 108-matrix test suite used in the
paper and produces the accuracy figures.

The 200-digit reference calculations require the optional Advanpix toolbox.

### Experiment 2

```matlab
run experiments/experiment2.m
```

This compares `phimv`, `bamphi`, `kiops`, and `phipm` on a Chebyshev
spectral discretization of the one-dimensional Laplacian for several values
of $t$.

The high-precision reference calculation requires the optional Advanpix
toolbox.

### Experiment 3

```matlab
run experiments/experiment3.m
```

This compares `phimv`, `bamphi`, and `kiops` on deterministic large,
low-rank matrix-free operators whose spectral and nonnormal behavior is
controlled by a small core matrix.

The reference is evaluated from the small core using `phi_funm`, so the
large matrix is never formed.

A 10-minute timeout is used for difficult solver calls.  With the Parallel
Computing Toolbox the timeout is enforced by a cancellable worker; otherwise
the script uses a serial soft timeout.

### Experiment 4

```matlab
run experiments/experiment4.m
```

This reproduces the advection-diffusion-reaction test using the exponential
Runge--Kutta scheme EXPRK4S6.  Each time step uses four grouped
$\varphi$-combination calls for each of `phimv`, `bamphi`, and `kiops`.

The spatial operator is stored as a sparse Kronecker sum.  The reference
solution is computed with `ode15s` using

```matlab
RelTol = 1e-12
AbsTol = 1e-12
```

and an exact sparse Jacobian.  The reference solution is cached in the
results directory so that it need not be recomputed on subsequent runs.

## Output files

The experiment scripts save numerical data and generated figures in

```text
experiments/results/
```

Execution times depend on the MATLAB release, operating system, hardware, and
current machine load.  Accuracy results should be reproducible up to normal
floating-point variation, while timings should not be expected to match the
paper exactly on a different system.

## Notes on external software

The third-party packages downloaded by `setup.m` remain the property of
their respective authors and are subject to their own licenses.

`setup.m` downloads the packages from their upstream repositories rather than
redistributing them with PHIMV.

## Citation

If you use this software in research, please cite:

```bibtex
@misc{almohy2025phimv,
  author       = {Awad H. Al-Mohy},
  title        = {Computing Linear Combinations of
                  {$\varphi$}-Function Actions for Exponential Integrators},
  year         = {2025},
  eprint       = {2509.26475},
  archivePrefix= {arXiv},
  primaryClass = {math.NA},
  url          = {https://arxiv.org/abs/2509.26475},
 note          = {To appear in IMA Journal of Numerical Analysis}
}
```

## Author

**Awad H. Al-Mohy**  
Department of Mathematics  
King Khalid University  
Abha, Saudi Arabia

GitHub: [aalmohy](https://github.com/aalmohy)
