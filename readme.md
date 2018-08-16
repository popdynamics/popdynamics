



# Popdynamics

Popdynamics is a library for building epidemiological models, or for that matter, any kind of compartmental model. It provides an easy way to switch between different integration methods

- direct integration of ODE
- integration of ODE using the Scipy library
- stochastic discrete-time
- stochastic continuous-time


## Installation

Popdynamics has been written to be compatible with Python 2 and Python 3. There are two main approaches to install the dependencies.

_1) Anaconda Python_

Install the [Anaconda Python](https://www.anaconda.com/download/#macos) environment. This package already contains all the dependencies required.

_2) Local Python_

If you have your own Python environment, you will have to install the dependencies. In the main Popdynamics directory, to install the dependencies:

```bash
> python setup.py develop
```

If you want to auto-generate flow-charts, you will also have to download [Graphviz](https://graphviz.gitlab.io/download/), which is a binary to generate well-formated flow-charts. Popydnamics will skip flow-charts if Graphviz is not found.

Requirements: _python 2/3; numpy, scipy, matplotlib, graphviz (python), graphviz (binary)_


## Quick start

Once the modules are installed, go the `examples` directory, and run:

```bash
> python sir_model.py
```

Then try the other examples.

- `sir_model.py` - a simple model showing a basic SIR model
- `seir_model.py` - basic SEIR model with a demography variation
- `stochastic_strains_model.py` - an SIR model using stochastic models to average over many runs
- `tb_model.py` - a basic TB model that represents the approach used by [AuTuMN](http://www.tb-modelling.com/home/index.php)
- `rmit_tb_model.py` - a reproduction of a TB model in the literature



## Why?

Epidemiological models are typically modeled as compartmental models where populations evolve through linear ordinary differential equations (ODE), or through stochastic models with probability transition matrices.

Approaches to modelling involve writing out the complete ODE for each compartment and typing these equations in matrix based programming enviromnents such as Fortran, Matlab and R.

In this approach, the management of these equations involve the use of poor mnemonics (single or two letter variables) and large chains of equation. This is an error-prone process.

The marjority of the terms in the ODE's is to connect coupled changes between compartments, with a smaller number of single entry/exit terms. Since the coupling is between different equations, it is exceedingly easy to lose track of the coupling.

A more modern approach is to track the coupling between compartments, and use meaningful descriptive variable names. The resultant ODE is then constructed from these couplings. This is essentially automating the double-entry book-keeping between coupled changes between compartments.

By approaching compartmental models in terms of coupling between compartments, it is then very easy to apply stochastic approaches to the integration, as well as to use the information to auto-generate flow-charts of the model.

The model has been carefully considered to use an underlying data structure that is as simple as possible, and allows the results to be easily exported for graphing.





