---
yaml: frontmatter
---

# basepop

a library for flexible epidemiological modelling

# Why?

- epidemiological models
  - most frequently consist of systems of differential equations
  - equations describe relationships between compartments
  - strict conservation of populations
- typical epidemiological modelling
  - inspired by fortran
  - today implemented in `matlab`, `python` or `R` 
  - implemented as differential equations
  - reliance on algebraic symbols
  - focus on the differential equation


# Modern programming approaches

- `Pythonic` approach
- the goal
  - focus on the epidemiology 
  - hide the maths
- alternative view
  - not equation based
  - focus on compartments
  - focus on flows
  - because we have constructed with flows then the flowchar can be built using the `grapgviz` library
- by intelligently splitting things up, allows:
  - creating automatic flow diagrams for debugging
  - allows constructing equations behind the scenes
  - thus freeing the researcher from general house-keeping and focusing on the model


  - dynamically create compartments
  - build up model through identify flows, let's the model build the vector representation and differentiation and derivatives
  - use strings to identify compartments and parameters
  - transfers between compartments built up using edges defined by the labels
  - auto generation of graphs of models
  - allows programmatic generation of compartments and transfers
  - allows easy scale-up functions


# General approach

- setting up compartments
- defining flows between compartments
- setting up parameters
- setting up dynamic transmissions
- collecting diagnostic variables
- saving trajectories


# Examples

- create a simple model
- create a model with dynamic transmission
- show how to extract out vectors and derivatives
- scale-up function
- programmatically generate compartments
- programmatically generate diagnostic variables
- flow chart
- produce graphs with matplotlib - javascript option

# Requirements

- python 2.7
- numpy
- scipy 
- graphviz.py 
- graphviz binary on path: http://www.graphviz.org/

