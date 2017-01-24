---
yaml: frontmatter
---

# basepop

a library for flexible epidemiological modelling

## Why?

- epidemiological models
  - systems of differential equations
  - equations connecting compartments
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
  - hide the math
- alternative view
  - not equation based
  - focus on compartments
  - and gflows
  - focus on compartments
  - focus on flows
  - because we have constructed with flows then the flowchar can be built using the `grapgviz` library
- By intellignently splitting things up, allows:
  - creating automatic flow diagrams for debugging
  - allows constructing equations behind the scenes
  - thus freeing the researcher from general house-keeping and focusing on the model


  - dynamically create compartments
  - build up model through identify flows, let's the model build the vector representation and differentation and derivatives
  - use strings to identify compartments and parameters
  - transfers between compartments built up using edges defined by the labels
  - auto generation of graphs of models
  - allows programmatic generation of compartments and transfers
  - allows easy scale-up functions


# General approach

- setting up compartments
- defining flows between compartments
- setting up parameters
- setting up dynamic ransmissions
- collecting diagnostic variables
- saving trajectories


# Examples

- Create a simple model
- Create a model with dynamic transmission
- Show how to extract out vectors and derivatives
- Scale-up function
- Programmatically generate compartments
- Programmtically generate diagnostic variables
- Flow chart
- Produce graphs with matplotlib - javascript option

## Requirements

- python 2.7
- numpy
- scipy 
- graphviz.py 
- graphviz binary on path: http://www.graphviz.org/

