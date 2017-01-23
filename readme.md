---
hello: dddd

---


# basepop

a library for flexible epidemiological modelling

## Why?

- epidemiological models are systems of differential equations
- reasonably simple equations
- compartments
- cascades
- disease vectors
- most come from a `matlab` or `R` - basically an engineering perspective
- write code from a `FORTRAN` like approach
  - think in terms of equations
  - reduce to algebraic symbols
- does not use modern programming techniques
- in `Python` there are lots of features 
  - modern programming
- in really: mapping concepts to algebraic symbols
- differential equations are written
- lots of paperwork to ensure consistency
- hard to change
- modern approach
- use descriptive string to build up model
- write code to autogenerate the differential equation
- using descriptive strings as proxies
- homoiconic the string themselves can be manipulated to be used index.



# Let Python build the model up

- programmatically build up the compartments

- set up links between compartments

- object-oriented programming -> inheritance of objects

- inheritable methods allows separation of concerns

- separate parameter values from generating relationships

- full strings
  - readable
  - no need for dereferencing
  - composable
  - searchable
  - holdover from early days 

- design code should be readable

- avoid all acronyms and abbreviations

- learn to construct compartments, cascades programtically

- isolated code is readable

- intermediate is a list and a callback derivative function on the list

- this allows farming out to any math library

- leverages numpy and scipy with hooks out

- because we have constructed with flows then the flowchar can be built using the `grapgviz` library

- use of dictionaries to keep track of items

- items

  - compartments (nodes)

  - flows (edges)

  - params

  - vars

  - saved arrays

  - reconstruct state

  - automatic graphing of nodes and edges

  - organized data in var

  - dynamic transmission

  - easily extensible

  - vectors, and derivatives

    ​


# Separation of concerns

One of the approaches in `basepop` is the separation of concerns.

- setting up compartments
- defining flows between compartments
- setting up parameters
- setting up dynamic ransmissions
- collecting diagnostic variables
- saving trajectories


By intellignently splitting things up, allows:

- creating automatic flow diagrams for debugging
- allows constructing equations behind the scenes
- thus freeing the researcher from hgeneral house-keeping and focusing on the model


# Concepts

- dynamically create compartments
- build up model through identify flows, let's the model build the vector representation and differentation and derivatives
- use strings to identify compartments and parameters
- transfers between compartments built up using edges defined by the labels
- auto generation of graphs of models
- allows programmatic generation of compartments and transfers
- allows easy scale-up functions


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


## Todo
- test_examples