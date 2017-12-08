# Popdynamics

Popdynamics is a library for building epidemiological models, or any kind of compartmental model.

It provides an easy way to switch between different integration methods - both deterministic ODE solutions, and stochastic discrete-time and continuous-time solutions.

# Why?

- epidemiological models
  - equations describe relationships between compartments
  - changes in compartments modelled with differential equations or stochastic models
  - changes are modelled with coupling to rates
  - dynamic transmission, the coupling is dependent on other populations
- typical epidemiological modelling
  - inspired by fortran, used in matlab and R
  - use matrices, indexing of compartments
  - long complex  differential equations
- as compartments are indexed, long chains of reasoning to force into matrices
- lots of manual housekeeping to ensure equations are properly entered
- changes are very difficult

# Modern programming approach

Basepop provides a modern approach to epidemiological modelling. It provides modularization, simple data-model, and lots of flexibility. 

In Basepop, the focus is not on compartments and their rate-of-change equation, but on individual flows between compartment. The rate-of-change equation is then constructed when the model is run. This means that the house-keeping to keep the equations balanced are done for you.

As well, Basepop can use the individual flows to construct equivalent stochastic models.

As much as possible, the objects are kept light and extensible. The relevant structures are stored in very simple data-structures: dictionaries of values and lists of tuples. The results are stored as numpy arrays.

Models can be built up with object inheritance or created dynamically upon instantiation.

Useful auxiliary functions include graphviz generation.

# General sequential approach

- setting up compartments
- defining flows between compartments
- setting up parameters
- setting up dynamic transmissions
- collecting diagnostic variables
- saving trajectories

# Examples

- create a model with dynamic transmission
- show how to extract out vectors and derivatives
- programmatically generate compartments 
- programmatically generate diagnostic variables
- saving variables for use later
- create a simple model @done
- scale-up function @done
- flow chart @done
- produce graphs with matplotlib - javascript option @done

# Requirements

- python 2.7
- numpy
- scipy 
- graphviz.py 
- graphviz binary on path: http://www.graphviz.org/

