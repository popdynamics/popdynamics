
# basepop

a library for flexible epidemiological modelling

## Preamble

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

# Concepts

- dynamically create compartments
- build up model through identify flows, let's the model build the vector representation and differentation and derivatives
- use strings to identify compartments and parameters
- transfers between compartments built up using edges defined by the labels
- auto generation of graphs of models
- allows programmatic generation of compartments and transfers
- allows easy scale-up functions



## Requirements

- python 2.7
- numpy
- scipy 
- graphviz.py 
- (graphviz)[http://www.graphviz.org/)


## Todo
- haven't initialized self.death_flow @done
- initialize birth_flows and death_flows to zero @done
- added birth_flows check @done
- generalized death rate/ background_death_rate @done
- add scale_ups @done
- examples of modular code: nest & brian simulators
- more documentation
- more test_examples