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
- increasing trend to more complex models
  - to simulate realistic dynamics, modern models may need stratification by:
        - compartment (as patient progresses through infection, activation, detection, treatment, recovery, death, etc.)
        - manifestation (organ involved, infectiousness)
        - drug resistance of infecting organism
        - population/demographic risk group
        - medical comorbidity
        - age
        - others
  - stratifications are often multiplicative, meaning that realistic models may require hundreds to thousands
        of compartments or more
  - it then becomes difficult to impossible to track and understand the behaviour of every compartment individuals
  - also difficulty/impossible to write down the system of ODEs and more prone to mistakes

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
  - build up model through identify flows,
        let's the model build the vector representation and differentiation and derivatives
  - use strings to identify compartments and parameters
  - transfers between compartments built up using edges defined by the labels
  - auto generation of graphs of models
  - allows programmatic generation of compartments and transfers
  - allows easy scale-up functions
  - avoidance of errors, such as:
    - transition flows that do not have either an entry or an exit compartment,
        or the entry value is not equal to the exit
    - entry (birth) flows that do not have an entry compartment or have an exit compartment
    - exit (death) flows that have an entry compartment or do not have an exit compartment
    - fixed flows that are not associated with a constant parameter or variable flows that are not associated with a
        time-variant parameter to be calculated at each model step
- by ensuring methods and functions are as independent as possible, adaptations to model can usually be achieved through
    changes to a single method/function
    - for example, further integration methods can be added to explicit and scipy by changing a single method only
    - similarly, additional output diagnostics can be added by changing a single method
        - e.g. those specific to a particular organism/disease

# General sequential approach

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

