"""
The following file creates and executes models with SIR structure. It
then graphs the results and presents the models" underlying
compartmental structure.

Although the structure of the model are simple and widely accepted,
the specific parameter values are taken from the following text: "An
Introduction to Infectious Disease Modelling" by Emilia Vynnycky and
Richard G White available at
http://www.anintroductiontoinfectiousdiseasemodelling.com/ with Excel-
based model solutions.

It uses methods from the BaseModel class in the basepop.py file from
this module (one directory above) to create the model objects for SIR
and SEIR models.

The purpose of this file is to present examples of how such models can
be built in Python within this popdynamics module. Specifically, the
user should note how inherited methods from BaseModel are used to
ensure processes such as compartment initiation and setting of flows
(entry, transfer and exit) are performed correctly.

Suggestion to get started:

- Adjust some parameters within the dictionaries of parameter values
in infection_param_dictionaries in line 317 and observe how model
outputs change.

- Try adapting the SEIR model without demography to an SEIS model, by
removing the recovered compartment and changing the recovery
transition to move patients from infectious to susceptible (rather
than recovered).
"""

from __future__ import print_function
from __future__ import division
from builtins import zip
from past.utils import old_div

import copy

import numpy

# # hack to allow basepop to be loaded from the parent directory
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
import basepop

import pylab
import matplotlib

# Create the SIR Model Object

class SirModel(basepop.BaseModel):
    """
    Based on the SIR models from Vynnycky and White Chapter 2
    and the corresponding on-line Excel difference equation-based
    models for measles and for flu.
    """

    def __init__(self, params={}):
        """
        In the initialization, params that are epidemiological
        meaningful are converted to parameters that can be
        expressed as coefficients in the resulting ODE.

        :param params: e.g. {
            "population": 500,
            "start_infectious": 1.,
            "r0": 2.,
            "duration_preinfectious": 2.,
            "duration_infectious": 2.,
            "life_expectancy": 70. * 365.
        }
        """
        basepop.BaseModel.__init__(self)

        default_params = {
            "population": 1000,
            "start_infectious": 1.,
            "r0": 12.,
            "duration_preinfectious": 8.,
            "duration_infectious": 7.,
            "life_expectancy": 70. * 365.
        }
        for key, value in default_params.items():
            self.set_param(key, value)
        for key, value in params.items():
            self.set_param(key, value)

        # define compartments and set their starting values
        self.set_compartment(
            "susceptible",
            self.params["population"] - self.params["start_infectious"])
        self.set_compartment(
            "infectious",
            self.params["start_infectious"])
        self.set_compartment(
            "immune", 0.)

        # set model parameters that can be refered to
        # by `self.set_fixed_transfer_rate_flow` and
        # by var calculations
        self.set_param(
            "infection_beta",
            self.params["r0"] /
            (self.params["duration_infectious"] * self.params["population"]))
        self.set_param(
            "infection_rate_recover",
            1. / self.params["duration_infectious"])

    def set_flows(self):
        """
        Connects up compartments in the flow diagram of disease
        progression. Each flow between two compartment refers
        to either a fixed param or a var. The values of a var
        in the calculation will be calculated in `self.calculate_vars`
        for each time-point.
        """

        # set variable infection transition flow, the rate refers
        # to values in self.vars, which are recalculated at every
        # time-step in self.calculate_vars.
        self.set_var_transfer_rate_flow(
            "susceptible", "infectious", "rate_force")

        # set fixed inter-compartmental flows, the rate refers
        # to values in self.params, which are fixed
        self.set_fixed_transfer_rate_flow(
            "infectious", "immune", "infection_rate_recover")

    def calculate_vars(self):
        """
        Calculates variables for every time-point in the simulation.
        These can be used as rates for dynamic transmission
        in var transfer flows.
        """

        # track total population size
        self.vars["population"] = sum(self.compartments.values())

        # force of infection, infection_beta is derived from R0 in __init__
        self.vars["rate_force"] = \
            self.params["infection_beta"] * \
                self.compartments["infectious"]

    def calculate_diagnostic_vars(self):
        """
        Calculates diagnostic variables at the end of the simulation.
        These are only calculated for the specified time-points, using
        cached values of the compartments from the simulation run.
        """

        # calculate incidence
        self.vars["incidence"] = 0.
        for from_label, to_label, rate in self.var_transfer_rate_flows:
            val = self.compartments[from_label] * self.vars[rate]
            if "infectious" in to_label:
                self.vars["incidence"] += old_div(val, self.vars["population"])

        # calculate prevalence
        self.vars["prevalence"] = \
            old_div(self.compartments["infectious"], self.vars["population"])



# Plotting functions for the model

def plot_overlays(times, solutions, ylabel, title, png):
    """
    :param times: list of [Float]
    :param solutions: list of ["key": Array(Float)]
    :param png: string
    """
    colors = []
    for name in "bgrykcm":
        rgb = matplotlib.colors.colorConverter.to_rgb(name)
        if len(solutions) > 1:
            rgba = list(rgb) + [0.1]
            colors.append(rgba)
        else:
            colors.append(rgb)

    pylab.clf()

    y_max = 0
    for i_soln, soln in enumerate(solutions):
        for i_key, key in enumerate(soln):
            y_vals = soln[key]
            color = colors[i_key % len(colors)]
            if i_soln == 0:
                # generate a fake dot so that legend can extract color/label
                pylab.plot([0], [0], label=key, color=color[:3])
            pylab.plot(times, y_vals, linewidth=2, color=color)
            y_max = max(max(y_vals), y_max)

    pylab.ylim([0, y_max * 1.1])
    pylab.legend()
    pylab.ylabel(ylabel)
    pylab.title(title)

    pylab.savefig(png)


# The main routine

out_dir = "mix_sir_graphs"
basepop.ensure_out_dir(out_dir)

model1 = SirModel({'population': 20})
model1.make_times(0, 10, 1)
model1.integrate_continuous_stochastic()

model2 = model1.clone()
model2.init_compartments = copy.deepcopy(model2.compartments)
model2.make_times(10, 50, 1)
model2.integrate()

solution = {}
times = None
for model in [model1, model2]:
    for key in model.compartments.keys():
        new_solution = numpy.copy(model.get_compartment_soln(key))
        if key not in solution:
            solution[key] = new_solution
        else:
            solution[key] = numpy.concatenate(
                (solution[key], new_solution[1:]), axis=0)
    if times is None:
        times = copy.deepcopy(model.times)
    else:
        times.extend(model.times[1:])

print(times)
plot_overlays(
    times,
    [solution],
    "Persons",
    "Compartments",
    os.path.join(out_dir, 'compartments'))

basepop.open_pngs_in_dir(out_dir)


