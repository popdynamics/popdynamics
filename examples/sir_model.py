"""

The following file creates and executes models with SIR and SEIR structures. It then graphs the results and presents the
models" underlying compartmental structure.

Although the structure of these models are simple and widely accepted, the specific parameter values are taken from the
following text:
"An Introduction to Infectious Disease Modelling" by Emilia Vynnycky and Richard G White
available at http://www.anintroductiontoinfectiousdiseasemodelling.com/ with Excel-based model solutions.

It uses methods from the basepop.BaseModel class in the basepop.py file from this module (one directory above) to create the
model objects for SIR and SEIR models.

The purpose of this file is to present examples of how such models can be built in Python within this popdynamics
module. Specifically, the user should note how inherited methods from basepop.BaseModel are used to ensure processes such as
compartment initiation and setting of flows (entry, transfer and exit) are performed correctly.

The first section of the code (to line 116) presents static functions for use in the master script.

The second section of the code (from line 119 to line 336) presents the creation of the model classes for SIR and
SEIR models.

The last section of the code (from line 336 to the end of the script) presents the execution of the example models and
calls the functions to graph their outputs.

Suggestion to get started:
- Adjust some parameters within the dictionaries of parameter values in infection_param_dictionaries in line 317 and
 observe how model outputs change.
- Try adapting the SEIR model without demography to an SEIS model, by removing the recovered compartment and changing
 the recovery transition to move patients from infectious to susceptible (rather than recovered).

"""
from __future__ import division
from builtins import zip
from past.utils import old_div

# hack to allow basepop to be loaded from the examples directory
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
import basepop

import pylab


def plot_epidemiological_indicators(model, indicators, out_dir, ylog=False):
    """
    Plot epidemiological outputs recorded in the model object.

    :param indicators: list of epidemiological indicators within the model"s var attribute
    """
    pylab.clf()
    for var_key in indicators:
        if not ylog:
            pylab.plot(model.times, model.get_var_soln(var_key), label=var_key)
        else:
            pylab.semilogy(model.times, model.get_var_soln(var_key), label=var_key)
    pylab.legend()
    pylab.ylabel("per day (except prevalence), per person")
    pylab.title("Indicators")
    if not ylog:
        pylab.savefig(os.path.join(out_dir, "indicators.png"))
    else:
        pylab.savefig(os.path.join(out_dir, "log_indicators.png"))


def plot_compartment_sizes(model, out_dir):
    pylab.clf()
    y_max = 0
    for compartment in model.compartments:
        soln = model.get_compartment_soln(compartment)
        pylab.plot(model.times, soln, label=compartment)
        y_max = max(soln.max(), y_max)
    pylab.ylim([0, y_max*1.1])
    pylab.legend()
    pylab.ylabel("persons")
    pylab.title("Populations")
    pylab.savefig(os.path.join(out_dir, "compartment_sizes.png"))


def plot_compartment_proportions(model, out_dir):
    pylab.clf()
    for compartment in model.compartments:
        sizes = model.get_compartment_soln(compartment)
        populations = model.get_var_soln("population")
        proportions = [old_div(i, j) for i, j in zip(sizes, populations)]
        pylab.plot(model.times, proportions, label=compartment)
    pylab.title("Compartment proportions")
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, "proportions.png"))


def plot_rn(model, r0, out_dir):
    pylab.clf()
    susceptibles = model.get_compartment_soln("susceptible")
    populations = model.get_var_soln("population")
    proportions = [old_div(i, j) for i, j in zip(susceptibles, populations)]
    r_n = [p * r0 for p in proportions]
    pylab.plot(model.times, r_n, label="Rn")
    pylab.title("Rn")
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, "rn.png"))


def generate_output(model, infection):
    # ensure creation of directory for output
    out_dir = infection + "_sir_graphs"
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # create the flow diagram of the model
    model.make_graph(os.path.join(out_dir, "flow_diagram"))

    indicators = ["incidence", "prevalence"]
    plot_epidemiological_indicators(model, indicators, out_dir)

    plot_compartment_sizes(model, out_dir)
    plot_compartment_proportions(model, out_dir)

    plot_rn(model, model.params["r0"], out_dir)

    basepop.open_pngs_in_dir(out_dir)



class SirModel(basepop.BaseModel):
    """
    Based on the SIR models from Vynnycky and White Chapter 2
    and the corresponding on-line Excel difference equation-based
    models for measles and for flu.
    """

    def __init__(self, params):
        """
        In the initialization, params that are epidemiological
        meaningful are converted to parameters that can be
        expressed as coefficients in the resulting ODE.

        :param params: {
            population: Total population size
            start_infectious: Number of infectious individuals
                at start of simulation
            r0: The R0 value for the infection
            duration_infectious: Number of days spent in the
                infectious compartment
        }
        """

        basepop.BaseModel.__init__(self)

        # define compartments and set their starting values
        self.set_compartment(
            "susceptible", params["population"] - params["start_infectious"])
        self.set_compartment(
            "infectious", params["start_infectious"])
        self.set_compartment(
            "immune", 0.)

        # set model parameters that can be refered to
        # by `self.set_fixed_transfer_rate_flow` and
        # by var calculations
        self.set_param("r0", params["r0"])
        self.set_param(
            "infection_beta",
            old_div(params["r0"], (params["duration_infectious"] * params["population"])))
        self.set_param(
            "infection_rate_recover",
            old_div(1., params["duration_infectious"]))

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


flu_params = {
    "population": 1e6,
    "start_infectious": 1.,
    "r0": 2.,
    "duration_preinfectious": 2.,
    "duration_infectious": 2.,
    "life_expectancy": 70. * 365.
}

measles_params = {
    "population": 1e6,
    "start_infectious": 1.,
    "r0": 12.,
    "duration_preinfectious": 8.,
    "duration_infectious": 7.,
    "life_expectancy": 70. * 365.
}


model = SirModel(measles_params)
model.make_times(0, 100, 1)
model.integrate()

# example graph generation
generate_output(model, "measles")


