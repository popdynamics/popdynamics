"""
The following file creates and executes models with SIR structure 
using both stochastic and differential models.

It then graphs the results and presents the models" underlying
compartmental structure.
"""

from __future__ import print_function
from __future__ import division
from past.utils import old_div

import copy
import numpy
import pylab
import matplotlib

# # hack to allow basepop to be loaded from the parent directory
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
import basepop


# Create the SIR Model Object

class SirModel(basepop.BaseModel):
    """
    Based on the SIR models from Vynnycky and White Chapter 2
    and the corresponding on-line Excel difference equation-based
    models for measles and for flu.
    """

    def __init__(self, params={}):
        """
        :params param: overridable params = 
        {
            "start_population": 500,
            "start_infectious": 1.,
            "r0": 2.,
            "duration_preinfectious": 2.,
            "duration_infectious": 2.,
            "life_expectancy": 70. * 365.
        }
        """
        basepop.BaseModel.__init__(self)

        default_params = {
            "start_population": 1000,
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
            self.params["start_population"] 
                - self.params["start_infectious"])
        self.set_compartment(
            "infectious",
            self.params["start_infectious"])
        self.set_compartment(
            "immune", 0.)

        self.set_param(
            "infection_beta",
            self.params["r0"] /
            (self.params["duration_infectious"] 
                * self.params["start_population"]))
        self.set_param(
            "infection_rate_recover",
            1. / self.params["duration_infectious"])

    def set_flows(self):
        self.set_var_transfer_rate_flow(
            "susceptible", "infectious", "rate_force")

        self.set_fixed_transfer_rate_flow(
            "infectious", "immune", "infection_rate_recover")

    def calculate_vars(self):
        self.vars["population"] = sum(self.compartments.values())

        self.vars["rate_force"] = \
            self.params["infection_beta"] * \
                self.compartments["infectious"]

    def calculate_diagnostic_vars(self):
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
    :param solutions: list of ["key": Array(Float)] to plot
    :param png: name of output file
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

n_day_stochastic = 10
n_day = 50

model = SirModel({'start_population': 30})

# run the model stochastically for n_day_stochastic
model.make_times(0, n_day_stochastic, 0.2)
model.integrate('continuous_time_stochastic')

model.make_times(n_day_stochastic, n_day, 1)
model.integrate(is_continue=True)

solution = {}
for key in model.compartments.keys():
    solution[key] = model.get_compartment_soln(key)

plot_overlays(
    model.times,
    [solution],
    "Persons",
    "Compartments",
    os.path.join(out_dir, 'compartments'))

basepop.open_pngs_in_dir(out_dir)


