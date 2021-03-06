"""
Michael Meehan's competitive-strain SIR model.

This demonstrates how to set up a stochastic continous-time model of
two competitive strain SIR models.

There are two infection strains - resident and invader.
They have different r0 - r0_resident and r0_invader.
The model starts off with 1 invader with a more effective r0.

@author: Bosco Ho, December 2017
"""

from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div

# hack to allow basepop to be loaded from the parent directory
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
import basepop

import math

import pylab
import matplotlib


class StrainsModel(basepop.BaseModel):
    """
    Michael Meehan's competitive strain models SIR model.

    There are two infection strains - resident and invader.
    They have different r0 - r0_resident and r0_invader.
    The model starts off with 1 invader with a more effective r0.

    :param input_params = {
        "r0_resident": 2,
        "r0_invader": 3,
        "rate_birth": 2, # lambda
        "rate_death": 1.0e-3, # mu
        "rate_recover": 1.0e-2, # a
        "rate_infection_death": 0, # phi
    }
    """

    def __init__(self, input_params=[]):
        basepop.BaseModel.__init__(self)

        default_params = {
            "r0_resident": 2,
            "r0_invader": 3,
            "rate_birth": 2,  # lambda
            "rate_death": 1.0e-3,  # mu
            "rate_recover": 1.0e-2,  # a
            "rate_infection_death": 0,  # phi
        }
        for key, value in default_params.items():
            self.params[key] = value
        for key, value in input_params:
            self.params[key] = value

        self.set_param(
            "s0", old_div(self.params["rate_birth"],
                          self.params["rate_death"]))
        self.set_param(
            "beta_resident", self.params["r0_resident"] *
            (self.params["rate_recover"] + self.params["rate_infection_death"]
             + self.params["rate_death"]) / self.params["s0"])
        self.set_param(
            "beta_invader", self.params["r0_invader"] *
            (self.params["rate_recover"] + self.params["rate_infection_death"]
             + self.params["rate_death"]) / self.params["s0"])

        # define compartments and set their starting values
        self.set_compartment(
            "susceptible",
            math.ceil(old_div(self.params["s0"], self.params["r0_resident"])))
        self.set_compartment(
            "infectious_resident",
            math.floor(self.params["rate_death"] / self.params["beta_resident"]
                       * (self.params["r0_resident"] - 1)))
        self.set_compartment("infectious_invader", 1)
        self.set_compartment("recovered_resident", 0)
        self.set_compartment("recovered_invader", 0)

    def calculate_vars(self):
        """
        Calculates variables for every time-point in the simulation.
        These can be used as rates for dynamic transmission
        in var transfer flows.
        """
        self.vars["population"] = sum(self.compartments.values())

        self.vars["rate_birth"] = self.params["rate_birth"]

        self.vars["rate_infection_invader"] = \
            self.compartments["infectious_invader"] * \
            self.params["beta_invader"]

        self.vars["rate_infection_resident"] = \
            self.compartments["infectious_resident"] * \
            self.params["beta_resident"]

    def set_flows(self):
        """
        Connects up compartments in the flow diagram of disease
        progression to either a param or var.
        """

        self.set_var_entry_rate_flow("susceptible", "rate_birth")

        self.set_background_death_rate("rate_death")

        self.set_infection_death_rate_flow("infectious_invader",
                                           "rate_infection_death")
        self.set_infection_death_rate_flow("infectious_resident",
                                           "rate_infection_death")

        self.set_var_transfer_rate_flow("susceptible", "infectious_invader",
                                        "rate_infection_invader")
        self.set_var_transfer_rate_flow("susceptible", "infectious_resident",
                                        "rate_infection_resident")

        self.set_fixed_transfer_rate_flow("infectious_resident",
                                          "recovered_resident", "rate_recover")
        self.set_fixed_transfer_rate_flow("infectious_invader",
                                          "recovered_invader", "rate_recover")

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
        self.vars["prevalence"] = 0.
        for label, val in list(self.compartments.items()):
            if "infectious" in label:
                self.vars["prevalence"] += old_div(val,
                                                   self.vars["population"])


def plot_populations(models, png):
    colors = []
    for name in "bgrykcm":
        rgb = matplotlib.colors.colorConverter.to_rgb(name)
        rgba = list(rgb) + [0.1]
        colors.append(rgba)

    pylab.clf()

    y_max = 0
    for i_model, model in enumerate(models):
        for i, compartment in enumerate(model.compartments):
            soln = model.get_compartment_soln(compartment)
            color = colors[i % len(colors)]
            if i_model == 0:
                pylab.plot([0], [0], label=compartment, color=color[:3])
            pylab.plot(model.target_times, soln, linewidth=2, color=color)
            y_max = max(soln.max(), y_max)

    pylab.ylim([0, y_max * 1.1])
    pylab.legend()
    pylab.ylabel("persons")
    pylab.title("Populations")
    pylab.savefig(png)


def plot_extinction(models, png):
    extinctions = []
    for model in models:
        invader_soln = model.get_compartment_soln("infectious_invader")
        is_extinct_invader = invader_soln[-1] == 0
        extinctions.append(is_extinct_invader)

    params = models[0].params
    asymptotic_prob_extinction = 1.0 * params["r0_resident"] / params["r0_invader"]
    title = "Comparing to r0_resident/r0_invader = %.3f" % asymptotic_prob_extinction

    n_replica = len(models)  # could be less for display purposes
    prob_extinction = []
    n_replica_range = list(range(1, n_replica))
    for n_replica in n_replica_range:
        n_extinct = extinctions[:n_replica].count(True)
        prob_extinction.append(1.0 * n_extinct / n_replica)

    pylab.clf()
    pylab.plot(n_replica_range, prob_extinction)
    pylab.ylabel("Probability of invader extinction")
    pylab.ylim([0, max(prob_extinction) * 1.1])
    pylab.xlabel("Number of replicas run")
    pylab.title(title)

    pylab.savefig(png)


# Run replicas of the model

n_replica = 200
models = []
for i_sim in range(n_replica):
    if i_sim % 10 == 0:
        print("Processed", i_sim, "replicas")
    model = StrainsModel()
    model.make_times(0, 1000, 1)
    model.integrate('continuous_time_stochastic')
    models.append(model)

# Generate output

out_dir = "strains_sir_graphs"
basepop.ensure_out_dir(out_dir)
plot_extinction(models, os.path.join(out_dir, "extinction.png"))
plot_populations(models, os.path.join(out_dir, "compartment_sizes.png"))
models[0].make_flow_diagram_png(os.path.join(out_dir, "flow_diagram"))
basepop.open_pngs_in_dir(out_dir)
