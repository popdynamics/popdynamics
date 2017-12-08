"""

Michael Meehan's competitive-strain SIR model.

@author: Bosco Ho, December 2017
"""

from __future__ import print_function
import platform
import os
import glob
import sys
import math

import pylab
import matplotlib

# hack to allow basepop to be loaded from the examples directory
parent_dir = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, parent_dir)
from basepop import BaseModel


class StrainsModel(BaseModel):
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

    def __init__(self, input_params=None):
        BaseModel.__init__(self)

        if input_params is None:
            input_params = {}

        default_params = {
            "r0_resident": 2,
            "r0_invader": 3,
            "rate_birth": 2,  # lambda
            "rate_death": 1.0e-3,  # mu
            "rate_recover": 1.0e-2,  # a
            "rate_infection_death": 0,  # phi
        }

        for key, value in default_params.items():
            if key not in input_params:
                input_params[key] = value

        for key, value in input_params.items():
            self.set_param(key, value)

        self.set_param(
            "s0",
            self.params["rate_birth"] / self.params["rate_death"])
        self.set_param(
            "beta_resident",
            self.params["r0_resident"] *
                (self.params["rate_recover"] +
                 self.params["rate_infection_death"] +
                 self.params["rate_death"]) /
            self.params["s0"])
        self.set_param(
            "beta_invader",
            self.params["r0_invader"] *
                (self.params["rate_recover"] +
                 self.params["rate_infection_death"] +
                 self.params["rate_death"]) /
            self.params["s0"] )

        # define compartments and set their starting values
        self.set_compartment(
            "susceptible",
            math.ceil(self.params["s0"] / self.params["r0_resident"]))
        self.set_compartment(
            "infectious_resident",
            math.floor(
                self.params["rate_death"] /
                self.params["beta_resident"] *
                (self.params["r0_resident"] - 1)))
        self.set_compartment(
            "infectious_invader", 1)
        self.set_compartment(
            "recovered_resident", 0)
        self.set_compartment(
            "recovered_invader", 0)

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

        self.set_infection_death_rate_flow(
            "infectious_invader", "rate_infection_death")
        self.set_infection_death_rate_flow(
            "infectious_resident", "rate_infection_death")

        self.set_var_transfer_rate_flow(
            "susceptible", "infectious_invader", "rate_infection_invader")
        self.set_var_transfer_rate_flow(
            "susceptible", "infectious_resident", "rate_infection_resident")

        self.set_fixed_transfer_rate_flow(
            "infectious_resident", "recovered_resident", "rate_recover")
        self.set_fixed_transfer_rate_flow(
            "infectious_invader", "recovered_invader", "rate_recover")

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
                self.vars["incidence"] += val / self.vars["population"]

        # calculate prevalence
        self.vars["prevalence"] = 0.
        for label, val in self.compartments.items():
            if "infectious" in label:
                self.vars["prevalence"] += val / self.vars["population"]


# Run replicas of the model

n_replica = 200
models = []
for i_sim in range(n_replica):
    if i_sim % 10 == 0:
        print("Processed", i_sim, "replicas")
    model = StrainsModel()
    model.make_times(0, 1000, 1)
    # model.integrate("explicit")
    model.integrate_continuous_stochastic()
    # model.integrate_discrete_time_stochastic()
    models.append(model)


# Set output directory

infection = "strains"
out_dir = infection + "_sir_graphs"
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)


# Create probability of extinction graph

extinction = []
for model in models:
    invader_soln = model.get_compartment_soln("infectious_invader")
    is_extinct_invader = invader_soln[-1] == 0
    extinction.append(is_extinct_invader)

params = models[0].params
asymptotic_prob_extinction = 1.0 * params["r0_resident"] / params["r0_invader"]

pylab.clf()

prob_extinction = []
n_replica_range = range(1, n_replica)
for n_replica in n_replica_range:
    n_extinct = extinction[:n_replica].count(True)
    prob_extinction.append(1.0 * n_extinct / n_replica)

pylab.plot(n_replica_range, prob_extinction)

pylab.ylabel("Probability of invader extinction")
pylab.ylim([0, max(prob_extinction)*1.1])
pylab.xlabel("Number of replicas run")
pylab.title("Comparing to r0_resident/r0_invader = %.3f" % asymptotic_prob_extinction)

pylab.savefig(os.path.join(out_dir, "extinction.png"))


# Create population overlay graphs

pylab.clf()
y_max = 0

print("Making png")
cm = matplotlib.cm.get_cmap("Pastel1")
colors = []
for name in "bgrykcm":
    rgb = matplotlib.colors.colorConverter.to_rgb(name)
    rgba = list(rgb) + [0.1]
    colors.append(rgba)

for i_model, model in enumerate(models):
    for i, compartment in enumerate(model.compartments):
        soln = model.get_compartment_soln(compartment)
        color = colors[i % len(colors)]
        if i_model == 0:
            pylab.plot([0], [0], label=compartment, color=color[:3])
        pylab.plot(model.times, soln, linewidth=2, color=color)
        y_max = max(soln.max(), y_max)

pylab.ylim([0, y_max*1.1])
pylab.legend()
pylab.ylabel("persons")
pylab.title("Populations")

pylab.savefig(os.path.join(out_dir, "compartment_sizes.png"))


# Create the flow diagram of the model

models[0].make_graph(os.path.join(out_dir, "flow_diagram"))


# Open output images for convenience

pngs = glob.glob(os.path.join(out_dir, "*png"))
operating_system = platform.system()
if "Windows" in operating_system:
    os.system("start " + " ".join(pngs))
elif "Darwin" in operating_system:
    os.system("open " + " ".join(pngs))
