""" 
The following file creates and executes models with SEIR structures. It
then graphs the results and presents the models' underlying compartmental
structure.

Although the structure of these models are simple and widely accepted, the
specific parameter values are taken from the following text: "An Introduction
to Infectious Disease Modelling" by Emilia Vynnycky and Richard G White
available at http://www.anintroductiontoinfectiousdiseasemodelling.com/ with
Excel-based model solutions.

It uses methods from the basepop.BaseModel class in the basepop.py file from
this module (one directory above) to create the model objects for SEIR models.

The purpose of this file is to present examples of how such models can be
built in Python within this popdynamics module. Specifically, the user should
note how inherited methods from basepop.BaseModel are used to ensure processes
such as compartment initiation and setting of flows (entry, transfer and exit)
are performed correctly.

Suggestion to get started: 
- Adjust some parameters within the dictionaries of
    parameter values in infection_param_dictionaries in line 317 and  observe how
    model outputs change. 
- Try adapting the SEIR model without demography to an
    SEIS model, by removing the recovered compartment and changing  the recovery
    transition to move patients from infectious to susceptible (rather than
    recovered). 
"""

from __future__ import print_function
from __future__ import division
from builtins import zip
from past.utils import old_div

import pylab

# hack to allow basepop to be loaded from the parent directory
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
import basepop



###################################
# Define model object classes
###################################


class SeirModel(basepop.BaseModel):
    """
    Based on the SEIR models from Vynnycky and White Chapters 2 and 3 and the
    corresponding online Excel difference equation-based models for measles
    and for flu. Nested inheritance from SirModel to use calculate_vars and
    calculate_diagnostic_vars
    """

    def __init__(self, params):
        """
        Takes a single dictionary of the parameter values to run the SEIR
        model in order that the rest of the SEIR model code can remain the
        same.

        :param seir_param_dictionary: {
            population: Total population size
            start_infectious:Number of infectious individuals at start of simulation
            r0: The R0 value for the infecion
            duration_preinfectious: Number of days spent in the preinfectious compartment
            duration_infectious: Number of days spent in the infectious compartment
        }
        """
        basepop.BaseModel.__init__(self)

        for key in params:
            self.params[key] = params[key]

        # set starting compartment values
        self.set_compartment(
            'susceptible',
            self.params['population'] - self.params['start_infectious'])
        self.set_compartment('preinfectious', 0.0)
        self.set_compartment('infectious', self.params['start_infectious'])
        self.set_compartment('immune', 0.0)

        # set model parameters
        self.set_param(
            'infection_beta',
            1.0 *
            self.params['r0'] /
            (self.params['duration_infectious'] *
             self.params['population']))
        self.set_param(
            'infection_rate_progress',
            1.0 / self.params['duration_preinfectious'])
        self.set_param(
            'infection_rate_recover',
            1.0 / self.params['duration_infectious'])

    def set_flows(self):
        # set variable infection transition flow
        self.set_var_transfer_rate_flow(
            'susceptible', 'preinfectious', 'rate_force')

        # set fixed inter-compartmental flows
        self.set_fixed_transfer_rate_flow(
            'preinfectious', 'infectious', 'infection_rate_progress')
        self.set_fixed_transfer_rate_flow(
            'infectious', 'immune', 'infection_rate_recover')

    def calculate_vars(self):
        # track total population size
        self.vars['population'] = sum(self.compartments.values())

        # calculate force of infection from beta (which was derived from R0 above)
        self.vars['rate_force'] = \
            self.params['infection_beta'] * self.compartments['infectious']

    def calculate_diagnostic_vars(self):
        # calculate incidence
        self.vars['incidence'] = 0.0
        for from_label, to_label, rate in self.fixed_transfer_rate_flows:
            val = self.compartments[from_label] * rate
            if 'infectious' in to_label:
                self.vars['incidence'] += val / self.vars['population']

        # calculate new infections
        self.vars['infections'] = 0.
        for from_label, to_label, rate in self.var_transfer_rate_flows:
            val = self.compartments[from_label] * self.vars[rate]
            if 'preinfectious' in to_label:
                self.vars['infections'] += val / self.vars['population']

        # calculate prevalence
        self.vars['prevalence'] = \
            self.compartments['infectious'] / self.vars['population']


class SeirDemographyModel(SeirModel):
    """
    This is model 3.2 from Vynnycky and White online material, although it is
    largely described in Chapter 4 of the textbook.
    """

    def __init__(self, param_dictionary):

        basepop.BaseModel.__init__(self)

        # set starting compartment values
        self.set_compartment(
            'susceptible',
            param_dictionary['population'] - param_dictionary['start_infectious'])
        self.set_compartment('preinfectious', 0.)
        self.set_compartment('infectious', param_dictionary['start_infectious'])
        self.set_compartment('immune', 0.)

        self.set_param(
            'r0', param_dictionary['r0'])

        # set model parameters
        self.set_param(
            'infection_beta',
            param_dictionary['r0'] / 
                (param_dictionary['duration_infectious'] * param_dictionary['population']))
        self.set_param(
            'infection_rate_progress', old_div(1., param_dictionary['duration_preinfectious']))
        self.set_param(
            'infection_rate_recover', old_div(1., param_dictionary['duration_infectious']))

        self.set_param(
            'demo_rate_death', old_div(1., param_dictionary['life_expectancy']))
        self.set_param(
            'demo_rate_birth', self.params['demo_rate_death'])  # closed population

    def set_flows(self):

        # set variable birth rate
        self.set_var_entry_rate_flow(
            'susceptible', 'rate_birth')

        # set variable infection transition flow
        self.set_var_transfer_rate_flow(
            'susceptible', 'preinfectious', 'rate_force')

        # set fixed inter-compartmental flows
        self.set_fixed_transfer_rate_flow(
            'preinfectious', 'infectious', 'infection_rate_progress')
        self.set_fixed_transfer_rate_flow(
            'infectious', 'immune', 'infection_rate_recover')

        # set background, population-wide death rate
        self.set_background_death_rate(
            'demo_rate_death')

    def calculate_vars(self):

        # track total population size
        self.vars['population'] = sum(self.compartments.values())

        # set birth rate
        self.vars['rate_birth'] = self.params['demo_rate_birth'] * self.vars['population']

        # calculate force of infection from beta (which was derived from R0 above)
        self.vars['rate_force'] = self.params['infection_beta'] * self.compartments['infectious']


###################################
# Plotting functions
###################################


def plot_epidemiological_indicators(model, indicators, png, ylog=False):
    """
    :param indicators: list of keys to model.vars
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
        pylab.savefig(png)
    else:
        pylab.savefig(png)


def plot_compartment_sizes(model, png):
    pylab.clf()
    for compartment in model.compartments:
        soln = model.get_compartment_soln(compartment)
        pylab.plot(model.times, soln, label=compartment)
    pylab.legend()
    pylab.ylabel("persons")
    pylab.title("Populations")
    pylab.savefig(png)


def plot_compartment_proportions(model, png):
    pylab.clf()
    compartment_props = {}
    for compartment in model.compartments:
        soln = model.get_compartment_soln(compartment)
        population = model.get_var_soln("population")
        compartment_props[compartment] = [i / j for i, j in zip(soln, population)]
        pylab.plot(model.times, compartment_props[compartment], label=compartment)
    r_n = [i * infection_params["flu"]["r0"] for i in compartment_props["susceptible"]]
    pylab.plot(model.times, r_n, label="Rn")
    pylab.title("Compartment proportions (and Rn)")
    pylab.legend()
    pylab.savefig(png)


def plot_cyclical_epidemics(model, out_dir):
    # find proportions and Rn for plotting in following sections
    props = {}
    for compartment in ['susceptible', 'immune']:
        soln = model.get_compartment_soln(compartment)
        populations = model.get_var_soln('population')
        props[compartment] = [i / j for i, j in zip(soln, populations)]
    r_n = [i * model.params['r0'] for i in props['susceptible']]

    # plot incidence and Rn (Figure 1)
    pylab.clf()
    fig = pylab.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(
        model.times[365 * 40:],
        model.get_var_soln('incidence')[365 * 40:],
        label='incidence',
        color='r')
    ax1.set_ylabel('incidence per person per day')
    ax1.set_ylim(bottom=0.)
    ax1.legend(loc=2)
    ax2 = ax1.twinx()
    ax2.plot(model.times[365 * 40:], r_n[365 * 40:], label='Rn')
    ax2.set_ylabel('Rn')
    ax2.set_ylim([0.5, 1.5])
    ax2.legend(loc=0)
    fig.savefig(os.path.join(out_dir, 'cyclical_r_n.png'))

    # plot incidence and susceptible proportion (second figure)
    pylab.clf()
    fig = pylab.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(
        model.times[365 * 40:],
        model.get_var_soln('incidence')[365 * 40:],
        label='incidence',
        color='r')
    ax1.set_ylabel('incidence per person per day')
    ax1.set_ylim(bottom=0.)
    ax1.legend(loc=2)
    ax2 = ax1.twinx()
    ax2.plot(
        model.times[365 * 40:],
        props['susceptible'][365 * 40:],
        label='susceptible')
    ax2.set_ylabel('proportion susceptible')
    ax2.set_ylim([0., 0.12])
    ax2.legend(loc=0)
    fig.savefig(os.path.join(out_dir, 'cyclical_prop_susceptible.png'))

    # plot incidence and immune proportion (Figure 3)
    pylab.clf()
    fig = pylab.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(
        model.times[365 * 40:],
        model.get_var_soln('incidence')[365 * 40:],
        label='incidence',
        color='r')
    ax1.set_ylabel('incidence per person per day')
    ax1.set_ylim(bottom=0.)
    ax1.legend(loc=2)
    ax2 = ax1.twinx()
    ax2.plot(
        model.times[365 * 40:],
        props['immune'][365 * 40:],
        label='immune')
    ax2.set_ylabel('proportion immune')
    ax2.set_ylim([0.88, 0.96])
    ax2.legend(loc=0)
    fig.savefig(os.path.join(out_dir, 'cyclical_prop_immune.png'))



# Main 


# SEIR model equivalent to that presented in spreadsheets "model 2.1",
# "model 2.1a" and "model 3.1" of the online materials from Vynnycky
# and White, although these output graphs separate compartment sizes,
# proportions and epidemiological rates from each other

# Define parameters
infection_params = {
    'measles':
        {'population': 1e6,
         'start_infectious': 1.,
         'r0': 13.,
         'duration_preinfectious': 8.,
         'duration_infectious': 7.,
         'life_expectancy': 70. * 365.},
    'flu':
        {'population': 1e6,
         'start_infectious': 1.,
         'r0': 2.,
         'duration_preinfectious': 2.,
         'duration_infectious': 2.,
         'life_expectancy': 70. * 365.}
}

out_dir = "seir_graphs"
basepop.ensure_out_dir(out_dir)

for infection in infection_params.keys():
    print("Running", infection, "SEIR Model for", 200, "days...")

    model = SeirModel(infection_params[infection])
    model.make_times(0, 200, 1)
    model.integrate('explicit')

    # plot results
    model.make_flow_diagram_png(
        os.path.join(out_dir, infection + '_flow_diagram'))

    plot_epidemiological_indicators(
        model,
        ['incidence', 'prevalence', 'infections'],
        os.path.join(out_dir, infection + "_indicators.png"))

    plot_compartment_sizes(
        model, 
        os.path.join(out_dir, infection + '_compartment_sizes.png'))

    if infection == 'flu':
        # equivalent to figure 3 from "model 4.1a" spreadsheet
        plot_epidemiological_indicators(
            model,
            ['incidence'],
            os.path.join(out_dir, infection + "_log_indicators.png"),
            ylog=True)

        # equivalent to figure 2 from "model 4.1a" spreadsheet
        plot_compartment_proportions(
            model,
            os.path.join(out_dir, infection + '_proportions.png'))

basepop.open_pngs_in_dir(out_dir)



# SEIR with demography will generate cyclical effects
# this model is equivalent to that from "model 3.2" spreadsheet

n_day = 365 * 50
infection = "measles"

print("Running SEIR Model with demography for", n_day, "days...")
model = SeirDemographyModel(infection_params[infection])
model.make_times(0, n_day, 1)
model.integrate('explicit')

# set output directory
out_dir = 'seir_demography_graphs'
basepop.ensure_out_dir(out_dir)

# plot results
model.make_flow_diagram_png(
    os.path.join(out_dir, infection + '_flow_diagram'))

plot_epidemiological_indicators(
    model,
    ['incidence', 'prevalence', 'infections'],
    os.path.join(out_dir, infection + "_indicators.png"))

plot_compartment_sizes(
    model,
    os.path.join(out_dir, infection + '_compartment_sizes.png'))

# illustrations of reasons for cyclical epidemics, from 'model 4.3a' spreadsheet
plot_cyclical_epidemics(model, out_dir)

# open output directory
basepop.open_pngs_in_dir(out_dir)



# SEIR with partial population immunity equivalent to figure
# from "model 4.2" spreadsheet in Vynnycky and White (note
# population sizes different and parameters slightly
# different)

partial_immunity_params = {
    'population': 5234.,
    'start_infectious': 2.,
    'r0': 2.1,
    'duration_preinfectious': 2.,
    'duration_infectious': 2.,
    'life_expectancy': 70. * 365.}

print("Running SEIR Model with partial immunity for", 120, "days...")

model = SeirModel(partial_immunity_params)

# redistribute 30% of susceptibles to the immune compartment
model.set_compartment(
    'immune', 0.3 * model.init_compartments['susceptible'])
model.set_compartment(
    'susceptible', 0.7 * model.init_compartments['susceptible'])

# integrate
model.make_times(0, 120, 1)
model.integrate('explicit')

# set output directory
out_dir = 'seir_partial_immunity_graphs'
basepop.ensure_out_dir(out_dir)

plot_compartment_sizes(
    model, os.path.join(out_dir, 'compartment_sizes.png'))

# open output directory
basepop.open_pngs_in_dir(out_dir)
