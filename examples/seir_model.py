
import platform
import os
import glob
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from basepop import BaseModel
import pylab


def check_out_dir_exists(out_dir):
    """
    Make sure the output directory exists and create if it doesn't.
    """
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)


def open_out_dir():
    """
    Open the output dir at the end of the model run.
    """

    pngs = glob.glob(os.path.join(out_dir, '*png'))
    operating_system = platform.system()
    if 'Windows' in operating_system:
        os.system("start " + " ".join(pngs))
    elif 'Darwin' in operating_system:
        os.system('open ' + " ".join(pngs))


def plot_epidemiological_indicators(model, infection, indicators, out_dir, ylog=False):
    """
    Plot epidemiological outputs recorded in the model object.

    Inputs:
        indicators: List of the strings that refer to the epidemiological indicators within the model's var attribute
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
        pylab.savefig(os.path.join(out_dir, infection + "_indicators.png"))
    else:
        pylab.savefig(os.path.join(out_dir, infection + "_log_indicators.png"))


def plot_compartment_sizes(model):
    """
    Plot compartment sizes over time.
    """

    pylab.clf()
    for compartment in model.compartments:
        pylab.plot(model.times, model.get_compartment_soln(compartment), label=compartment)
    pylab.legend()
    pylab.ylabel("persons")
    pylab.title("Populations")
    pylab.savefig(os.path.join(out_dir, 'compartment_sizes.png'))


def plot_compartment_proportions(model, infection):
    """
    Plot compartment proportions over time.
    """

    pylab.clf()
    compartment_props = {}
    for compartment in model.compartments:
        compartment_props[compartment] \
            = [i / j for i, j in zip(model.get_compartment_soln(compartment), model.get_var_soln("population"))]
        pylab.plot(model.times, compartment_props[compartment], label=compartment)
    r_n = [i * infection_param_dictionaries["flu"]["r0"] for i in compartment_props["susceptible"]]
    pylab.plot(model.times, r_n, label="Rn")
    pylab.title("Compartment proportions (and Rn)")
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, infection + '_proportions.png'))


###################################
### Define model object classes ###
###################################


class SirModel(BaseModel):
    """
    Based on the SIR models from Vynnycky and White Chapter 2
    and the corresponding on-line Excel difference equation-based models for measles and for flu.
    """

    def __init__(self, param_dictionary):
        """
        Takes a single dictionary of the parameter values to run the SIR model in order that the rest of the SIR
        model code can remain the same.

        Inputs:
            param sir_param_dictionary: Dictionary with keys as follows:
                population:             Total population size
                start_infectious:       Number of infectious individuals at start of simulation
                r0:                     The R0 value for the infecion
                duration_infectious:    Number of days spent in the infectious compartment
        """

        BaseModel.__init__(self)

        # set starting compartment values
        self.set_compartment("susceptible",
                             param_dictionary["population"] - param_dictionary["start_infectious"])
        self.set_compartment("infectious", param_dictionary["start_infectious"])
        self.set_compartment("immune", 0.)

        # set model parameters
        self.set_param("infection_beta",
                       param_dictionary["r0"]
                       / (param_dictionary["duration_infectious"] * param_dictionary["population"]))
        self.set_param("infection_rate_recover", 1. / param_dictionary["duration_infectious"])

    def calculate_vars(self):

        # track total population size
        self.vars["population"] = sum(self.compartments.values())

        # calculate force of infection from beta (which was derived from R0 above)
        self.vars["rate_force"] = self.params["infection_beta"] * self.compartments["infectious"]

    def set_flows(self):

        # set variable infection transition flow
        self.set_var_transfer_rate_flow("susceptible", "infectious", "rate_force")

        # set fixed inter-compartmental flows
        self.set_fixed_transfer_rate_flow("infectious", "immune", "infection_rate_recover")

    def calculate_diagnostic_vars(self):

        # calculate incidence
        self.vars["incidence"] = 0.
        for from_label, to_label, rate in self.var_transfer_rate_flows:
            val = self.compartments[from_label] * self.vars[rate]
            if "infectious" in to_label:
                self.vars["incidence"] += val / self.vars["population"]

        # calculate prevalence
        self.vars["prevalence"] = self.compartments["infectious"] / self.vars["population"]


class SeirModel(SirModel):
    """
    Based on the SEIR models from Vynnycky and White Chapters 2 and 3
    and the corresponding on-line Excel difference equation-based models for measles and for flu.
    Nested inheritance from SirModel to use calculate_vars and calculate_diagnostic_vars
    """

    def __init__(self, param_dictionary):
        """
        Takes a single dictionary of the parameter values to run the SEIR model in order that the rest of the SEIR
        model code can remain the same.

        Inputs:
            param seir_param_dictionary: Dictionary with keys as follows:
                population:             Total population size
                start_infectious:       Number of infectious individuals at start of simulation
                r0:                     The R0 value for the infecion
                duration_preinfectious: Number of days spent in the preinfectious compartment
                duration_infectious:    Number of days spent in the infectious compartment
        """

        BaseModel.__init__(self)

        # set starting compartment values
        self.set_compartment("susceptible",
                             param_dictionary["population"] - param_dictionary["start_infectious"])
        self.set_compartment("preinfectious", 0.)
        self.set_compartment("infectious", param_dictionary["start_infectious"])
        self.set_compartment("immune", 0.)

        # set model parameters
        self.set_param("infection_beta",
                       param_dictionary["r0"]
                       / (param_dictionary["duration_infectious"] * param_dictionary["population"]))
        self.set_param("infection_rate_progress", 1. / param_dictionary["duration_preinfectious"])
        self.set_param("infection_rate_recover", 1. / param_dictionary["duration_infectious"])

    def set_flows(self):

        # set variable infection transition flow
        self.set_var_transfer_rate_flow("susceptible", "preinfectious", "rate_force")

        # set fixed inter-compartmental flows
        self.set_fixed_transfer_rate_flow("preinfectious", "infectious", "infection_rate_progress")
        self.set_fixed_transfer_rate_flow("infectious", "immune", "infection_rate_recover")

    def calculate_vars(self):

        # track total population size
        self.vars["population"] = sum(self.compartments.values())

        # calculate force of infection from beta (which was derived from R0 above)
        self.vars["rate_force"] = self.params["infection_beta"] * self.compartments["infectious"]

    def calculate_diagnostic_vars(self):

        # calculate incidence
        self.vars["incidence"] = 0.
        for from_label, to_label, rate in self.fixed_transfer_rate_flows:
            val = self.compartments[from_label] * rate
            if "infectious" in to_label:
                self.vars["incidence"] += val / self.vars["population"]

        # calculate new infections
        self.vars["infections"] = 0.
        for from_label, to_label, rate in self.var_transfer_rate_flows:
            val = self.compartments[from_label] * self.vars[rate]
            if "preinfectious" in to_label:
                self.vars["infections"] += val / self.vars["population"]

        # calculate prevalence
        self.vars["prevalence"] = self.compartments["infectious"] / self.vars["population"]


class SeirDemographyModel(SeirModel):

    """
    This is model 3.2 from Vynnycky and White online material, although it is largely described in Chapter 4 of the
    textbook.
    """

    def __init__(self, param_dictionary):

        BaseModel.__init__(self)

        # set starting compartment values
        self.set_compartment("susceptible",
                             param_dictionary["population"] - param_dictionary["start_infectious"])
        self.set_compartment("preinfectious", 0.)
        self.set_compartment("infectious", param_dictionary["start_infectious"])
        self.set_compartment("immune", 0.)

        # set model parameters
        self.set_param("infection_beta",
                       param_dictionary["r0"]
                       / (param_dictionary["duration_infectious"] * param_dictionary["population"]))
        self.set_param("infection_rate_progress", 1. / param_dictionary["duration_preinfectious"])
        self.set_param("infection_rate_recover", 1. / param_dictionary["duration_infectious"])
        self.set_param("demo_rate_death", 1. / param_dictionary["life_expectancy"])
        self.set_param("demo_rate_birth", self.params["demo_rate_death"])  # closed population

    def set_flows(self):

        # set variable birth rate
        self.set_var_entry_rate_flow("susceptible", "rate_birth")

        # set variable infection transition flow
        self.set_var_transfer_rate_flow("susceptible", "preinfectious", "rate_force")

        # set fixed inter-compartmental flows
        self.set_fixed_transfer_rate_flow("preinfectious", "infectious", "infection_rate_progress")
        self.set_fixed_transfer_rate_flow("infectious", "immune", "infection_rate_recover")

        # set background, population-wide death rate
        self.set_background_death_rate("demo_rate_death")

    def calculate_vars(self):

        # track total population size
        self.vars["population"] = sum(self.compartments.values())

        # set birth rate
        self.vars["rate_birth"] = self.params["demo_rate_birth"] * self.vars["population"]

        # calculate force of infection from beta (which was derived from R0 above)
        self.vars["rate_force"] = self.params["infection_beta"] * self.compartments["infectious"]


######################################
### Create and run SIR/SEIR models ###
######################################

# define parameter values for two SEIR infections - measles and influenza
infection_param_dictionaries = {
    "measles":
        {"population": 1e6,
         "start_infectious": 1.,
         "r0": 13.,
         "duration_preinfectious": 8.,
         "duration_infectious": 7.,
         "life_expectancy": 70. * 365.},
    "flu":
        {"population": 1e6,
         "start_infectious": 1.,
         "r0": 2.,
         "duration_preinfectious": 2.,
         "duration_infectious": 2.,
         "life_expectancy": 70. * 365.}
}

#############################
### Create and run models ###
#############################

# loop over infection types
for infection in ["flu", "measles"]:

    # SIR
    model = SirModel(infection_param_dictionaries[infection])
    model.make_times(0, 200, 1)
    model.integrate("explicit")

    # set output directory
    out_dir = infection + '_sir_graphs'
    check_out_dir_exists(out_dir)

    # plot results
    model.make_graph(os.path.join(out_dir, infection + '_flow_diagram'))
    plot_epidemiological_indicators(model, infection, ["incidence", "prevalence"], out_dir)
    plot_compartment_sizes(model)

    # SEIR
    # this model is equivalent to that presented in spreadsheets "model 2.1", "model 2.1a" and "model 3.1" of the online materials
    # from Vynnycky and White, although output graphs separate compartment sizes, proportions and epidemiological rates
    model = SeirModel(infection_param_dictionaries[infection])
    model.make_times(0, 200, 1)
    model.integrate("explicit")

    # set output directory
    out_dir = infection + '_seir_graphs'
    check_out_dir_exists(out_dir)

    # plot results
    model.make_graph(os.path.join(out_dir, infection + '_flow_diagram'))
    plot_epidemiological_indicators(model, infection, ["incidence", "prevalence", "infections"], out_dir)
    plot_compartment_sizes(model)

    if infection == "flu":
        # equivalent to figure 3 from "model 4.1a" spreadsheet in Vynnycky and White online materials
        plot_epidemiological_indicators(model, "flu", ["incidence"], out_dir, ylog=True)

        # equivalent to figure 2 from "model 4.1a" spreadsheet in Vynnycky and White online materials
        plot_compartment_proportions(model, infection)

    # open output directory
    open_out_dir()

# SEIR demography
# this model is equivalent to that from "model 3.2" spreadsheet of the Vynnycky and White online materials
model = SeirDemographyModel(infection_param_dictionaries["measles"])
model.make_times(0, 36500, 1)
model.integrate("explicit")

# set output directory
out_dir = "measles_seir_demography_graphs"
check_out_dir_exists(out_dir)

# plot results
model.make_graph(os.path.join(out_dir, "measles_flow_diagram"))
plot_epidemiological_indicators(model, "measles", ["incidence", "prevalence", "infections"], out_dir)
plot_compartment_sizes(model)

# illustrations of reasons for cyclical epidemics, from "model 4.3a" spreadsheet
model.make_times(0, 365 * 50, 1)
model.integrate("explicit")

# find proportions
compartment_props = {}
for compartment in ["susceptible", "immune"]:
    compartment_props[compartment] \
        = [i / j for i, j in zip(model.get_compartment_soln(compartment), model.get_var_soln("population"))]
r_n = [i * infection_param_dictionaries["measles"]["r0"] for i in compartment_props["susceptible"]]

# plot incidence and susceptible proportion (second figure)
pylab.clf()
fig = pylab.figure()
ax1 = fig.add_subplot(111)
ax1.plot(model.times[365 * 40:], model.get_var_soln("incidence")[365 * 40:], label="incidence", color='r')
ax1.set_ylabel("incidence per person per day")
ax1.set_ylim(bottom=0.)
ax1.legend(loc=2)
ax2 = ax1.twinx()
ax2.plot(model.times[365 * 40:], compartment_props["susceptible"][365 * 40:], label="susceptible")
ax2.set_ylabel("proportion susceptible")
ax2.set_ylim([0., 0.12])
ax2.legend(loc=0)
fig.savefig(os.path.join(out_dir, "cyclical_prop_susceptible.png"))

# plot incidence and immune proportion (Figure 3)
pylab.clf()
fig = pylab.figure()
ax1 = fig.add_subplot(111)
ax1.plot(model.times[365 * 40:], model.get_var_soln("incidence")[365 * 40:], label="incidence", color='r')
ax1.set_ylabel("incidence per person per day")
ax1.set_ylim(bottom=0.)
ax1.legend(loc=2)
ax2 = ax1.twinx()
ax2.plot(model.times[365 * 40:], compartment_props["immune"][365 * 40:], label="immune")
ax2.set_ylabel("proportion immune")
ax2.set_ylim([0.88, 0.96])
ax2.legend(loc=0)
fig.savefig(os.path.join(out_dir, "cyclical_prop_immune.png"))

# plot incidence and Rn (Figure 1)
pylab.clf()
fig = pylab.figure()
ax1 = fig.add_subplot(111)
ax1.plot(model.times[365 * 40:], model.get_var_soln("incidence")[365 * 40:], label="incidence", color='r')
ax1.set_ylabel("incidence per person per day")
ax1.set_ylim(bottom=0.)
ax1.legend(loc=2)
ax2 = ax1.twinx()
ax2.plot(model.times[365 * 40:], r_n[365 * 40:], label="Rn")
ax2.set_ylabel("Rn")
ax2.set_ylim([0.5, 1.5])
ax2.legend(loc=0)
fig.savefig(os.path.join(out_dir, "cyclical_r_n.png"))

# SEIR with partial population immunity from start
# equivalent to figure from "model 4.2" spreadsheet in Vynnycky and White
partial_immunity_params = {"population": 5234.,
                           "start_infectious": 2.,
                           "r0": 2.1,
                           "duration_preinfectious": 2.,
                           "duration_infectious": 2.,
                           "life_expectancy": 70. * 365.}
model = SeirModel(partial_immunity_params)

# redistribute 30% of susceptibles to the immune compartment
model.set_compartment("immune", 0.3 * model.init_compartments["susceptible"])
model.set_compartment("susceptible", 0.7 * model.init_compartments["susceptible"])

# integrate
model.make_times(0, 120, 1)
model.integrate("explicit")

# set output directory
out_dir = 'partial_immunity_graphs'
check_out_dir_exists(out_dir)
plot_compartment_sizes(model)
