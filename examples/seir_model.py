
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


def plot_epidemiological_indicators(model, infection, indicators, out_dir):
    """
    Plot epidemiological outputs recorded in the model object.

    Inputs:
        indicators: List of the strings that refer to the epidemiological indicators within the model's var attribute
    """

    pylab.clf()
    for var_key in indicators:
        pylab.plot(model.times, model.get_var_soln(var_key), label=var_key)
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, infection + '_indicators.png'))


def plot_compartment_sizes(model):
    """
    Plot compartment sizes over time.
    """

    pylab.clf()
    for compartment in model.compartments:
        pylab.plot(model.times, model.get_compartment_soln(compartment), label=compartment)
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, 'compartment_sizes.png'))


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
        for from_label, to_label, rate in self.fixed_transfer_rate_flows:
            val = self.compartments[from_label] * rate
            if "infectious" in to_label:
                self.vars["incidence"] += val / self.vars["population"] * 1E5

        # calculate prevalence
        self.vars["prevalence"] = self.compartments["infectious"] / self.vars["population"] * 1E5


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
    model = SeirModel(infection_param_dictionaries[infection])
    model.make_times(0, 200, 1)
    model.integrate("explicit")

    # set output directory
    out_dir = infection + '_seir_graphs'
    check_out_dir_exists(out_dir)

    # plot results
    model.make_graph(os.path.join(out_dir, infection + '_flow_diagram'))
    plot_epidemiological_indicators(model, infection, ["incidence", "prevalence"], out_dir)
    plot_compartment_sizes(model)

    # open output directory
    open_out_dir()

# SEIR demography
model = SeirDemographyModel(infection_param_dictionaries["measles"])
model.make_times(0, 36500, 1)
model.integrate("explicit")

# set output directory
out_dir = "measles_seir_demography_graphs"
check_out_dir_exists(out_dir)

# plot results
model.make_graph(os.path.join(out_dir, "measles_flow_diagram"))
plot_epidemiological_indicators(model, "measles", ["incidence", "prevalence"], out_dir)
plot_compartment_sizes(model)


