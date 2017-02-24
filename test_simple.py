
from basepop import BaseModel, make_sigmoidal_curve
import platform
import os
import glob
import pylab


# Some static functions to avoid code repetition
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


def make_two_step_curve(y_low, y_med, y_high, x_start, x_med, x_end):

    curve1 = make_sigmoidal_curve(
        y_high=y_med, y_low=y_low,
        x_start=x_start, x_inflect=(x_med-x_start)*0.5 + x_start,
        multiplier=4)

    curve2 = make_sigmoidal_curve(
        y_high=y_high, y_low=y_med,
        x_start=x_med, x_inflect=(x_end-x_med)*0.5 + x_med,
        multiplier=4)

    def curve(x):
        if x < x_start:
            return y_low
        if x < x_med:
            return curve1(x)
        if x < x_end:
            return curve2(x)
        return y_high

    return curve


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


def plot_scaleup_functions(model, infection):
    """
    Plot scale-up functions (time-variant parameters) in model.
    """

    pylab.clf()
    y_vals = [model.scaleup_fns["program_rate_detect"](t) for t in model.times]
    pylab.plot(model.times, y_vals)
    pylab.title(infection + " scaleup test curve")
    pylab.savefig(os.path.join(out_dir, infection + "_scaleup.png"))


def plot_compartment_sizes(model):
    """
    Plot compartment sizes over time.
    """

    pylab.clf()
    for compartment in model.compartments:
        pylab.plot(model.times, model.get_compartment_soln(compartment), label=compartment)
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, 'compartment_sizes.png'))


# define model objects based on BaseModel from basepop.py
class SeirModel(BaseModel):
    """
    Based on the SEIR models from Vynnycky and White Chapter 3 and the corresponding on-line Excel
    difference equation-based models for measles and for flu.
    """

    def __init__(self, seir_param_dictionary):
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
                             seir_param_dictionary["population"] - seir_param_dictionary["start_infectious"])
        self.set_compartment("preinfectious", 0.)
        self.set_compartment("infectious", seir_param_dictionary["start_infectious"])
        self.set_compartment("immune", 0.)

        # set model parameters
        self.set_param("infection_beta",
                       seir_param_dictionary["r0"]
                       / (seir_param_dictionary["duration_infectious"] * seir_param_dictionary["population"]))
        self.set_param("infection_rate_progress", 1. / seir_param_dictionary["duration_preinfectious"])
        self.set_param("infection_rate_recover", 1. / seir_param_dictionary["duration_infectious"])

    def calculate_vars(self):

        # track total population size
        self.vars["population"] = sum(self.compartments.values())

        # calculate force of infection from beta (which was derived from R0 above)
        self.vars["rate_force"] = self.params["infection_beta"] * self.compartments["infectious"]

    def set_flows(self):

        # set variable infection transition flow
        self.set_var_transfer_rate_flow("susceptible", "preinfectious", "rate_force")

        # set fixed inter-compartmental flows
        self.set_fixed_transfer_rate_flow("preinfectious", "infectious", "infection_rate_progress")
        self.set_fixed_transfer_rate_flow("infectious", "immune", "infection_rate_recover")

    def calculate_diagnostic_vars(self):

        # calculate incidence
        self.vars["incidence"] = 0.
        for from_label, to_label, rate in self.fixed_transfer_rate_flows:
            val = self.compartments[from_label] * rate
            if "infectious" in to_label:
                self.vars["incidence"] += val / self.vars["population"] * 1E5

        # calculate prevalence
        self.vars["prevalence"] = \
            self.compartments["infectious"] \
            / self.vars["population"] * 1E5


class TbModel(BaseModel):
    """
    Initial Autumn model designed by James
    """

    def __init__(self):

        BaseModel.__init__(self)

        # set compartment sizes
        self.set_compartment("susceptible", 1e6)
        self.set_compartment("latent_early", 0.)
        self.set_compartment("latent_late", 0.)
        self.set_compartment("active", 1.)
        self.set_compartment("treatment_infect", 0.)
        self.set_compartment("treatment_noninfect", 0.)

        # set demographic rates
        self.set_param("demo_rate_birth", 20. / 1e3)
        self.set_param("demo_rate_death", 1. / 65)

        # set disease-specific parameters
        self.set_param("tb_n_contact", 40.)
        self.set_param("tb_rate_earlyprogress", .1 / .5)
        self.set_param("tb_rate_lateprogress", .1 / 100.)
        self.set_param("tb_rate_stabilise", .9 / .5)
        self.set_param("tb_rate_recover", .6 / 3.)
        self.set_param("tb_rate_death", .4 / 3.)

        # set programmatic parameters
        self.set_param("program_rate_detect", 1.)
        time_treatment = .5
        self.set_param("program_time_treatment", time_treatment)
        self.set_param("program_rate_completion_infect", .9 / time_treatment)
        self.set_param("program_rate_default_infect", .05 / time_treatment)
        self.set_param("program_rate_death_infect", .05 / time_treatment)
        self.set_param("program_rate_completion_noninfect", .9 / time_treatment)
        self.set_param("program_rate_default_noninfect", .05 / time_treatment)
        self.set_param("program_rate_death_noninfect", .05 / time_treatment)

        # set time-variant programmatic parameter (for case detection rate)
        curve1 = make_sigmoidal_curve(y_high=2, y_low=0, x_start=1950, x_inflect=1970, multiplier=4)
        curve2 = make_sigmoidal_curve(y_high=4, y_low=2, x_start=1995, x_inflect=2003, multiplier=3)
        test_curve = lambda x: curve1(x) if x < 1990 else curve2(x)
        self.set_scaleup_fn("program_rate_detect", test_curve)

    def calculate_vars(self):

        # demographic variables
        self.vars["population"] = sum(self.compartments.values())
        self.vars["rate_birth"] = \
            self.params["demo_rate_birth"] * self.vars["population"]

        # infectious population for force of infection calculations
        self.vars["infectious_population"] = 0.
        for label in self.labels:
            if 'active' in label or '_infect' in label:
                self.vars["infectious_population"] += \
                    self.compartments[label]

        # force of infection for direct use as a flow rate
        self.vars["rate_force"] = \
            self.params["tb_n_contact"] \
            * self.vars["infectious_population"] \
            / self.vars["population"]

    def set_flows(self):

        # demographic
        self.set_var_entry_rate_flow("susceptible", "rate_birth")
        self.set_population_death_rate()

        # infection
        self.set_var_transfer_rate_flow("susceptible", "latent_early", "rate_force")

        # disease progression
        self.set_fixed_transfer_rate_flow("latent_early", "active", "tb_rate_earlyprogress")
        self.set_fixed_transfer_rate_flow("latent_early", "latent_late", "tb_rate_stabilise")
        self.set_fixed_transfer_rate_flow("latent_late", "active", "tb_rate_lateprogress")
        self.set_fixed_transfer_rate_flow("active", "latent_late", "tb_rate_recover")
        self.set_infection_death_rate_flow("active", "tb_rate_death")
        self.set_infection_death_rate_flow("treatment_infect", "program_rate_death_infect")
        self.set_infection_death_rate_flow("treatment_noninfect", "program_rate_death_noninfect")

        # programmatically variable rate(s)
        self.set_var_transfer_rate_flow("active", "treatment_infect", "program_rate_detect")

        # progression through treatment
        self.set_fixed_transfer_rate_flow("treatment_infect", "treatment_noninfect", "program_rate_completion_infect")
        self.set_fixed_transfer_rate_flow("treatment_infect", "active", "program_rate_default_infect")
        self.set_fixed_transfer_rate_flow("treatment_noninfect", "susceptible", "program_rate_completion_noninfect")
        self.set_fixed_transfer_rate_flow("treatment_noninfect", "active", "program_rate_default_noninfect")

    def calculate_diagnostic_vars(self):

        # progressively sum incidence
        rate_incidence = 0.
        for from_label, to_label, rate in self.fixed_transfer_rate_flows:
            val = self.compartments[from_label] * rate
            if 'latent' in from_label and 'active' in to_label:
                rate_incidence += val

        # calculate main epidemiological indicators
        self.vars["prevalence"] \
            = self.vars["infectious_population"] / self.vars["population"] * 1E5
        self.vars["incidence"] \
            = rate_incidence / self.vars["population"] * 1E5
        self.vars["mortality"] \
            = self.vars["rate_infection_death"] / self.vars["population"] * 1E5

        # calculate proportion of population latently infected (including early and late)
        self.vars["latent"] = 0.
        for label in self.labels:
            if "latent" in label:
                self.vars["latent"] += self.compartments[label] / self.vars["population"] * 1E5


# run models
if __name__ == "__main__":

    ###################
    ### SEIR models ###
    ###################

    # define parameter values for two SEIR infections - measles and influenza
    seir_param_dictionary = {
        "measles":
            {"population": 1e6,
             "start_infectious": 1.,
             "r0": 13.,
             "duration_preinfectious": 8.,
             "duration_infectious": 7.},
        "flu":
            {"population": 1e6,
             "start_infectious": 1.,
             "r0": 2.,
             "duration_preinfectious": 2.,
             "duration_infectious": 2.}
    }

    # loop over SEIR infections
    for infection in ["flu", "measles"]:

        # set output directory
        out_dir = infection + '_graphs'
        check_out_dir_exists(out_dir)

        # instantiate model object
        model = SeirModel(seir_param_dictionary[infection])
        model.make_times(0, 200, 1)
        model.integrate_explicit()

        # plot results
        model.make_graph(os.path.join(out_dir, infection + '_flow_diagram'))
        plot_epidemiological_indicators(model, infection, ["incidence", "prevalence"], out_dir)
        plot_compartment_sizes(model)

        # open output directory
        open_out_dir()

    #################
    ### TB models ###
    #################

    # make directory for TB outputs
    out_dir = "tb_graphs"
    check_out_dir_exists(out_dir)

    # create and run TB model
    model = TbModel()
    model.make_times(1900, 2050, 0.05)
    model.integrate_explicit()

    # plot results
    model.make_graph(os.path.join(out_dir, "tb_flow_diagram"))
    plot_epidemiological_indicators(model, "tb", ["incidence", "prevalence", "mortality"], out_dir)
    plot_scaleup_functions(model, "tb")

    # open output directory
    open_out_dir()



