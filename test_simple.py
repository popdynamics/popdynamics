
from basepop import BaseModel, make_sigmoidal_curve

import platform
import os
import glob

import pylab


def make_two_step_curve(
    y_low, y_med, y_high, x_start, x_med, x_end):

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


class SimplifiedModel(BaseModel):
    """
    Initial Autumn model designed by James
    """

    def __init__(self):

        BaseModel.__init__(self)

        self.set_compartment("susceptible", 1e6)
        self.set_compartment("latent_early", 0.)
        self.set_compartment("latent_late", 0.)
        self.set_compartment("active", 1.)
        self.set_compartment("treatment_infect", 0.)
        self.set_compartment("treatment_noninfect", 0.)

        self.set_param("demo_rate_birth", 20. / 1e3)
        self.set_param("demo_rate_death", 1. / 65)

        self.set_param("tb_n_contact", 40.)
        self.set_param("tb_rate_earlyprogress", .1 / .5)
        self.set_param("tb_rate_lateprogress", .1 / 100.)
        self.set_param("tb_rate_stabilise", .9 / .5)
        self.set_param("tb_rate_recover", .6 / 3.)
        self.set_param("tb_rate_death", .4 / 3.)

        self.set_param("program_rate_detect", 1.)
        time_treatment = .5
        self.set_param("program_time_treatment", time_treatment)
        self.set_param("program_rate_completion_infect", .9 / time_treatment)
        self.set_param("program_rate_default_infect", .05 / time_treatment)
        self.set_param("program_rate_death_infect", .05 / time_treatment)
        self.set_param("program_rate_completion_noninfect", .9 / time_treatment)
        self.set_param("program_rate_default_noninfect", .05 / time_treatment)
        self.set_param("program_rate_death_noninfect", .05 / time_treatment)

        curve1 = make_sigmoidal_curve(y_high=2, y_low=0, x_start=1950, x_inflect=1970, multiplier=4)
        curve2 = make_sigmoidal_curve(y_high=4, y_low=2, x_start=1995, x_inflect=2003, multiplier=3)
        test_curve = lambda x: curve1(x) if x < 1990 else curve2(x)
        self.set_scaleup_fn("program_rate_detect", test_curve)


    def calculate_vars(self):
        self.vars["population"] = sum(self.compartments.values())
        self.vars["rate_birth"] = \
            self.params["demo_rate_birth"] * self.vars["population"]

        self.vars["infectious_population"] = 0.0
        for label in self.labels:
            if 'active' in label or '_infect' in label:
                self.vars["infectious_population"] += \
                    self.compartments[label]

        self.vars["rate_force"] = \
            self.params["tb_n_contact"] \
            * self.vars["infectious_population"] \
            / self.vars["population"]

    def set_flows(self):
        self.set_var_entry_rate_flow("susceptible", "rate_birth")

        self.set_var_transfer_rate_flow(
            "susceptible", "latent_early", "rate_force")

        self.set_fixed_transfer_rate_flow(
            "latent_early", "active", "tb_rate_earlyprogress")
        self.set_fixed_transfer_rate_flow(
            "latent_early", "latent_late", "tb_rate_stabilise")

        self.set_fixed_transfer_rate_flow(
            "latent_late", "active", "tb_rate_lateprogress")

        self.set_fixed_transfer_rate_flow(
            "active", "latent_late", "tb_rate_recover")

        self.set_var_transfer_rate_flow(
            "active", "treatment_infect", "program_rate_detect")

        self.set_fixed_transfer_rate_flow(
            "treatment_infect", "treatment_noninfect", "program_rate_completion_infect")
        self.set_fixed_transfer_rate_flow(
            "treatment_infect", "active", "program_rate_default_infect")

        self.set_fixed_transfer_rate_flow(
            "treatment_noninfect", "susceptible", "program_rate_completion_noninfect")
        self.set_fixed_transfer_rate_flow(
            "treatment_noninfect", "active", "program_rate_default_noninfect")

        self.set_background_death_rate("demo_rate_death")
        self.set_infection_death_rate_flow(
            "active", "tb_rate_death")
        self.set_infection_death_rate_flow(
            "treatment_infect", "program_rate_death_infect")
        self.set_infection_death_rate_flow(
            "treatment_noninfect", "program_rate_death_noninfect")

    def calculate_diagnostic_vars(self):

        rate_incidence = 0.0
        for from_label, to_label, rate in self.fixed_transfer_rate_flows:
            val = self.compartments[from_label] * rate
            if 'latent' in from_label and 'active' in to_label:
                rate_incidence += val

        # Main epidemiological indicators - note that denominator is not individuals
        self.vars["prevalence"] = \
            self.vars["infectious_population"] \
            / self.vars["population"] * 1E5

        self.vars["incidence"] = \
            rate_incidence \
            / self.vars["population"] * 1E5

        self.vars["mortality"] = \
            self.vars["rate_infection_death"] \
            / self.vars["population"] * 1E5

        self.vars["latent"] = 0.0
        for label in self.labels:
            if "latent" in label:
                self.vars["latent"] += (
                    self.compartments[label]
                    / self.vars["population"] * 1E5)


out_dir = 'test_graphs'
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

model = SimplifiedModel()
model.make_times(1900, 2050, 0.05)
model.integrate_explicit()

model.make_graph(os.path.join(out_dir, 'workflow'))

y_vals = [model.scaleup_fns["program_rate_detect"](t) for t in model.times]
pylab.plot(model.times, y_vals)
pylab.title('scaleup test curve')
pylab.savefig(os.path.join(out_dir, 'scaleup.png'))

pylab.clf()
for var_key in ['mortality', 'incidence', 'prevalence']:
    pylab.plot(model.times, model.get_var_soln(var_key), label=var_key)
pylab.legend()
pylab.savefig(os.path.join(out_dir, 'fraction.png'))

pngs = glob.glob(os.path.join(out_dir, '*png'))
operating_system = platform.system()
if 'Windows' in operating_system:
    os.system("start " + " ".join(pngs))
elif 'Darwin' in operating_system:
    os.system('open ' + " ".join(pngs))


