import sys
import platform
import os
import glob
import copy

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
from basepop import BaseModel, make_sigmoidal_curve


class SimplifiedModel(BaseModel):
    """
    Initial TB by James Trauer
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
        two_step_curve = lambda x: curve1(x) if x < 1990 else curve2(x)
        self.set_scaleup_fn("program_rate_detect", two_step_curve)


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


def ensure_dir(out_dir):
    try: 
        os.makedirs(out_dir)
    except OSError:
        if not os.path.isdir(out_dir):
            raise


def make_plots(model, out_dir):
    import pylab

    scaleup_fn = model.scaleup_fns["program_rate_detect"]
    y_vals = map(scaleup_fn, model.times)
    pylab.plot(model.times, y_vals)
    pylab.title('scaleup test curve')
    pylab.savefig(os.path.join(out_dir, 'scaleup.png'))

    pylab.clf()
    for var_key in ['mortality', 'incidence', 'prevalence']:
        soln = model.get_var_soln(var_key)
        pylab.plot(model.times, soln, label=var_key)
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, 'fraction.png'))


def show_pngs(out_dir):
    pngs = glob.glob(os.path.join(out_dir, '*png'))
    operating_system = platform.system()
    if 'Windows' in operating_system:
        os.system("start " + " ".join(pngs))
    elif 'Darwin' in operating_system:
        os.system('open ' + " ".join(pngs))


if __name__ == "__main__":

    model = SimplifiedModel()
    model.make_times(1900, 2050, 0.05)
    model.integrate(method="explicit")

    out_dir = 'tb_graphs'
    ensure_dir(out_dir)
    model.make_graph(os.path.join(out_dir, 'workflow'))
    make_plots(model, out_dir)
    show_pngs(out_dir)



