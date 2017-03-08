import sys
import platform
import os
import glob

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
from basepop import BaseModel


"""

Simple example of SIR (Susceptible-Infected-Recovered) model

@author: sjarvis    21 july 2016

"""


class SIRModel(BaseModel):
    
    # class defaults
    initial_populations = {
        'susceptible': 1e6,
        'infected': 1e2,
        'recovered': 0
    }

    default_params = {
        'beta': 0.2,
        'gamma': 0.01,
        'zero_rate': 0.,
        'rate_birth': 0.00,
    }

    def __init__(self):
        
        BaseModel.__init__(self)
        
        self.name = 'sir'

        for p, v in self.initial_populations.items():
            self.set_compartment(p, v)

        for k, v in self.default_params.items():
            self.set_param(k, v)

    def calculate_vars(self):

        self.vars["population"] = sum(self.compartments.values())
        self.vars["birth"] = self.vars["population"] * self.params["rate_birth"]

    def set_flows(self):

        self.set_var_entry_rate_flow("susceptible", "birth")

        self.set_background_death_rate("zero_rate")
        
        self.set_fixed_transfer_rate_flow("susceptible", "infected", "beta")
        self.set_fixed_transfer_rate_flow("infected", "recovered", "gamma")
        
    def calculate_diagnostic_vars(self):
        
        rate_incidence = 0.0
        for from_label, to_label, rate in self.fixed_transfer_rate_flows:
            val = self.compartments[from_label] * rate
            if 'susceptible' in from_label and 'infected' in to_label:
                rate_incidence += val

        self.vars["incidence"] = rate_incidence / self.vars["population"]

        self.vars["infectious_population"] = self.compartments['infected']
        self.vars["prevalence"] = \
            self.vars["infectious_population"] \
                / self.vars["population"]
        
  
def ensure_dir(out_dir):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)


def make_plots(model, out_dir):
    import pylab

    pylab.clf()
    for var_key in ['incidence', 'prevalence']:
        soln = model.get_var_soln(var_key)
        pylab.plot(model.times, soln, label=var_key)
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, 'fraction.png'))

    pylab.clf()
    for population in model.initial_populations.keys():
        soln = model.population_soln[population]
        pylab.plot(model.times, soln, label=population)
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, 'populations.png'))


def show_pngs(out_dir):
    pngs = glob.glob(os.path.join(out_dir, '*png'))
    operating_system = platform.system()
    if 'Windows' in operating_system:
        os.system("start " + " ".join(pngs))
    elif 'Darwin' in operating_system:
        os.system('open ' + " ".join(pngs))


model = SIRModel()      
model.make_times(1900, 2150, 0.05)
model.integrate(method="explicit")

out_dir = 'sir_graphs'
ensure_dir(out_dir)
model.make_graph(os.path.join(out_dir, 'workflow'))
make_plots(model, out_dir)
show_pngs(out_dir)
