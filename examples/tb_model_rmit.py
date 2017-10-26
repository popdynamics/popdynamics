
import sys
import platform
import os
import glob

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
from basepop import BaseModel, make_sigmoidal_curve, make_constant_function
import tool_kit


''' Static functions for graphing and managing files '''


def make_plots(model, out_dir):
    """
    Make some basic graphs of the scaling time-variant parameters and the basic model outputs.

    Args:
        model: Instance of the model object to be interrogated
        out_dir: The directory to put the graphs in
    """

    import pylab

    # main epidemiological outputs
    pylab.clf()
    for var_key in ['prevalence']:
        soln = model.get_var_soln(var_key)
        pylab.plot(model.times, soln, label=var_key)
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, 'fraction.png'))


''' Define model object class '''


class SimpleTbModel(BaseModel):
    """
    Initial TB model by James Trauer.
    Nested inheritance from BaseModel, which applies to any infectious disease generally.
    """

    def __init__(self, parameters):
        """
        Inputs:
            interventions: List of interventions to be simulated in the run of the model
        """

        BaseModel.__init__(self)

        # define all compartments, initialise as empty and then populate
        model_compartments = ['S', 'L1', 'L2', 'P', 'I', 'R']
        for each_compartment in model_compartments: self.set_compartment(each_compartment, 0.)
        self.set_compartment('S', 1e6)
        self.set_compartment('I', 1.)

        # parameter setting
        for parameter, value in parameters.items(): self.set_param(parameter, value)

        # parameter processing
        self.set_param('f_phi', self.params['f'] * self.params['phi'])
        self.set_param('1-f_phi', 1. - self.params['f'] * self.params['phi'])

    def calculate_vars(self):
        """
        Calculate values that change with time over the course of model integration.
        """

        # demographic
        self.vars['population'] = sum(self.compartments.values())
        self.vars['lambda'] = self.params['demo_rate_birth'] * self.vars['population']

        # infection
        self.vars['beta_I'] = self.params['beta'] * self.compartments['I'] / self.vars['population']
        for i in range(1, 4):
            self.vars['sigma' + str(i) + '_beta_I'] = self.vars['beta_I'] * self.params['sigma' + str(i)]

    def set_flows(self):
        """
        Set inter-compartmental flows, whether time-variant or constant
        """

        # demographic
        self.set_var_entry_rate_flow('S', 'lambda')
        self.set_background_death_rate('demo_rate_death')

        # infection
        self.set_var_transfer_rate_flow('S', 'L1', 'beta_I')
        self.set_var_transfer_rate_flow('L2', 'L1', 'sigma1_beta_I')
        self.set_var_transfer_rate_flow('P', 'L1', 'sigma2_beta_I')
        self.set_var_transfer_rate_flow('R', 'L1', 'sigma3_beta_I')

        # latency
        self.set_fixed_transfer_rate_flow('L1', 'I', 'f_phi')
        self.set_fixed_transfer_rate_flow('L1', 'L2', '1-f_phi')
        self.set_fixed_transfer_rate_flow('L2', 'I', 'eta')

        # recovery
        self.set_fixed_transfer_rate_flow('I', 'R', 'program_rate_detect')


        self.set_fixed_transfer_rate_flow('L2', 'P', 'rho')
        self.set_fixed_transfer_rate_flow('L1', 'P', 'theta')
        self.set_fixed_transfer_rate_flow('I', 'L2', 'tb_rate_recover')
        self.set_infection_death_rate_flow('I', 'tb_rate_death')
        self.set_fixed_transfer_rate_flow('R', 'I', 'omega')
        self.set_fixed_transfer_rate_flow('R', 'S', 'program_rate_completion')
        self.set_infection_death_rate_flow('R', 'program_rate_death')

    def calculate_diagnostic_vars(self):
        """
        Calculate output variables from model quantities.
        """

        # main epidemiological indicators
        self.vars['prevalence'] = self.compartments['I'] / self.vars['population'] * 1e5


''' Run models '''


if __name__ == '__main__':
    """
    Run and graph a simple TB model with time-variant case detection rate, then run the same model with an intervention
    (BCG vaccination) applied.

    Create a simple TB model without any interventions and a single scaling parameter for case detection rate
    (as shown in the instantiation of the TB model object).
    """

    time_treatment = .5
    fixed_parameters = {
        'demo_rate_birth': 20. / 1e3,
        'demo_rate_death': 1. / 65,
        'f': .1,
        'theta': 0.,
        'beta': 40.,
        'tb_rate_earlyprogress': .1 / .5,
        'eta': .1 / 100.,
        'tb_rate_stabilise': .9 / .5,
        'tb_rate_recover': .6 / 3.,
        'tb_rate_death': .4 / 3.,
        'program_rate_completion': .9 / time_treatment,
        'omega': .05 / time_treatment,
        'program_rate_death': .05 / time_treatment,
        'rho': 0.,
        'sigma1': .25,
        'sigma2': .25,
        'sigma3': .25,
        'program_rate_detect': .7,
        'phi': 12.
    }

    model = SimpleTbModel(fixed_parameters)
    model.make_times(1950, 2000, .05)
    model.integrate(method='explicit')

    # graph outputs
    out_dir = 'tb_graphs'
    tool_kit.ensure_out_dir(out_dir)
    model.make_graph(os.path.join(out_dir, 'workflow'))
    make_plots(model, out_dir)
    tool_kit.open_out_dir(out_dir)

