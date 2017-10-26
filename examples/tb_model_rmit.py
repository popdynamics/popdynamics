
import os
from basepop import BaseModel
import tool_kit


''' static functions for graphing and managing files '''


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
    pylab.savefig(os.path.join(out_dir, 'prevalence_time.png'))

    # main epidemiological outputs
    pylab.clf()
    for var_key in ['population']:
        soln = model.get_var_soln(var_key)
        pylab.plot(model.times, soln, label=var_key)
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, 'population_time.png'))


''' define model object class '''


class RmitTbModel(BaseModel):
    """
    TB model by Isaac Mwangi.
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
        self.set_param('1-f_phi', (1. - self.params['f']) * self.params['phi'])
        self.set_param('tau_alpha', self.params['tau'] + self.params['alpha'])

    def calculate_vars(self):
        """
        Calculate values that change with time over the course of model integration.
        """

        # demographic
        self.vars['population'] = sum(self.compartments.values())
        self.vars['lambda'] = self.params['mu'] * self.vars['population'] + self.params['d'] * self.compartments['I']

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
        self.set_background_death_rate('mu')

        # infection
        self.set_var_transfer_rate_flow('S', 'L1', 'beta_I')
        self.set_var_transfer_rate_flow('L2', 'L1', 'sigma1_beta_I')
        self.set_var_transfer_rate_flow('P', 'L1', 'sigma2_beta_I')
        self.set_var_transfer_rate_flow('R', 'L1', 'sigma3_beta_I')

        # latency
        self.set_fixed_transfer_rate_flow('L1', 'I', 'f_phi')
        self.set_fixed_transfer_rate_flow('L1', 'L2', '1-f_phi')
        self.set_fixed_transfer_rate_flow('L2', 'I', 'eta')

        # latent treatment
        self.set_fixed_transfer_rate_flow('L1', 'P', 'theta')
        self.set_fixed_transfer_rate_flow('L2', 'P', 'rho')

        # recovery
        self.set_fixed_transfer_rate_flow('I', 'R', 'tau_alpha')
        self.set_infection_death_rate_flow('I', 'd')

        # relapse
        self.set_fixed_transfer_rate_flow('R', 'I', 'omega')

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

    fixed_parameters = {
        'beta': 40.,
        'mu': 1. / 70.,
        'd': .1,
        'phi': 12.,
        'f': .1,
        'eta': 2e-4,
        'tau': 1. / 2.,
        'alpha': 2. / 9.,
        'theta': 0.,
        'rho': 0.,
        'omega': 2e-5,
        'sigma1': .25,
        'sigma2': .5,
        'sigma3': .5
    }

    model = RmitTbModel(fixed_parameters)
    model.make_times(1950, 2000, .05)
    model.integrate(method='explicit')

    # graph outputs
    out_dir = 'tb_graphs'
    tool_kit.ensure_out_dir(out_dir)
    model.make_graph(os.path.join(out_dir, 'workflow'))
    make_plots(model, out_dir)
    tool_kit.open_out_dir(out_dir)

