
import os
from basepop import BaseModel
import numpy
import pylab
import tool_kit


class RmitTbModel(BaseModel):
    """
    TB model by Isaac Mwangi.
    Nested inheritance from BaseModel, which applies to any infectious disease generally.
    """

    def __init__(self, parameters):
        BaseModel.__init__(self)

        # define all compartments, initialise as empty and then populate
        model_compartments = ['S', 'L1', 'L2', 'P', 'I', 'R']
        for each_compartment in model_compartments: self.set_compartment(each_compartment, 0.)
        self.set_compartment('S', 1.)
        self.set_compartment('I', 1e-3)

        # parameter setting
        for parameter, value in parameters.items(): self.set_param(parameter, value)

        # parameter processing
        self.set_param('f_phi', self.params['f'] * self.params['phi'])
        self.set_param('1-f_phi', (1. - self.params['f']) * self.params['phi'])
        self.set_param('tau_alpha', self.params['tau'] + self.params['alpha'])

    def calculate_vars(self):
        """
        Calculate values that change with time over the course of model integration - only forces of infection.
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

        self.vars['proportion'] = self.compartments['I'] / self.vars['population']


''' Run model '''


if __name__ == '__main__':
    fixed_parameters = {
        'mu': 1. / 70.,
        'd': .1,
        'phi': 1.5,  # not sure whether 1.5 or 12 or value in between?
        'f': .05,  # not sure whether 0.05 or 0.1?
        'eta': 2e-4,
        'tau': 1. / 2.,
        'alpha': 2. / 9.,
        'theta': 1.,  # not sure whether 0. or 1.?
        'rho': 0.,
        'omega': 2e-5,
        'sigma1': 0.,
        'sigma2': 0.,
        'sigma3': 0.
    }

    # figure 2
    betas = list(numpy.linspace(1., 99., 99))
    betas += list(numpy.linspace(100., 500., 51))
    proportions = []
    for beta in betas:
        model = RmitTbModel(fixed_parameters)
        model.set_param('beta', float(beta))
        model.make_times(0, 500., 1.)
        model.integrate(method='explicit')
        proportions.append(model.vars['proportion'])
    out_dir = 'tb_graphs'
    print(betas)
    pylab.clf()
    pylab.semilogy(betas, proportions)
    pylab.text(400., 1e-2, r'$\sigma_1$=0' + '\n' + r'$\sigma_2$=0' + '\n' + r'$\sigma_3$=0')
    pylab.ylabel('proportion of infectives, I')
    pylab.xlabel(r'transmission coefficient, $\beta$')
    pylab.ylim([9e-7, 1.])
    pylab.legend()
    pylab.savefig(os.path.join(out_dir, 'proportion_beta.png'))
    tool_kit.open_out_dir(out_dir)


