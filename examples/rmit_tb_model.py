"""
TB model by Isaac Mwangi.
Nested inheritance from basepop.BaseModel, which applies to any infectious disease generally.
"""
from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div

# hack to allow basepop to be loaded from the parent directory
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
import basepop

import numpy
import pylab


class RmitTbModel(basepop.BaseModel):
    def __init__(self, parameters):
        basepop.BaseModel.__init__(self)

        # define all compartments, initialise as empty and then populate
        model_compartments = ['S', 'L1', 'L2', 'P', 'I', 'R']
        for each_compartment in model_compartments: self.set_compartment(each_compartment, 0.)
        self.set_compartment('S', 1.)
        self.set_compartment('I', 1e-3)

        # parameter setting
        for parameter, value in list(parameters.items()):
            self.set_param(parameter, value)

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

        self.vars['proportion'] = old_div(self.compartments['I'], self.vars['population'])


''' run model '''


# setup directory

out_dir = 'rmit_tb_graphs'
basepop.ensure_out_dir(out_dir)


# setup parameters

fixed_parameters = {
    'mu': old_div(1., 70.),
    'd': .1,
    'phi': 12,
    'f': .05,  # not sure whether 0.05 or 0.5?
    'eta': 2e-4,
    'tau': old_div(1., 2.),
    'alpha': old_div(2., 9.),
    'theta': 1.,
    'rho': 0.1,
    'omega': 2e-5,
}
sigmas = [0., 0., 0.]
sigmas_dict = {}
for sigma in range(len(sigmas)):
    sigmas_dict['sigma' + str(sigma + 1)] = sigmas[sigma]
fixed_parameters.update(sigmas_dict)


# figure 2

betas = list(numpy.linspace(1., 99., 99))
betas += list(numpy.linspace(100., 500., 51))
proportions = []
for beta in betas:
    print('beta', float(beta))
    model = RmitTbModel(fixed_parameters)
    model.set_param('beta', float(beta))
    model.make_times(0, 500., 1.)
    model.integrate(method='explicit')
    proportions.append(model.vars['proportion'])

pylab.clf()
pylab.semilogy(betas, proportions)
text = r''
for sigma in range(len(sigmas)):
    text += r'$\sigma_' + str(sigma + 1) + '$=%.0f\n' % sigmas_dict['sigma' + str(sigma + 1)]
pylab.text(400., 1e-2, text)
pylab.ylabel('proportion of infectives, I')
pylab.xlabel(r'transmission coefficient, $\beta$')
pylab.ylim([9e-7, 1.])
pylab.legend()
pylab.savefig(os.path.join(out_dir, 'figure_2.png'))


# figure 8

sets_of_sigmas = {'a': [0.25, 0.125, 0.125], 'b': [0.25, 0.5, 0.5]}
for fig_letter in sets_of_sigmas:
    pylab.clf()
    sigmas = sets_of_sigmas[fig_letter]
    sigmas_dict = {}
    for sigma in range(len(sigmas)):
        sigmas_dict['sigma' + str(sigma + 1)] = sigmas[sigma]
    fixed_parameters.update(sigmas_dict)
    fixed_parameters['tau'] = 2.
    for rho in [0., 0.5, 1., 10.]:
        proportions = []
        for beta in betas:
            print('rho', float(rho), 'beta', float(beta))
            model = RmitTbModel(fixed_parameters)
            model.set_param('beta', float(beta))
            model.set_param('rho', float(rho))
            model.make_times(0, 500., 1.)
            model.integrate(method='explicit')
            proportions.append(model.vars['proportion'])
        pylab.semilogy(betas, proportions, label=r'$\rho$=' + '{0:.1g}'.format(rho))
    text = r''
    for sigma in range(len(sigmas)):
        text += r'$\sigma_' + str(sigma + 1) + '$=%.2f\n' % sigmas_dict['sigma' + str(sigma + 1)]
    pylab.text(400., 1e-4, text)
    pylab.ylabel('proportion of infectives, I')
    pylab.xlabel(r'transmission coefficient, $\beta$')
    pylab.ylim([9e-7, 1.])
    pylab.legend(frameon=False, loc=8)
    pylab.savefig(os.path.join(out_dir, 'figure_8' + fig_letter + '.png'))


# Open images

basepop.open_pngs_in_dir(out_dir)

