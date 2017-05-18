# -*- coding: utf-8 -*-


"""
Base Population Model to handle different type of models
"""

import os
from scipy.integrate import odeint
import numpy
from math import exp


def add_unique_tuple_to_list(a_list, a_tuple):
    """
    Adds or modifies a list of tuples, compares only the items before the last in the tuples,
    the last value in the tuple is assumed to be a value.
    """

    for i, test_tuple in enumerate(a_list):
        if test_tuple[:-1] == a_tuple[:-1]:
            a_list[i] = a_tuple
            break
    else:
        a_list.append(a_tuple)


def label_intersects_tags(label, tags):
    """
    Determine whether a string is contained within a list of strings for use in functions such as calculation of the
    force of infection, where we might want a list of all the compartments that contain a particular string
    (such as 'active' or 'infectious').

    Args:
        label: The string we're searching for
        tags: List for comparison
    """

    for tag in tags:
        if tag in label:
            return True
    return False


def make_sigmoidal_curve(y_low=0, y_high=1., x_start=0, x_inflect=0.5, multiplier=1.):
    """
    Function to make a sigmoidal curve for smooth scaling of time-variant parameter values

    Args:
        y_low: lowest y value
        y_high: highest y value
        x_inflect: inflection point of graph along the x-axis
        multiplier: if 1, slope at x_inflect goes to (0, y_low), larger
                    values makes it steeper

    Returns:
        function that increases sigmoidally from 0 y_low to y_high
        the halfway point is at x_inflect on the x-axis and the slope
        at x_inflect goes to (0, y_low) if the multiplier is 1.
    """

    amplitude = y_high - y_low
    if amplitude == 0:
        def curve(x):
            return y_low

        return curve

    x_delta = x_inflect - x_start
    slope_at_inflection = multiplier * 0.5 * amplitude / x_delta
    b = 4. * slope_at_inflection / amplitude

    def curve(x):
        arg = b * (x_inflect - x)
        # check for large values that will blow out exp
        if arg > 10.:
            return y_low
        return amplitude / (1. + exp(arg)) + y_low

    return curve


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


class BaseModel:
    """
    Basic concepts
      var - values that are calculated at every time step
      param - values that remain fixed throughout the model run
    """

    def __init__(self):

        # list of labels for all compartments
        self.labels = []

        # stores the initial value for all compartments
        self.init_compartments = {}

        # stores the values of all parameters
        # there should be no hard-coded values except as contained in this structure
        self.params = {}

        # stored list of time points
        self.times = None

        # scale-up functions, generally used for time-variant parameters,
        # whose values change in a way that is predictable before the model has been run
        self.scaleup_fns = {}

        # stores any auxillary variables used to calculate dynamic effects (such as transmission) at each time-step
        self.vars = {}

        # total flow of each compartment
        self.flows = {}

        # variable entry (birth) rate flow(s)
        # list of 2-tuple (label, var_label)
        #   - label: name of compartment
        #   - var_label: name of var that holds the entry rate
        self.var_entry_rate_flow = []

        # fixed transfer rates (often for progression between infected states after infection has occurred)
        # list of 3-tuple (from_label, to_label, param_label)
        #   - from_label: name of compartment that loses population
        #   - to_label: name of compartment that gains population
        #   - param_label: name of param that holds the rate
        self.fixed_transfer_rate_flows = []

        # variable transfer rates (can be used either for flows that vary due to predictable effects
        # - scale-up functions - or for flows that vary due to the dynamics of the infection - force of infection)
        # list of 3-tuple (from_label, to_label, var_label)
        #   - from_label: name of compartment that loses population
        #   - to_label: name of compartment that gains population
        #   - var_label: name of var that holds the rate
        self.var_transfer_rate_flows = []

        # fixed infection death rate flows (Should only be applied to persons with active disease as a consequence of
        # the infection, while latent or preinfectious individuals should be assumed to take the population-wide
        # death rates only. Also note that for actively diseased persons, this rate is added to the population rate.)
        # list of 2-tuple (label, rate)
        #   - label: name of compartment
        #   - rate: disease specific death rate
        self.infection_death_rate_flows = []

        # the generalised death rate of all compartments
        self.background_death_rate = 0.

        # other universally required model attributes
        self.soln_array = None
        self.var_labels = None
        self.var_array = None
        self.flow_array = None

    def make_times(self, start, end, delta):
        """
        Make a list of times for integration to be performed at. Units are arbitrary.

        Args:
            start: Numerical time to start at (can be a calendar year, or zero for epidemic time, or other)
            end: Numerical end time
            delta: Interval between times
        Creates:
            self.times: List of numerical times for integration to be assessed at
        """

        assert type(start) is float or type(start) is int, 'Start time not specified with float'
        assert type(end) is float or type(end) is int, 'End time not specified with a number'
        assert type(delta) is float or type(delta) is int, 'Time increment not specified with a number'
        assert end >= start, 'End time is before start time'
        self.times = []
        step = start
        while step <= end:
            self.times.append(step)
            step += delta

    def make_times_with_n_step(self, start, end, n):
        """
        Alternative to make_times. For this one, the third argument is the number of time points required, rather than
        the step between time points.

        Args:
            start: As for make_times
            end: As for make_times
            n: Number of time points that are needed
        """

        self.times = []
        step = start
        delta = (end - start) / float(n)
        while step <= end:
            self.times.append(step)
            step += delta

    def set_compartment(self, label, init_val=0.):
        """
        Create a population compartment for the model.

        Args:
            label: String to describe the compartment
            init_val: Starting value for that compartment (default behaviour is to start from empty compartment)
        """

        assert type(label) is str, 'Compartment label for initial setting not string'
        assert type(init_val) is float, 'Value to start % compartment from not string' % label
        assert init_val >= 0., 'Start with negative compartment not permitted'
        if label not in self.labels:
            self.labels.append(label)
        self.init_compartments[label] = init_val

    def set_param(self, label, val):
        """
        Add a parameter value to the dictionary of parameter values

        Args:
            label: String name of the compartment
            val: Value (generally float) for the parameter
        """

        assert type(label) is str, 'Parameter name is not string'
        assert type(val) is float or type(val) is int, 'Fixed parameter value is not numeric for %' % label
        self.params[label] = val

    def convert_list_to_compartments(self, vec):
        """
        Distribute the list of compartment values to create a dictionary (reverse of convert_compartments_to_list).

        Args:
            vec: List of compartment values ordered according to the labels list
        Returns:
            Dictionary with keys from labels attribute and values from vec
        """

        return {l: vec[i] for i, l in enumerate(self.labels)}

    def convert_compartments_to_list(self, compartments):
        """
        Distribute compartments dictionary to a list for making derivative function.

        Args:
            compartments: Dictionary of compartments
        Returns:
            List of compartment values with equivalent ordering to labels
        """

        return [compartments[l] for l in self.labels]

    def get_init_list(self):
        """
        Convert starting compartment size dictionary to list for integration.

        Returns:
            List of values for starting compartments
        """

        return self.convert_compartments_to_list(self.init_compartments)

    def set_scaleup_fn(self, label, fn):
        """
        Simple method to add a scale-up function to the dictionary of scale-ups.

        Args:
            label: String for name of function
            fn: The function to be added
        """

        assert type(label) is str, 'Name of scale-up function is not string'
        self.scaleup_fns[label] = fn

    def set_background_death_rate(self, param_label='demo_rate_death'):
        """
        Sets the population death rate to be applied to all compartments.

        Args:
            param_label: String for the population death rate
        """

        assert type(self.background_death_rate) is float, 'Background death rate is not float'
        assert self.background_death_rate >= 0., 'Background death rate is negative'
        self.background_death_rate = self.params[param_label]

    def set_infection_death_rate_flow(self, label, param_label):
        """
        Set an additional death rate for those with active infection.

        Args:
            label: Compartment to apply the additional death rate to
            param_label: String of the parameter to be used for setting the death rate
        """

        assert type(label) is str, 'Compartment label not string when setting infection death rate'
        assert type(param_label) is str, 'Parameter label not string when setting infection death rate'
        add_unique_tuple_to_list(self.infection_death_rate_flows, (label, self.params[param_label]))

    def set_fixed_transfer_rate_flow(self, from_label, to_label, param_label):
        """
        Set constant inter-compartmental flow, such as those related to progression through infection states after
        infection.

        Args:
            from_label: Compartment that this flow comes from
            to_label: Compartment that this flow goes into
            param_label: String of the parameter to be used for setting this transition rate
        """

        assert type(from_label) is str, 'Origin compartment label not string for setting fixed transfer rate'
        assert type(to_label) is str, 'Destination compartment label not string for setting fixed transfer rate'
        add_unique_tuple_to_list(self.fixed_transfer_rate_flows, (from_label, to_label, self.params[param_label]))

    def set_var_transfer_rate_flow(self, from_label, to_label, var_label):
        """
        Set variable inter-compartmental flow - as can be used for setting predictable time-variant inter-compartmental
        flows (such as scale-up functions) or for those that have to be calculated from the model (such as the force
        of infection).

        Args:
            from_label: Compartment that this flow comes from
            to_label: Compartment that this flow goes into
            param_label: String of the parameter to be used for setting this transition rate
        """

        assert type(from_label) is str, 'Origin compartment label not string for setting variable transfer rate'
        assert type(to_label) is str, 'Destination compartment label not string for setting variable transfer rate'
        assert type(var_label) is str, 'Function label not string for setting fixed transfer rate'
        add_unique_tuple_to_list(self.var_transfer_rate_flows, (from_label, to_label, var_label))

    def set_var_entry_rate_flow(self, label, var_label):
        """
        Set variable entry rate to model, to be used for new births into the population (and migration if required).

        Args:
            label: Compartment that this flow goes in to
            var_label: String of the var to be used for setting this entry rate
        """

        assert type(label) is str, 'Destination compartment label not string for setting variable transfer rate'
        assert type(var_label) is str, 'Function label not string for setting variable entry rate'
        add_unique_tuple_to_list(self.var_entry_rate_flow, (label, var_label))

    def calculate_scaleup_vars(self):
        """
        Find the values of the scale-up functions at a specific point in time
        (to be called within the integration process).
        """

        for label, fn in self.scaleup_fns.iteritems(): 
            self.vars[label] = fn(self.time)

    def calculate_vars(self):
        """
        Calculate self.vars that only depend on compartment values
        (i.e. those not pre-specified, such as force of infection)
        """
        pass

    def calculate_flows(self):
        """
        Main method to calculate changes in all compartment sizes based on flow rate,
        which should only depend on compartment values and self.vars already calculated in self.calculate_vars.
        """

        # start from zero for all compartments
        for label in self.labels:
            self.flows[label] = 0.

        # birth flows
        for label, var_label in self.var_entry_rate_flow:
            self.flows[label] += self.vars[var_label]

        # dynamic transmission flows
        for from_label, to_label, var_label in self.var_transfer_rate_flows:
            val = self.compartments[from_label] * self.vars[var_label]
            self.flows[from_label] -= val
            self.flows[to_label] += val

        # fixed-rate flows
        for from_label, to_label, rate in self.fixed_transfer_rate_flows:
            val = self.compartments[from_label] * rate
            self.flows[from_label] -= val
            self.flows[to_label] += val

        # background death flows
        self.vars['rate_death'] = 0.
        for label in self.labels:
            val = self.compartments[label] * self.background_death_rate
            self.flows[label] -= val
            self.vars['rate_death'] += val

        # extra infection-related death flows
        self.vars['rate_infection_death'] = 0.
        for label, rate in self.infection_death_rate_flows:
            val = self.compartments[label] * rate
            self.flows[label] -= val
            self.vars['rate_infection_death'] += val

    def prepare_vars_and_flows(self):
        """
        Calculate all scale-up function, vars and flows during integration
        """

        # clear previously populated vars dictionary
        self.vars.clear()

        # calculate vars and flows sequentially
        self.calculate_scaleup_vars()
        self.calculate_vars()
        self.calculate_flows()

    def make_derivate_fn(self):
        """
        Create the derivative function for integration
        """

        def derivative_fn(y, t):
            self.time = t
            self.compartments = self.convert_list_to_compartments(y)
            self.prepare_vars_and_flows()
            flow_vector = self.convert_compartments_to_list(self.flows)
            self.checks()
            return flow_vector

        return derivative_fn

    def init_run(self):
        """
        Starting method to integration processes
        """

        self.set_flows()
        self.var_labels = None
        self.soln_array = None
        self.var_array = None
        self.flow_array = None

        # check that each compartment has at least one entry flow or exit flow
        # (i.e. that all compartments are connected to the model)
        labels_with_entry = [v[0] for v in self.var_entry_rate_flow]
        labels_with_var_out = [v[0] for v in self.var_transfer_rate_flows]
        labels_with_var_in = [v[1] for v in self.var_transfer_rate_flows]
        labels_with_fixed_out = [v[0] for v in self.fixed_transfer_rate_flows]
        labels_with_fixed_in = [v[1] for v in self.fixed_transfer_rate_flows]

        connected_compartments \
            = labels_with_entry + labels_with_var_out + labels_with_var_in \
              + labels_with_fixed_out + labels_with_fixed_in

        for label in self.labels:
            msg = 'Compartment "%s" doesn\'t have any entry or transfer flows' % label
            assert label in connected_compartments, msg

    def integrate(self, method='explicit'):
        """
        Master integration function to prepare for integration run and then call the process to actually run the
        integration process according to the process selected in the arguments to this method.

        Args:
            method: Either 'explicit' or 'scipy' to select the integration process required
        """

        self.init_run()
        assert self.times is not None, 'Haven\'t set times yet'
        derivative = self.make_derivate_fn()
        y = self.get_init_list()
        if method == 'explicit':
            self.integrate_explicit(y, derivative)
        elif method == 'scipy':
            self.soln_array = odeint(derivative, y, self.times)
        self.calculate_diagnostics()

    def integrate_explicit(self, y, derivative, min_dt=0.05):
        """
        Integrate with Euler explicit method

        Args:
            min_dt:  Minimum time change allowed
        """

        n_component = len(y)
        n_time = len(self.times)
        self.soln_array = numpy.zeros((n_time, n_component))

        time = self.times[0]
        self.soln_array[0, :] = y
        for i_time, new_time in enumerate(self.times):
            while time < new_time:
                f = derivative(y, time)
                old_time = time
                time = time + min_dt
                dt = min_dt
                if time > new_time:
                    dt = new_time - old_time
                    time = new_time
                for i in range(n_component):
                    y[i] = y[i] + dt * f[i]
                    # hack to avoid errors due to time-step
                    if y[i] < 0.:
                        y[i] = 0.
            if i_time < n_time - 1:
                self.soln_array[i_time + 1, :] = y

    def calculate_diagnostic_vars(self):
        """
        Calculate diagnostic vars that can depend on self.flows as well as self.vars calculated in calculate_vars
        """

        pass

    def calculate_diagnostics(self):
        """
        Run diagnostic functions to get outcomes at end of integration
        """

        # create dictionary of lists of compartment values
        self.population_soln = {}
        for label in self.labels:
            if label in self.population_soln:
                continue
            self.population_soln[label] = self.get_compartment_soln(label)

        n_time = len(self.times)
        for i in range(n_time):

            self.time = self.times[i]

            for label in self.labels:
                self.compartments[label] = self.population_soln[label][i]

            self.calculate_vars()
            self.calculate_flows()
            self.calculate_diagnostic_vars()

            # only set after self.calculate_diagnostic_vars has been run so that we have all var_labels,
            # including the ones in calculate_diagnostic_vars
            if self.var_labels is None:
                self.var_labels = self.vars.keys()
                self.var_array = numpy.zeros((n_time, len(self.var_labels)))
                self.flow_array = numpy.zeros((n_time, len(self.labels)))

            for i_label, label in enumerate(self.var_labels):
                self.var_array[i, i_label] = self.vars[label]
            for i_label, label in enumerate(self.labels):
                self.flow_array[i, i_label] = self.flows[label]

        # calculate compartment sizes as fractions of population
        self.fraction_soln = {}
        for label in self.labels:
            self.fraction_soln[label] = [
                v / t
                for v, t
                in zip(
                    self.population_soln[label],
                    self.get_var_soln('population')
                )
            ]

    def get_compartment_soln(self, label):
        """
        Find the values of a named compartment over the time steps of the integration

        Args:
            label: String for the name of the compartment of interest
        Returns:
            List of values for this compartment's values over time
        """

        assert self.soln_array is not None, 'calculate_diagnostics has not been run'
        i_label = self.labels.index(label)
        return self.soln_array[:, i_label]

    def get_var_soln(self, label):
        """
        Find the values of a var over the time steps of the integration

        Args:
            label: String for the name of the var of interest
        Returns:
            List of values for this var's values over time
        """

        assert self.var_array is not None, 'calculate_diagnostics has not been run'
        i_label = self.var_labels.index(label)
        return self.var_array[:, i_label]

    def get_flow_soln(self, label):
        """
        Find the values of a compartment's flows over the time steps of the integration

        Args:
            label: String for the name of the compartment of interest
        Returns:
            List of values for this compartment's flows over time
        """

        assert self.flow_array is not None, 'calculate_diagnostics has not been run'
        i_label = self.labels.index(label)
        return self.flow_array[:, i_label]

    def load_state(self, i_time):
        """
        Reload the compartmental state of a model from a previous time point
        (for running alternative scenarios if epidemiological changes have been made from a point in time)

        Args:
            i_time: The index of the list of times to load from
        """

        self.time = self.times[i_time]
        for i_label, label in enumerate(self.labels):
            self.compartments[label] = \
                self.soln_array[i_time, i_label]
        self.calculate_vars()

    def checks(self, error_margin=0.1):
        """
        Assertions run during the simulation, should be over-ridden for each model, but given as an example

        Args:
            error_margin: acceptable difference between target invariants
        """

        # Check all compartments are positive
        for label in self.labels:
            assert self.compartments[label] >= 0.

        # Check population is conserved across compartments
        # population_change = \
        #     self.vars['rate_birth'] \
        #     - self.vars['rate_death'] \
        #     - self.vars['rate_infection_death']
        # assert abs(sum(self.flows.values()) - population_change) < error_margin

    def make_graph(self, png):
        """
        Use external module (graphviz) to create a flow diagram of the compartmental structure of the model

        Args:
            png: String for the filename of the file for the diagram to be stored in
        """

        from graphviz import Digraph

        styles = {
            'graph': {
                'label': 'Dynamic Transmission Model',
                'fontsize': '16',
            },
            'nodes': {
                'fontname': 'Helvetica',
                'shape': 'box',
                'style': 'filled',
                'fillcolor': '#CCDDFF',
            },
            'edges': {
                'style': 'dotted',
                'arrowhead': 'open',
                'fontname': 'Courier',
                'fontsize': '10',
            }
        }

        def apply_styles(graph, styles):
            graph.graph_attr.update(
                ('graph' in styles and styles['graph']) or {}
            )
            graph.node_attr.update(
                ('nodes' in styles and styles['nodes']) or {}
            )
            graph.edge_attr.update(
                ('edges' in styles and styles['edges']) or {}
            )
            return graph

        def num_str(f):

            return '%.3g' % f

            # alternative code to uncomment if you prefer more natural language, rather than scientific notation
            # abs_f = abs(f)
            # if abs_f > 1e9:
            #     return '%.1f billion' % (f / 1e9)
            # if abs_f > 1e6:
            #     return '%.1f million' % (f / 1e6)
            # if abs_f > 1e3:
            #     return '%.1f thousand' % (f / 1e3)
            # if abs_f > 1e2:
            #     return '%.0f' % f
            # if abs_f > 0.5:
            #     return '%.1f' % f
            # if abs_f > 5e-2:
            #     return '%.2f' % f
            # if abs_f > 5e-4:
            #     return '%.4f' % f
            # if abs_f > 5e-6:
            #     return '%.6f' % f
            # return str(f)

        self.graph = Digraph(format='png')
        for label in self.labels:
            self.graph.node(label)
        if len(self.infection_death_rate_flows) > 0:
            self.graph.node('infection_death')
        for from_label, to_label, var_label in self.var_transfer_rate_flows:
            self.graph.edge(from_label, to_label, label=var_label)
        for from_label, to_label, rate in self.fixed_transfer_rate_flows:
            self.graph.edge(from_label, to_label, label=num_str(rate))
        if len(self.infection_death_rate_flows) > 0:
            for label, rate in self.infection_death_rate_flows:
                self.graph.edge(label, 'infection_death', label=num_str(rate))
        base, ext = os.path.splitext(png)
        if ext.lower() != '.png':
            base = png

        self.graph = apply_styles(self.graph, styles)

        self.graph.render(base)

    def check_converged_compartment_fraction(self, label, equil_time, test_fraction_diff):
        """
        Numerically determine whether a compartment's proportion has converged to reach an approximate equilibrium

        Args:
            label: Compartment to check
            equil_time: Point in time to check at
            test_fraction_diff: Difference in compartment value below which convergence can be considered
        Returns:
            Boolean for whether convergence has occurred
        """

        self.calculate_diagnostics()
        times = self.times
        fraction = self.fraction_soln[label]
        i = -2
        max_fraction_diff = 0
        time_diff = 0
        while time_diff < equil_time:
            i -= 1
            if -i >= len(times):
                break
            time_diff = abs(times[-1] - times[i])
            frac_diff = (fraction[-1] - fraction[i])
            if abs(frac_diff) > max_fraction_diff:
                max_fraction_diff = frac_diff
            if abs(frac_diff) > test_fraction_diff:
                return False
        return True

