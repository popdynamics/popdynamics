# -*- coding: utf-8 -*-


"""

Base Population Model to handle different type of models.

Implicit time unit: years

"""

import os
from scipy.integrate import odeint
import numpy
from math import exp


def add_unique_tuple_to_list(a_list, a_tuple):
    """
    Adds or modifies a list of tuples, compares only the items
    before the last in the tuples, the last value in the tuple
    is assumed to be a value.
    """
    for i, test_tuple in enumerate(a_list):
        if test_tuple[:-1] == a_tuple[:-1]:
            a_list[i] = a_tuple
            break
    else:
        a_list.append(a_tuple)


def label_intersects_tags(label, tags):
    for tag in tags:
        if tag in label:
            return True
    return False


class BaseModel():
    """
    Basic concept
      - var - values that are calculated at every time step
      - param - values that are set at the beginning 

    """
    def __init__(self):

        # labels for all compartments
        self.labels = []

        # stores the initial value for all compartments
        self.init_compartments = {}

        # stores the values of all parameters, there 
        # should be no hard-coded values except as contained
        # in this structure
        self.params = {}

        # stores list of time points
        self.times = None

        self.scaleup_fns = {}

        # stores any auxillary variables used to calculate
        # dynamic transmission at every time-step
        self.vars = {}

        # total flow of each compartment
        self.flows = {}

        # list of 2-tuple (label, var_label)
        #   - label: name of compartment
        #   - var_label: name of var that holds the entry rate
        self.var_entry_rate_flow = []

        # list of 3-tuple (from_label, to_label, param_label)
        #   - from_label: name of compartment that loses population
        #   - to_label: name of compartment that gains population
        #   - param_label: name of param that holds the rate
        self.fixed_transfer_rate_flows = []

        # list of 3-tuple (from_label, to_label, var_label)
        #   - from_label: name of compartment that loses population
        #   - to_label: name of compartment that gains population
        #   - var_label: name of var that holds the rate
        self.var_transfer_rate_flows = []

        # list of 2-tuple (label, rate)
        #   - label: name of compartment
        #   - rate: disease specific death rate
        self.infection_death_rate_flows = []

        # the generalized death rate of all compartments
        self.background_death_rate = 0.0

        self.soln_array = None
        self.var_labels = None
        self.var_array = None
        self.flow_array = None

    def make_times(self, start, end, delta):
        "Return steps with n or delta"
        self.times = []
        step = start
        while step <= end:
            self.times.append(step)
            step += delta

    def make_times_with_n_step(self, start, end, n):
        "Return steps with n or delta"
        self.times = []
        step = start
        delta = (end - start) / float(n)
        while step <= end:
            self.times.append(step)
            step += delta

    def set_compartment(self, label, init_val=0.0):
        if label not in self.labels:
            self.labels.append(label)
        self.init_compartments[label] = init_val
        assert init_val >= 0, 'Start with negative compartment not permitted'

    def set_param(self, label, val):
        self.params[label] = val

    def convert_list_to_compartments(self, vec):
        return {l: vec[i] for i, l in enumerate(self.labels)}

    def convert_compartments_to_list(self, compartments):
        return [compartments[l] for l in self.labels]

    def get_init_list(self):
        return self.convert_compartments_to_list(self.init_compartments)

    def set_scaleup_fn(self, label, fn):

        """
        Simple method to add a scale-up function to the dictionary of scale-ups.

        Args:
            label: String for name of function.
            fn: The function to be added.
        """

        self.scaleup_fns[label] = fn

    def set_background_death_rate(self, param_label):
        self.background_death_rate = self.params[param_label]

    def set_infection_death_rate_flow(self, label, param_label):
        add_unique_tuple_to_list(
            self.infection_death_rate_flows,
            (label, self.params[param_label]))

    def set_fixed_transfer_rate_flow(self, from_label, to_label, param_label):
        add_unique_tuple_to_list(
            self.fixed_transfer_rate_flows,
            (from_label, to_label, self.params[param_label]))

    def set_var_transfer_rate_flow(self, from_label, to_label, var_label):
        add_unique_tuple_to_list(
            self.var_transfer_rate_flows,
            (from_label, to_label, var_label))

    def set_var_entry_rate_flow(self, label, var_label):
        add_unique_tuple_to_list(
            self.var_entry_rate_flow,
            (label, var_label))

    def calculate_scaleup_vars(self):

        """
        Find the values of the scale-up functions at a specific point in time. Called within the integration process.
        """

        for label, fn in self.scaleup_fns.iteritems(): 
            self.vars[label] = fn(self.time)

    def calculate_vars(self):
        """
        Calculate self.vars that only depend on compartment values
        """
        pass

    def calculate_flows(self):
        """
        Calculate flows, which should only depend on compartment values
        and self.vars calculated in self.calculate_vars.
        """
        for label in self.labels:
            self.flows[label] = 0.0

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

        # normal death flows
        self.vars["rate_death"] = 0.0
        for label in self.labels:
            val = self.compartments[label] * self.background_death_rate
            self.flows[label] -= val
            self.vars['rate_death'] += val

        # extra death flows
        self.vars["rate_infection_death"] = 0.0
        for label, rate in self.infection_death_rate_flows:
            val = self.compartments[label] * rate
            self.flows[label] -= val
            self.vars["rate_infection_death"] += val


    def prepare_vars_and_flows(self):
        """
        This function collects some other functions that previously 
        led to a bug because not all of them were called
        in the diagnostics round.
        """

        # Clear previously populated vars dictionary
        self.vars.clear()

        # Calculate vars and flows sequentially
        self.calculate_scaleup_vars()
        self.calculate_vars()
        self.calculate_flows()


    def make_derivate_fn(self):

        def derivative_fn(y, t):
            self.time = t
            self.compartments = self.convert_list_to_compartments(y)
            self.prepare_vars_and_flows()
            flow_vector = self.convert_compartments_to_list(self.flows)
            self.checks()
            return flow_vector

        return derivative_fn

    def init_run(self):
        self.set_flows()
        self.var_labels = None
        self.soln_array = None
        self.var_array = None
        self.flow_array = None

        # check that each compartment has an entry flow
        labels_with_entry = [v[0] for v in self.var_entry_rate_flow]
        labels_with_var = [v[1] for v in self.var_transfer_rate_flows]
        labels_with_fixed = [v[1] for v in self.fixed_transfer_rate_flows]
        labels = labels_with_entry + labels_with_var + labels_with_fixed
        for label in self.labels:
            msg = "Compartment '%s' doesn't have any entry or transfer flows" % label
            assert label in labels, msg

    def integrate_scipy(self):
        self.init_run()
        assert not self.times is None, "Haven't set times yet"
        init_y = self.get_init_list()
        derivative = self.make_derivate_fn()
        self.soln_array = odeint(derivative, init_y, self.times)

        self.calculate_diagnostics()

    def integrate_explicit(self, min_dt=0.05):
        self.init_run()
        assert not self.times is None, "Haven't set times yet"
        y = self.get_init_list()
        n_component = len(y)
        n_time = len(self.times)
        self.soln_array = numpy.zeros((n_time, n_component))

        derivative = self.make_derivate_fn()
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
                    if y[i] < 0.0:
                        y[i] = 0.0
            if i_time < n_time - 1:
                self.soln_array[i_time + 1, :] = y

        self.calculate_diagnostics()

    def calculate_diagnostic_vars(self):
        """
        Calculate diagnostic vars that can depend on self.flows as
        well as self.vars calculated in calculate_vars
        """
        pass

    def calculate_diagnostics(self):
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

            # only set after self.calculate_diagnostic_vars is
            # run so that we have all var_labels, including
            # the ones in calculate_diagnostic_vars
            if self.var_labels is None:
                self.var_labels = self.vars.keys()
                self.var_array = numpy.zeros((n_time, len(self.var_labels)))
                self.flow_array = numpy.zeros((n_time, len(self.labels)))

            for i_label, label in enumerate(self.var_labels):
                self.var_array[i, i_label] = self.vars[label]
            for i_label, label in enumerate(self.labels):
                self.flow_array[i, i_label] = self.flows[label]

        self.fraction_soln = {}
        for label in self.labels:
            self.fraction_soln[label] = [
                v / t
                for v, t
                in zip(
                    self.population_soln[label],
                    self.get_var_soln("population")
                )
            ]

    def get_compartment_soln(self, label):
        assert self.soln_array is not None, "calculate_diagnostics has not been run"
        i_label = self.labels.index(label)
        return self.soln_array[:, i_label]

    def get_var_soln(self, label):
        assert self.var_array is not None, "calculate_diagnostics has not been run"
        i_label = self.var_labels.index(label)
        return self.var_array[:, i_label]

    def get_flow_soln(self, label):
        assert self.flow_array is not None, "calculate_diagnostics has not been run"
        i_label = self.labels.index(label)
        return self.flow_array[:, i_label]

    def load_state(self, i_time):
        self.time = self.times[i_time]
        for i_label, label in enumerate(self.labels):
            self.compartments[label] = \
                self.soln_array[i_time, i_label]
        self.calculate_vars()

    def checks(self, error_margin=0.1):
        """
        Assertion run during the simulation, should be overriden
        for each model.

        Args:
            error_margin: acceptable difference between target invariants

        Returns:

        """
        # # Check all compartments are positive
        # for label in self.labels:
        #     assert self.compartments[label] >= 0.0
        # Check population is conserved across compartments
        population_change = \
            self.vars['rate_birth'] \
            - self.vars['rate_death'] \
            - self.vars['rate_infection_death']
        assert abs(sum(self.flows.values()) - population_change) < error_margin

    def make_graph(self, png):
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
            abs_f = abs(f)
            if abs_f > 1E9:
                return "%.1fB" % (f / 1E9)
            if abs_f > 1E6:
                return "%.1fM" % (f / 1E6)
            if abs_f > 1E3:
                return "%.1fK" % (f / 1E3)
            if abs_f > 100:
                return "%.0f" % f
            if abs_f > 0.5:
                return "%.1f" % f
            if abs_f > 0.05:
                return "%.2f" % f
            if abs_f > 0.0005:
                return "%.4f" % f
            if abs_f > 0.000005:
                return "%.6f" % f
            return str(f)

        self.graph = Digraph(format='png')
        for label in self.labels:
            self.graph.node(label)
        self.graph.node("infection_death")
        for from_label, to_label, var_label in self.var_transfer_rate_flows:
            self.graph.edge(from_label, to_label, label=var_label)
        for from_label, to_label, rate in self.fixed_transfer_rate_flows:
            self.graph.edge(from_label, to_label, label=num_str(rate))
        for label, rate in self.infection_death_rate_flows:
            self.graph.edge(label, "infection_death", label=num_str(rate))
        base, ext = os.path.splitext(png)
        if ext.lower() != '.png':
            base = png

        self.graph = apply_styles(self.graph, styles)

        self.graph.render(base)

    def check_converged_compartment_fraction(
            self, label, equil_time, test_fraction_diff):
        labels = self.labels
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


def make_sigmoidal_curve(y_low=0, y_high=1.0, x_start=0, x_inflect=0.5, multiplier=1.):
    """
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
        arg = b * ( x_inflect - x )
        # check for large values that will blow out exp
        if arg > 10.0:
            return y_low
        return amplitude / ( 1. + exp( arg ) ) + y_low

    return curve
