#!/usr/bin/env python2

import sys
import os

root_path = os.path.join(os.path.dirname(os.path.realpath( __file__ )),os.path.pardir)
sys.path.insert(0, root_path)

import phase_space_generators as PS
import madgraph.core.base_objects as base_objects
import vectors as vectors
import madgraph.various.cluster as cluster
import madgraph.various.misc as misc
import models.import_ufo as import_ufo
import math

import integrands
import vegas3_integrator as vegas3


import matplotlib.pyplot as plt
import numpy as np

from madgraph import MG5DIR

pjoin = os.path.join

import copy
import math
import random

FINAL=True
INITIAL=False

DUMMY=99

def get_model():
    """Instantiate a model,
    which will be useful for the non-flat phase-space generator test.
    """
    
    model_with_params_set = import_ufo.import_model(
            pjoin(MG5DIR,'models','sm'), prefix=True,
            complex_mass_scheme = False )
    model_with_params_set.pass_particles_name_in_mg_default()
    model_with_params_set.set_parameters_and_couplings(
            param_card = pjoin(MG5DIR,'models','sm','restrict_default.dat'),
            complex_mass_scheme=False)
    return model_with_params_set


def get_topology(channel_name):

        # A specific sets of s- and t-channels for this example:

        ##############################################################################
        # a) e+(1) e-(2) > mu+(3) mu-(4) a(5) with (p3+p4)^2 sampled as a z-resonance
        ##############################################################################
        topology = None

        if channel_name=='2_to_3_s34':

                topology = (
                    # s-channels first:
                    base_objects.VertexList([
                        base_objects.Vertex({
                            'id': DUMMY, # Irrelevant
                            'legs': base_objects.LegList([
                                base_objects.Leg({
                                    'id': -12,
                                    'number': 3,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 12,
                                    'number': 4,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 23,
                                    'number': -1,
                                    'state': FINAL,
                                })
                            ])
                        }),
                    ]),
                    # t-channels then:
                    base_objects.VertexList([
                        base_objects.Vertex({
                            'id': DUMMY, # Irrelevant
                            'legs': base_objects.LegList([
                                base_objects.Leg({
                                    'id': 11,
                                    'number': 1,
                                    'state': INITIAL,
                                }),
                                base_objects.Leg({
                                    'id': 22,
                                    'number': 5,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 11,
                                    'number': -2,
                                    'state': FINAL,
                                })
                            ])
                        }),
                        # The dummy vertex below is just to close the diagram and connect
                        # with the number -3 which is to be understood here as the initial state #2.
                        base_objects.Vertex({
                            'id': DUMMY, # Irrelevant
                            'legs': base_objects.LegList([
                                base_objects.Leg({
                                    'id': 11,
                                    'number': -2,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 23,
                                    'number': -1,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': -11,
                                    'number': -3,
                                    'state': INITIAL,
                                })
                            ])
                        }),
                    ])
                )
            
        else:
            print("Channel name '%s' not implemented."%channel_name)

        return topology

class IntegrandForTest(integrands.VirtualIntegrand):
    """An integrand for this phase-space volume test."""
    def __init__(self, phase_space_generator):
        super(IntegrandForTest, self).__init__(phase_space_generator.get_dimensions())
        self.phase_space_generator = phase_space_generator
        self.counter               = 0
        #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
        #    self.my_random_path = self.phase_space_generator.generate_random_path()
    
    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
        #    PS_point, wgt, x1, x2 = self.phase_space_generator.get_PS_point(continuous_inputs,self.my_random_path)
        #else:
        self.counter += 1
        PS_point, wgt, x1, x2 = self.phase_space_generator.get_PS_point(continuous_inputs)
        return wgt

#=========================================================================================
# Standalone main for running the example
#=========================================================================================
if __name__ == '__main__':

    model = get_model()
    topology = get_topology('2_to_3_s34')
    
    # Modify below the mass of the resonance
    # It wil be centered here at 150 GeV
    model.get('parameter_dict')['mdl_MZ'] = 150.0
    # And will have a width of 60 GeV
    model.get('parameter_dict')['mdl_WZ'] = 60.0

    E_cm = 1000.0

    my_SCPS_generator = PS.SingleChannelPhasespace([0.]*2, [0.]*3,
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0), 
              model=model, topology=topology, path = [[0,],[]])
    random_variables = [random.random() for _ in range(my_SCPS_generator.nDimPhaseSpace())]
    print(random_variables)
    
    print "Considering the following topology:"
    print "-"*10
    print my_SCPS_generator.get_topology_string(my_SCPS_generator.topology, path_to_print=my_SCPS_generator.path)
    print "-"*10

    # Now generate a point
    momenta, wgt, _, _ = my_SCPS_generator.get_PS_point(random_variables)
    
    print "Kinematic configuration generated for random variables %s:"%str(random_variables)
    print momenta
    print "Jacobian weight: %.16f"%wgt

    print ''
    print '-'*10
    print ''

    my_integrand = IntegrandForTest(my_SCPS_generator)
    # Increase verbosity to 2 in order to see the evolution of the integration
    verbosity = 1
    runner = cluster.MultiCore(8)
    # Use the line below instead to run on one core only
    #runner = cluster.onecore
    my_integrator = vegas3.Vegas3Integrator(my_integrand, 
            n_points_survey=400, n_points_refine=400, accuracy_target=None,
            verbosity=verbosity, cluster=runner
    )
    # Finally integrate
    print "Now integrating the phase-space volume with the Single Channel Phase-Space parametrisation"
    result = my_integrator.integrate()
    print 'Phase-space volume using single channel phase-space: %.4e +/- %.2e'%result

    my_flat_PS_generator = PS.FlatInvertiblePhasespace([0.]*2, [0.]*3, beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0))
    my_integrand = IntegrandForTest(my_flat_PS_generator)
    my_integrator = vegas3.Vegas3Integrator(my_integrand, 
            n_points_survey=400, n_points_refine=400, accuracy_target=None,
            verbosity=verbosity, cluster=runner
    )
    print "Now integrating the phase-space volume with the flat phase-space generator"
    result = my_integrator.integrate()

    print 'Phase-space volume using flat phase-space parametrisation: %.4e +/- %.2e'%result
    print ''
    print '>>>>> QUESTION 1: Explain why integrating the PS volume witht he SCPS parametrisation is less efficient.'

    print ''
    print '-'*100
    print ''
    print 'Now plot the distributions of the invariant mass (p3+p4)^2'
    n_samples = 10000
    flat_data = []
    scps_data = []
    flat_data_weighted = []
    scps_data_weighted = []
    for i in range(n_samples):
        random_variables = [random.random() for _ in range(my_SCPS_generator.nDimPhaseSpace())]
        flat_momenta, flat_wgt, _, _ = my_flat_PS_generator.get_PS_point(random_variables)
        scps_momenta, scps_wgt, _, _ = my_SCPS_generator.get_PS_point(random_variables)
        flat_momenta_s34 = math.sqrt((flat_momenta[2]+flat_momenta[3]).square())
        scps_momenta_s34 = math.sqrt((scps_momenta[2]+scps_momenta[3]).square())
        flat_data.append((flat_momenta_s34, 1.0))
        scps_data.append((scps_momenta_s34, 1.0))
        flat_data_weighted.append((flat_momenta_s34, flat_wgt))
        scps_data_weighted.append((scps_momenta_s34, scps_wgt))
         
    fig, axs = plt.subplots(2, 2)
    values, weights = zip(*flat_data)    
    axs[0,0].hist(values, weights=weights, density=1.0, bins=100)
    axs[0,0].title.set_text('Flat dw/ds34')
    values, weights = zip(*flat_data_weighted)    
    axs[0,1].hist(values, weights=weights, density=1.0, bins=100)
    axs[0,1].title.set_text('Flat dw/ds34 weighted')
    values, weights = zip(*scps_data)    
    axs[1,0].hist(values, weights=weights, density=1.0, bins=100)
    axs[1,0].title.set_text('SCPS dw/ds34')
    values, weights = zip(*scps_data_weighted)    
    axs[1,1].hist(values, weights=weights, density=1.0, bins=100)
    axs[1,1].title.set_text('SCPS dw/ds34 weighted')

    print ''
    print '>>>>> QUESTION 2: Explain what you see in the four plots, and in particular what is the difference for the weighted histogram.'
    print ''
    plt.show()

