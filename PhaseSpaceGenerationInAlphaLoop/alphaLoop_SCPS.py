#!/usr/bin/env python3

import sys
import os
from multiprocessing import Value 
import time
from pprint import pformat

WORKSPACE = os.path.abspath(os.path.dirname(os.path.realpath( __file__ )))
MGPATH="/Users/valentin/Documents/MG5/3.0.2.py3"
if not os.path.exists(MGPATH):
    print("ERROR: Please adjust the MGPATH variable in the file alphaLoop_SCPS.py.")
root_path = MGPATH
sys.path.insert(0, root_path)

import progressbar
import yaml
import shutil
import random
from scipy import optimize
from argparse import ArgumentParser
import phase_space_generators as PS
import madgraph.core.base_objects as base_objects
import vectors as vectors
import madgraph.various.cluster as cluster
import madgraph.various.misc as misc
import models.import_ufo as import_ufo
import re
import math
from vectors import LorentzVector, LorentzVectorList

import integrands
import vegas3_integrator as vegas3


import matplotlib.pyplot as plt
import numpy as np

import logging
#logging.basicConfig(format='%(asctime)s %(message)s')

from madgraph import MG5DIR

pjoin = os.path.join

AL_PATH=pjoin(pjoin(root_path,'PLUGIN','alphaloop'))
sys.path.insert(0, AL_PATH)

import copy
import math
import random

FINAL=True
INITIAL=False

DUMMY=99

CORES=8

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
            print(("Channel name '%s' not implemented."%channel_name))

        return topology


def set_hyperparameters(process_name, sg_name, workspace=None, min_jets=0, multi_settings={}):
    '''Create the correct hyperparameter file'''

    path = pjoin(AL_PATH, 'LTD', 'hyperparameters.yaml')
    hyperparameters = yaml.load(open(path, 'r'), Loader=yaml.Loader)

    try:
        min_jpt = multi_settings['Selectors']['jet']['min_jpt']
    except KeyError:
        min_jpt = hyperparameters['Selectors']['jet']['min_jpt']

    if min_jets == 0:
        min_jpt = 'nocut'
    # Set some default values for the hyperparameter
    # Set values for General
    hyperparameters['General']['topology'] = "%s_%s" % (sg_name, min_jpt)
    hyperparameters['General']['res_file_prefix'] = "%s_" % process_name
    hyperparameters['General']['log_file_prefix'] =\
        "stats/%s_%s_%s_" % (process_name, sg_name, min_jpt)
    hyperparameters['General']['partial_fractioning_multiloop'] = True
    hyperparameters['General']['partial_fractioning_threshold'] = -1
    # Set values for Integrator
    hyperparameters['Integrator']['dashboard'] = False
    hyperparameters['Integrator']['integrator'] = 'vegas'
    hyperparameters['Integrator']['internal_parallelization'] = True
    hyperparameters['Integrator']['n_vec'] = int(1e4)
    # Set values for CrossSection
    hyperparameters['CrossSection']['numerator_source'] = 'FORM_integrand'
    hyperparameters['CrossSection']['picobarns'] = True

    # General settings can be parsed though the multi_settings valiable
    # It must have the same structure as hyperparameters
    # It's not necessary that it contains all the elements
    for name, settings in multi_settings.items():
        for setting, value in settings.items():
            if name not in hyperparameters or setting not in hyperparameters[name]:
                print("Unknown setting hyperparameters[{}][{}] does not exits!".format(
                    name, setting))
                sys.exit(1)
            if isinstance(hyperparameters[name][setting], dict):
                hyperparameters[name][setting].update(value)
            else:
                hyperparameters[name][setting] = value

    # Create the hyperpamater file
    output_path = pjoin(workspace, 'hyperparameters',
                        '%s_%s_%s.yaml' % (process_name, sg_name, min_jpt))
    with open(output_path, "w") as stream:
        yaml.dump(hyperparameters, stream)  # default_flow_style=None)
    return output_path


def get_rust_worker(
        proc_path,
        diag_yaml_path,
        hyperparameters_path):

    try:
        # Import the rust bindings
        from ltd import CrossSection
    except ImportError:
        print("ERROR: Could not import the rust back-end 'ltd' module. Compile it first with:"
            " ./make_lib from within the pyNLoop directory.")
        raise
    
    os.environ['MG_NUMERATOR_PATH'] = proc_path if proc_path.endswith('/') else '%s/'%proc_path

    if not os.path.isdir(pjoin(WORKSPACE,'stats')):
        os.makedirs(pjoin(WORKSPACE,'stats'))

    return CrossSection( diag_yaml_path, hyperparameters_path )

class TestHFuncIntegrand(integrands.VirtualIntegrand):
    """An integrand for this phase-space volume test."""

    def __init__(self, dimensions, h_function, debug=0):
    
        super(TestHFuncIntegrand, self).__init__( dimensions )

        self.h_function = h_function
        self.debug = debug

        self.n_evals = Value('i', 0)

        #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
        #    self.my_random_path = self.phase_space_generator.generate_random_path()
    
    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
        #    PS_point, wgt, x1, x2 = self.phase_space_generator.get_PS_point(continuous_inputs,self.my_random_path)
        #else:
        self.n_evals.value += 1

        x = list(continuous_inputs)[0]
        t, wgt = self.h_function(x)


        final_res = (1./(1+t**2))*(1./(math.pi/2.)) * wgt

        if self.debug: print("Final wgt returned to integrator: ",final_res)
        if self.debug > 2: time.sleep(self.debug)

        return final_res

class DefaultALIntegrand(integrands.VirtualIntegrand):
    """An integrand for this phase-space volume test."""

    def __init__(self, rust_worker, generator, debug=0, phase='real'):

        self.generator = generator
        self.dimensions = generator.dimensions
        self.debug = debug
    
        super(DefaultALIntegrand, self).__init__( self.dimensions )

        self.rust_worker = rust_worker
        self.phase = phase
        self.n_evals = Value('i', 0)

        #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
        #    self.my_random_path = self.phase_space_generator.generate_random_path()
    
        self.first_final_res = None

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        #if type(self.phase_space_generator).__name__ == 'SingleChannelPhasespace':
        #    PS_point, wgt, x1, x2 = self.phase_space_generator.get_PS_point(continuous_inputs,self.my_random_path)
        #else:
        self.n_evals.value += 1

        xs, wgt = self.generator(continuous_inputs)
        re, im = self.rust_worker.evaluate_integrand( xs )

        if self.phase=='real':
            final_res = re*wgt
        else:
            final_res = im*wgt
        
        if self.debug: print("Final wgt returned to integrator: %.16e"%final_res)
        if hasattr(self.generator,"test_inverse_jacobian") and self.generator.test_inverse_jacobian:
            if final_res == 0.0:
                self.generator.i_count = 0
            else:
                if self.first_final_res is None:
                    self.first_final_res = final_res
                else:
                    if self.debug: print("Diff w.r.t previous final wgt   : %.16e"%(self.first_final_res-final_res))
        if self.debug > 2: time.sleep(self.debug/10.0)

        return final_res

class HFunction(object):
    """
    Implements sampling from the noramlised distibution:

            PDF(t) = frac { t^\sigma } { 1 + t^2\sigma } frac { 2 \sigma } { \pi } Cos[ frac {\pi} {2 \sigma} ]
 
    The resulting CDF is expressed for generic \sigma in terms of a 2_F_1, which is not very practical, so we 
    implement it here only for \sigma=2 and \sigma=3
    """

    max_iteration = 500

    def __init__(self, sigma, debug=0):

        self.debug = debug

        self.sigma = sigma

        self.normalisation = ( (2*self.sigma) / math.pi ) * math.cos( math.pi / (2 * self.sigma) )

        self.PDF = (
            lambda t : (t**self.sigma / (1 + t**(2*self.sigma))) * self.normalisation
        )

        if self.sigma == 2:
            self.CDF = lambda t, x : (
                (1 / (2*math.pi))*(
                      math.log( 1 - ((2.*math.sqrt(2)*t)/(t*(t+math.sqrt(2))+1)) ) 
                    - 2.*math.atan(1-math.sqrt(2)*t) 
                    + 2.*math.atan(math.sqrt(2)*t+1)
                ) - x
            )
            self.CDF_prime = lambda t,x : (
                (2.*math.sqrt(2)*(t**2)) / (math.pi*((t**4)+1))
            )
            self.CDF_double_prime = lambda t,x : (
                - ( ( 4. * math.sqrt(2) * t * ( t**4 - 1 ) ) / ( math.pi * (t**4 + 1)**2 ) )
            )
        elif self.sigma == 3:
            self.CDF = lambda t,x : (
                (1 / (4*math.pi))*(
                    - 2*math.sqrt(3)*math.log(t**2+1)
                    + math.sqrt(3)*math.log(t**4-t**2+1)
                    - 6*math.atan(math.sqrt(3)-2*t)
                    - 6*math.atan(2*t+math.sqrt(3))
                    + 4*math.pi
                ) - x
            )
            self.CDF_prime = lambda t,x : (
                (3.*math.sqrt(3)*(t**3)) / (math.pi*(((t**6)+1)**2))
            )
            self.CDF_double_prime = lambda t,x : (
                - ( ( 9. * math.sqrt(3) * (t**2) * ( t**6 - 1 ) ) / ( math.pi * (t**6 + 1)**2 ) )
            )
        else:
            raise BaseException("Currently the H function is only implemented for sigma in [2,3].")

    def __call__(self, x):

        # Find the inverse point of the CDF
        bracket = [0.,1000000.0]
        if self.CDF(bracket[0],x) > 0.:
            print("ERROR: For x=%.16f, lower value %.16f for Netwons method invalid as it yields a positive value %.16f"%(x,bracket[0],self.CDF(bracket[0],x)))
            raise Exception("Incorrect seeds for Newtons method.")
        
        count=0
        while self.CDF(bracket[1],x) < 0. and count<10:
            bracket[1] *= 100
            count += 1        
        if self.CDF(bracket[1],x) < 0.:
            print("ERROR: For x=%.16f, upper value %.16f for Netwons method invalid as it yields a negative value %.16f"%(x,bracket[1],self.CDF(bracket[1],x)))
            raise Exception("Incorrect seeds for Newtons method.")

        root_result = optimize.root_scalar(self.CDF,
#            method="newton",
#            x0=0.01,x1=100.0,
            bracket=bracket,
            args=tuple([x,]),
            fprime=self.CDF_prime,
            fprime2=self.CDF_double_prime,
            maxiter=self.max_iteration)

    
        if self.debug or not root_result.converged:
            print("--------------------------------------------------------------")
            print("Numerical solver called with x=%.16e"%x)
            print("CDF(bracket[0]),CDF(bracket[1])=%.16e,%.16e"%(self.CDF(bracket[0],x),self.CDF(bracket[1],x)))
            print("CDF_prime(bracket[0]),CDF_prime(bracket[1])=%.16e,%.16e"%(self.CDF_prime(bracket[0],x),self.CDF_prime(bracket[1],x)))
            print("CDF_double_prime(bracket[0]),CDF_double_prime(bracket[1])=%.16e,%.16e"%(self.CDF_double_prime(bracket[0],x),self.CDF_double_prime(bracket[1],x)))
            print("Numerical solution found root=%.16e"%root_result.root)
            print("Number of iterations required=%d"%root_result.iterations)
            print("Number of function calls required=%d"%root_result.function_calls)
            print("Did converge? = %s"%root_result.converged)
            print("Algorithm termination cause: %s"%root_result.flag)
            print("PDF of for that t solution: %.16e"%self.PDF(root_result.root))
            print("--------------------------------------------------------------")
        
        return root_result.root, 1./self.PDF(root_result.root)

class CustomGenerator(object):
    pass

class generator_aL(CustomGenerator):

    def __init__(self, dimensions, rust_worker, SG_info, model, h_function, hyperparameters, debug=0, **opts):

        self.rust_worker = rust_worker
        self.model = model
        self.SG_info = SG_info
        self.hyperparameters = hyperparameters
        self.dimensions = dimensions
        self.debug = debug

    def __call__(self, random_variables, **opts):
        """ 
        Generates a point using a sampling reflecting the topolgy of SG_QG0.
        It will return the point directly in x-space as well as the final weight to combine
        the result with.
        """
        return list(random_variables), 1.0

class generator_spherical(CustomGenerator):

    def __init__(self, dimensions, rust_worker, SG_info, model, h_function, hyperparameters, debug=0, **opts):

        self.debug = debug
        self.rust_worker = rust_worker
        self.model = model
        self.SG_info = SG_info
        self.hyperparameters = hyperparameters
        self.dimensions = dimensions
        self.h_function = h_function

        E_cm = math.sqrt(max(sum(LorentzVector(vec) for vec in hyperparameters['CrossSection']['incoming_momenta']).square(),0.0))
        self.E_cm = E_cm
        if debug: print("E_cm detected=%.16e"%self.E_cm)

        # This generator is actually only used for its 'self.generator.dim_name_to_position' attribute
        self.generator = PS.FlatInvertiblePhasespace(
            [0.]*2, [0.]*(SG_info['topo']['n_loops']+1), 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
            dimensions = self.dimensions,
        )

    def __call__(self, random_variables):
        """ 
        Generates a point using a sampling reflecting the topolgy of SG_QG0.
        It will return the point directly in x-space as well as the final weight to combine
        the result with.
        """

#        self.debug=True
        x_t = random_variables[self.generator.dim_name_to_position['t']]
        if self.debug: print('x_t=',x_t)

        # We must now promote the point to full 3^L using the causal flow.
        # First transform t in [0,1] to t in [0,\inf]
#        rescaling_t, wgt_t = self.h_function(x_t)
        rescaling_t = self.E_cm * ( 1. / ( 1. - x_t ) - 1. / x_t )
        # Not clear where the factor 1./2. come from.
        wgt_t = (1./2.) * self.E_cm * ( 1. / x_t**2 + 1. / ( ( 1.-x_t ) * ( 1.-x_t ) ) )

        if self.debug: print('t, wgt_t=',rescaling_t, wgt_t)


        theta = random_variables[self.generator.dim_name_to_position['x_1']]*math.pi
        phi = random_variables[self.generator.dim_name_to_position['x_2']]*2*math.pi

        rescaled_PS_point_LMB = [ [
            rescaling_t*math.sin(theta)*math.cos(phi) , 
            rescaling_t*math.sin(theta)*math.sin(phi) , 
            rescaling_t*math.cos(theta)
        ], ]

        wgt_param = rescaling_t**2 * math.sin(theta) * (2*math.pi) * math.pi

        if self.debug: print("before: %s"%rescaled_PS_point_LMB[0]) 

        # Now we inverse parameterise this point using alphaLoop (non-multichanneled internal param)
        aL_xs = []
        inv_aL_jac = 1.0
        for i_v, v in enumerate(rescaled_PS_point_LMB):
            kx, ky, kz, inv_jac = self.rust_worker.inv_parameterize(v, i_v, self.E_cm**2)
            inv_aL_jac *= inv_jac
            aL_xs.extend([kx, ky, kz])
        
#        if self.debug: print("inv parameterise: %s, %.16e"%(str([kx, ky, kz]), inv_aL_jac**-1))
#        kkx, kky, kkz, kjac = self.rust_worker.parameterize([kx, ky, kz], 0, self.E_cm**2)
#        if self.debug: print("parameterise: %s, %.16e"%(str([kkx, kky, kkz]), kjac))

        # The final jacobian must then be our param. jac together with that of t divided by the one from alphaloop.
        final_jacobian = wgt_param * wgt_t * inv_aL_jac

        if self.debug: print('xs=',aL_xs)
        if self.debug: print('jacobian=',final_jacobian)

        return aL_xs, final_jacobian


class generator_flat(CustomGenerator):

    def __init__(self, dimensions, rust_worker, SG_info, model, h_function, hyperparameters, debug=0, **opts):

        self.debug = debug
        self.rust_worker = rust_worker
        self.model = model
        self.SG_info = SG_info
        self.hyperparameters = hyperparameters
        self.dimensions = dimensions
        self.h_function = h_function

        E_cm = math.sqrt(max(sum(LorentzVector(vec) for vec in hyperparameters['CrossSection']['incoming_momenta']).square(),0.0))
        self.E_cm = E_cm
        if debug: print("E_cm detected=%.16e"%self.E_cm)

        self.generator = PS.FlatInvertiblePhasespace(
            [0.]*2, [0.]*(SG_info['topo']['n_loops']+1), 
            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
            dimensions = self.dimensions,
        )

        # Variables for testing that the inverse jacobian works well.
        # Sett test_inverse_jacobian to True in order to debug the inverse jacobian.
        self.test_inverse_jacobian = False
        if self.test_inverse_jacobian:
            print("WARNING: Running in debug mode 'test_inverse_jacobian'.")
            time.sleep(1.0)
        self.rv = None
        self.i_count = 0

    def __call__(self, random_variables):
        """ 
        Generates a point using a sampling reflecting the topolgy of SG_QG0.
        It will return the point directly in x-space as well as the final weight to combine
        the result with.
        """
        #random_variables=[4.27617846e-01, 5.55958133e-01, 1.57593910e-01, 1.13340867e-01, 2.74746432e-01, 2.16284116e-01, 1.99368173e-04, 8.08726361e-02, 1.40842072e-01]

        if self.test_inverse_jacobian:
            if self.i_count==0:
                self.rv = copy.deepcopy(random_variables)
            else:
                random_variables = self.rv
            self.i_count += 1

        x_t = random_variables[self.generator.dim_name_to_position['t']]

        if self.debug: print('x_t=',x_t)

        # We must now promote the point to full 3^L using the causal flow.
        # First transform t in [0,1] to t in [0,\inf]
        rescaling_t, wgt_t = self.h_function(x_t)

        if self.test_inverse_jacobian:
            rescaling_t /= float(self.i_count)

        if self.debug: print('t, wgt_t=',rescaling_t, wgt_t)

        PS_point, wgt_param, x_1, x_2 = self.generator.get_PS_point(random_variables)
        if self.debug: print("orig wgt_param:",wgt_param)

        # Correct for the 1/(2E) of each cut propagator that alphaLoop includes but which is already included in the normal PS parameterisation
        for v in PS_point[2:]:
            wgt_param *= 2*v[0]
        if self.debug: print("after wgt_param:",wgt_param)

        if self.debug: print("PS point from param:\n", LorentzVectorList(PS_point))

        # Now rotate this point to the lMB of the topology.
        # In this case, since we generate flat, we don't care and simply take the first :-1 vectors
        # (and remove the initial external states)
        # and strip off the energies

#       START for ./run_ddx_NLO_SG_QG2.sh        
#        PS_point[2]=-PS_point[2]
#       END for ./run_ddx_NLO_SG_QG2.sh        

#       START for run_jjj_NLO_SG_QG36.sh
        PS_point[2]=-PS_point[2]
#       END for run_jjj_NLO_SG_QG36.sh

        PS_point_LMB = [ list(v[1:]) for v in PS_point[2:-1] ]
        if self.debug: print("PS point in LMB:\n", pformat(PS_point_LMB))

        # Apply the rescaling
        rescaled_PS_point_LMB = [ [ rescaling_t*ki for ki in k ] for k in PS_point_LMB]

        # Exclude normalising func when testing the inverse jacobian.
        if not self.test_inverse_jacobian:
            normalising_func = self.h_function.PDF(rescaling_t)
        else:
            normalising_func = 1.0
        if self.debug: print('normalising_func=',normalising_func)

        # Now we inverse parameterise this point using alphaLoop (non-multichanneled internal param)
        aL_xs = []
        # The rescaling jacobian comes both the rescaling change of variable and solving the delta function
        # TODO adjust so as to align wiith the particular LMB chosen for this topology.
        # TODO massive case not supported, issue a crash and warning in this case.
        dependent_momentum = -sum(vectors.Vector(v) for v in PS_point_LMB)

#       START for ./run_ddx_NLO_SG_QG2.sh
#        dependent_momentum = vectors.Vector(PS_point_LMB[0])-vectors.Vector(PS_point_LMB[1])
#       END for ./run_ddx_NLO_SG_QG2.sh        

#       START for run_jjj_NLO_SG_QG36.sh
        dependent_momentum = -vectors.Vector(PS_point_LMB[0])+vectors.Vector(PS_point_LMB[1])+vectors.Vector(PS_point_LMB[2])
#       END for run_jjj_NLO_SG_QG36.sh

        k_norms = [ math.sqrt(sum([ ki**2 for ki in k ])) for k in (PS_point_LMB+[list(dependent_momentum),]) ]

        if self.debug: print('h(t)=',self.h_function.PDF(rescaling_t))
        if self.debug: print('1/t=',1./rescaling_t)
        if self.debug: print('h(1/t)=',self.h_function.PDF(1./rescaling_t))
        if self.debug: print('rescaling_t**(3*len(PS_point_LMB))=',rescaling_t**(3*len(PS_point_LMB)))
        if self.debug: print('sum(k_norms)=',sum(k_norms))
        if self.debug: print('1./(rescaling_t*sum(k_norms))=',1./(rescaling_t*sum(k_norms)))

        inv_aL_jac = ( 1.0 / self.h_function.PDF(1./rescaling_t) ) * (rescaling_t**(len(PS_point_LMB)*3)) * ( rescaling_t * sum(k_norms) )

        for i_v, v in enumerate(rescaled_PS_point_LMB):
            kx, ky, kz, inv_jac = self.rust_worker.inv_parameterize(v, i_v, self.E_cm**2)
            inv_aL_jac *= inv_jac
            aL_xs.extend([kx, ky, kz])

        if self.debug: print('inv_aL_jac=',inv_aL_jac)

        # The final jacobian must then be our param. jac together with that of t divided by the one from alphaloop.
        final_jacobian = wgt_param * wgt_t * normalising_func * inv_aL_jac

        if self.debug: print('xs=',aL_xs)
        if self.debug: print('jacobian=',final_jacobian)

        return aL_xs, final_jacobian

class generator_epem_a_ddx_SG_QG0(CustomGenerator):

    def __init__(self, dimensions, rust_worker, SG_info, model, h_function, hyperparameters, debug=0, **opts):

        self.rust_worker = rust_worker
        self.model = model
        self.SG_info = SG_info
        self.dimensions = dimensions
        self.h_function = h_function
        self.debug = debug

        self.topology = (
                    # s-channels first:
                    base_objects.VertexList([]),
                    # t-channels then:
                    base_objects.VertexList([
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

        E_cm = math.sqrt(max(sum(LorentzVector(vec) for vec in hyperparameters['CrossSection']['incoming_momenta']).square(),0.0))


        self.generator = PS.SingleChannelPhasespace([0.]*2, [0.]*3,
                beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0), 
                model=model, topology=self.topology, path = [[0,],[]],
                dimensions = self.dimensions
            )
        
        print("Considering the following topology:")
        print("-"*10)
        print(self.generator.get_topology_string(self.generator.topology, path_to_print=self.generator.path))
        print("-"*10)

    def __call__(self, random_variables):
        """ 
        Generates a point using a sampling reflecting the topolgy of SG_QG0.
        It will return the point directly in x-space as well as the final weight to combine
        the result with.
        """

        return [1.0,2.0,2.0], 1.0

class generator_epem_a_jjjj_SG_QG36(CustomGenerator):

    def __init__(self, dimensions, rust_worker, SG_info, model, h_function, hyperparameters, debug=0, **opts):

        self.rust_worker = rust_worker
        self.model = model
        self.SG_info = SG_info
        self.dimensions = dimensions
        self.h_function = h_function
        self.debug = debug

        # The process is e+(1) e-(2) > a > d(3) d~(4) g(5) g(6)
        # The SG_QG36 topoology assigns externals as follows:
        # LMB vec #1 = pq3 = -d(3)
        # LMB vec #2 = pq6 = +d~(4)
        # LMB vec #3 = pq10 = +g(5)
        # dependnet momentum = pq9 = -(-#1 +#2 +#3 +Q) 

        # Topology is 
        # g(5) g(6) -> -1
        # -1 d~(4) -> -2
        # -2 d(3) -> incoming "decaying" particle

        self.topology = (
                    # s-channels first:
                    base_objects.VertexList([
                        base_objects.Vertex({
                            'id': DUMMY, # Irrelevant
                            'legs': base_objects.LegList([
                                base_objects.Leg({
                                    'id': 21,
                                    'number': 5,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 21,
                                    'number': 6,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 21,
                                    'number': -1,
                                    'state': FINAL,
                                })
                            ])
                        }),
                        base_objects.Vertex({
                            'id': DUMMY, # Irrelevant
                            'legs': base_objects.LegList([
                                base_objects.Leg({
                                    'id': 21,
                                    'number': -1,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 1,
                                    'number': 4,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 1,
                                    'number': -2,
                                    'state': FINAL,
                                })
                            ])
                        }),
                        base_objects.Vertex({
                            'id': DUMMY, # Irrelevant
                            'legs': base_objects.LegList([
                                base_objects.Leg({
                                    'id': 1,
                                    'number': -2,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 1,
                                    'number': 3,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 22,
                                    'number': -3,
                                    'state': FINAL,
                                })
                            ])
                        }),
                    ]),
                    # t-channels then:
                    base_objects.VertexList([
                        # The dummy vertex below is just to close the diagram and connect
                        # with the number -3 which is to be understood here as the initial state #2.
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
                                    'number': -3,
                                    'state': FINAL,
                                }),
                                base_objects.Leg({
                                    'id': 11,
                                    'number': -4,
                                    'state': INITIAL,
                                })
                            ])
                        }),
                    ])
                )

        E_cm = math.sqrt(max(sum(LorentzVector(vec) for vec in hyperparameters['CrossSection']['incoming_momenta']).square(),0.0))
        self.E_cm = E_cm
        if debug: print("E_cm detected=%.16e"%self.E_cm)

        self.generator = PS.SingleChannelPhasespace([0.]*2, [0.]*(SG_info['topo']['n_loops']+1),
                beam_Es =(self.E_cm/2.,self.E_cm/2.), beam_types=(0,0), 
                model=model, topology=self.topology, #path = [[0,1], []],
                dimensions = self.dimensions
        )

        print("Considering the following topology:")
        print("-"*10)
        print(self.generator.get_topology_string(self.generator.topology, path_to_print=self.generator.path))
        print("-"*10)
        print("With the following path = ",self.generator.path)

#       For testing we can substitute here a flat PS generator.
#        self.generator = PS.FlatInvertiblePhasespace(
#            [0.]*2, [0.]*(SG_info['topo']['n_loops']+1), 
#            beam_Es =(E_cm/2.,E_cm/2.), beam_types=(0,0),
#            dimensions = self.dimensions
#        )

        # Variables for testing that the inverse jacobian works well.
        # Sett test_inverse_jacobian to True in order to debug the inverse jacobian.
        self.test_inverse_jacobian = False
        if self.test_inverse_jacobian:
            print("WARNING: Running in debug mode 'test_inverse_jacobian'.")
            time.sleep(1.0)
        self.rv = None
        self.i_count = 0

    def __call__(self, random_variables):
        """ 
        Generates a point using a sampling reflecting the topolgy of SG_QG36.
        It will return the point directly in x-space as well as the final weight to combine
        the result with.
        """

        #random_variables=[4.27617846e-01, 5.55958133e-01, 1.57593910e-01, 1.13340867e-01, 2.74746432e-01, 2.16284116e-01, 1.99368173e-04, 8.08726361e-02, 1.40842072e-01]


        if self.test_inverse_jacobian:
            if self.i_count==0:
                self.rv = copy.deepcopy(random_variables)
            else:
                random_variables = self.rv
            self.i_count += 1

        x_t = random_variables[self.generator.dim_name_to_position['t']]

        if self.debug: print('x_t=',x_t)

        # We must now promote the point to full 3^L using the causal flow.
        # First transform t in [0,1] to t in [0,\inf]
        rescaling_t, wgt_t = self.h_function(x_t)

        if self.test_inverse_jacobian:
            rescaling_t /= float(self.i_count)

        if self.debug: print('t, wgt_t=',rescaling_t, wgt_t)

        PS_point, wgt_param, x_1, x_2 = self.generator.get_PS_point(random_variables)

        if self.debug: print("PS point from param:\n", LorentzVectorList(PS_point))
        
        if self.debug: print("orig wgt_param:",wgt_param)
        # Correct for the 1/(2E) of each cut propagator that alphaLoop includes but which is already included in the normal PS parameterisation
        for v in PS_point[2:]:
            wgt_param *= 2*v[0]
        if self.debug: print("after wgt_param:",wgt_param)

        # Now rotate this point to the lMB of the topology.
        # In this case, since we generate flat, we don't care and simply take the first :-1 vectors
        # (and remove the initial external states)
        # and strip off the energies

        # The process is e+(1) e-(2) > a > d(3) d~(4) g(5) g(6)
        # The SG_QG36 topoology assigns externals as follows:
        # LMB vec #1 = pq3 = -d(3)
        # LMB vec #2 = pq6 = +d~(4)
        # LMB vec #3 = pq10 = +g(5)
        # dependnet momentum = pq9 = -(-#1 +#2 +#3 +Q) 
        PS_point_LMB = vectors.LorentzVectorList([])
        PS_point_LMB.append(-PS_point[2])
        PS_point_LMB.append(PS_point[3])
        PS_point_LMB.append(PS_point[4])
        PS_point_LMB.append(-PS_point[2]-PS_point[3]-PS_point[4]-vectors.LorentzVector([self.E_cm,0.,0.,0.]))

        PS_point_LMB_vec = [ list(v[1:]) for v in PS_point_LMB ]
        if self.debug: print("PS point in LMB:\n", pformat(PS_point_LMB_vec))

        # Apply the rescaling
        rescaled_PS_point_LMB_vec = [ [ rescaling_t*ki for ki in k ] for k in PS_point_LMB_vec]

        k_norms = [ math.sqrt(sum([ ki**2 for ki in k ])) for k in PS_point_LMB_vec ]

        # Exclude normalising func when testing the inverse jacobian.
        if not self.test_inverse_jacobian:
            normalising_func = self.h_function.PDF(rescaling_t)
        else:
            normalising_func = 1.0
        
        if self.debug: print('normalising_func=',normalising_func)

        # Now we inverse parameterise this point using alphaLoop (non-multichanneled internal param)
        aL_xs = []

        if self.debug: print('h(t)=',self.h_function.PDF(rescaling_t))
        if self.debug: print('1/t=',1./rescaling_t)
        if self.debug: print('h(1/t)=',self.h_function.PDF(1./rescaling_t))
        if self.debug: print('rescaling_t**(3*(len(PS_point_LMB)-1))=',rescaling_t**(3*(len(PS_point_LMB)-1)))
        if self.debug: print('sum(k_norms)=',sum(k_norms))
        if self.debug: print('1./(rescaling_t*sum(k_norms))=',1./(rescaling_t*sum(k_norms)))

        inv_aL_jac = ( 1.0 / self.h_function.PDF(1./rescaling_t) ) * (rescaling_t**((len(PS_point_LMB)-1)*3)) * ( rescaling_t * sum(k_norms) )

        for i_v, v in enumerate(rescaled_PS_point_LMB_vec[:-1]):
            kx, ky, kz, inv_jac = self.rust_worker.inv_parameterize(v, i_v, self.E_cm**2)
            inv_aL_jac *= inv_jac
            aL_xs.extend([kx, ky, kz])

        if self.debug: print('inv_aL_jac=',inv_aL_jac)

        # The final jacobian must then be our param. jac together with that of t divided by the one from alphaloop.
        final_jacobian = wgt_param * wgt_t * normalising_func * inv_aL_jac

        if self.debug: print('xs=',aL_xs)
        if self.debug: print('jacobian=',final_jacobian)

        return aL_xs, final_jacobian

#=========================================================================================
# Standalone main for running the example
#=========================================================================================
if __name__ == '__main__':


    # Parse arguments
    parser = ArgumentParser()
    parser.add_argument("--process", dest="process", 
                        help="define either a process_name form the card/ subfolder or a path to a alphaLoop@MG5 folder")
    parser.add_argument("--phase", dest="phase",default="real",
                        help="which phase to integral")
    parser.add_argument("--deformation", dest="deformation", default="none",
                        help="select deformation strategy")
    parser.add_argument("--min_jets", dest="min_jets", default=0, type=int,
                        help="minimum number of jets in the dinal state")
    parser.add_argument("--min_jpt",  dest="min_jpt", default=100.0, type=float,
                        help="Jet cutoff")
    parser.add_argument("--dr",  dest="dr", default=0.4, type=float,
                        help="Jet delta R")
    parser.add_argument("--n_max", dest="n_max", default=int(1e9), type=int,
                        help="max number of evaluation with Vegas")
    parser.add_argument("--n_new", dest="n_new", default=int(1e5), type=int,
                        help="n_new Vegas")
    parser.add_argument("--n_start", dest="n_start", default=int(1e5), type=int,
                        help="n_start Vegas")
    parser.add_argument("--n_increase", dest="n_increase", default=int(1e5), type=int,
                        help="n_increase Vegas")
    parser.add_argument("--multi_channeling", action="store_true", dest="multi_channeling", default=True,
                        help="enable multi-channeling during integration")
    parser.add_argument("--no_multi_channeling", action="store_false", dest="multi_channeling", default=True,
                        help="Disable multi-channeling during integration")
    parser.add_argument("--cores", dest="cores", default=CORES, type=int,
                        help="number of cores used during integration")
    parser.add_argument("--diag_name", dest="diag_name", default=None, type=str,
                        help="Integrate single SG: selecting by name")
    parser.add_argument("--seed", dest="seed", default=1, type=int,
                        help="Fix random seed")
    parser.add_argument("--h_function", dest="h_function", default="HFunction(2)", type=str,
                        help="Fix random seed")
    parser.add_argument("--sigma", dest="sigma", default=2, type=int,
                        help="H function spread to specify to alphaloop")
    parser.add_argument("--verbosity", dest="verbosity", default=1, type=int,
                        help="Fix verbosity ( set to 2 in order to see the evolution of the integration )")
    parser.add_argument("--n_points_survey", dest="n_points_survey", default=3, type=int,
                        help="Fix the number of points for the survey")
    parser.add_argument("--n_points_refine", dest="n_points_refine", default=3, type=int,
                        help="Fix the number of points for the refine")
    parser.add_argument("--run_mode", dest="run_mode", default="al", type=str,
                        help="Choose to run nomal multi-channeled alphaLoop (al) or the Single Channel PS generator (scps) ")
    parser.add_argument("--debug", dest="debug", default=0, type=int,
                        help="debug level")

    args = parser.parse_args()
    process_name = None
    if args.process is None:
        if args.run_mode!="test_h_function":
            raise ValueError(
                "Missing process definition: %s --process=epem_a_ttx --run" % sys.argv[0])
    else:            
        process_name = os.path.basename(re.sub(r'/$','',args.process))


    model = get_model()
    
    # If we had resonannces we could modify the model parameters below as follows
    # An BreitWigner centered here at MZ set here to 150 GeV
#    model.get('parameter_dict')['mdl_MZ'] = 150.0
    #  An BreitWigner of width here WZ set to  60 GeV
#    model.get('parameter_dict')['mdl_WZ'] = 60.0

    E_cm = 1000.0

    random.seed(args.seed)

#   Now set hyperparameters

    # Extra settings
    hyper_settings = {'General': {},
                      'CrossSection': {},
                      'Deformation': {},
                      'Integrator': {},
                      'Observables': {},
                      'Selectors': {},
                      'Parameterization': {},
                      }

    if args.n_new < 20:
        args.n_new = int(10**args.n_new)
    if args.n_max < 20:
        args.n_max = int(10**args.n_max)
    if args.n_points_survey < 20:
        args.n_points_survey = int(10**args.n_points_survey)
    if args.n_points_refine < 20:
        args.n_points_refine = int(10**args.n_points_refine)
    if args.n_increase < 20:
        if args.n_increase == 0:
            args.n_increase = 0
        else:
            args.n_increase = int(10**args.n_increase)
    if args.n_start < 20:
        args.n_start = int(10**args.n_start)

    selected_h_function = None
    try:
        selected_h_function = eval(args.h_function)
        selected_h_function.debug = args.debug
    except Exception as e:
        print("Cannot instantiate specified H-function: '%s'"%args.h_function)
        raise

    CORES = args.cores
    runner = cluster.MultiCore(CORES)

    if args.run_mode == 'test_h_function':

        print("Dummy integrand for testing h function:")

        my_integrand = TestHFuncIntegrand(
            integrands.DimensionList([ integrands.ContinuousDimension('x_1',lower_bound=0.0, upper_bound=1.0), ])
            , selected_h_function, debug=args.debug
        )

        # Use the line below instead to run on one core only
        #runner = cluster.onecore
        my_integrator = vegas3.Vegas3Integrator(my_integrand, 
                n_points_survey=args.n_points_survey, n_points_refine=args.n_points_refine, accuracy_target=None,
                verbosity=args.verbosity, cluster=runner
        )

        result = my_integrator.integrate()

        print('')
        print("Result of the test integration using h function '%s' with %d function calls (target is 1.0): %.4e +/- %.2e"%(
            args.h_function, my_integrator.tot_func_evals, result[0], result[1]))
        print('')
        sys.exit(0)

    # Integrator Settings
    hyper_settings['Integrator']['integrated_phase'] = args.phase
    hyper_settings['General']['deformation_strategy'] = args.deformation
    hyper_settings['General']['debug'] = args.debug
    hyper_settings['Integrator']['state_filename_prefix'] = "%s_%s_" % (
        process_name, "nocut" if args.min_jets == 0 else args.min_jpt)
    hyper_settings['Integrator']['keep_state_file'] = True
    hyper_settings['Integrator']['load_from_state_file'] = True
    hyper_settings['Integrator']['n_max'] = args.n_max
    hyper_settings['Integrator']['n_new'] = args.n_new
    hyper_settings['Integrator']['n_start'] = args.n_start
    hyper_settings['Integrator']['n_increase'] = args.n_increase
    # Selector Settings
    hyper_settings['Selectors']['active_selectors'] = [] if args.min_jets == 0 else ['jet']
    hyper_settings['Selectors']['jet'] = {
        'min_jets': args.min_jets,
        'dR': args.dr,
        'min_jpt': args.min_jpt}
    # Set multi-channeling strategy
    hyper_settings['General']['multi_channeling'] = args.multi_channeling
    hyper_settings['General']['multi_channeling_including_massive_propagators'] = args.multi_channeling
    # Cross section settings
    hyper_settings['CrossSection']['incoming_momenta'] = [[500, 0, 0 ,500],[500, 0, 0, -500]]
#    hyper_settings['CrossSection']['incoming_momenta'] = [[0.5, 0, 0 ,0.5],[0.5, 0, 0, -0.5]]
    hyper_settings['CrossSection']['m_uv_sq'] = (75.0)**2
    hyper_settings['CrossSection']['mu_r_sq'] = (91.188)**2
#    hyper_settings['CrossSection']['m_uv_sq'] = 1.0**2
#    hyper_settings['CrossSection']['mu_r_sq'] = 1.0**2

    hyper_settings['CrossSection']['mu_r_sq'] = (91.188)**2

    hyper_settings['CrossSection']['gs'] = 1.2177157847767195

    hyper_settings['CrossSection']['NormalisingFunction'] = {
        'name'                              :   'left_right_polynomial',
        'center'                            :   1.0,
        'spread'                            :   float(args.sigma)
    }

    diag_yaml_path = pjoin(args.process,'FORM','Rust_inputs','%s.yaml'%args.diag_name)
    if not os.path.isfile(diag_yaml_path):
        diag_yaml_path = pjoin(args.process,'Rust_inputs','%s.yaml'%args.diag_name)
        if not os.path.isfile(diag_yaml_path):
            print("Supergraph yaml file '%s.yaml' could not be found."%args.diag_name)
            raise BaseException("Could not find yaml file '%s.yaml'."%args.diag_name)

    # SG yaml info
    SG_info = yaml.load(open(diag_yaml_path, 'r'), Loader=yaml.Loader)
    n_loops = SG_info['topo']['n_loops']

    if args.run_mode != 'aL':
        # We need to disable the multichanneling when not running default aL
        hyper_settings['General']['multi_channeling'] = False

    hyperparameters_path = set_hyperparameters(process_name, args.diag_name, workspace=WORKSPACE, min_jets=args.min_jets, multi_settings=hyper_settings)
    hyperparameters = yaml.load(open(hyperparameters_path, 'r'), Loader=yaml.Loader)
    rust_worker = get_rust_worker(
        args.process,
        diag_yaml_path,
        hyperparameters_path
    )

    dimensions = integrands.DimensionList(
            [ integrands.ContinuousDimension('x_%d'%i,lower_bound=0.0, upper_bound=1.0) for i in range(1, n_loops*3) ]+
            [ integrands.ContinuousDimension('t',lower_bound=0.0, upper_bound=1.0), ] 
            )
    my_PS_generator = None
    try:
        my_PS_generator = eval("generator_%s(dimensions, rust_worker, SG_info, model, selected_h_function, hyperparameters, debug=args.debug)"%args.run_mode)
    except Exception as e:
        print("Unreckognized run mode: '%s'"%args.run_mode)
        raise

    print("Integrating '%s' with the parameterisation %s:"%(args.diag_name,args.run_mode))

    my_integrand = DefaultALIntegrand(rust_worker, my_PS_generator, debug=args.debug)

    # Use the line below instead to run on one core only
    #runner = cluster.onecore
    my_integrator = vegas3.Vegas3Integrator(my_integrand, 
            n_points_survey=args.n_points_survey, n_points_refine=args.n_points_refine, accuracy_target=None,
            verbosity=args.verbosity, cluster=runner
    )

    result = my_integrator.integrate()

    print('')
    print("Result of integration using the '%s' parameterisation with %d function calls: %.7e +/- %.2e"%(
        args.run_mode, my_integrator.tot_func_evals, result[0], result[1]))
    print('')
