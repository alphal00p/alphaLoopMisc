#!/usr/bin/env python
import os
import sys
import argparse
from pprint import pprint

_VALID_OBSERVABLES = ('ptj1','ptj2','ptj3','xsec')
_VALID_ORDERS = ['NLO','LO']

parser = argparse.ArgumentParser(description='Compute semi differential cross-sections.')
requiredNamed = parser.add_argument_group('required named arguments')

# Required options
requiredNamed.add_argument('--min_value', '--min', type=float,
                    help='Minimum allowed value of the observable', required=True)
requiredNamed.add_argument('--max_value', '--max', type=float,
                    help='Maximum allowed value of the observable', required=True)
requiredNamed.add_argument('--HwU-path', '--HwU', type=str,
                    help='Path of the raw source of the differential distribution histograms.', required=True)

# Optional options
parser.add_argument('--observable', '-o', default=_VALID_OBSERVABLES[0], choices=_VALID_OBSERVABLES,
                    help='Selected observable to compute the semi-inclusive cross-section from.')
parser.add_argument('--MG5aMC-path', '--MGpath', type=str, default='/Users/valentin/Documents/MG5/2.6.6',
                    help='Path of your MadGraph5_aMC@NLO distribution.')
parser.add_argument('--perturbative-order', '--po', type=str, default=_VALID_ORDERS[0],choices=_VALID_ORDERS,
                    help='Perturbative order considered')

args = parser.parse_args()

if args.observable=='xsec':
    args.min_value = 0.5
    args.max_value = 1.5

# EDIT below the local path of your MadGraph5_aMC@NLO distribution
# Alternatively set the environment variable MG5AMC when running this script
MG5aMC_root_path = os.getenv('MG5AMC',args.MG5aMC_path)
sys.path.insert(0,MG5aMC_root_path)
try:
    import madgraph.various.histograms as HwU 
except ZeroDivisionError:
    print "WARNING: Could not import the histograms package from MadGraph5_aMC@NLO.\n"+\
          "Make sure that the specified root path of your MG5aMC installation is correct: %s\n"%MG5aMC_root_path+\
          "Note: you can specify this path either by editing this script or through the env. variable MG5AMC or with the option --MGpath."
    sys.exit(1)


observable_to_histo_title = {'ptj1' : 'j1 pT',
                             'ptj2' : 'j2 pT',
                             'ptj3' : 'j3 pT',
                             'xsec' : 'total rate %s'%args.perturbative_order
                            }

# Find the histogram with the specified observable
all_histos = HwU.HwUList(args.HwU_path)
selected_histo = None
for histo in all_histos:
    #print "title: %s vs %s"%(histo.title, observable_to_histo_title[args.observable])
    if histo.title != observable_to_histo_title[args.observable]:
        continue
    # We do no need to check the type when interested in the overall xsec
    if args.observable=='xsec':
        selected_histo = histo
        break
    #print "type: %s vs %s"%(histo.type, args.perturbative_order)    
    if histo.type != args.perturbative_order:
        continue
    selected_histo = histo
    break
if selected_histo is None:
    print "Could not find the histogram for the desired observable '%s' and perturbative order '%s'."%(
        args.observable, args.perturbative_order
    )
    sys.exit(1)

# Sanity check
if selected_histo.bins[0].boundaries[0] > args.min_value:
    print "Minimal value specified %.4e is smaller than the minimum value of the '%s' histogram (%.4e)"%(
        args.min_value, args.observable, selected_histo.bins[0].boundaries[0]
    )
    sys.exit(1)
if selected_histo.bins[-1].boundaries[1] < args.max_value:
    print "Maximal value specified %.4e is smaller than the minimum value of the '%s' histogram (%.4e)"%(
        args.max_value, args.observable, selected_histo.bins[-1].boundaries[1]
    )
    sys.exit(1)

# Now aggregate the results
cum_xsec = 0.
cum_error = 0.
for abin in selected_histo.bins:
    start_x = None
    if abin.boundaries[1]>args.min_value:
        start_x = max(abin.boundaries[0],args.min_value)
    end_x = None
    if abin.boundaries[0]<args.max_value:
        end_x = min(abin.boundaries[1],args.max_value)
    if start_x is None or end_x is None:
        continue
    if start_x != abin.boundaries[0]:
        print "WARNING: your specified minimal value %.4f does not match exactly the lower bound of a bin. Use %.4f instead to avoid interpolation."%(
            args.min_value, abin.boundaries[0]
        )
    if end_x != abin.boundaries[1]:
        print "WARNING: your specified maximal value %.4f does not match exactly the lower bound of a bin. Use %.4f instead to avoid interpolation."%(
            args.max_value, abin.boundaries[1]
        )
    #print("%s, %s"%(abin.boundaries, abin.wgts))
    factor = abs((end_x-start_x)/(abin.boundaries[1]-abin.boundaries[0]))
    cum_xsec += factor*abin.get_weight('central')
    cum_error += factor*abin.get_weight('stat_error')

print ""
print "Semi-inclusive xsec for %.4e < %s < %.4e:"%(args.min_value, args.observable, args.max_value)
print ""
print " %.8e +/- %.8e [pb] (%.2g%%)"%(cum_xsec, cum_error, abs(cum_error/cum_xsec)*100.0)
print ""
