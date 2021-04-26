import os 
import sys

# ADJUST THESE PATH's
if True:
    sys.path.append('/home/armin/my_programs/MG5_py_3_alphaloop_master/PLUGIN/alphaloop/LTD/')
    runcard = "/home/armin/my_programs/alphaLoopMisc/amps_in_standalone_rewrite/1l_pentagon_subtracted/runcard_pentagon.yaml"
    from amplitudes_rewrite.amp_exporter import AmpExporter


os.system('rm -r -f pentagon_zeno')
exporter = AmpExporter.from_runcard(runcard)
exporter.export()