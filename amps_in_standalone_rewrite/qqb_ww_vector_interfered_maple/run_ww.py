
import os 
import sys

# ADJUST THESE PATH's. ALSO: adjust PATH in runcard.json
if True:
    sys.path.append('/home/armin/my_programs/MG5_py_3_alphaloop_master/PLUGIN/alphaloop/LTD/')
    runcard = "/home/armin/my_programs/alphaLoopMisc/amps_in_standalone_rewrite/qqb_ww_vector_interfered_maple/WW_runcard.json"
    from amplitudes_rewrite.amp_exporter import AmpExporter


os.system('rm -r -f WW')

exporter = AmpExporter.from_runcard(runcard)
exporter.export()
