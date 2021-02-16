import model aL_sm-no_widths

set_alphaLoop_option n_rust_inputs_to_generate -1
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False
set_alphaLoop_option consider_edge_orientation_in_graph_isomorphism False
set_alphaLoop_option consider_vertex_id_in_graph_isomorphism False 
set_alphaLoop_option apply_graph_isomorphisms False


# Memory limit (max 1T)
set_alphaLoop_option virtual_memory 100G
set_alphaLoop_option FORM_compile_cores 10

# FORM options
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_optimization 3 
set_alphaLoop_option FORM_compile_arg all 

# SET PATH TO AMPLITUDE RUN-CARD
set_alphaLoop_option AMPLITUDE_runcard_path PATHTOHERE/pentagon_runcard.yaml

# Process
process_amplitude psiS psiS > psiS psiS psiS  []
!rm -rf TEST_pentagon
output amplitude TEST_pentagon
