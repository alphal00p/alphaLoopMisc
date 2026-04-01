python3 evaluator.py dot_files/ScalarGravity_2L_processed.dot \
 -ve \
 --n_loops 2 \
 --n_edges 11 \
 --n_graphs 2 \
 --build_dot_products_form \
 --build_spenso_parametric_form \
 --n_core_for_evaluator_building 100 \
 --n_core_for_dot_processing 100 \
 --do-check-emr-energies-linearity

python3 evaluator.py dot_files/ScalarGravity_3L_processed.dot \
 -ve \
 --n_loops 3 \
 --n_edges 14 \
 --n_graphs 100000 \
 --build_dot_products_form \
 --n_core_for_evaluator_building 100 \
 --n_core_for_dot_processing 100 \
 --do-check-emr-energies-linearity

python3 evaluator.py dot_files/ScalarGravity_4L_processed.dot \
 -ve \
 --n_loops 4 \
 --n_edges 17 \
 --n_graphs 100000 \
 --build_dot_products_form \
 --n_core_for_evaluator_building 100 \
 --n_core_for_dot_processing 100 \
 --do-check-emr-energies-linearity
