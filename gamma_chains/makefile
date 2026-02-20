F2PY = f2py
FC = gfortran

all: compute_chain gamma_chain.so

gamma_chain.so: compute_chain.f90  
	touch __init__.py
	$(F2PY) -m gamma_chain -c compute_chain.f90 --fcompiler=$(FC)  -I ..

compute_chain_processed.f90: compute_chain.f90
	gfortran -E compute_chain.f90 -cpp  -o compute_chain_processed.f90

compute_chain:	compute_chain.f90
	$(FC) -o compute_chain compute_chain.f90 -O5 -I ..

clean:
	rm -f *.o compute_chain gamma_chain.so
