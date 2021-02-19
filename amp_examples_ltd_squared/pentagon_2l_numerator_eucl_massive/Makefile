.PHONY: process cards all clean

MG5Path = /home/armin/my_programs/MG5_py_3_alphaloop_master


process:
	python3 $(MG5Path)/bin/mg5_aMC --mode alphaloop $(CURDIR)/pentagon_amp.aL;

cards:
	sed  "s+PATHTOHERE+$(CURDIR)+g" pentagon_amp_template.aL > pentagon_amp.aL;
	sed  "s+PATHTOHERE+$(CURDIR)+g" pentagon_runcard_template.yaml > pentagon_runcard.yaml;

	

all: cards process
	


clean:
	rm -rf pentagon_amp.aL
	rm -rf pentagon_runcard.yaml;
	rm -rf TEST_pentagon