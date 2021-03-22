.PHONY: process cards all clean

MG5Path = /home/armin/my_programs/MG5_py_3_alphaloop_master


process:
	python3 $(MG5Path)/bin/mg5_aMC --mode alphaloop $(CURDIR)/pentagon_amp.aL;
	python3 $(MG5Path)/bin/mg5_aMC --mode alphaloop $(CURDIR)/pentagon_amp_new.aL;


cards:
	sed  "s+PATHTOHERE+$(CURDIR)+g" pentagon_amp_template.aL > pentagon_amp.aL;
	sed  "s+PATHTOHERE+$(CURDIR)+g" pentagon_runcard_template.yaml > pentagon_runcard.yaml;
	sed  "s+PATHTOHERE+$(CURDIR)+g" pentagon_amp_template_new.aL > pentagon_amp_new.aL;
	sed  "s+PATHTOHERE+$(CURDIR)+g" pentagon_runcard_template_new.yaml > pentagon_runcard_new.yaml;


	

all: cards process
	


clean:
	rm -rf pentagon_amp.aL pentagon_amp_new.aL
	rm -rf pentagon_runcard.yaml pentagon_runcard_new.yaml;
	rm -rf TEST_pentagon*
