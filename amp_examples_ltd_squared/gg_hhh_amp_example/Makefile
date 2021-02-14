.PHONY: process cards all clean

MG5Path = /home/armin/my_programs/MG5_py_3_alphaloop_master


process:
	python3 $(MG5Path)/bin/mg5_aMC --mode alphaloop $(CURDIR)/gg_hhh_amp.aL;

cards:
	sed  "s+PATHTOHERE+$(CURDIR)+g" gg_hhh_amp_template.aL > gg_hhh_amp.aL;
	sed  "s+PATHTOHERE+$(CURDIR)+g" gg_hhh_runcard_template.yaml > gg_hhh_runcard.yaml;

	

all: cards process
	


clean:
	rm -rf gg_hhh_amp.aL
	rm -rf gg_hhh_runcard.yaml
	rm -rf TEST_AMP_gg_hhh