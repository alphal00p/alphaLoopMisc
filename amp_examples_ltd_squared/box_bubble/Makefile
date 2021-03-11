.PHONY: process cards all clean

MG5Path = /scratch/zcapatti/MG5_aMC_v3_0_2_py3/
CURDIR = /scratch/zcapatti/alphaLoopMisc/amp_examples_ltd_squared/box_bubble

process: cards
	python3 $(MG5Path)/bin/mg5_aMC --mode alphaloop $(CURDIR)/bubble_box_amp.aL;

cards:
	sed  "s+PATHTOHERE+$(CURDIR)+g" bubble_box_amp_template.aL > bubble_box_amp.aL;
	sed  "s+PATHTOHERE+$(CURDIR)+g" bubble_box_runcard_template.yaml > bubble_box_runcard.yaml;

	

all: cards process
	


clean:
	rm -rf bubble_box_amp.aL
	rm -rf bubble_box_runcard.yaml;
	rm -rf TEST_bubble_box
