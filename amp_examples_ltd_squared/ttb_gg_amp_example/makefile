.PHONY: process cards all clean

MG5Path = /home/armin/my_programs/MG5_py_3_alphaloop_master

CC = gcc
CFLAGS = -Wall
LFLAGS = $(CURDIR)/TEST_AMP_ttb_gg_LO/color_struc_0/lib/libFORM_numerators.so -lm -lquadmath
TARGET = ./numerical_eval/eval_num


process:
	python3 $(MG5Path)/bin/mg5_aMC --mode alphaloop $(CURDIR)/ttb_gg_LO.aL;

cards:
	sed  "s+PATHTOHERE+$(CURDIR)+g" ttb_gg_LO_amp_template.aL > ttb_gg_LO.aL;
	sed  "s+PATHTOHERE+$(CURDIR)+g" ttb_gg_LO_runcard_template.yaml > ttb_gg_LO_runcard.yaml;

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) $(TARGET).c $(LFLAGS) -o $(TARGET)
	

all: cards process $(TARGET)
	


clean:
	rm -rf ttb_gg_LO.aL
	rm -rf ttb_gg_LO_runcard.yaml
	rm -rf TEST_AMP_ttb_gg_LO
	rm -rf $(TARGET)