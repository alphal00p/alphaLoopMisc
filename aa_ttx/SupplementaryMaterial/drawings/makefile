OUTPUT = feynman_diagrams.pdf
N_ROWS_IN_MERGED_OUTPUT = 3
N_COLUMNS_IN_MERGED_OUTPUT = 3
COMPILE_LOG = compilation.log

ALLGRAPHS=all_graphs.pdf
GRAPHS_SRC=$(wildcard *.dot)
# GRAPHS_PDF=$(sort $(patsubst %.tex,%.pdf,$(wildcard *.tex)))
GRAPHS_PDF= GL0.pdf GL1.pdf GL4.pdf GL5.pdf GL8.pdf GL9.pdf GL12.pdf GL13.pdf GL16.pdf GL18.pdf GL24.pdf GL26.pdf GL28.pdf GL30.pdf GL32.pdf GL34.pdf GL36.pdf GL38.pdf GL40.pdf GL42.pdf GL44.pdf GL46.pdf GL80.pdf GL82.pdf GL84.pdf GL86.pdf GL88.pdf GL94.pdf GL96.pdf GL98.pdf GL100.pdf GL102.pdf GL104.pdf GL106.pdf GL112.pdf GL114.pdf GL116.pdf GL118.pdf GL120.pdf GL122.pdf GL124.pdf GL130.pdf GL132.pdf GL134.pdf GL136.pdf GL138.pdf GL144.pdf GL146.pdf GL148.pdf GL150.pdf GL152.pdf GL154.pdf GL156.pdf GL158.pdf GL160.pdf GL162.pdf GL164.pdf GL166.pdf GL168.pdf GL172.pdf GL176.pdf GL178.pdf GL180.pdf GL182.pdf GL188.pdf GL190.pdf GL196.pdf GL198.pdf GL200.pdf GL202.pdf GL204.pdf GL206.pdf GL208.pdf GL212.pdf GL218.pdf GL220.pdf GL222.pdf GL224.pdf GL226.pdf GL228.pdf GL230.pdf GL232.pdf GL238.pdf GL244.pdf GL246.pdf GL256.pdf GL258.pdf GL260.pdf GL262.pdf GL268.pdf GL270.pdf GL272.pdf GL274.pdf GL276.pdf GL278.pdf GL280.pdf GL282.pdf GL284.pdf GL286.pdf GL288.pdf GL290.pdf GL292.pdf GL294.pdf GL296.pdf GL298.pdf GL300.pdf GL362.pdf GL364.pdf GL366.pdf GL368.pdf GL370.pdf GL372.pdf GL374.pdf GL376.pdf GL378.pdf GL380.pdf GL382.pdf GL384.pdf GL386.pdf GL388.pdf GL394.pdf GL396.pdf GL402.pdf GL404.pdf GL406.pdf GL408.pdf GL416.pdf GL418.pdf GL424.pdf GL426.pdf GL428.pdf GL430.pdf GL440.pdf GL442.pdf GL448.pdf GL450.pdf GL452.pdf GL454.pdf

all: $(OUTPUT)

%.pdf: %.dot
	@echo "\nCompiling "$*.pdf" ...\n"
	dot -Tpdf $*.dot -o $*.pdf 2>&1 >> $(COMPILE_LOG)
	sleep 0.1
	pdfcrop $*.pdf $*.pdf 2>&1 >> $(COMPILE_LOG)

$(ALLGRAPHS): $(GRAPHS_PDF)
	pdfunite $^ $@

$(OUTPUT): $(ALLGRAPHS)
	python3 ./combine_pages.py $(ALLGRAPHS) $(OUTPUT) --margins 0.15 0.15 --grid_shape $(N_ROWS_IN_MERGED_OUTPUT) $(N_COLUMNS_IN_MERGED_OUTPUT)
	@echo "\n> Feynman diagrams all written to file '"$(OUTPUT)"'\n"

clean:
	rm -f *.t1 *.log *.dvi *.mp *.log *.1 *.ps *.aux

clean_pdf:
	rm *.pdf
