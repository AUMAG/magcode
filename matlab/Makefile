
NAME = magnetforces
SRC = $(NAME).web
SCPT = $(NAME).m
DRV = magnetocode-matlab
DOC = $(DRV).pdf
DIR = /tmp/mweb
MWEB = $(DIR)/mweb.tex $(DIR)/swebbind.sty $(DIR)/webfiles.sty $(DIR)/matlabweb.sty

all: $(SCPT) $(DOC)
code: $(SCPT)

$(SCPT):  $(SRC)
	mtangle $(SRC)

$(DOC): $(DIR)/$(DOC)
	cp -f $< $@

$(DIR)/$(DOC): $(DIR)/$(SRC) $(DIR)/$(DRV).tex $(MWEB)
	cd $(DIR)/; mweave $(SRC)
	cd $(DIR)/; pdflatex -draftmode $(DRV)
	cd $(DIR)/; bibtex $(DRV)
	cd $(DIR)/; pdflatex $(DRV)

$(DIR)/$(SRC): $(SRC)
	mkdir -p /tmp/mweb/
	cp -f $< $@

$(DIR)/$(DRV).tex: $(DRV).tex
	cp -f $< $@

$(DIR)/mweb.tex: mweb/mweb.tex
	cp -f $< $@

$(DIR)/%.sty: mweb/%.sty
	cp -f $< $@

