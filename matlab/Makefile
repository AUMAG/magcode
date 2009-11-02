
NAME = magnetforces
SRC = $(NAME).web
SCPT = $(NAME).m
DOC = $(NAME).pdf
DRV = mweb-driver
DIR = /tmp/mweb
MWEB = $(DIR)/mweb.tex $(DIR)/swebbind.sty $(DIR)/webfiles.sty

all: $(SCPT) $(DOC)
tangle: $(SCPT)

$(SCPT):  $(SRC)
	mtangle $(SRC)

$(DOC): $(DIR)/$(DRV).pdf
	cp -f $(DIR)/$(DRV).pdf $(DOC)

$(DIR)/$(DRV).pdf: $(DIR)/$(SRC) $(DIR)/$(DRV).tex $(MWEB)
	cd $(DIR)/; mweave $(SRC)
	cd $(DIR)/; pdflatex "\def\WEBNAME{$(NAME)}\input{$(DRV)}"
	cd $(DIR)/; bibtex $(DRV)
	cd $(DIR)/; pdflatex "\def\WEBNAME{$(NAME)}\input{$(DRV)}"
	cd $(DIR)/; pdflatex "\def\WEBNAME{$(NAME)}\input{$(DRV)}"

$(DIR)/$(SRC): $(SRC)
	mkdir -p /tmp/mweb/
	cp $(SRC) $(DIR)/$(SRC)

$(DIR)/$(DRV).tex: $(DRV).tex
	cp $(DRV).tex $(DIR)/$(DRV).tex

$(DIR)/mweb.tex: mweb/mweb.tex
	cp -f $< $@

$(DIR)/swebbind.sty: mweb/swebbind.sty
	cp -f $< $@

$(DIR)/webfiles.sty: mweb/webfiles.sty
	cp -f $<  $@
