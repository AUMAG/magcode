
help:
	@echo 'This is the MAGCODE makefile.'
	@echo 'The following targets are available to generate respectively:'
	@echo '  '
	@echo '     code - Matlab .m files and test suite'
	@echo '      doc - documentation in PDF format'
	@echo '      all - both of the above'
	@echo '  '
	@echo 'E.g., type `make all` to generate everything.'

Name = magcode-matlab
Source  = $(Name).nw
DocName = $(Name)
Doc  = $(DocName).pdf

Dir = /tmp/$(Name)
TestDir = testsuite/

Outputs = magnetforces.m multipoleforces.m

TestNames = testall.m \
        test001a.m test001b.m test001c.m test001d.m test001e.m \
        test002a.m test002b.m test002c.m test002d.m \
        test003a.m
Tests = $(addprefix $(TestDir),$(TestNames))

all: $(Outputs) $(Tests) $(Doc)
code: $(Outputs) $(Tests)
doc: $(Doc)

.PHONY: all code doc

$(Outputs):	$(Source)
	notangle -R' $@ ' $< > $@

$(Tests):	$(Source)
	@mkdir -p testsuite
	notangle -R' $(shell basename $@) ' $< > $@

$(Doc):	$(Dir)/$(Source)
	cd $(Dir); \
	noweave -delay -index $(Source) > $(DocName).tex \
		&& pdflatex $(DocName) \
		&& bibtex $(DocName) \
		&& pdflatex $(DocName)
	cp $(Dir)/$(Doc) $@

$(Dir)/$(Source): $(Source)
	@mkdir -p $(Dir)
	cp $< $@