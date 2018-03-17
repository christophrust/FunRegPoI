# Makefile for FunRegPoI
all: doc pkg
.PHONY: doc pkg

doc:
	R -e 'devtools::document()'

pkg:
	cd .. && \
	  R CMD build FunRegPoI && \
	  R CMD check FunRegPoI_1.2.tar.gz
	rm R/*.R~
