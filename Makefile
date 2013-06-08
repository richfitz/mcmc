TARGETS = mcmc.pdf

all: ${TARGETS}

%.pdf: %.md
	pandoc $<  --include-in-header=include.tex -V linkcolor:black -o $@
%.html: %.md
	pandoc $< -o $@ --standalone

clean:
	rm -f ${TARGETS}
