
SRC= formal.tex 

all: main.pdf

main.bbl: main.tex
	pdflatex -file-line-error main.tex

main.pdf: $(SRC) main.bbl
	pdflatex -file-line-error main.tex

clean:
	rm *.syn* *.aux *.log *~ *.bbl *.bst *.out
