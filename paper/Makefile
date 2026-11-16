# The paper is generated from the Verso source in verso/; see verso/README.md.
.PHONY: all submission clean

all:
	./verso/build.sh

# A self-contained TeX upload (no BibTeX run needed), in verso/_out/submission/.
submission: all
	python3 package.py

clean:
	rm -rf verso/_out build
