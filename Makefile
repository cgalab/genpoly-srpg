srpg: srpg.c
	cc -O2 -Wall -o "$@" $<

clean:
	rm -f srpg

.PHONY: clean
