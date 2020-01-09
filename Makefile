srpg: srpg.c
	$(CC) -O2 -Wall -o "$@" $<

clean:
	rm -f srpg

.PHONY: clean
