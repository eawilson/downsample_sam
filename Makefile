CC          = gcc
CFLAGS      = -Wall -O2
LDFLAGS     = -lz
prefix      = /usr/local
exec_prefix = $(prefix)/bin

src = $(wildcard *.c)
obj = $(src:.c=.o)

downsample_sam: $(obj)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) downsample_sam

.PHONY: install
install:
	cp downsample_sam $(exec_prefix)







