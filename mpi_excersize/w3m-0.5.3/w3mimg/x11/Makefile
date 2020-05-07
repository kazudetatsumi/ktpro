#
# w3mimg/x11/Makefile
#
#

SHELL=/bin/sh
srcdir = .
top_srcdir = ../..
VPATH = .
CFLAGS=$(OPTS) -I../.. -I$(top_srcdir) -I$(srcdir) -g -O2 -I$(srcdir)/libwc  -I/usr/local/opt/openssl/include -I/opt/X11/include -DHAVE_CONFIG_H $(IMGCFLAGS)
RM=rm

IMGCFLAGS=-I/usr/local/Cellar/imlib2/1.5.1/include
OBJS=x11_w3mimg.o

all: x11_w3mimg.o

x11_w3mimg.o: x11_w3mimg.c
	$(CC) $(CFLAGS) -c $<

clean:
	@-$(RM) -f *.o

distclean: clean
	-$(RM) -f Makefile

#


