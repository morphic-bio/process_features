#add library for openmp
#add library for math functions

CC=gcc  
CFLAGS=-Ofast -Wall $(INCLUDES) 
LDFLAGS=-lz -lglib-2.0 -fopenmp -lm
ifdef DEBUG
	CFLAGS=-ggdb -Wunused-function -fsanitize=address $(INCLUDES)
endif
ifdef NO_HEATMAP
	CFLAGS+=-DNO_HEATMAP
else
	LDFLAGS+=$(shell pkg-config --libs cairo)
	CFLAGS+=$(shell pkg-config --cflags cairo)
endif

INCLUDES=-I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include -I/usr/lib64/glib-2.0/include -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/hdf5/serial
#INCLUDES=-I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include -I/usr/lib64/glib-2.0/include -I/usr/lib/x86_64-linux-gnu/glib-2.0/include 

default: clean assignBarcodes
all: clean heatmap assignBarcodes

assignBarcodes: assignBarcodes.c queue.c
	$(CC) $(CFLAGS) -o assignBarcodes assignBarcodes.c queue.c $(LDFLAGS) $(INCLUDES)
heatmap: heatmap.c
	$(CC) $(CFLAGS) -o heatmap heatmap.c $(LDFLAGS) $(INCLUDES)
clean:
	rm -f heatmap assignBarcodes

