CC=gcc

# Directories
SRCDIR = src
INCDIR = include
VPATH = $(SRCDIR)

# Library flags (no cairo - heatmaps now use Plotly HTML+JSON)
HTS_LIBS := -lhts

DEFINES :=

# Source files (basenames only)
# Note: heatmap.c now uses Plotly (HTML+JSON), no external dependencies
SRCS_BASE=main.c assignBarcodes.c queue.c globals.c utils.c memory.c io.c plot_histogram.c barcode_match.c heatmap.c
SRCS = $(SRCS_BASE)

# Compiler flags: optimise by default, but if DEBUG=1 build with symbols and no optimisations
ifeq ($(DEBUG),1)
    CFLAGS=-g -ggdb -Wall -O0 -I$(INCDIR) -fopenmp $(DEFINES)
else
    CFLAGS=-g -Wall -O3 -I$(INCDIR) -fopenmp $(DEFINES)
endif
LDFLAGS=-lm -lpthread -lz -fopenmp $(HTS_LIBS)

# Object files
OBJS=$(SRCS:.c=.o)

# Executable names
TARGET=assignBarcodes
DEMUX_TARGET=demux_fastq

# demux_fastq sources/objects
DEMUX_SRCS=demux_fastq.c io.c utils.c globals.c barcode_match.c
DEMUX_OBJS=$(DEMUX_SRCS:.c=.o)

DEMUX_BAM_TARGET=demux_bam
DEMUX_BAM_SRCS=demux_bam.c utils.c globals.c barcode_match.c io.c memory.c
DEMUX_BAM_OBJS=$(DEMUX_BAM_SRCS:.c=.o)

.PHONY: all clean

all: $(TARGET) $(DEMUX_TARGET) $(DEMUX_BAM_TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(DEMUX_TARGET): $(DEMUX_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(DEMUX_BAM_TARGET): $(DEMUX_BAM_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

# The implicit rule will now use VPATH to find the .c files in src/
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET) $(DEMUX_OBJS) $(DEMUX_TARGET) $(DEMUX_BAM_OBJS) $(DEMUX_BAM_TARGET)
