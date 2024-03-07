CC = g++
CFLAGS = -Wall -Wextra -Isrc/include `pkg-config opencv4 --cflags`
DEBUG = 1
LDFLAGS = `pkg-config opencv4 --libs`
EXEC = slic
SRCDIR = ./src
OBJDIR = ./obj
VPATH = $(SRCDIR)
OBJ_NAMES = SLIC.o 
OBJS = $(addprefix $(OBJDIR)/,$(OBJ_NAMES))

ifeq ($(DEBUG), 1)
CFLAGS += -g -DDEBUG
else
CFLAGS += -O2
endif

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
    
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) -o $@ -c $(CFLAGS) $^
	
test: $(OBJDIR)/test.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/test.o: $(SRCDIR)/testCountNeighbours.cpp
	$(CC) -o $@ -c $(CFLAGS) $^
	
.PHONY: clean cleanest init

clean:
	rm $(EXEC)
cleanest: clean
	rm $(OBJDIR)/*.o
init:
	mkdir -p $(OBJDIR) $(SRCDIR)
    
