CC = g++

CFLAGS = `root-config --cflags`
LFLAGS = `root-config --libs` -lMinuit2
LFLAGS_DATA = -Levent -l:Dict.so

SOURCES = $(subst ./src/,,$(shell find ./src -mindepth 2 -name '*.cc' -type f))

HEADER_DIR = ./include
HEADERS = $(subst ./src/,$(HEADER_DIR)/,$(shell find ./src -name '*.hh' -type f))

OBJ_DIR = ./build
OBJECTS = $(SOURCES:%.cc=$(OBJ_DIR)/%.o)

INCLUDES = -Iinclude

EXE_DIR = ./bin
EXE_SOURCES = $(subst ./src/,,$(shell find ./src -maxdepth 1 -name '*.cc' -type f))
EXE = $(EXE_SOURCES:%.cc=$(EXE_DIR)/%)

.PHONY: all includes

all: build_dirs $(HEADERS) $(OBJECTS) $(EXE)

clean:
	-$(RM) build/*.o $(EXE)
	-$(RM) event/Dict.*
	-$(RM) include/*

build_dirs:
	@mkdir -p build
	@mkdir -p include
	@mkdir -p include/event
	@mkdir -p bin

event/Dict.cc:
	rootcling -f event/Dict.cc -c event/Hit.h event/LinkDef.h
	cp event/*.h include/event/.

event/Dict.so: event/Dict.cc
	$(CC) -shared -fPIC -o $@ $< $(CFLAGS) $(INCLUDES) $(LFLAGS) -I`root-config --incdir`

$(HEADER_DIR)/%.hh: ./src/%.hh
	cp $^ $@

$(OBJ_DIR)/%.o: ./src/%.cc $(HEADERS) event/Dict.so
	$(CC) -c -o $@ $< -g $(CFLAGS) $(INCLUDES)

$(EXE_DIR)/%: $(OBJECTS) $(HEADERS) src/%.cc event/Dict.so
	$(CC) -o $@ $^ -g $(CFLAGS) $(INCLUDES) $(LFLAGS) $(LFLAGS_DATA)
