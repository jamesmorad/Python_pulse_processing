CC=g++


CFLAGS=-Wall -g -std=c++0x -O3 -I$(BOOST_INC) -I$(PYTHON_INC)
LDFLAGS=-Wall -g -std=c++0x -O3 -L$(BOOST_LIB) -L$(PYTHON_LIB) -lboost_python -lpython2.7


OBJECTS  = obj/PulseFinder.o
OBJECTS += obj/SplitGaussian.o

SHAREDLIB=PulseFinder_BlackBox

################################################################################

all: $(SHAREDLIB)

clean:
	rm -f $(OBJECTS) $(SHAREDLIB).so


################################################################################
$(OBJECTS): obj/%.o : %.cxx %.h
obj/%.o: %.cxx %.h
	$(CC) $< -c $(CFLAGS) -o $@

$(SHAREDLIB): %: %.cpp $(OBJECTS)
	$(CC) $< -c $(CFLAGS) -o obj/$@.o
	$(CC) -shared -dynamiclib obj/$@.o $(OBJECTS) $(LDFLAGS) -o $@.so
