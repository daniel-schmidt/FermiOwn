EIGENDIR=/work/libs/eigen3.2.7/
CXXFLAGS =	-std=c++11 -O2 -g -Wall -fmessage-length=0 -I$(EIGENDIR)

OBJS =		Lattice.o Action.o
INCL = Lattice.h FieldScalar.h Action.h
LIBS =

TARGET =	FermiOwn

%.o: %.cpp $(INCL)
	$(CXX) -c $< $(CXXFLAGS)

$(TARGET):	$(OBJS) FermiOwn.o
	$(CXX) FermiOwn.o $(OBJS) $(LIBS) -o $(TARGET)

all:	$(TARGET)

test:	$(OBJS) FieldScalar_test.o
	$(CXX) FieldScalar_test.o $(OBJS) -o test $(CXXFLAGS)
clean:
	rm -f $(OBJS) $(TARGET) test FieldScalar_test.o FermiOwn.o