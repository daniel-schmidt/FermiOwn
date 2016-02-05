EIGENDIR=/work/libs/eigen3.2.7/
CXXFLAGS =	-std=c++11 -O2 -g -Wall -fmessage-length=0 -I$(EIGENDIR)

OBJS =		Lattice.o Action.o

LIBS =

TARGET =	FermiOwn

$(TARGET):	$(OBJS) FermiOwn.o
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

test:	$(OBJS) FieldScalar_test.o
	$(CXX) FieldScalar_test.o $(OBJS) -o test $(CXXFLAGS)
clean:
	rm -f $(OBJS) $(TARGET) test FieldScalar_test.o FermiOwn.o