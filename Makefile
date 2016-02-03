CXXFLAGS =	-std=c++11 -O2 -g -Wall -fmessage-length=0

OBJS =		FermiOwn.o Lattice.o

LIBS =

TARGET =	FermiOwn

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
