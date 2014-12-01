ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif



EIGENDIR = /home/nick/opt/eigen-eigen-1306d75b4a21/

# include headers the object include directory
CPPFLAGS += -std=gnu++0x -I./  -I${EIGENDIR}/

CXX = g++
EXEC = hf
CXXTESTSRC = $(EXEC).cc
CXXTESTOBJ = $(CXXTESTSRC:%.cc=%.o)


$(EXEC): $(CXXTESTOBJ) 
	$(CXX) ${CPPFLAGS} -o $@  $^ 

clean::
	-rm -rf $(TEST) *.o *.d

distclean:: realclean

realclean:: clean

targetclean:: clean
