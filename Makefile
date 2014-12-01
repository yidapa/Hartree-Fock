ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif


#$(CXX) -I[libint_prefix]/include -I[libint_prefix]/include/libint2 -L[libint_prefix]/lib -o test_eri ./test_eri.cc -lint2

EIGENDIR = /home/nick/opt/eigen-eigen-1306d75b4a21/
LIBINTDIR = /home/nick/opt/libint/

# include headers the object include directory
CPPFLAGS += -std=c++0x -I./ -I${LIBINTDIR}/include/ -I${LIBINTDIR}/include/libint2  -L${LIBINTDIR}/lib/ -I${EIGENDIR}/
CXXPOSTFLAGS = -lint2

CXX = g++
EXEC = hf
CXXTESTSRC = $(EXEC).cc
CXXTESTOBJ = $(CXXTESTSRC:%.cc=%.o)


$(EXEC): $(CXXTESTOBJ) 
	$(CXX) ${CPPFLAGS} -o $@  $^ ${CXXPOSTFLAGS}

clean::
	-rm -rf $(TEST) *.o *.d

distclean:: realclean

realclean:: clean

targetclean:: clean
