include ../../mitool/Makefile.in
include ../../mitool/milib/Makefile.in
include ../../mitool/miblaslib/Makefile.in
include ../../mitool/mifitslib/Makefile.in
include ../../mitool/mirootlib/Makefile.in
include ../../mitool/mifunclib/Makefile.in

INCLUDE := $(MIFUNCLIB_INC) $(MIROOTLIB_INC) $(MIFITSLIB_INC) $(MILIB_INC) \
	$(CFITSIO_INC) $(ROOT_INC) $(MT_INC)

TARGET := libmvdethghlib.a

OBJ := mvdethghlib.o

all: $(TARGET) 

libmvdethghlib.a: ${OBJ}
	$(AR) -r $@ $(OBJ)

libmvdethghlib.so: ${OBJ}
	$(CC) -shared -g -O2 -o $@ ${OBJ}

.SUFFIXES:
.SUFFIXES: .o .c .f .F .cc

.c.o:
	$(CCN) -c $<

.F.o : 
	$(FCN) -c $<

.f.o:
	$(FCN) -c $<

.cc.o:
	${CXXN} -c $<

clean:
	-rm -f *.o *.a *.so

cleaner: clean
	-rm -f ${TARGET} *% *~ core
