include ../../mitool/Makefile.in
include ../../mitool/milib/Makefile.in
include ../../mitool/miblaslib/Makefile.in
include ../../mitool/mifitslib/Makefile.in
include ../../mitool/mirootlib/Makefile.in
include ../../mitool/mifunclib/Makefile.in
include ../Makefile.in

INCLUDE = $(MVDETHGHLIB_INC) \
	$(MIFUNCLIB_INC) $(MIROOTLIB_INC) $(MIFITSLIB_INC) $(MILIB_INC) \
	$(CFITSIO_INC) $(ROOT_INC) $(MT_INC)
LDLIBS  = $(MVDETHGHLIB_LIB) \
	$(MIFUNCLIB_LIB) $(MIROOTLIB_LIB) $(MIFITSLIB_LIB) $(MILIB_LIB) \
	$(CFITSIO_LIB) $(ROOT_ALLLIBS) $(MT_LIB) $(LIB)

TARGET = mkflatxyt

all: $(TARGET) 

mkflatxyt: mkflatxyt.o  arg_mkflatxyt.o
	$(CXXN) $@.o arg_$@.o ${LDLIBS} -o $@

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

install: all
	if [ ! -d $(INST_BIN_DIR) ]; then mkdir $(INST_BIN_DIR); fi
	install ${TARGET} ${INST_BIN_DIR}

clean:
	-rm -f *.o

cleaner: clean
	-rm -f ${TARGET} *% *~ core
