LIBS = lib  
BINS = mkflat  mkpsf  mvdethgh  preproc  mkdetimg  mvdetsum  mvdethghxyt

SUBDIR = $(LIBS) $(BINS)  doc  script

ALL = $(LIBS) $(BINS)

all:
	for dir in $(ALL); do \
	(cd $$dir; ${MAKE} all); \
	done
test:
	for dir in $(LIBS); do \
	(cd $$dir; ${MAKE} test); \
	done

clean:
	for dir in $(SUBDIR); do \
	(cd $$dir; ${MAKE} clean); \
	done

cleaner:
	-rm -f *~
	for dir in $(SUBDIR); do \
	(cd $$dir; ${MAKE} cleaner); \
	done
