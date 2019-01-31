ifndef ROOTDIR
ROOTDIR=.
endif

REALMAKEFILE=../Makefile.2

default: lib

all: lib

lib: objdir serialobjs
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) libs)

clean:
	rm libwannier.a
	cd $(ROOTDIR) && rm -f *~
	cd $(ROOTDIR) && rm -f src/*~
	@( cd $(ROOTDIR) && if [ -d src/obj ] ; \
		then cd src/obj && \
		$(MAKE) -f $(REALMAKEFILE) clean && \
		cd ../ && rm -rf obj ; \
	fi )
	@( cd $(ROOTDIR) && if [ -d src/objp ] ; \
		then cd src/objp && \
		$(MAKE) -f $(REALMAKEFILE) clean && \
		cd ../ && rm -rf objp ; \
	fi )


objdir:
	@( cd $(ROOTDIR) && if [ ! -d src/obj ] ; \
	    then mkdir src/obj ; \
	fi ) ;


serialobjs: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) serialobjs)
