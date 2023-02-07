SRC_lib = external/primme 
SRC_exe = \
		cygutz/GINIT  \
		cygutz/GBANDS \
		cygutz/GERR  \
		cygutz/GFINAL \
		cygutz/GSPCI

all:
	mkdir -p lib include
	for dir in ${SRC_lib} - ; do \
    if [ "$$dir" = "-" ] ; then break ; fi ; \
    (cd $$dir && ${MAKE} ; ) ; done
	for dir in ${SRC_exe} - ; do \
    if [ "$$dir" = "-" ] ; then break ; fi ; \
    (cd $$dir/real && ${MAKE} && cd ../cmplx && \
	${MAKE}; ) ; done

real:
	for dir in ${SRC_exe} - ; do \
    if [ "$$dir" = "-" ] ; then break ; fi ; \
    (cd $$dir/real && ${MAKE} && cd ../cmplx && \
	${MAKE}; ) ; done

clean:
	for dir in ${SRC_lib} ${SRC_exe} - ; do \
    if [ "$$dir" = "-" ] ; then break ; fi ; \
    (cd $$dir/real && ${MAKE} clean_all && cd ../cmplx && \
	${MAKE} clean_all; ) ; done

install :
	for dir in ${SRC_exe} - ; do \
    if [ "$$dir" = "-" ] ; then break ; fi ; \
    (cd $$dir/real && ${MAKE} install && \
	cd ../cmplx && ${MAKE} install; ) ; done

