include Makefile.in

install:
	#--------------------------------------------------------
	mkdir -p ${DESTDIR} 
	# Build Gutzwiller solver 
	cd src && make && make install && cd ..
	cp -r ./pygtool/*  ${DESTDIR}

pygtools:
	cp -r ./pygtool/*  ${DESTDIR}

clean:
	cd src && make clean && cd ..
	rm -f ${DESTDIR}/Cy* ${DESTDIR}/exe_*
