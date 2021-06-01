include arch.mk

default: all

all: env com_wannier90 com_comgw com_comwann

env:
	mkdir -p bin
	if ! grep -Fq "COMRISB_BIN=" ~/.bashrc; then \
		echo -e '\nexport COMRISB_BIN=${PWD}/bin' >> ~/.bashrc; \
	fi
com_wannier90:  	
	cd wannier90_2.1 && $(MAKE) all && cd ../
com_comgw:
	cd gw && $(MAKE) && cd ../
com_comwann:
	cd ComWann && $(MAKE) && cd ../  


clean: clean_wannier90 clean_comgw clean_comwann 

clean_comgw:
	cd gw && $(MAKE) clean && cd ../
clean_comwann:
	cd ComWann && $(MAKE) clean && cd ../
clean_wannier90:  	
	cd wannier90_2.1 && $(MAKE) clean && cd ../
clean_bin:
	rm -rf ./bin/*
