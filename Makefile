.PHONY:
all: build

.PHONY:
configure:
	./liblast/configure
	./libbbrc/configure
	./fminer/configure

.PHONY: 
build:
	make -C liblast/
	make -C libbbrc/
	make -C fminer/
	@echo
	@echo "Please start 'fminer' in fminer/ subdirectory in case of successful compilation!"
	@echo

.PHONY:
clean:
	make -C liblast/ clean
	make -C libbbrc/ clean
	make -C fminer/ clean

