.PHONY: depend clean

all:    dirs torque
	@echo  executables stored in bin/

dirs:
	mkdir -p bin

torque: force_look
	cd src/ ; $(MAKE) all ; mv torque ../bin

clean:
	(if [ -d bin ]; then /bin/rm -r bin ; fi)
	(cd src/ ; $(MAKE) clean)

force_look:
	true
