SHELL := /bin/bash # Use bash as shell
TARGET = binodal

# List sets for multirun
HADRON_MULTIRUN_SETS = eNJL2mSigmaRho1 eNJL3SigmaRho1

QUARK_MULTIRUN_SETS = Buballa_1 \
		      PCP-0.1

.PHONY: all debug run graph tests tgraph clean

all:
	@echo "[Compiling...]"
	@cd src; make
	@echo "[done.]"
debug:
	@echo "[Compiling with debug symbols...]"
	@cd src; make debug
	@echo "[done]"
run:
	@./$(TARGET) -d $(ARGS)
drun:
	@./$(TARGET) -d -! $(ARGS)
arun:
	@./$(TARGET) -d -e $(ARGS)
graph:
	@echo "[Plotting...]"
	@cd output; \
	gnuplot gnuplot.gpi; \
	cd ..;
	@echo "[done.]"
multirun:
	@echo "[Running for multiple parameterizations...]"
	@for hadron_set in $(HADRON_MULTIRUN_SETS); do \
		for quark_set in $(QUARK_MULTIRUN_SETS); do \
			echo ""; \
			echo "[$$hadron_set-$$quark_set]"; \
			./$(TARGET) -d -h "$$hadron_set" -q "$$quark_set" $(ARGS); \
			if [ -d multioutput/"$$quark_set-$$hadron_set" ]; then \
				rm -r multioutput/"$$quark_set-$$hadron_set"; \
			fi; \
			cp -r output multioutput/"$$quark_set-$$hadron_set"; \
		done; \
	done
	@echo "[done.]"
multidrun:
	@echo "[Running for multiple parameterizations...]"
	@for hadron_set in $(HADRON_MULTIRUN_SETS); do \
		for quark_set in $(QUARK_MULTIRUN_SETS); do \
			echo ""; \
			echo "[$$hadron_set-$$quark_set]"; \
			./$(TARGET) -d  -e -h "$$hadron_set" -q "$$quark_set" $(ARGS); \
			if [ -d multioutput/"$$quark_set-$$hadron_set" ]; then \
				rm -r multioutput/"$$quark_set-$$hadron_set"; \
			fi; \
			cp -r output multioutput/"$$quark_set-$$hadron_set"; \
		done; \
	done
	@echo "[done.]"
mgraph:
	@echo "[Plotting for multiple parameterizations...]"
	@cd multioutput; \
	for dir in `echo */`; do \
		echo "$$dir"; \
		cd "$$dir"; \
		if [ -e gnuplot.gpi ]; \
		then \
			gnuplot gnuplot.gpi; \
		fi; \
		cd ..; \
	done;
	@echo "[done.]"
tests:
	@cd tests; make
tgraph:
	@cd tests; make graph
clean:
	@echo "[Cleaning...]"
	@-rm -f $(TARGET)
	@cd src; make clean
	@cd tests; make clean
	@find output -name "*.dat" -type f -delete
	@find output -name "*.log" -type f -delete
	@find output -name "*.png" -type f -delete
	@cd multioutput; rm -rf */
	@echo "[done.]"
