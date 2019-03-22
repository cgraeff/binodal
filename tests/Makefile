SHELL := /bin/bash # Use bash as shell
TARGET = binodal

all:
	@echo "[Running tests...]"
	@cd ..; ./$(TARGET) -a $(ARGS)
	@echo "[done.]"
graph:
	@echo "[Plotting tests ...]"
	for dir in `echo */`; do \
		cd "$$dir"; \
		echo "$$dir"; \
		if [ -e gnuplot.gpi ]; \
		then \
			gnuplot gnuplot.gpi; \
		fi; \
		if [ -e script.sh ]; \
		then \
			bash script.sh; \
		fi; \
		if [ -e script.py ]; \
		then \
			python3 script.py; \
		fi; \
		cd ../; \
	done;
	@echo "[done.]"
clean:
	@echo "[Cleaning...]"
	@-rm -f $(TARGET)
	@cd binodal_point_graph/data; rm -rf */ 
	@cd binodal_point_graph/graph; rm -rf */
	@find . -name "*.dat" -type f -delete
	@find . -name "*.log" -type f -delete
	@find . -name "*.png" -type f -delete
	@echo "[done.]"
