CXX = g++
CPPFLAGS = -std=c++2a
CFLAGS = 

HPC_sim: run_sim.o Quarantine_outputs.o Quarantine_breach_model.o 
	$(MAKE) -C ../vaccine_abm
	$(CXX) $(CPPFLAGS) -o runsim run_sim.o Quarantine_outputs.o Quarantine_breach_model.o ../vaccine_abm/*.o

main_GI: main_GISITOST.o
	$(MAKE) -C ../vaccine_abm
	$(CXX) $(CPPFLAGS) -o Run_GI main_GISITOST.o ../vaccine_abm/*.o

main: main.o 
	$(MAKE) -C ../vaccine_abm
	$(CXX) $(CPPFLAGS) -o Run main.o ../vaccine_abm/*.o
	
run_sim.o: run_sim.cpp
Quarantine_breach_model.o: Quarantine_breach_model.cpp
Quarantine_outputs.o: Quarantine_outputs.cpp
main.o: main.cpp
main_GISITOST.o: main_GISITOST.cpp

clean: 
	rm Run  Run_GI runsim *.o
	$(MAKE) clean -C ../vaccine_abm
