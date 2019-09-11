################################################
## RVAR makefile ###############################
################################################

# Include machine specific variables
HOSTNAME = $(shell hostname)
-include opts/Makefile_$(HOSTNAME) 

# Compiler
CXX = g++

# Figure out all the objects
raw_OBJS = lin_maths.o observation_function.o rvar.o location.o txt_files.o geo_observations.o readings.o wang_p_model_correction.o db_tables.o accessors.o settings_interpreter.o climate_data_reader.o main.o
OBJS = $(addprefix $(OD),$(raw_OBJS)) 

# Directories that we'll use
OD = out/
SD = src/
LD = lib/

# Shared compiler and linker flags
CXXSFLAGS=-fPIC
CXXSLFLAGS=-shared

# Final executable name
RUN = rmin

# Standard and shared compilation
$(OD)%.o : $(SD)%.cpp
	@mkdir -p $(OD)
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

$(OD)%.o.sh : $(SD)%.cpp
	@mkdir -p $(OD)
	$(CXX) $(CXXSFLAGS) $(CXXFLAGS)  -c $< -o $@


.PHONY: all clean run configure install lib uninstall

# Make everything
all : $(RUN) #lib

# Make main executable
$(RUN) : $(OBJS) $(MAKEFILE) 
	$(CXX) $(OBJS) $(CXXLFLAGS) -o $(RUN) 

run : 
	./$(RUN)

# Make the library
lib : $(LD)libpmodel_corrector.so

# Make the machine specific files
configure:
	@if [ ! -f "opts/Makefile_$(HOSTNAME)" ]; then\
		cp opts/Makefile_hostname opts/Makefile_$(HOSTNAME);\
	fi	

	@if [ ! -f "config/$(HOSTNAME).cfg" ]; then\
		cp config/example.cfg config/$(HOSTNAME).cfg;\
	fi	

# Clean out the executable and object files
clean :
	rm rmin out/*

# Install the library
install : $(LD)libpmodel_corrector.so
	sudo ln -s `pwd`/$< /usr/local/lib/libpmodel_corrector.so

uninstall :
	sudo rm /usr/local/lib/libpmodel_corrector.so

# Perform memory checking and profiling
mem:
	valgrind --tool=memcheck --leak-check=yes ./$(RUN) > /tmp/mem.txt 2>&1

prof:
	valgrind --tool=callgrind --callgrind-out-file=/tmp/prof.txt ./$(RUN)

prof_view:
	kcachegrind /tmp/prof.txt 


## Dependencies #############################################

# Regular dependencies 

$(OD)lin_maths.o : $(SD)lin_maths.cpp $(SD)lin_maths.h $(SD)lin_maths_types.h $(SD)lin_maths_types.cpp
$(OD)observation_function.o : $(SD)observation_function.cpp $(SD)observation_function.h 
$(OD)rvar.o : $(SD)rvar.cpp $(SD)rvar.h $(SD)lin_maths.h $(SD)observation_function.h 
$(OD)location.o : $(SD)location.cpp $(SD)location.h
$(OD)asa047.o : $(SD)asa047.cpp $(SD)asa047.h
$(OD)readings.o : $(SD)readings.cpp $(SD)readings.h $(SD)lin_maths.h $(SD)location.h
$(OD)settings_interpreter.o : $(SD)settings_interpreter.cpp $(SD)settings_interpreter.h
$(OD)climate_data_reader.o : $(SD)climate_data_reader.cpp $(SD)climate_data_reader.h $(SD)geo_observations.h $(SD)accessors.h $(SD)lin_maths.h $(SD)settings_interpreter.h
$(OD)wang_p_model_correction.o : $(SD)wang_p_model_correction.cpp $(SD)wang_p_model_correction.h $(SD)orbital.h $(SD)location.h
$(OD)geo_observations.o : $(SD)geo_observations.cpp $(SD)geo_observations.h $(SD)wang_p_model_correction.h $(SD)readings.h $(SD)observation_function.h $(SD)lin_maths.h $(SD)location.h
$(OD)accessor.o : $(SD)accessor.cpp $(SD)accessor.h $(SD)db_tables.h $(SD)txt_files.h $(SD)readings.h $(SD)db_tables.h $(SD)geo_observations.h
$(OD)db_tables.o : $(SD)db_tables.cpp $(SD)db_tables.h $(SD)lin_maths.h $(SD)location.h
$(OD)txt_files.o : $(SD)txt_files.cpp $(SD)txt_files.h $(SD)readings.h
$(OD)main.o : $(SD)main.cpp $(SD)main.h $(SD)climate_data_reader.h	

# Library dependencies

$(LD)libpmodel_corrector.so : $(OD)wang_p_model_correction.o.sh $(OD)lin_maths.o.sh $(OD)location.o.sh
	$(CXX) $(CXXSLFLAGS) $(CXXLFLAGS) -o $@ $^
	
$(OD)location.o.sh : $(SD)location.cpp $(SD)location.h 
$(OD)lin_maths.o.sh : $(SD)lin_maths.cpp $(SD)lin_maths.h $(SD)location.h
$(OD)wang_p_model_correction.o.sh : $(SD)wang_p_model_correction.cpp $(SD)wang_p_model_correction.h $(SD)orbital.h $(SD)lin_maths.h $(SD)location.h
