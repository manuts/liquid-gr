# paths
srcdir			= .

# programs
CC				:= gcc
CXX				:= g++
MV				:= mv -f
RM				:= rm -f

# flags
INCLUDE_CFLAGS	= -I./ -I/target/include -I./include
CFLAGS			= 	$(INCLUDE_CFLAGS) -g -O2 -Wall -fPIC
LDFLAGS			= 	-L/target/lib
LDFLAGS			+= 	-luhd -lliquid -lm -lc -lboost_system-mt -lboost_thread-mt -lliquidgr
LDFLAGS			+= 	-lboost_program_options-mt -lpthread -lvolk -lfftw3f
env				 	= 	LD_LIBRARY_PATH="/target/lib"

sup_src			:=									\

sup_obj					= $(patsubst %.cc, %.o, $(sup_src))

all							: siso.exe

siso.exe				: $(sup_obj) siso.o
	$(env) $(CXX) $^ -o $@ $(LDFLAGS)

siso.o					: siso.cc
	$(CXX) $(CFLAGS) -c $< -o $@

$(sup_obj)			: %.o :	%.cc
	$(CXX) $(CFLAGS) -c $< -o $@

clean						:
	$(RM) siso.exe siso.o
	$(RM) $(sup_obj)
