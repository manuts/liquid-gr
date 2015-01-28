# paths
srcdir			= .
prefix			= /target

# programs
CC				:= gcc
CXX				:= g++
AR				:= ar
RANLIB			:= ranlib
MV				:= mv -f
RM				:= rm -f

# flags
INCLUDE_CFLAGS	= -I./ -I./include/ -I/target/include/
CONFIG_CFLAGS	= -g -O2 -march=core2
CFLAGS			+= $(INCLUDE_CFLAGS) -Wall -fPIC $(CONFIG_CFLAGS)
LDFLAGS			= -L/target/lib -L/target/lib64
LDFLAGS			+= -lm -lc -lliquid
ARFLAGS			= r

library_src		:=							\
  src/test/copy.cc						\
  src/mimo/framegen.cc				\
  src/mimo/framesync.cc				\
  src/alamouti/framegen.cc		\
  src/alamouti/framesync.cc		\
  src/alamouti/channel_estimator.cc		\
  src/alamouti/delay.cc				\
	src/simo/framegen.cc				\
  src/simo/framesync.cc				\
  src/math/liquid_math.cc			\

library_hdr		:=							\
  include/test.h							\
  include/mimo.h							\
  include/alamouti.h					\
	include/simo.h							\
  include/liquid_math.h				\

library_objs	= $(patsubst %.cc, %.o, $(library_src))

all				: libliquidgr.a libliquidgr.so

libliquidgr.a	: $(library_objs)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

libliquidgr.so	: $(library_objs)
	$(CXX) -shared -Xlinker -soname=$@ -o $@ -Wl,-whole-archive $^ -Wl,-no-whole-archive $(LDFLAGS)

$(library_objs)	: %.o : %.cc
	$(CXX) $(CFLAGS) -c $< -o $@

clean			:
	$(RM) $(library_objs)
	$(RM) libliquidgr.a
	$(RM) libliquidgr.so

install			:
	@echo "installing..."
	mkdir -p $(prefix)/lib
	install -m 644 -p libliquidgr.so libliquidgr.a $(prefix)/lib
	mkdir -p $(prefix)/include
	mkdir -p $(prefix)/include/liquid
	install -m 644 -p $(library_hdr) $(prefix)/include/liquid
	@echo ""
	@echo "---------------------------------------------------------"
	@echo "  liquid-gr was successfully installed.     "
	@echo ""
	@echo "  On some machines (e.g. Linux) you should rebind your"
	@echo "  libraries by running 'ldconfig' to make the shared"
	@echo "  object available.  You might also need to modify your"
	@echo "  LD_LIBRARY_PATH environment variable to include the"
	@echo "  directory $(prefix)"
	@echo "---------------------------------------------------------"
	@echo ""

## 
## TARGET : uninstall - uninstalls the libraries and header files in the host system
##

uninstall:
	@echo "uninstalling..."
	$(RM) $(addprefix $(prefix)/include/liquid/, $(library_hdr))
	$(RM) $(prefix)/lib/libliquidgr.a
	$(RM) $(prefix)/lib/libliquidgr.so
	@echo "done."
