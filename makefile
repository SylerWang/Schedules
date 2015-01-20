### Makefile for my custom Radiance-based tools
### This assumes that in there is a "src" directory, with the "ray" Radiance's directory and
### the "view3d-3.5" view directory.

# You have to unpack the libraries from Radiance (lrtrad). I do it by running ./makeall install.

## Germán Molina

####  Mac OS X 
MACH= -DBSD -DNOSTEREO -Dfreebsd #-I/usr/X11R6/include -L/usr/X11R6/lib
OPT= -O0 

 
SPECIAL= ogl
CC = cc
RAD_COMP = $(CC) $(CFLAGS)  -I./src/ray/src/common -I./src/ray/src/rt -L./src/ray/src/lib 

####  LINUX
#mach= -Dlinux -D_FILE_OFFSET_BITS=64 -L/usr/X11R6/lib -I/usr/include/X11 -DNOSTEREO
#opt= -O0
#arch= IBMPC
#compat= erf.o
#extras= CC=gcc


CFLAGS = $(OPT) $(MACH) 
MLIB = -lm

programs = mkSchedule uSchedule epw2daymtx 

all:	$(programs)

redo: clean
	rm -f $(programs)

clean: 
	rm -f *.o *.c *.h *.gch




## uSchedule
###############

uSchedule:	uSchedule.o ./src/ray/src/gen/sun.o
	$(RAD_COMP) -llua -o uSchedule uSchedule.o ./src/ray/src/gen/sun.o -lrtrad $(MLIB)

uSchdule_depend = ./src/ray/src/common/standard.h \
./src/ray/src/common/rtmisc.h ./src/ray/src/common/rtio.h \
./src/ray/src/common/rtmath.h ./src/ray/src/common/mat4.h ./src/ray/src/common/fvect.h \
./src/ray/src/common/rterror.h ./src/ray/src/common/platform.h ./src/ray/src/common/paths.h \
./src/ray/src/common/color.h ./src/ray/src/common/bsdf.h ./src/ray/src/common/bsdf_m.h ./src/ray/src/common/resolu.h \
/usr/local/include/lua.h /usr/local/include/luaconf.h /usr/local/include/lualib.h \
/usr/local/include/lauxlib.h

uSchedule.o:	$(uSchedule_depend)
	$(RAD_COMP) -Wall -c $(uSchedule_depend) ./src/uSchedule.c



## mkSchedule
###############

mkSchedule:	mkSchedule.o ./src/ray/src/gen/sun.o
	$(RAD_COMP) -llua -o mkSchedule mkSchedule.o ./src/ray/src/gen/sun.o -lrtrad $(MLIB)

mkSchdule_depend = ./src/ray/src/common/standard.h \
./src/ray/src/common/rtmisc.h ./src/ray/src/common/rtio.h \
./src/ray/src/common/rtmath.h ./src/ray/src/common/mat4.h ./src/ray/src/common/fvect.h \
./src/ray/src/common/rterror.h ./src/ray/src/common/platform.h ./src/ray/src/common/paths.h \
./src/ray/src/common/color.h ./src/ray/src/common/bsdf.h ./src/ray/src/common/bsdf_m.h ./src/ray/src/common/resolu.h \
/usr/local/include/lua.h /usr/local/include/luaconf.h /usr/local/include/lualib.h \
/usr/local/include/lauxlib.h

mkSchedule.o:	$(mkSchedule_depend)
	$(RAD_COMP) -Wall -c $(mkSchedule_depend) ./src/mkSchedule.c

## epw2daymtx
###############

epw2daymtx:	epw2daymtx.o ./src/ray/src/gen/sun.o
	$(RAD_COMP) -o epw2daymtx epw2daymtx.o ./src/ray/src/gen/sun.o -lrtrad $(MLIB)	

epw2daymtx.o: ./src/ray/src/common/rtmath.h ./src/ray/src/common/color.h
	$(RAD_COMP) -c ./src/ray/src/common/rtmath.h ./src/ray/src/common/color.h ./src/epw2daymtx.c	





