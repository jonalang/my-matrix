CFLAGS = -Wall -g 
mylibdir = ../my-lib
mysortdir = ../my-sort
Ldirs = -L$(mylibdir) -L$(mysortdir)
headerdirs = -iquote$(mylibdir) -iquote$(mysortdir)
linkerdirs = -Wl,-rpath,$(mylibdir),-rpath,$(mysortdir)

# Option explanations: 
# -L tells gcc to add a directory to library search path
# -Wl following are comma separated instructions to the linker
# -rpath Add a directory to the RUNTIME library search path
# -c create object files, do not link
# fPIC create position independent code (characteristic of a shared library)
# -lmylib a library to include when building
# -shared create a shared library
all: 
	gcc $(CFLAGS)  -c $(headerdirs) -fPIC  mymatrix.c
	gcc -shared $(Ldirs) $(linkerdirs) -o libmymatrix.so mymatrix.o -lmylib -lmysort

clean: rm libmymatrix*


