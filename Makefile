CFLAGS = -Wall -g 
mylibdir = ../my-lib

# Option explanations: 
# -L tells gcc to add a directory to library search path
# -Wl following are comma separated instructions to the linker
# -rpath Add a directory to the RUNTIME library search path
# -c create object files, do not link
# fPIC create position independent code (characteristic of a shared library)
# -lmylib a library to include when building
# -shared create a shared library
all: 
	gcc $(CFLAGS) -L$(mylibdir) -Wl,-rpath,$(mylibdir) -iquote$(mylibdir) -c -fPIC  mymatrix.c -lmylib
	gcc -shared -o libmymatrix.so mymatrix.o

clean: rm libmymatrix*


