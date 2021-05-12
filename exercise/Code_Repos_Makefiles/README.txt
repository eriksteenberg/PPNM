THe git part was done as you can see from the fact that you can read this.
Question "make":
1.)
$@ name of the target file.
$< name of the the first requirement/prerequisite
$^ the names of all the requirements(prerequisites)
CFALGS is extra stuff for the C compiler.
LDLFLAGS stuff to give to the compiler when it has to invok stuff.
LDLIBS is similar to LDFLAGS but for libaries instead.
2.)
echo CFLAGS
CFLAGS
echo FFLAGS
FFLAGS
echo FFLAGS
FFLAGS
echo -Wall -Ofast -std=clx
-Wall -Ofast -std=clx
3.)
As far as I was able to determine, only number 6 would not be able to link. 
