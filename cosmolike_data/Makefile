
# If you already have your own Makefile you can 
# replace all of this, but you need to keep this line
# at the top:
include ${COSMOSIS_SRC_DIR}/config/compilers.mk


#For python, leave this as is.  For C/C++/Fortran, 
#comment the first line and uncomment the second:
all: 
#all load_cosmolike_covariance.so

# ... and then provide the command to compile it, perhaps
# like this for a simple example in C:
#load_cosmolike_covariance.so: my_files.c
#	$(CC) -shared -o load_cosmolike_covariance.so my_files.c


#Replace this to put in some kind of test command for
#your module, if you have one.
test:
	@echo "Alas, 6x2pt_test/load_cosmolike_covariance module has no tests"

#Add anything else to clean here
clean:
	rm -f *.o *.mod

