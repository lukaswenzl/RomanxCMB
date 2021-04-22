include ${COSMOSIS_SRC_DIR}/config/compilers.mk


OBJ = \
	precision.o \
	constants.o \
	physics.o \
	basic_operations.o \
	array_operations.o \
	random_numbers.o \
	file_info.o \
	table_integer.o \
	special_functions.o \
	interpolate.o \
	solve_equations.o \
	string_operations.o \
	calculus_table.o \
	camb_stuff.o \
	sorting.o \
	statistics.o \
	calculus.o \
	minimization.o \
	multidark_stuff.o \
	cosmology_functions.o \
	hmx.o 

# Build directory
#BUILD_DIR = build
# Fortran library directory
MOD_DIR = library/src
# This is a rule to make directories
make_dirs = @mkdir -p $(@D)

#OBJ = $(addprefix $(BUILD_DIR)/,$(_OBJ))

# Rule to make object files
%.o: $(MOD_DIR)/%.f90
	$(make_dirs)
	$(FC) -c -o $@ $< $(LDFLAGS) $(FFLAGS)

all: mead_interface.so

mead_interface.so: $(OBJ) mead_interface.f90
	#$(FC) $(FFLAGS) -shared -o $@ $+ -L. -lhmcode $(LDFLAGS) -lcosmosis_fortran -lcosmosis
	$(FC) $(FFLAGS) -shared -o $@ $+ -L. $(LDFLAGS) -lcosmosis_fortran -lcosmosis



#test: test.f90 libhmcode.so
#	$(FC) $(FFLAGS) -o $@ $< -L. -lhmcode $(LDFLAGS)


clean:
	#rm -f hmcode
	rm -f libhmcode.so
	rm -f hmcode.o
	rm -f mead_interface.so
	#rm -rf test
	rm -rf power.dat
	rm -rf *.dSYM/
	rm -rf *.mod
	rm -rf *.o

#libhmcode.so: HMcode.f90
#	$(FC) $(FFLAGS) -c  $+ $(LDFLAGS)
#	$(FC) -shared -o $(CURDIR)/$@ HMcode.o