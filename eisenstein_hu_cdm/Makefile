include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: eisenstein_hu_cdm_module.so

test:
	@echo "Alas, CRL ST has no tests"

clean:
	rm -f eisenstein_hu_cdm_module.so *.o *.mod
	rm -rf eisenstein_hu_cdm_module.so.dSYM/

eisenstein_hu_cdm_module.so:   spline.f90 interface_tools.f90 eisenstein_hu_cdm.f90  eisenstein_hu_cdm_module.f90 ${COSMOSIS_SRC_DIR}/cosmosis/datablock/libcosmosis_fortran.so
	$(FC) $(FFLAGS) -shared -o $@ $+
