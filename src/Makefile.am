AUTOMAKE_OPTIONS = foreign

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = IMEXSfloW
bindir=@bindir@

# list of sources for the binary
IMEXSfloW_SOURCES = \
		complexify.f90 \
		parameters.f90 \
		geometry.f90 \
		constitutive.f90 \
		solver.f90 \
		init.f90 \
		inpout.f90 \
		IMEX_SfloW.f90


FCFLAGS         = -g -O3 -ffast-math
FCDEBUGFLAGS    = -Wall -g

LIBS            = -lm -llapack

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran

