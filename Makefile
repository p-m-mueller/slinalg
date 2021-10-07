all: GMRES-driver schurGMRES-driver 

GMRES-driver: utils.o ksp.o gmres_driver.o
	gfortran -o $@ $^

schurGMRES-driver: utils.o ksp.o schur_driver.o
	gfortran -o $@ $^

%.o %.mod: %.f90
	gfortran -g -O0 -fbounds-check -c -o $@ $<

clean:
	@rm -vf *.o *.mod run *.mat GMRES-driver schurGMRES-driver
