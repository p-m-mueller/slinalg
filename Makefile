all: GMRES-driver schurGMRES-driver 

GMRES-driver: utils.o linalg.o gmres_driver.o
	gfortran -o $@ $^

schurGMRES-driver: utils.o linalg.o schur_driver.o
	gfortran -o $@ $^

%.o %.mod: %.f90
	gfortran -g -O0 -fbounds-check -c -o $@ $<

clean:
	@rm -vf *.o *.mod run *.mat GMRES-driver schurGMRES-driver
