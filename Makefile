all: GMRES-driver schurGMRES-driver  lu-driver

GMRES-driver: timing.o utils.o ksp.o gmres_driver.o
	gfortran -o $@ $^

schurGMRES-driver: utils.o ksp.o schur_driver.o
	gfortran -o $@ $^

lu-driver: timing.o utils.o lu.o lu_driver.o
	gfortran -o $@ $^

%.o %.mod: %.f90
	gfortran -g -O0 -fbounds-check -c -o $@ $<

%.o: %.c
	gcc -g -O0 -c -o $@ $<

clean:
	@rm -vf *.o *.mod run *.mat GMRES-driver schurGMRES-driver lu-driver
