Modified Krylov subspace methods for special application
===

# Requirements
For compiling the following is needed
- Fortran compiler (e.g. GNU gfortran or Intel ifort)
- make

# Krylov subspace methods

| Method | Description |
| --- | --- |
| GMRES | Algorithm from Saad and Schultz 1986 with 2-times Gram-Schmidt |
| schurGMRES | Modified GMRES for solving a saddlepoint problem with Schur complement |

# Helper Routines

| Routine | Descritpion |
| --- | --- |
| givens | Compute and apply the Givens rotation for a given vector |
| printMatrix | Print an array of dimension 2 to stdout |
| writeMatrix | Save an array of dimension 2 to a file |
| printVector | Print an array of dimension 1 to stdout |
| writeVector | Save an array of dimension 1 to a file |
