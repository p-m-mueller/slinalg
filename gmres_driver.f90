program gmres_driver

  use linalg
  use utils

  implicit none

  double precision, allocatable :: A(:,:), x(:), b(:)
  double precision :: tol
  integer :: i, n, maxIter

  n = 5

  allocate(A(n,n))
  allocate(b(n))
  allocate(x(n))

  A = 0.d0
  do i = 2, n-1
    A(i,i-1) = -1.d0
    A(i,i) = 2.d0
    A(i,i+1) = -1.d0
  end do
  A(1,1) = 2.d0
  A(1,2) = -1.d0
  A(n,n) = 2.d0
  A(n,n-1) = -1.d0

  write(*,*)
  write(*,'(A)') "Operator A:"
  call printMatrix(A)

  x = 0.d0
  b = matmul(A,x)

  write(*,*)
  write(*,'(A)') "RHS b:"
  call printVector(b)
  
  write(*,*)
  write(*,'(A)') "Exact solution x:"
  call printVector(x)

  maxIter = 4
  tol = 1.0e-16

  x = 0.d0
  call GMRES(A, b, maxIter, tol, x)

  write(*,*)
  write(*,'(A)') "Solution x:"
  call printVector(x)

  deallocate(x)
  deallocate(b)
  deallocate(A)

end program gmres_driver
