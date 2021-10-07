program schurGMRES_driver

  use utils
  use linalg

  implicit none 

  integer :: na, n
  double precision, allocatable :: A(:,:), B(:,:)
  double precision, allocatable :: f(:), g(:), x(:), y(:)
  double precision :: res, tol
  integer :: maxIter

  integer :: i, j

  na = 12
  n = 4

  allocate(A(na,na))
  allocate(B(na,n))
  allocate(f(na))
  allocate(g(n))
  allocate(x(na))
  allocate(y(n))

  do i = 2, na-1
    A(i,i-1) = -1.d0; A(i,i) = 2.d0; A(i,i+1) = -1.d0
  end do
  A(1,1) = 2.d0; A(1,2) = -1.d0
  A(na,na-1) = -1.d0; A(na,na) = 2.d0;

  do i = 1, n
    B(i:na,i) = 1.d0;
  end do

  ! Exact solution
  x = 1.d0
  y = 1.d0

  ! Create right hand side
  f = matmul(A,x) + matmul(B,y)
  g = matmul(transpose(B),x)


  ! Reset initial guess
  x = 0.d0
  y = 0.d0

  maxIter = 4
  tol = 1.0e-16

  verbose=.true.
  call schurGMRES(A, B, f, g, maxIter, tol, x, y)
  write(*,*) "schurGMRES:"
  write(*,*) "Iteration: ", maxIter 
  write(*,*) "Residual:  ", tol
  write(*,*)
  write(*,*) "x: "
  call printVector(x)
  
  write(*,*) "y: "
  call printVector(y)
  
  deallocate(y)
  deallocate(x)
  deallocate(g)
  deallocate(f)
  deallocate(B)
  deallocate(A)

end program schurGMRES_driver
