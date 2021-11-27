program gmres_driver

  use ksp 
  use utils

  implicit none

  double precision, allocatable :: A(:,:), x_ext(:), x(:), b(:), resVec(:)
  double precision :: tol, resNorm, t_start, t_end
  integer :: i, n, maxIter

  n = 500

  allocate(A(n,n))
  allocate(b(n))
  allocate(x_ext(n))
  allocate(x(n))
  allocate(resVec(n))

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

  x_ext = 1.d0
  b = matmul(A,x_ext)

  maxIter = n
  tol = 1.0e-16

  x = 0.d0
  call get_wtime(t_start)
  call GMRES(A, b, maxIter, tol, x)
  call get_wtime(t_end)

  resVec = abs(x_ext - x)
  resNorm = norm2(resVec)

  write(*,'(A20,E16.7,A)') 'Elapsed time: ', t_end - t_start, ' seconds'
  write(*,'(A20,E16.7)') 'Error norm: ', resNorm
  write(*,'(A20,I16)') 'Iterations: ', maxIter

  deallocate(resVec)
  deallocate(x)
  deallocate(x_ext)
  deallocate(b)
  deallocate(A)

end program gmres_driver
