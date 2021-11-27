program lu_driver

  use lu
  use utils

  implicit none

  double precision, allocatable :: A(:,:), x(:), b(:,:), resVec(:)
  integer, allocatable :: ipiv(:)
  integer :: i, n

  double precision :: resNorm
  double precision :: t_start, t_end

  n = 500

  allocate(A(n,n))
  allocate(b(n,1))
  allocate(x(n))
  allocate(ipiv(n))
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

  x = 1.d0
  b(:,1) = matmul(A,x)


  call get_wtime(t_start) 
  call lu_solve_lgs(n,A,1,b,ipiv)
  call get_wtime(t_end)
  
  resVec = abs(b(:,1) - x)
  resNorm = norm2(resVec)

  write(*,'(A,E16.7,A)') 'Elapsed time: ', t_end - t_start, ' seconds'
  write(*,'(A,E16.7)') 'Error norm: ', resNorm

  deallocate(resVec)
  deallocate(ipiv)
  deallocate(x)
  deallocate(b)
  deallocate(A)

end program lu_driver
