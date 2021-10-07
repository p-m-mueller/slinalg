program matmul_test

  implicit none 

  integer :: n
  double precision, allocatable :: A(:,:), x(:), b(:)

  integer :: i, j

  n = 5

  allocate(A(n,n))
  allocate(x(n))
  allocate(b(n))

  do i = 2, n-1
      A(i-1,i) = -1.d0; A(i,i) = 2.d0; A(i+1,i) = -1.d0
  end do
  A(1,1) = 2.d0
  A(n-1,n) = -1.d0; A(n,n) = 2.d0

  x = 1.d0
  b = 0.d0

  write(*,*) "intrinsic matmul"
  b = matmul(A,x)
  call printVector(b)

  write(*,*) "own matmul"
  b = AXY(n,n,A,x)
  call printVector(b)

  deallocate(b)
  deallocate(x)
  deallocate(A)

  contains

    subroutine printVector(x)
      double precision, intent(in) :: x(:)
      integer :: i, n
  
      n = size(x) 
      do i = 1, n 
        write(*,'(A,E10.3,A)') "[", x(i), "]"
      end do
    end subroutine printVector

    function AXY(n,m,A,y)
      double precision, intent(in) :: A(n,m), y(m)
      integer, intent(in) :: n, m
      double precision :: AXY(n)
      integer :: i, j

      AXY = 0.d0
      do j = 1, m
        do i = 1, n
          AXY(i) = AXY(i) + A(i,j) * y(j)
        end do
      end do

    end function

end program
