module utils

  implicit none

  contains
    
    subroutine givens(x,cs,sn)
      double precision, intent(in out) :: x(2)
      double precision, intent(out) :: cs, sn
      double precision :: r

      r = sqrt(x(1)*x(1) + x(2)*x(2))

      cs = x(1) / r
      sn = x(2) / r

      x(1) = r
      x(2) = 0.d0

    end subroutine givens

    subroutine printMatrix(A)
      double precision, intent(in) :: A(:,:)
      integer :: i, j, n, m

      n = size(A,1)
      m = size(A,2)

      do i = 1, n
        write(*,'(A)', advance='no') "["
        do j = 1, m
          write(*,'(1X,E10.2)', advance='no') A(i,j)
        end do
        write(*,'(A)') "]"
      end do

    end subroutine printMatrix

    subroutine writeMatrix(A, filename)
      double precision, intent(in) :: A(:,:)
      character(len=*), intent(in) :: filename 

      integer :: i, j, n, m
      integer :: file_unit, ierr

      open(newunit=file_unit, file=filename, action='write', iostat=ierr)
      if (ierr /= 0) then 
        write(*,*) "Error in writeMatrix(): Could not open file "//filename
        return
      end if 

      n = size(A,1)
      m = size(A,2)

      do i = 1, n
        do j = 1, m
          write(file_unit,'(1X,F23.16)',advance='no') A(i,j)
        end do
        write(file_unit,*)
      end do

      close(file_unit, iostat=ierr)
      if (ierr /= 0) then 
        write(*,*) "Error in writeMatrix(): Could not close file "//filename
        return
      end if 

    end subroutine writeMatrix

    subroutine printVector(v)
      double precision, intent(in) :: v(:)
      integer :: i, n

      n = size(v)

      do i = 1, n
        write(*,'(A,E10.3,A)') "[", v(i), "]"
      end do

    end subroutine printVector
    
    subroutine writeVector(v, filename)
      double precision, intent(in) :: v(:)
      character(len=*), intent(in) :: filename 

      integer :: i, n
      integer :: file_unit, ierr

      open(newunit=file_unit, file=filename, action='write', iostat=ierr)
      if (ierr /= 0) then 
        write(*,*) "Error in writeVector(): Could not open file "//filename
        return
      end if 

      n = size(v)

      do i = 1, n
        write(file_unit,'(1X,F23.16)',advance='no') v(i)
      end do

      close(file_unit, iostat=ierr)
      if (ierr /= 0) then 
        write(*,*) "Error in writeVector(): Could not close file "//filename
        return
      end if 

    end subroutine writeVector

end module utils
