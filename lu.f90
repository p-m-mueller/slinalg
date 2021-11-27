module lu

  implicit none 

  double precision :: small = 1.0e-25
  logical :: verbose = .false. 


  contains 

    subroutine lu_decomp(n,A,ipiv)

      integer, intent(in) :: n
      integer, intent(in out) :: ipiv(:)
      double precision, intent(in out) :: A(:,:)
  
      double precision, allocatable :: tmpRow(:)
      double precision :: tmp, maxColumnEntry
      integer :: i, j, k, imax

      allocate(tmpRow(n))

      do i = 1, n
        ipiv(i) = i
      end do

      do i = 1, n-1

        ! Find largest entry in column i
        maxColumnEntry = 0.d0 
        imax = i 
        do k = i, n
          tmp = abs(A(k,i))
          if (tmp > maxColumnEntry) then 
            maxColumnEntry = tmp
            imax = k
          end if 
        end do

        if (maxColumnEntry < small) then 
          write(0,*) 'Error in lu::lu_decomp(): Matrix close to singular'
          stop 1
        end if

        ! Pivoting
        if (imax /= i) then 
          ! Swap indices
          j = ipiv(i)
          ipiv(i) = ipiv(imax)
          ipiv(imax) = j

          ! Swap rows 
          tmpRow = A(i,:)
          A(i,:) = A(imax,:)
          A(imax,:) = tmpRow
        end if

        do  j = i+1, n
          ! Entries of L
          A(j,i) = A(j,i) / A(i,i)

          do k = i+1, n
            ! Enties of R
            A(j,k) = A(j,k) - A(j,i) * A(i,k)
          end do
        end do
      end do

      deallocate(tmpRow)

    end subroutine lu_decomp

    subroutine lu_solve_lgs(n,A,nRhs,B,ipiv)

      integer, intent(in) :: n, nRhs
      integer, intent(in out) :: ipiv(:)
      double precision, intent(in out) :: A(:,:), B(:,:)

      integer :: i, k

      call lu_decomp(n, A, ipiv)

      B = B(ipiv,:)

      ! Forward substitution
      do i = 1, n
        do k = 1, i-1
          B(i,:) = B(i,:) - A(i,k) * B(k,:)
        end do
      end do

      ! Backwards substitution
      do i = n, 1, -1
        do k = i+1, n
          B(i,:) = B(i,:) - A(i,k) * B(k,:)
        end do

        B(i,:) = B(i,:) / A(i,i)
      end do

    end subroutine lu_solve_lgs


end module lu
