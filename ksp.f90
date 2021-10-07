module ksp

  use utils 

  implicit none 

  logical :: verbose = .false.

  contains

    subroutine GMRES(A, b, maxIter, tol, x)
      double precision, intent(in) :: A(:,:)
      double precision, intent(in) :: b(:)
      double precision, intent(in out) :: x(:)
      integer, intent(in) :: maxIter
      double precision, intent(in) :: tol

      double precision, allocatable :: Q(:,:), H(:,:)
      double precision, allocatable :: r(:), c(:), beta(:), y(:)
      double precision, allocatable :: cs(:), sn(:)
      double precision :: reltol, res, tmp
      integer :: i, j, k, n, m

      n = size(A,1)
      if (n /= size(b) .or. n /= size(x)) then 
        write(*,'(A)') "Error in gmres(): Size mismatch"
        return
      end if

      m = maxIter
      reltol = tol * norm2(b)

      if (verbose) then
        write(*,*)
        write(*,'(A)')       '---------------------------------------'
        write(*,'(A,E10.3)') 'atol:     ', tol 
        write(*,'(A,E10.3)') 'rtol:     ', reltol 
        write(*,'(A,I0)')    'Max iter: ', maxIter
        write(*,'(A)')       '---------------------------------------'
      end if

      allocate(Q(n,m+1))
      allocate(H(m+1,m))
      allocate(beta(m+1))
      allocate(cs(m), sn(m))
      allocate(r(n), c(n))
     
      beta = 0.d0
      Q = 0.d0; H = 0.d0

      r = b - matmul(A,x)
      beta(1) = norm2(r)
     
      res = beta(1) 
      if (res > reltol) then
        Q(:,1) = r / beta(1)
   
        if (verbose) then  
          write(*,'(A10,1X,A10)') "Iter", "Residual"
          write(*,'(A10,1X,A10)') "----", "--------"
          write(*,'(I10,1X,E10.3)') 0, beta(1)
        end if 

        j = 1
        do !j = 1, m
          
          ! Arnoldi
          r = matmul(A,Q(:,j))
          do i = 1, 2                                 ! Gram-Schmidt 2-times
            do k = 1, j
              c(k) = dot_product(Q(:,k),r)
            end do
            r = r - matmul(Q(:,1:j),c(1:j))
            H(1:j,j) = H(1:j,j) + c(1:j)
          end do
          H(j+1,j) = norm2(r)
          Q(:,j+1) = r / H(j+1,j)
  
          ! Givens
          do i = 1, j-1
            tmp = cs(i) * H(i,j) + sn(i) * H(i+1,j)
            H(i+1,j) = - sn(i) * H(i,j) + cs(i) * H(i+1,j)
            H(i,j) = tmp
          end do
          call givens(H(j:j+1,j),cs(j),sn(j))
  
          beta(j+1) = -sn(j) * beta(j)
          beta(j)   =  cs(j) * beta(j)
          res = abs(beta(j+1))
  
          if (verbose) write(*,'(I10,1X,E10.3)') j, res
          if (res <= reltol .or. j >= maxIter) exit
          j = j + 1
        end do
       
        if (verbose) then  
          write(*,'(A)')       '---------------------------------------'
          write(*,'(A,I0,A)')  "GMRES exit after ", j, " iterations"
          write(*,'(A)')       '---------------------------------------'
          write(*,*)
        end if
        
        allocate(y(j))
        do i = j, 1, -1
          y(i) = ( beta(i) - dot_product(H(i,i+1:j), y(i+1:j)) ) / H(i,i)
        end do  
        x = x + matmul(Q(:,1:j),y)
        deallocate(y) 
      end if

      deallocate(r, c)
      deallocate(cs, sn)
      deallocate(beta) 
      deallocate(H)
      deallocate(Q)
    end subroutine GMRES


    ! Solves a block structured system 
    !
    !   [ A  B ] [x]   [f]
    !   [ B' O ] [y] = [g]
    ! 
    ! by solving the Schur complement reduction
    !
    !  - B'inv(A)By = g - B'inv(A)f
    !            Ax = f - By
    !
    ! with GMRES.
    subroutine schurGMRES(A, B, f, g, maxIter, tol, x, y)
      double precision, intent(in) :: A(:,:), B(:,:)
      double precision, intent(in) :: f(:), g(:)
      double precision, intent(in out) :: x(:), y(:)
      integer, intent(in out) :: maxIter
      double precision, intent(in out) :: tol

      double precision, allocatable :: gg(:), Q(:,:), H(:,:), beta(:), u(:), v(:), w(:), z(:)
      double precision, allocatable :: r(:), c(:)
      double precision, allocatable :: cs(:), sn(:)
      double precision, allocatable :: ff(:)
      double precision :: reltol, res, tmp
      integer :: i, j, k, nra, nca, nrb, ncb, n, na, m

      nra = size(A,1)
      nca = size(A,2)
      nrb = size(B,1)
      ncb = size(B,2)
      if (nra /= nrb) then 
        write(*,'(A)') "Error in gmres(): Size mismatch. Rows of A and B are different."
        return
      end if
      if (nca /= nrb) then 
        write(*,'(A)') "Error in gmres(): Size mismatch. Cols of A and transpose(B) are different."
        return
      end if
      
      n = ncb
      na = nra

      m = maxIter
      reltol = tol * norm2(g)

      if (verbose) then
        write(*,*)
        write(*,'(A)')       '---------------------------------------'
        write(*,'(A,E10.3)') 'atol:     ', tol 
        write(*,'(A,E10.3)') 'rtol:     ', reltol 
        write(*,'(A,I0)')    'Max iter: ', maxIter
        write(*,'(A)')       '---------------------------------------'
      end if

      allocate(Q(n,m+1))
      allocate(H(m+1,m))
      allocate(beta(m+1))
      allocate(cs(m), sn(m))
      allocate(c(n), gg(n), r(n))
      allocate(u(na), v(na), w(na))
     
      beta = 0.d0
      Q = 0.d0; H = 0.d0

      w = 0.d0
      call gmres(A, f, na, tol, w)
      gg = g - matmul(transpose(B),w);

      u = matmul(B, y)
      v = 0.d0
      call gmres(A, u, na, tol, v)
      r = gg + matmul(transpose(B),v)
      beta(1) = norm2(r)
      
      res = beta(1) 
      if (res > reltol) then
        Q(:,1) = r / beta(1)
 
        if (verbose) then
          write(*,'(A10,1X,A10)') "Iter", "Residual"
          write(*,'(A10,1X,A10)') "----", "--------"
          write(*,'(I10,1X,E10.3)') 0, beta(1)
        end if
  
        j = 1
        do !j = 1, m
          
          ! Arnoldi
          u = matmul(B,Q(:,j))
          v = 0.d0;
          call gmres(A, u, na, tol, v)
          r = - matmul(transpose(B), v)
          do i = 1, 2                                 ! Gram-Schmidt 2-times
            do k = 1, j
              c(k) = dot_product(Q(:,k),r)
            end do
            r = r - matmul(Q(:,1:j),c(1:j))
            H(1:j,j) = H(1:j,j) + c(1:j)
          end do
          H(j+1,j) = norm2(r)
          Q(:,j+1) = r / H(j+1,j)
  
          ! Givens
          do i = 1, j-1
            tmp = cs(i) * H(i,j) + sn(i) * H(i+1,j)
            H(i+1,j) = - sn(i) * H(i,j) + cs(i) * H(i+1,j)
            H(i,j) = tmp
          end do
          call givens(H(j:j+1,j),cs(j),sn(j))
  
          beta(j+1) = -sn(j) * beta(j)
          beta(j) = cs(j) * beta(j)
          res = abs(beta(j+1))
  
          if (verbose) write(*,'(I10,1X,E10.3)') j, res
          if (res <= reltol .or. j >= maxIter) exit
          j = j + 1
        end do
      
        if (verbose) then   
          write(*,'(A)')       '---------------------------------------'
          write(*,'(A,I0,A)')  "GMRES exit after ", j, " iterations"
          write(*,'(A)')       '---------------------------------------'
          write(*,*)
        end if
        
        call writeMatrix(matmul(transpose(Q(:,1:j)),Q(:,1:j)), 'eye.mat')
  
        allocate(z(j))
        do i = j, 1, -1
          z(i) = ( beta(i) - dot_product(H(i,i+1:j), z(i+1:j)) ) / H(i,i)
        end do  
        y = y + matmul(Q(:,1:j),z)
        deallocate(z) 
      end if

      deallocate(gg, r, u, v, w)
      deallocate(cs, sn)
      deallocate(beta)
      deallocate(H)
      deallocate(Q)

      allocate(ff(na))
      ff = f - matmul(B,y)
      call gmres(A, ff, na, tol, x)
      deallocate(ff)

      maxIter = j
      tol = res
    end subroutine schurGMRES
    

end module ksp
