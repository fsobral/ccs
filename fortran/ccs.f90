MODULE ccs

  implicit none

CONTAINS

  subroutine perm(n,p)

    ! ARGUMENTS
    integer :: n
    integer :: p(n)

    ! LOCAL ARRAYS
    integer :: t(n)

    ! LOCAL SCALARS
    integer :: sizet,i,r
    real(8) :: rd

    do i = 1,n
       t(i) = i
    end do

    sizet = n

    do i = 1,n
       call RANDOM_NUMBER(rd)
       r = int(1 + sizet * rd)
       p(i) = t(r)
       t(r) = t(sizet)
       sizet = sizet - 1
    end do

  end subroutine perm

  subroutine serialccs(m1,n,M,b,x,r,iter,eps,maxit)

    ! SCALAR ARGUMENTS
    integer :: n,iter,m1,maxit
    real(8) :: eps,r

    ! ARRAY ARGUMENTS
    real(8) :: M(n,n),b(n),x(n)

    intent(in   ) :: M,b,n,eps,maxit
    intent(out  ) :: r,iter
    intent(inout) :: x

    ! LOCAL SCALARS
    integer :: i,j,k,curr,delay
    real(8) :: s,rd

    ! LOCAL ARRAYS
    integer :: p(m1)
    real(8) :: back(n,m1)

    r = 0
    do i = 1,n
       s = b(i)
       do j = 1,n
          s = s - M(i,j) * x(j)
       end do
       r = max(r,abs(s))
    end do

    iter = 0

    do j = 1,m1
       do i = 1,n
          back(i,j) = x(i)
       end do
    end do

    do while (r .gt. eps)

       call perm(n,p)

       do i = n+1,m1
          call RANDOM_NUMBER(rd)
	  p(i) = INT(1 + n * rd)
       end do

       ! UPDATING X

       do k = 1,m1

          curr = 1 + mod(iter,m1)

          do i = 1,n
             back(i,curr) = x(i)
          end do

          call RANDOM_NUMBER(rd)

          delay = int(1 + m1 * rd)

          call updatexh(p(k),p(k),n,M,b,back(:,1+mod(m1 + curr - delay,m1)),x)

          ! STOPPING CRITERIUM

          r = 0
          do i = 1,n
             s = b(i)
             do j = 1,n
                s = s - M(i,j) * x(j)
             end do

             r = max(r,abs(s))
             if (r > eps) then
                exit
             end if
          end do

          iter = iter + 1

          if (r .le. eps) then
             exit
          end if

       end do

    end do

  end subroutine serialccs

  subroutine updatexh(start,end,n,M,b,x,xout)

    ! This subroutine applies a Jacobi step in a subset of the rows
    ! and variables of the system. The input vector considered is 'x'
    ! and the output vector is 'xout'.

    ! ARGUMENTS
    integer :: start,end,n
    real(8) :: M(n,n),b(n),x(n),xout(n),tmp(start:end)

    intent(in   ) :: start,end,n,M,b
    intent(inout) :: xout

    ! LOCAL SCALARS
    integer :: i,j
    real(8) :: s

    do i = start,end
       s = b(i)
       do j = 1,n
          if (i .ne. j) then
             s = s - M(i,j) * x(j)
          end if
       end do
       tmp(i) = s / M(i,i)
    end do

    do i = start,end
       xout(i) = tmp(i)
    end do

  end subroutine updatexh

END MODULE ccs
