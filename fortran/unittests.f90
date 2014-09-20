program unittests

  implicit none

  integer :: info

  write(*,*) 'Test 1 - Random DP matrix generation'
  call testrandDP(info)
  if (info .eq. 0) then
     write(*,*) 'OK'
  else
     write(*,*) 'NOT OK'
  end if

!!$  write(*,*) 'Test 2 - Random Lower Chol generation'
!!$  call testrandchol(info)
!!$  if (info .eq. 0) then
!!$     write(*,*) 'OK'
!!$  else
!!$     write(*,*) 'NOT OK'
!!$  end if

contains

  subroutine testrandDP(info)

    use gentests

    implicit none
    
    real(8),allocatable :: A(:,:)
    integer :: i,j,n,info

    interface
       subroutine dpotrf (UPLO, N, A, LDA, INFO)
         character :: UPLO
         integer   :: N,LDA,INFO
         real(8)   :: A(LDA,*)

         intent(   in) :: UPLO,N,LDA
         intent(inout) :: A
         intent(  out) :: INFO
       end subroutine dpotrf
    end interface

    do n = 1,50,1

       allocate(A(n,n))
    
       call genRandDP(n,A,5.0D-1)

       CALL dpotrf('L',n,A,n,info)

       if (info .ne. 0) then
          write(*,*) 'Error in the SDP creation: Matrix is not DP. n=', n
          exit
       end if

       deallocate(A)

    end do

  end subroutine testrandDP

  subroutine testrandchol(info)

    use gentests

    implicit none
    
    real(8),allocatable :: L(:,:)
    integer :: i,j,n,info
    real(8) :: tmp

    info = 0

    do n = 10,1010,100

       allocate(L(n,n))
    
       call genRandCholL(n,L,5.0D-1)

       tmp = 1.0D+20
       do i = 1,n
          tmp = min(tmp,L(i,i))
       end do

       if (tmp .le. 1.0D-5) then
          write(*,*) 'Error in the lower Chol. creation. Negative diag. n=', &
               n
          info = -1
          exit
       end if

       tmp = 0.0D0
       do j = 1,n
          do i = 1,j - 1
             tmp = max(tmp,L(i,j))
          end do
       end do

       if (tmp .gt. 0.0D0) then
          write(*,*) 'Error in the lower Chol. creation. Nonzero upper. n=', &
               n
          info = -1
          exit
       end if

       deallocate(L)

    end do

  end subroutine testrandchol

end program unittests
