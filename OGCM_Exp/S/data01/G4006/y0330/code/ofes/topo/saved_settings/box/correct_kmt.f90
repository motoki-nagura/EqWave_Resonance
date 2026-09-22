!**************************************************
!      program correct_kmt
!      
!
!        read kmt.dta (ASCCI 722 x 218)
!             kmt.dta.orig
!        write kmt.dta (ASCCI 722 x 218)
!**************************************************
!
  parameter (imt=722,jmt=218)
  parameter (imt1=imt-2,jmt1=jmt-2)
!
  integer, dimension(imt,jmt) ::  kmt_orig,kmt_old,kmt, tmp_kmt
  real    rkmt(imt,jmt)
  character*20 cell_type
  character*80 flname

  integer :: iounit
!
!==== READ original kmt.dta =======

  do iounit = 21, 22

    read(iounit,'(a)') flname
    write (*,*) ' IN << '//trim(flname)

    open(unit=10,FILE=flname,STATUS='OLD',FORM='FORMATTED')
    read (10,'(i6,i6,i6,a20)') imax, jmax, kmax ,cell_type
    write(*,*) 'READ kmt.dta file '
    write(*,*) 'imax,jmax,kmax:',imax,jmax,kmax,'  ',cell_type

    select case ( iounit )
      case ( 21 )  ;  inc = 40
      case ( 22 )  ;  inc = 30
    end select

    do l=0,imax,inc
      incr = min(inc,imax-l)

      select case ( iounit )
        case ( 21 ) ;  read (10,1001) (ii,i=1,incr,2)
        case ( 22 ) ;  read (10,1002) (ii,i=1,incr,2)
      end select

!      write (*,*) 'l,incr,ii = ',l,incr,ii

      do jrow=jmax,1,-1
        select case ( iounit )
          case ( 21 ) ;  read (10,2001) jrowin,(tmp_kmt(l+i,jrow),i=1,incr)
          case ( 22 ) ;  read (10,2002) jrowin,(tmp_kmt(l+i,jrow),i=1,incr)
        end select
      enddo
    enddo
    close(10)

    select case ( iounit )
      case ( 21 )  ;  kmt_old  = tmp_kmt
      case ( 22 )  ;  kmt_orig = tmp_kmt
    end select

  end do

  1001 format(/,/,/,2x,20i6/)
  1002 format(/,/,/,1x,15i8)

  2001 format(1x,i4,1x,40i3)
  2002 format(i4,1x,30i4)

  write (*,*)
  write (*,*) 'kmt_old  = ',minval(kmt_old), maxval(kmt_old)
  write (*,*) 'kmt_orig = ',minval(kmt_orig),maxval(kmt_orig)
  write (*,*)
!
!-----------------------------------------------------------------
!     correct kmt
!-----------------------------------------------------------------
!
  kmt = kmt_orig

  where ( kmt_old == 0 )
    kmt = 0             !  same land grid as mine_0-360e
  end where

  !----------

  do j=2,jmt-1
  do i=2,imt-1

    if( kmt(i,j) > 0 ) then

!******* 1 row path
      if ( (kmt(i-1,j) == 0) .and. (kmt(i+1,j) == 0) ) then
        kmt(i,j) = 0
      end if

      if ( (kmt(i,j-1) == 0) .and. (kmt(i,j+1) == 0) ) then
        kmt(i,j) = 0
      end if
!
!******* isolated hole 1
      if( (kmt(i,j) > kmt(i+1,j)) .and. &
          (kmt(i,j) > kmt(i-1,j)) .and. &
          (kmt(i,j) > kmt(i,j+1)) .and. &
          (kmt(i,j) > kmt(i,j-1)) )       then
        kmt(i,j) = max(kmt(i+1,j), kmt(i-1,j), kmt(i,j+1), kmt(i,j-1))
      endif

!**** isolated hole 2
      kmt1=min(kmt(i+1,j),kmt(i,j+1),kmt(i+1,j+1))
      !             . x x
      !             . . x
      !             . . .
      kmt2=min(kmt(i+1,j),kmt(i,j-1),kmt(i+1,j-1))
      !             . . .
      !             . . x
      !             . x x
      kmt3=min(kmt(i-1,j),kmt(i,j-1),kmt(i-1,j-1))
      !             . . .
      !             x . .
      !             x x .
      kmt4=min(kmt(i-1,j),kmt(i,j+1),kmt(i-1,j+1))
      !             x x .
      !             x . .
      !             . . .
      if( (kmt(i,j) > kmt1)  .and. &
          (kmt(i,j) > kmt2)  .and. &
          (kmt(i,j) > kmt3)  .and. &
          (kmt(i,j) > kmt4) )        then
          kmt(i,j) = max(kmt1, kmt2, kmt3, kmt4)
      endif
!
    end if
!
  end do
  end do

!
!-----------------------------------------------------------------
!     write kmt
!-----------------------------------------------------------------
!
  iounit = 41

  read(iounit,'(a)') flname
  write (*,*) ' OUT >> '//trim(flname)

  open(unit=10,FILE=flname,STATUS='NEW',FORM='FORMATTED')
  write (10,'(i6,i6,i6,a20)') imax, jmax, kmax ,cell_type
  write(*,*) 'WRITE kmt.dta file '

  inc = 30

  do l=0,imax,inc
    incr = min(inc,imax-l)
    write (10,1002) (i+l,i=1,incr,2)
    do jrow=jmax,1,-1
      write (10,2002) jrow,(kmt(l+i,jrow),i=1,incr)
    enddo
  enddo
  close(10)

  stop
end

