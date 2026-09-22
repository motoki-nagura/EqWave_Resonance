c**************************************************
c      program read_kmt
c      
c
c        read kmt.dta(ASCCI 3602x1502)
c        make kmt.direct(binary 3602x1502)
c**************************************************
c
      parameter (imt=622,jmt=402)
      parameter (imt1=imt-2,jmt1=jmt-2)
c
      integer kmt(imt,jmt)
      real    rkmt(imt,jmt)
      character*20 cell_type
      character*80 flname
c
c==== READ kmt.dta =======
      flname='kmt.dta'
      open(unit=10,FILE=flname,STATUS='OLD',FORM='FORMATTED')
      read (10,'(i6,i6,i6,a20)') imax, jmax, kmax ,cell_type
      write(*,*) 'READ kmt.dta file '
      write(*,*) 'imax,jmax,kmax:',imax,jmax,kmax,'  ',cell_type

      inc = 30
      !inc = 40

      do l=0,imax,inc
        incr = min(inc,imax-l)
        !read (10,'(/,/,/,2x,20i6/)') (ii,i=1,incr,2)
        read (10,1002) (ii,i=1,incr,2)
        do jrow=jmax,1,-1
          !read (10,'(1x,i4,1x,40i3)') jrowin,(kmt(l+i,jrow),i=1,incr)
          read (10,2002) jrowin,(kmt(l+i,jrow),i=1,incr)
          write(*,*) jrow
        enddo
      enddo
      close(10)

 1002 format(/,/,/,1x,15i8)
 2002 format(i4,1x,30i4)
c
c-----------------------------------------------------------------
c     Write kmt.direct
c-----------------------------------------------------------------
c
      do i=1,imt
      do j=1,jmt
        rkmt(i,j)=real(kmt(i,j))
      enddo
      enddo
      write(6,*) (rkmt(i,jmt/2),i=imt/2,imt/2+10)
c
      flname='kmt.dta.bin'
      open(81,file=flname,form='unformatted',access='direct'
     &,    recl=imt1*jmt1*4)
        write(81,rec=1) ((rkmt(i,j),i=2,imt-1),j=2,jmt-1)
c        write(81,rec=1) ((rkmt(i,j),i=1,imt),j=1,jmt)
      close(81)
c
c***** Check Error *****
      write(*,*)
      write(*,*) '===  Error check begins  ==='
c
      do j=2,jmt-1
      do i=2,imt-1
c******* 1 row path
        if( kmt(i,j) .gt. 0 ) then
          if( kmt(i-1,j).eq.0 .and. kmt(i+1,j).eq.0 )
     #      write(*,*)'Error. 1 row path :  i,j = ',i,j
          if( kmt(i,j-1).eq.0 .and. kmt(i,j+1).eq.0 )
     #      write(*,*)'Error. 1 row path :  i,j = ',i,j
c
c******* isolated kmt
          if( kmt(i,j).gt.kmt(i+1,j) .and.
     #        kmt(i,j).gt.kmt(i-1,j) .and.
     #        kmt(i,j).gt.kmt(i,j+1) .and.
     #        kmt(i,j).gt.kmt(i,j-1) )
     #      write(*,*)'Error. Isolated hole (1) :  i,j = ',i,j
c
c*******
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
          if( kmt(i,j).gt.kmt1  .and.
     #        kmt(i,j).gt.kmt2  .and.
     #        kmt(i,j).gt.kmt3  .and.
     #        kmt(i,j).gt.kmt4 )
     #      write(*,*)'Error. Isolated hole (2) :  i,j = ',i,j
c
        end if
c
      end do
      end do
c
      stop
      end

