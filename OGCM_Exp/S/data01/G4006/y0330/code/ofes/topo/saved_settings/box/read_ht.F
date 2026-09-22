
      parameter (imt=802,jmt=302,km=100)
      parameter (imt1=imt-2, jmt1=jmt-2)

      real*8 ht(imt,jmt)
      integer*8 imax,jmax,kmax
      real rht(imt,jmt)
      character*80 flname
c
c-----------------------------------------------------------------
c     Read ht data
c-----------------------------------------------------------------
c
      io=10
c
      flname='ht.dta'
      open(io,file=flname,
     &     form='unformatted')
      read (io)
      read (io) imax, jmax, kmax
       write(6,*) imax,jmax,kmax
       if(imax.ne.imt .or. jmax.ne.jmt .or. kmax.ne.km)
     &   stop 'Different parameters'
      read (io)
      read (io) ht
      close(io)
      write(6,*) (ht(i,jmt/2),i=imt/2,imt/2+10)
c
c-----------------------------------------------------------------
c     Write ht data
c-----------------------------------------------------------------
c
      do i=1,imt
      do j=1,jmt
        rht(i,j)=real(ht(i,j))/100.
      enddo
      enddo
      write(6,*) (rht(i,jmt/2),i=imt/2,imt/2+10)
c
      flname='ht.dta.bin'
      open(81,file=flname,form='unformatted',access='direct'
     &,    recl=imt1*jmt1*4)
        write(81,rec=1) ((rht(i,j),i=2,imt-1),j=2,jmt-1)
      close(81)
c
      stop
      end

