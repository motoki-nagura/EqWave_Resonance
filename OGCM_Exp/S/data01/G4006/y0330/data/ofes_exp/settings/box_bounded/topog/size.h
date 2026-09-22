!====================== include file "size.h" ==========================
!                                                                       
!-----------------------------------------------------------------------
!     USER INPUT: 
!-----------------------------------------------------------------------
!
!     imt    = number of grid cells in the longitudinal direction
!              (calculated cells are from 2 through imt-1. end cells
!               are boundaries)          
!
!     jmt    = number of grid cells (latitude rows) in the latitudinal
!              direction (calculated cells are from 2 through jmt-1.
!              end cells are boundaries)   
!                                
!     km     = number of grid cells in the vertical direction          
!              (calculated cells are from 1 through km. no boundary
!               cells)
!
!     nt     = number of tracers (temperature, salinity, ...)
!
!     mnisle = maximum number of islands (unconnected land masses)
!
!     maxipp = maximum number of all island perimeter points
!-----------------------------------------------------------------------
!                                                                       
#include "derived_options.h"
#ifdef bbl_ag
      parameter (imt=122, jmt=66, km=16)
#else
!MODIFY--------------------------------------
!+Modified by Sasaki
!      parameter (imt=3602, jmt=1502, km= 54)
! for OIFES of the north Pacific
!      parameter (imt=762, jmt=354, km= 54)
! for OIFES of the north Pacific with 0.1 deg
!      parameter (imt=1902, jmt=882, km= 54)
! for OIFES of the north Pacific with 0.5 deg
!      parameter (imt=382, jmt=178, km= 54)

!==========================================
!=======   USER INPUT  from here  =========
!==========================================

      parameter (imt= 622, jmt= 402, km= 100)

!==========================================
!=======   USER INPUT  until here  ========
!==========================================


!!! for idealized domain
!!      parameter (imt=102, jmt=82, km= 54)


!--------------------------------------------
#endif
      parameter (nt=2)
!+Modified by H.Sasaki 
      parameter (mnisle=5000, maxipp=500000)
!ssk      parameter (mnisle=500, maxipp=500000)
#ifdef bbl_ag
      parameter (kbot = km-1)
#else
      parameter (kbot = km)
#endif
!
#ifdef obctest
!
!     open boundary tests for Stevens (1990)
!                                                                       
# if defined obc_south || defined obc_north
      parameter (imt= 21, jmt= 24, km= 6)
# else
      parameter (imt= 21, jmt= 41, km= 6)
# endif
#endif
#ifdef obctest2
# ifdef obc_south
      parameter (imt= 21, jmt= 11, km= 6)
# else
      parameter (imt= 21, jmt= 21, km= 6)
# endif
#endif
!
!-----------------------------------------------------------------------
!     END USER INPUT
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     "jmw"  = number of latitude rows within the MW (memory window).
!     "jbuf" = number of northermost and southernmost buffer rows.
!              The MW is symmetric because it has the same number of
!              northermost and southernmost buffer rows "jbuf".
!              total number of buffer rows = 2*jbuf
!-----------------------------------------------------------------------
!
#if defined second_order_window
      parameter (jbuf=1)
#endif
#if defined fourth_order_window
      parameter (jbuf=2)
#endif
#if defined sixth_order_window
      parameter (jbuf=3)
#endif
#if defined max_window
!   note: jmw is a function of the number of processors
!         jmw is set appropriately for 1 processor below
      parameter (jmw=(jmt-2) + 2*jbuf)
#else
# if defined xlandmix
      parameter (jmw=3 + 2*jbuf)
# else
!MODIFY---------------------------------
!+Modified by Masumoto
!      parameter (jmw=70 + 2*jbuf)
!+Modified by Sasaki
      parameter (jmw=75 + 2*jbuf)
!-Original------------------------------
!     parameter (jmw=1 + 2*jbuf)
!---------------------------------------
# endif                                                                       
#endif                                                                       
!
!-----------------------------------------------------------------------
!      num_rows_to_copy = number of northernmost rows to copy into
!                         the southernmost rows of the MW when moving
!                         the MW northward
!-----------------------------------------------------------------------
!
      parameter (num_rows_to_copy = 2*jbuf)
