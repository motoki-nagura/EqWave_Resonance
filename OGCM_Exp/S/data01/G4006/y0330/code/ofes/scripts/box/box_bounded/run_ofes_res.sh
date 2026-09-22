#!/bin/sh

#======================================================================
#   set namelist variables
#======================================================================
#
# running length = 3.0 month
# &contrl init=.true., runlen=4.0, rununits='months',
# &contrl init=.true., runlen=$days, rununits='days',
#days=1.0
diag=1.0
cat > namelist.contrl << ENDNAMELIST
 &contrl init=.false., runlen=12.0, rununits='months',
         restrt=.true., initpt=.true.
 /
ENDNAMELIST

cat > namelist.mbcin << ENDNAMELIST
 &mbcin   /
ENDNAMELIST

#======================================================================
# time steps
#======================================================================
cat > namelist.tsteps << ENDNAMELIST
 &tsteps dtts=200., dtuv=200., dtsf=4.,
 /
ENDNAMELIST

cat > namelist.riglid << ENDNAMELIST
 &riglid mxscan=2000, tolrsf=1.0e8, tolrsp=1.0e-4,
         tolrfs=1.0e-4,
 /
ENDNAMELIST

# &mixing am=0.0e0, ah=0.0e0, ambi=1.e22, ahbi=1.e15,
cat > namelist.mixing << ENDNAMELIST
 &mixing am=0.0e0, ah=0.0e0, ambi=2.7e18, ahbi=9.e17,
         kappa_m=10.0, kappa_h=1.0, aidif=1.0,
         nmix=17, eb=.true.,  ncon=1, cdbot=2.5e-3,
         acor=0.0, dampts=10.0, 15.0
 /
ENDNAMELIST

cat > namelist.isopyc << ENDNAMELIST
 &isopyc ahisop=1.e6, slmx=0.01, athkdf=1.e6,   
         ahsteep=1.e6,
 /
ENDNAMELIST

cat > namelist.bbl << ENDNAMELIST
 &bbl raydrag=2.0e-5, cdbotbbl=3.0e-3, kvbbltop=200.0, entrain_kbbl=0.0,
 /
ENDNAMELIST

cat > namelist.blmix << ENDNAMELIST
 &blmix   
 /
ENDNAMELIST

cat > namelist.ncdiff << ENDNAMELIST
 &ncdiff dbot=2000.e2, dtop=100.e2, 
         diffmin=1.e5, diffmax=1.e8,  
         vmhs_alpha=0.015, vmhs_rate2_limit=1.96e-12,
         ijvmhs_maxlen=10,diffint=4.0, diffsnapint=1.000
 /
ENDNAMELIST

cat > namelist.kppmix << ENDNAMELIST
 &kppmix
  lri=.true., ldd=.false.,
  visc_con_limit=50.0, diff_con_limit=50.0,
  visc_cbu_iw=1.0, diff_cbt_iw=0.0,
 /
ENDNAMELIST

cat > namelist.tcvmix << ENDNAMELIST
 &tcvmix1
  diff_cbt_back=0.01,  visc_cbu_back=0.1,
  visc_cbu_limit=1.e4, diff_cbt_limit=1.e4, tune_noh_a=30., tune_noh_b=0.0,
 /
ENDNAMELIST

cat > namelist.smagnl << ENDNAMELIST
 &smagnl k_smag=4, diff_c_back=0.0, visc_c_back=0.0, prandtl=1.0
 /
ENDNAMELIST

#  ----------> From here (added, nagura, 2025/2/6)
#  ----------   wndmix:  10.0 --> 0.0  (2025/2/6)
cat > namelist.ppmix << ENDNAMELIST
 &ppmix wndmix=0.0,  fricmx=50.0, diff_cbt_back=0.01,  visc_cbu_back=0.1,
        visc_cbu_limit=1.e4, diff_cbt_limit=1.e4
 /
ENDNAMELIST
#  <---------- Until here (added, nagura, 2025/2/6)

#  ----------> From here (commented out, nagura, 2025/2/6)
#cat > namelist.ppmix << ENDNAMELIST
# &ppmix wndmix=10.0,  
# /
#ENDNAMELIST
#  <---------- Until here (commented out, nagura, 2025/2/6)

cat > namelist.diagn << ENDNAMELIST
 &diagn  tsiint=0.0, 
         tavgint=$diag,  itavg=.true.,
         tmbint=$diag,   itmb=.true., tmbper=$diag,
         trmbint=$diag,  itrmb=.true., 
         gyreint=$diag,  igyre=.true., 
	 snapint=$diag,
	 timavgint=1.000, timavgper=1.000,
         glenint=$diag,  vmsfint=$diag, stabint=$diag, zmbcint=$diag,
         extint=$idag,   prxzint=$diag, dspint=$diag, dspper=$diag,
         trajint=$idag,  xbtint=$diag, xbtper=$diag, exconvint=$diag, cmixint=$diag,
         crossint=$diag, pressint=$diag, fctint=$diag, tyzint=$diag, rhoint=$diag
 /
ENDNAMELIST

cat > namelist.io << ENDNAMELIST
 &io     iotavg=-2, iotmb=-2,   iotrmb=-2, iozmbc=-2,
         ioglen=-2, iovmsf=-2, iogyre=-2,
         ioprxz=-2, ioext=-2,   iodsp=-2,
         iotsi=-2,   iotraj=-2, ioxbt=-2,
         iosnap=-2,
         restart_in_SINGLE=.false.,
         restart_out_SINGLE=.false.,
         init3d_file_SINGLE=.false.,
         obc_file_SINGLE=.false.,
 /         
ENDNAMELIST

cat > namelist.ictime << ENDNAMELIST
 &ictime eqyear=.true., eqmon=.false., refinit=.true., 
         year0=1990, month0=1, day0=1, hour0=0, min0=0, sec0=0,
 /         
ENDNAMELIST

cat > namelist.ice << ENDNAMELIST
 &ice    pstar=10.0e3,sref=0.0348,cdw=5.5,
         evp_sub_steps=10,
         heff_min=0.01, gamma_io=0.00965,
         diff1=0.004, diff2=0.004,
 /         
ENDNAMELIST

cat > namelist.parNPZD << ENDNAMELIST
 &parNPZD alphabio = 0.025,    abio  = 0.6,    bbio  = 1.066,
          cbio  = 1.0,         gbio  = 2.0,    epsbio= 1.0,
          a_npz = 0.75,        d_npz = 0.03,   parbio= 0.43,
          rk1bio= 0.5,         phiphy= 0.05,   phiphyq = 0.05,
          phizoo= 0.20,
          remina= 0.05,        w_detr= 5.0,   dkcbio= 0.03,
          dbio  = 25.,         kmeuph=28,      d_mrtn= 10000.0,
          c_mrtn = 0.858,
 /
ENDNAMELIST

### END OF NAMELIST OUTPUT ###

./oifes.LM

exit
