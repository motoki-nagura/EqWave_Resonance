dset ^kmt.dta.bin
undef 0.0
*
xdef 720 linear   0.25 0.5
ydef 216 linear -64.75 0.5
zdef    1 levels     2.50
tdef    1 linear 01jan1979 1yr
*
options big_endian
vars  1
kmt    1 0 Number of levels
endvars

