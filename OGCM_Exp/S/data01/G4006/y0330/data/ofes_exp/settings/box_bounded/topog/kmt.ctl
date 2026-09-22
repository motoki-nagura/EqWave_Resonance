dset ^kmt.dta.bin
undef 0.0
*
xdef 620 linear   0.05 0.1
ydef 400 linear -19.95 0.1
zdef   1 levels     2.50
tdef   1 linear 01jan1979 1yr
*
options big_endian
vars  1
kmt    1 0 Number of levels
endvars

