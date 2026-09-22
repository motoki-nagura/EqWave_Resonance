import numpy as np

def reorder(k, ik_lookup):
    """ Sort in terms of ``k'' """
    i = np.argsort(k)
    ik_lookup = ik_lookup[i,:][:,i]
    k = k[i]
    return k, ik_lookup


def generate_k3(ik3_lookup,k):
    k3 = np.full(ik3_lookup.shape, np.nan)
    for ik1 in range(len(ik3_lookup[:,0])):
        for ik2 in range(len(ik3_lookup[0,:])):
            if ik3_lookup[ik1,ik2] >= 0:
                k3[ik1,ik2] = k[ik3_lookup[ik1,ik2]]
    return k3


def check_ω1ω2ω3(txt, ω1s=None,ω2s=None,ω3s=None):
    iflag = 0

    for i1, ω1 in enumerate(ω1s):
        for i2, ω2 in enumerate(ω2s):
            ω3 = ω3s[i1,i2]
            if np.logical_not( np.isnan(ω3) ):
                if np.isclose(ω1+ω2,ω3,rtol=1e-5):
                    continue
                else:
                    print(txt+f"1={ω1:.3e}, "+txt+f"2={ω2:.3e}, "+txt+f"3={ω3:.3e}, "+
                          txt+"1+"+txt+f"2={ω1+ω2:.3e}")
                    iflag = 1

    if iflag == 0:
        print(txt+"1+"+txt+"2="+txt+"3 is satisfied")
