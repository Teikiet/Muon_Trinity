# Derive the basis from the ACTUAL callback: alpha(h) = (tau(h+25)-tau(h-25))/(h_hi-h_lo)
# with h_lo=max(0,h-25), h_hi=min(50000,h+25) and a max(0,.) clamp.
p = '/uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat'
tok = open(p).read().split(); it = iter(tok)
nl = int(next(it)); nh = int(next(it))
wl = [float(next(it)) for _ in range(nl)]; hh = [float(next(it)) for _ in range(nh)]
ec = [[float(next(it)) for _ in range(nh)] for _ in range(nl)]
HMIN, HMAX = hh[0], hh[-1]
def tau(lam, h):
    lam = max(wl[0], min(wl[-1], lam)); h = max(HMIN, min(HMAX, h))
    iw = 0
    while iw < nl-2 and wl[iw+1] < lam: iw += 1
    ih = 0
    while ih < nh-2 and hh[ih+1] < h: ih += 1
    fw = 0.0 if wl[iw+1]==wl[iw] else (lam-wl[iw])/(wl[iw+1]-wl[iw])
    fh = 0.0 if hh[ih+1]==hh[ih] else (h-hh[ih])/(hh[ih+1]-hh[ih])
    return (ec[iw][ih]*(1-fw)+ec[iw+1][ih]*fw)*(1-fh) + (ec[iw][ih+1]*(1-fw)+ec[iw+1][ih+1]*fw)*fh
def alphaCB(lam, h):
    lo = max(HMIN, h-25.0); hi = min(HMAX, h+25.0)
    a = hi - lo
    if a <= 0.0: return 0.0
    return max(0.0, (tau(lam, hi) - tau(lam, lo)) / a)
def slope(k, lam):
    return (tau(lam, hh[k+1]) - tau(lam, hh[k])) / (hh[k+1] - hh[k])
def phi(k, h):
    a, b, w = hh[k], hh[k+1], 25.0
    if h <= a-w or h >= b+w: return 0.0
    if h < a+w: return (h-(a-w))/(2*w)
    if h > b-w: return ((b+w)-h)/(2*w)
    return 1.0
lam = 400.0
print('POINTWISE reconstruction, lam=%.0f nm' % lam)
print('  alpha_cb(h) vs sum_k s_k phi_k(h)  [interior ramp basis]')
print('%10s %16s %16s %14s' % ('h [m]', 'alpha_cb', 'basis', 'diff'))
for h in (0.0, 5.0, 12.5, 24.9, 25.1, 500.0, 987.5, 1000.0, 1012.5, 1025.0, 5000.0, 49975.0, 49987.5, 49999.0, 50000.0):
    b = sum(slope(k, lam)*phi(k, h) for k in range(nh-1))
    a = alphaCB(lam, h)
    print('%10.1f %16.9e %16.9e %14.3e' % (h, a, b, abs(a-b)))
print()
print('ratio basis/alpha_cb near the ground (should be 1.0 if the basis is right):')
for h in (0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 24.0, 25.0, 30.0):
    b = sum(slope(k, lam)*phi(k, h) for k in range(nh-1))
    a = alphaCB(lam, h)
    print('   h=%6.1f  alpha_cb=%.6e  basis=%.6e  ratio=%.6f' % (h, a, b, b/a if a>0 else float('nan')))
