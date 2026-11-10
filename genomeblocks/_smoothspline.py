"""R's ``smooth.spline(x, y, cv = TRUE)``, ported to numpy.

FitHiChIP fits its distance and bias-regression curves with R's cubic smoothing
spline, whose smoothing parameter is chosen by leave-one-out cross-validation.
scipy's smoothers use other knots, another penalty scale and GCV, so they give
different curves; this module follows R step by step instead:
``smooth.spline`` (ties, knot choice ``.nknots.smspl``), ``sbart.c`` (the
penalty ratio, the Brent search for ``spar`` in [-1.5, 1.5]), ``sgram.f``,
``stxwx.f``, ``sslvrg.f`` and ``sinerp.f`` (Gram matrix, normal equations,
banded Cholesky solve, leverages, CV criterion), and
``predict.smooth.spline`` (linear extrapolation), keeping the Fortran
summation order. As in R, the returned fit is the one of the last criterion
evaluation of the search. On FitHiChIP's curves predictions agree with R
(4.1–4.5) to 1e-9 relative or better; where the CV criterion is nearly flat
(near-interpolating fits, spar < -1) the search can stop at another spar.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


def _interv(xt, lxt, x):
    """R ``interv`` (rightmost_closed = all_inside = 0): 1-based ``ileft`` with
    ``xt[ileft] <= x < xt[ileft+1]``, and ``mflag`` (-1 left of, 1 right of the knots)."""
    i = int(np.searchsorted(xt[:lxt], x, side="right"))
    if i == 0:
        return 1, -1
    if i >= lxt:
        return lxt, 1
    return i, 0


def _bsplvb(t, x, left, k=4):
    """{order m: values of the m B-splines of order m nonzero on [t[left], t[left+1]]} (de Boor)."""
    out = {1: [1.0]}
    biatx = [1.0]
    deltar, deltal = [], []
    for j in range(1, k):
        deltar.append(t[left + j] - x)
        deltal.append(x - t[left + 1 - j])
        saved = 0.0
        new = []
        for i in range(j):
            term = biatx[i] / (deltar[i] + deltal[j - 1 - i])
            new.append(saved + deltar[i] * term)
            saved = deltal[j - 1 - i] * term
        new.append(saved)
        biatx = new
        out[j + 1] = biatx
    return out


def _bsplvd(t, x, left1, nderiv, k=4):
    """de Boor's ``bsplvd``: column ``m`` = (m)-th derivative of the k B-splines nonzero at ``x`` on the interval
    with 1-based index ``left1``, in Fortran's arithmetic order."""
    db = np.zeros((k, max(nderiv, 1)))
    st = {"j": 1, "dl": [0.0] * 20, "dr": [0.0] * 20}

    def bsplvb(jhigh, index):                 # writes into db[:, 0], keeps its state between calls (Fortran SAVE)
        b = db[:, 0]
        if index != 2:
            st["j"] = 1
            b[0] = 1.0
            if st["j"] >= jhigh:
                return
        while True:
            j = st["j"]
            st["dr"][j - 1] = t[left1 + j - 1] - x
            st["dl"][j - 1] = x - t[left1 - j]
            saved = 0.0
            for i in range(1, j + 1):
                term = b[i - 1] / (st["dr"][i - 1] + st["dl"][j - i])
                b[i - 1] = saved + st["dr"][i - 1] * term
                saved = st["dl"][j - i] * term
            b[j] = saved
            st["j"] = j + 1
            if st["j"] >= jhigh:
                return

    mhigh = max(min(nderiv, k), 1)
    kp1 = k + 1
    bsplvb(kp1 - mhigh, 1)
    if mhigh == 1:
        return db
    ideriv = mhigh
    for _ in range(2, mhigh + 1):
        jp1mid = 1
        for j in range(ideriv, k + 1):
            db[j - 1, ideriv - 1] = db[jp1mid - 1, 0]
            jp1mid += 1
        ideriv -= 1
        bsplvb(kp1 - ideriv, 2)
    a = np.zeros((k, k))
    for i in range(k):
        a[i, i] = 1.0
    for m in range(2, mhigh + 1):
        kp1mm = kp1 - m
        fkp1mm = float(kp1mm)
        il, i = left1, k
        for _ in range(kp1mm):
            factor = fkp1mm / (t[il + kp1mm - 1] - t[il - 1])
            for j in range(1, i + 1):
                a[i - 1, j - 1] = (a[i - 1, j - 1] - a[i - 2, j - 1]) * factor
            il -= 1
            i -= 1
        for i in range(1, k + 1):
            s = 0.0
            for j in range(max(i, m), k + 1):
                s = a[j - 1, i - 1] * db[j - 1, m - 1] + s
            db[i - 1, m - 1] = s
    return db


def _bvalue(t, bcoef, n, x, jderiv=0, k=4):
    """de Boor's ``bvalue`` (R's version: the right end of the knots belongs to the last interval), vectorised
    over ``x`` in [t[k-1], t[n]]: value or ``jderiv``-th derivative of the spline at each x."""
    x = np.asarray(x, float)
    i = np.empty(len(x), np.int64)
    for r, xv in enumerate(x):
        if xv == t[n] and t[n] == t[n + k - 1]:
            i[r] = n
        else:
            i[r], mflag = _interv(t, n + k, xv)
            if mflag != 0:
                raise ValueError("bvalue: x outside the knots")
    km1 = k - 1
    imk = i - k
    if (imk < 0).any() or (n - i < 0).any():
        raise ValueError("bvalue: x outside the spline's support")
    dm = [None] + [x - t[i + 1 - j - 1] for j in range(1, k)]
    dp = [None] + [t[i + j - 1] - x for j in range(1, k)]
    aj = [None] + [bcoef[imk + jc - 1] for jc in range(1, k + 1)]
    for j in range(1, jderiv + 1):
        kmj = k - j
        ilo = kmj
        for jj in range(1, kmj + 1):
            aj[jj] = ((aj[jj + 1] - aj[jj]) / (dm[ilo] + dp[jj])) * float(kmj)
            ilo -= 1
    for j in range(jderiv + 1, km1 + 1):
        kmj = k - j
        ilo = kmj
        for jj in range(1, kmj + 1):
            aj[jj] = (aj[jj + 1] * dm[ilo] + aj[jj] * dp[jj]) / (dm[ilo] + dp[jj])
            ilo -= 1
    return aj[1]


def _sgram(tb, nb):
    """The four nonzero diagonals of SIGMA[i, j] = integral B_i'' B_j'' (``sgram.f``, its 0.3330 included)."""
    sg = np.zeros((4, nb))
    for i in range(nb):
        ileft, _ = _interv(tb, nb + 1, tb[i])
        yw1 = _bsplvd(tb, tb[i], ileft, 3)[:, 2].copy()
        yw2 = _bsplvd(tb, tb[i + 1], ileft, 3)[:, 2] - yw1
        wpt = tb[i + 1] - tb[i]
        for ii in range(4):
            for dd in range(4 - ii):
                jj = ii + dd
                sg[dd, ileft - 4 + ii] += wpt * (yw1[ii] * yw1[jj] + (yw2[ii] * yw1[jj] + yw2[jj] * yw1[ii]) * 0.5
                                                 + yw2[ii] * yw2[jj] * 0.3330)
    return sg


def _stxwx(x, z, w, knot, n):
    """X'W z and the four nonzero diagonals of X'W X (``stxwx.f``; ``w`` are square-rooted weights)."""
    y = np.zeros(n)
    hs = np.zeros((4, n))
    for i in range(len(x)):
        ileft, mflag = _interv(knot, n + 1, x[i])
        if mflag == 1:
            if x[i] <= knot[ileft - 1] + 1e-10:
                ileft -= 1
            else:
                return y, hs
        v = _bsplvd(knot, x[i], ileft, 1)[:, 0]
        j = ileft - 4
        ww = w[i] ** 2
        y[j:j + 4] += ww * z[i] * v
        for r in range(4):
            hs[0, j + r] += ww * v[r] ** 2
            if r < 3:
                hs[1, j + r] += ww * v[r] * v[r + 1]
            if r < 2:
                hs[2, j + r] += ww * v[r] * v[r + 2]
        hs[3, j] += ww * v[0] * v[3]
    return y, hs


def _sinerp(abd, nk):
    """Banded inner products of the columns of R^-1 (``sinerp.f``, pass 1), for the leverages."""
    p1ip = np.zeros((4, nk))
    wjm3, wjm2, wjm1 = [0.0, 0.0, 0.0], [0.0, 0.0], [0.0]
    for i in range(nk):
        j = nk - 1 - i
        c0 = 1.0 / abd[3, j]
        c1 = abd[0, j + 3] * c0 if j <= nk - 4 else 0.0
        c2 = abd[1, j + 2] * c0 if j <= nk - 3 else 0.0
        c3 = abd[2, j + 1] * c0 if j <= nk - 2 else 0.0
        p1ip[0, j] = 0.0 - (c1 * wjm3[0] + c2 * wjm3[1] + c3 * wjm3[2])
        p1ip[1, j] = 0.0 - (c1 * wjm3[1] + c2 * wjm2[0] + c3 * wjm2[1])
        p1ip[2, j] = 0.0 - (c1 * wjm3[2] + c2 * wjm2[1] + c3 * wjm1[0])
        p1ip[3, j] = (c0 ** 2 + c1 ** 2 * wjm3[0] + 2.0 * c1 * c2 * wjm3[1] + 2.0 * c1 * c3 * wjm3[2]
                      + c2 ** 2 * wjm2[0] + 2.0 * c2 * c3 * wjm2[1] + c3 ** 2 * wjm1[0])
        wjm3 = [wjm2[0], wjm2[1], p1ip[1, j]]
        wjm2 = [wjm1[0], p1ip[2, j]]
        wjm1 = [p1ip[3, j]]
    return p1ip


def _dpbfa(abd, n, m=3):
    """LINPACK ``dpbfa``: Cholesky R'R of a banded SPD matrix in place (upper band storage, row ``m`` the
    diagonal), with LINPACK's summation order. Returns 0, or the 1-based column where it is not positive."""
    for j in range(n):
        s = 0.0
        ik = m
        jk = max(j - m, 0)
        mu = max(m + 1 - (j + 1), 0)
        for k in range(mu, m):
            dot = 0.0
            for r in range(k - mu):
                dot += abd[ik + r, jk] * abd[mu + r, j]
            t = (abd[k, j] - dot) / abd[m, jk]
            abd[k, j] = t
            s += t * t
            ik -= 1
            jk += 1
        s = abd[m, j] - s
        if s <= 0.0:
            return j + 1
        abd[m, j] = np.sqrt(s)
    return 0


def _dpbsl(abd, n, b, m=3):
    """LINPACK ``dpbsl``: solve R'R x = b with the factor of ``_dpbfa`` (in place on ``b``)."""
    for k in range(n):
        lm = min(k, m)
        la, lb = m - lm, k - lm
        t = 0.0
        for r in range(lm):
            t += abd[la + r, k] * b[lb + r]
        b[k] = (b[k] - t) / abd[m, k]
    for k in range(n - 1, -1, -1):
        lm = min(k, m)
        la, lb = m - lm, k - lm
        b[k] = b[k] / abd[m, k]
        t = -b[k]
        for r in range(lm):
            b[lb + r] += t * abd[la + r, k]
    return b


def _nknots(n):
    """R ``.nknots.smspl``."""
    if n < 50:
        return n
    a1, a2, a3, a4 = np.log2(50), np.log2(100), np.log2(140), np.log2(200)
    if n < 200:
        v = 2 ** (a1 + (a2 - a1) * (n - 50) / 150)
    elif n < 800:
        v = 2 ** (a2 + (a3 - a2) * (n - 200) / 600)
    elif n < 3200:
        v = 2 ** (a3 + (a4 - a3) * (n - 800) / 2400)
    else:
        v = 200 + (n - 3200) ** 0.2
    return int(np.trunc(v))


@dataclass
class SmoothSpline:
    """A fitted cubic smoothing spline: B-spline coefficients on knots over [0, 1] (x rescaled)."""
    knot: np.ndarray
    coef: np.ndarray
    xmin: float
    xrange: float
    spar: float
    df: float
    crit: float

    def predict(self, x) -> np.ndarray:
        """Values at ``x``, linear beyond the data range (R ``predict.smooth.spline``)."""
        x = np.asarray(x, float)
        nk = len(self.coef)
        xs = (x - self.xmin) / self.xrange
        y = xs.copy()
        left, right = xs < 0, xs > 1
        inside = ~(left | right)
        if inside.any():
            u, inv = np.unique(xs[inside], return_inverse=True)
            y[inside] = _bvalue(self.knot, self.coef, nk, u)[inv]
        if (~inside).any():
            ends = (np.array([self.xmin, self.xmin + self.xrange]) - self.xmin) / self.xrange
            val = _bvalue(self.knot, self.coef, nk, ends)
            slope = (_bvalue(self.knot, self.coef, nk, ends, 1) / self.xrange) * self.xrange
            y[left] = val[0] + slope[0] * (xs[left] - 0)
            y[right] = val[1] + slope[1] * (xs[right] - 1)
        return y


def smooth_spline(x, y, *, spar: Optional[float] = None) -> SmoothSpline:
    """R ``smooth.spline(x, y, cv = TRUE)``: unweighted, default knots; ``spar`` is
    searched in [-1.5, 1.5] by leave-one-out CV unless given. Tied x values are
    averaged with weight = their count, as R does."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("missing or infinite values in inputs are not allowed")
    n = len(x)
    q1, q3 = np.percentile(x, [25, 75])
    xx = np.round((x - x.mean()) / (1e-6 * (q3 - q1)))
    uxx, first = np.unique(xx, return_index=True)
    ux = np.sort(x[first])
    nx = len(ux)
    if nx <= 3:
        raise ValueError("need at least four unique 'x' values")
    if nx == n:
        wbar = np.ones(n)
        ybar = y[np.argsort(x, kind="stable")]
    else:
        ox = np.searchsorted(uxx, xx)
        wbar = np.bincount(ox, minlength=nx).astype(float)
        ysum = np.zeros(nx, np.longdouble)                    # R's sum() accumulates in long double
        np.add.at(ysum, ox, y.astype(np.longdouble))
        ybar = ysum.astype(float) / np.where(wbar > 0, wbar, 1)
    r_ux = ux[-1] - ux[0]
    xbar = (ux - ux[0]) / r_ux
    nknots = _nknots(nx)
    sel = 1 + np.arange(nknots) * ((nx - 1) / (nknots - 1)) if nknots > 1 else np.array([1.0])
    sel[-1] = nx
    sel = np.trunc(sel).astype(int) - 1
    knot = np.concatenate([np.repeat(xbar[0], 3), xbar[sel], np.repeat(xbar[-1], 3)])
    nk = nknots + 2

    ws = np.where(wbar > 0, np.sqrt(wbar), wbar)
    sg = _sgram(knot, nk)
    xwy, hs = _stxwx(xbar, ybar, ws, knot, nk)
    ratio = np.cumsum(hs[0, 2:nk - 3])[-1] / np.cumsum(sg[0, 2:nk - 3])[-1]

    # sslvrg(): B-spline values at each x for the leverages, the ends nudged inside by 1e-11
    jl = np.empty(nx, np.int64)
    Bl = np.empty((nx, 4))
    for i, xv in enumerate(xbar):
        ileft, mflag = _interv(knot, nk + 1, xv)
        if mflag == -1:
            ileft, xv = 4, knot[3] + 1e-11
        elif mflag == 1:
            ileft, xv = nk, knot[nk] - 1e-11
        jl[i] = ileft - 4
        Bl[i] = _bsplvd(knot, xv, ileft, 1)[:, 0]
    b0, b1, b2, b3 = Bl.T
    state = {"ier": 0}

    def sslvrg(sp):
        lam = ratio * 16.0 ** (sp * 6.0 - 2.0)
        abd = np.zeros((4, nk))
        abd[3] = hs[0] + lam * sg[0]
        abd[2, 1:] = hs[1, :-1] + lam * sg[1, :-1]
        abd[1, 2:] = hs[2, :-2] + lam * sg[2, :-2]
        abd[0, 3:] = hs[3, :-3] + lam * sg[3, :-3]
        if _dpbfa(abd, nk):
            state["ier"] = 1
            return np.inf
        coef = _dpbsl(abd, nk, xwy.copy())
        sz = _bvalue(knot, coef, nk, xbar)
        p = _sinerp(abd, nk)
        j = jl
        lev = (p[3, j] * b0 ** 2 + 2.0 * p[2, j] * b0 * b1 + 2.0 * p[1, j] * b0 * b2 + 2.0 * p[0, j] * b0 * b3
               + p[3, j + 1] * b1 ** 2 + 2.0 * p[2, j + 1] * b1 * b2 + 2.0 * p[1, j + 1] * b1 * b3
               + p[3, j + 2] * b2 ** 2 + 2.0 * p[2, j + 2] * b2 * b3 + p[3, j + 3] * b3 ** 2) * ws ** 2
        crit = np.cumsum((((ybar - sz) * ws) / (1 - lev)) ** 2)[-1] / nx      # Fortran's left-to-right sum
        state.update(coef=coef, lev=lev)
        return crit

    if spar is not None:
        crit = sslvrg(float(spar))
        xk = float(spar)
    else:
        xk, crit = _brent(sslvrg, state)
    if state["ier"]:
        raise ValueError("smoothing parameter value too small or too large")
    return SmoothSpline(knot=knot, coef=state["coef"], xmin=float(ux[0]), xrange=float(r_ux), spar=float(xk),
                        df=float(state["lev"].sum()), crit=float(crit))


def _brent(f, state, a=-1.5, b=1.5, tol=1e-4, eps=2e-8, maxit=500):
    """``sbart.c``'s golden-section / parabolic minimisation of the criterion over spar."""
    c_gold = 0.381966011250105151795413165634
    big = 1e100
    v = w = x = a + c_gold * (b - a)
    e = d = 0.0
    fx = fv = fw = f(x)
    it = 0
    while state["ier"] == 0:
        xm = (a + b) * 0.5
        tol1 = eps * abs(x) + tol / 3.0
        tol2 = tol1 * 2.0
        it += 1
        if abs(x - xm) <= tol2 - (b - a) * 0.5 or it > maxit:
            break
        golden = abs(e) <= tol1 or fx >= big or fv >= big or fw >= big
        if not golden:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = (q - r) * 2.0
            if q > 0.0:
                p = -p
            q = abs(q)
            r = e
            e = d
            if abs(p) >= abs(0.5 * q * r) or q == 0.0 or p <= q * (a - x) or p >= q * (b - x):
                golden = True
            else:
                d = p / q
                u = x + d
                if u - a < tol2 or b - u < tol2:
                    d = np.copysign(tol1, xm - x)
        if golden:
            e = (a - x) if x >= xm else (b - x)
            d = c_gold * e
        u = x + (d if abs(d) >= tol1 else np.copysign(tol1, d))
        fu = f(u)
        if not np.isfinite(fu):
            fu = 2.0 * big
        if fu <= fx:
            if u >= x:
                a = x
            else:
                b = x
            v, fv = w, fw
            w, fw = x, fx
            x, fx = u, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v, fv = w, fw
                w, fw = u, fu
            elif fu <= fv or v == x or v == w:
                v, fv = u, fu
    return x, fx
