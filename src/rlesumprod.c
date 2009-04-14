
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP 
rle_sum_prod(SEXP px1, SEXP pn1, SEXP ps1, 
	     SEXP px2, SEXP pn2, SEXP ps2, 
	     SEXP plen)
{
    double 
	ans,
	*x1 = REAL(px1),
	*x2 = REAL(px2);
    int 
	i1, i2, k, next_k,
	*n1 = INTEGER(pn1),
	*n2 = INTEGER(pn2),
	*s1 = INTEGER(ps1),
	*s2 = INTEGER(ps2),
	len = asInteger(plen);
    i1 = i2 = k = 0; /* i: RLE index, k: virtual index for underlying seq */
    ans = 0;
    while (k < len) {
        if (x1[i1] == 0 || x2[i2] == 0) { /* either 0, no contribution */
            i1 = i1++;
            i2 = i2++;
            /* move lagging pointer ahead to location of the one ahead */
            while (s1[i1-1] < s2[i2-1]) i1 = i1++;
            while (s1[i1-1] > s2[i2-1]) i2 = i2++;
            k = 1 + imax2(s1[i1-1], s2[i2-1]);
        }
        else {
	    next_k = 1 + imin2(s1[i1], s2[i2]);
	    ans += (next_k - k) * x1[i1] * x2[i2];
	    if (s1[i1] == (next_k - 1)) i1 = i1++;
	    if (s2[i2] == (next_k - 1)) i2 = i2++;
	    k = next_k;
        }
    }
    return ScalarReal(ans);
}



SEXP 
rle_sum_any(SEXP px1, SEXP pn1, SEXP ps1, 
	    SEXP px2, SEXP pn2, SEXP ps2, 
	    SEXP plen)
{
    int 
	i1, i2, k, next_k, ans,
	*x1 = INTEGER(px1),
	*x2 = INTEGER(px2),
	*n1 = INTEGER(pn1),
	*n2 = INTEGER(pn2),
	*s1 = INTEGER(ps1),
	*s2 = INTEGER(ps2),
	len = asInteger(plen);
    i1 = i2 = k = 0; /* i: RLE index, k: virtual index for underlying seq */
    ans = 0;
    while (k < len) {
	next_k = 1 + imin2(s1[i1], s2[i2]);
	if (x1[i1] || x2[i2]) ans += (next_k - k);
	if (s1[i1] == (next_k - 1)) i1 = i1++;
	if (s2[i2] == (next_k - 1)) i2 = i2++;
	k = next_k;
    }
    return ScalarInteger(ans);
}



SEXP 
do_naive_density(SEXP x, SEXP dk, SEXP width)
{
    SEXP ans;
    double *pans;
    int 
	i, j, xi,
	n = length(x), 
	m = length(dk),
	x1 = INTEGER(x)[0], 
	xn = INTEGER(x)[n-1];
    ans = PROTECT(allocVector(REALSXP, xn - x1 + 2 + 2 * asInteger(width)));
    pans = REAL(ans);
    for (i = 0; i < length(ans); i++) pans[i] = 0;
    /* for each point xi in x, add dk suitably shifted */
    for (i = 0; i < n; i++) { 	
	xi = INTEGER(x)[i] - x1 + 1;
	for (j = 0; j < m; j++) {
	    pans[xi + j] += REAL(dk)[j];
	}
    }
    UNPROTECT(1);
    return ans;
}

