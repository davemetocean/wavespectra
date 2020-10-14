/* specpart.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__600 = 600;
static integer c__3 = 3;
static integer c__200 = 200;

/* -*- f90 -* */
doublereal minval_(real *arr, integer *n)
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --arr;

    /* Function Body */
    ret_val = arr[1];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (arr[i__] < ret_val) {
	    ret_val = arr[i__];
	}
    }
    return ret_val;
} /* minval_ */

doublereal maxval_(real *arr, integer *n)
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --arr;

    /* Function Body */
    ret_val = arr[1];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (arr[i__] > ret_val) {
	    ret_val = arr[i__];
	}
    }
    return ret_val;
} /* maxval_ */

integer iminval_(integer *arr, integer *n)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --arr;

    /* Function Body */
    ret_val = arr[1];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (arr[i__] < ret_val) {
	    ret_val = arr[i__];
	}
    }
    return ret_val;
} /* iminval_ */

/* Subroutine */ int partition_(real *spec, integer *ipart, integer *nk, 
	integer *nth)
{
    /* System generated locals */
    integer spec_dim1, spec_offset, ipart_dim1, ipart_offset, i__1, i__2, 
	    i__3, i__4, i__5;
    real r__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer i_nint(real *);

    /* Local variables */
    static integer ik, mk;
    static real zp[600];
    static integer ind[600], imi[600], imo[600], mth;
    static real fact;
    static integer iang;
    static real zmin, zmax;
    static integer neigh[5400]	/* was [9][600] */, nspec, npart;
    extern /* Subroutine */ int pt_fld__(integer *, real *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), ptnghb_(
	    integer *, integer *, integer *, integer *);
    extern doublereal minval_(real *, integer *), maxval_(real *, integer *);
    extern /* Subroutine */ int ptsort_(integer *, integer *, integer *, 
	    integer *);

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 6, 0, 0, 0 };


/*     ---------------------------------------------------------------- */
/*       imi     i.a.   i   input discretized spectrum. */
/*       ind     i.a.   i   sorted addresses. */
/*       imo     i.a.   o   output partitioned spectrum. */
/*       zp      r.a.   i   spectral array. */
/*       npart   int.   o   number of partitions found. */
/*     ---------------------------------------------------------------- */
    /* Parameter adjustments */
    ipart_dim1 = *nk;
    ipart_offset = 1 + ipart_dim1;
    ipart -= ipart_offset;
    spec_dim1 = *nk;
    spec_offset = 1 + spec_dim1;
    spec -= spec_offset;

    /* Function Body */
    npart = 0;
    nspec = *nk * *nth;
    if (nspec > 600) {
	s_wsle(&io___6);
	do_lio(&c__9, &c__1, "Error: Spectrum size exceeds maximum of ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&c__600, (ftnlen)sizeof(integer));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    ptnghb_(neigh, &nspec, nk, nth);
    i__1 = mth;
    for (iang = 1; iang <= i__1; ++iang) {
	i__2 = mk;
	for (ik = 1; ik <= i__2; ++ik) {
	    zp[(iang - 1) * mk + ik - 1] = spec[ik + iang * spec_dim1];
	}
    }
    zmin = minval_(zp, &nspec);
    zmax = maxval_(zp, &nspec);
    if (zmax - zmin < 1e-9f) {
	i__1 = mth;
	for (iang = 1; iang <= i__1; ++iang) {
	    i__2 = *nk;
	    for (ik = 1; ik <= i__2; ++ik) {
		ipart[ik + iang * ipart_dim1] = 0;
	    }
	}
	npart = 0;
	return 0;
    }
    fact = 199.f / (zmax - zmin);
    i__1 = nspec;
    for (ik = 1; ik <= i__1; ++ik) {
	zp[ik - 1] = zmax - zp[ik - 1];
/* Computing MAX */
/* Computing MIN */
	r__1 = zp[ik - 1] * fact + 1.f;
	i__4 = 200, i__5 = i_nint(&r__1);
	i__2 = 1, i__3 = min(i__4,i__5);
	imi[ik - 1] = max(i__2,i__3);
    }
    ptsort_(ind, imi, &nspec, &c__200);
    pt_fld__(neigh, zp, imo, imi, ind, &nspec, &npart, &c__200);
    i__1 = mth;
    for (iang = 1; iang <= i__1; ++iang) {
	i__2 = *nk;
	for (ik = 1; ik <= i__2; ++ik) {
	    ipart[ik + iang * ipart_dim1] = imo[(iang - 1) * mk + 1 + ik - 1];
	}
    }
    return 0;
} /* partition_ */

/* Subroutine */ int ptsort_(integer *ind, integer *imi, integer *nspec, 
	integer *ihmax)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, in, iv, numv[200], iaddr[200], iorder[600];


/* -------------------------------------------------------------------- / */
/* 1.  occurences per height */

    /* Parameter adjustments */
    --imi;
    --ind;

    /* Function Body */
    i__1 = *ihmax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	numv[i__ - 1] = 0;
    }
    i__1 = *nspec;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++numv[imi[i__] - 1];
    }

/* -------------------------------------------------------------------- / */
/* 2.  starting address per height */

    iaddr[0] = 1;
    i__1 = *ihmax - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iaddr[i__] = iaddr[i__ - 1] + numv[i__ - 1];
    }

/* -------------------------------------------------------------------- / */
/* 3.  order points */

    i__1 = *nspec;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iv = imi[i__];
	in = iaddr[iv - 1];
	iorder[i__ - 1] = in;
	iaddr[iv - 1] = in + 1;
    }

/* -------------------------------------------------------------------- / */
/* 4.  sort points */

    i__1 = *nspec;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ind[iorder[i__ - 1]] = i__;
    }

    return 0;
/* / */
/* / end of ptsort ----------------------------------------------------- / */
/* / */
} /* ptsort_ */

/* / ------------------------------------------------------------------- / */
/* Subroutine */ int ptnghb_(integer *neigh, integer *nspec, integer *mk, 
	integer *mth)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, n;

/* -------------------------------------------------------------------- / */
/* 2.  build map */

    /* Parameter adjustments */
    neigh -= 10;

    /* Function Body */
    i__1 = *nspec;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 9; ++j) {
	    neigh[j + i__ * 9] = 0;
	}
    }

/* ... base loop */

    i__1 = *nspec;
    for (n = 1; n <= i__1; ++n) {

	j = (n - 1) / *mk + 1;
	i__ = n - (j - 1) * *mk;
	k = 0;

/* ... point at the left(1) */

	if (i__ != 1) {
	    ++k;
	    neigh[k + n * 9] = n - 1;
	}

/* ... point at the right (2) */

	if (i__ != *mk) {
	    ++k;
	    neigh[k + n * 9] = n + 1;
	}

/* ... point at the bottom(3) */

	if (j != 1) {
	    ++k;
	    neigh[k + n * 9] = n - *mk;
	}

/* ... add point at bottom_wrap to top */

	if (j == 1) {
	    ++k;
	    neigh[k + n * 9] = *nspec - (*mk - i__);
	}

/* ... point at the top(4) */

	if (j != *mth) {
	    ++k;
	    neigh[k + n * 9] = n + *mk;
	}

/* ... add point to top_wrap to bottom */

	if (j == *mth) {
	    ++k;
	    neigh[k + n * 9] = n - (*mth - 1) * *mk;
	}

/* ... point at the bottom, left(5) */

	if (i__ != 1 && j != 1) {
	    ++k;
	    neigh[k + n * 9] = n - *mk - 1;
	}

/* ... point at the bottom, left with wrap. */

	if (i__ != 1 && j == 1) {
	    ++k;
	    neigh[k + n * 9] = n - 1 + *mk * (*mth - 1);
	}

/* ... point at the bottom, right(6) */

	if (i__ != *mk && j != 1) {
	    ++k;
	    neigh[k + n * 9] = n - *mk + 1;
	}

/* ... point at the bottom, right with wrap */

	if (i__ != *mk && j == 1) {
	    ++k;
	    neigh[k + n * 9] = n + 1 + *mk * (*mth - 1);
	}

/* ... point at the top, left(7) */

	if (i__ != 1 && j != *mth) {
	    ++k;
	    neigh[k + n * 9] = n + *mk - 1;
	}

/* ... point at the top, left with wrap */

	if (i__ != 1 && j == *mth) {
	    ++k;
	    neigh[k + n * 9] = n - 1 - *mk * (*mth - 1);
	}

/* ... point at the top, right(8) */

	if (i__ != *mk && j != *mth) {
	    ++k;
	    neigh[k + n * 9] = n + *mk + 1;
	}

/* ... point at top, right with wrap */


	if (i__ != *mk && j == *mth) {
	    ++k;
	    neigh[k + n * 9] = n + 1 - *mk * (*mth - 1);
	}

	neigh[n * 9 + 9] = k;

    }

    return 0;
/* / */
/* / end of ptnghb ----------------------------------------------------- / */
/* / */
} /* ptnghb_ */

/* / ------------------------------------------------------------------- / */
/* Subroutine */ int pt_fld__(integer *neigh, real *zp, integer *imo, integer 
	*imi, integer *ind, integer *nspec, integer *npart, integer *ihmax)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Local variables */
    static integer ic_label__;
    extern /* Subroutine */ int fifo_add__(integer *, integer *, integer *, 
	    integer *);
    static integer iq_start__, i__, j, m, ih, jl, jn, ip, iq[300];
    static real ep1;
    extern /* Subroutine */ int fifo_first__(integer *, integer *, integer *, 
	    integer *), fifo_empty__(integer *, integer *, integer *);
    static integer imd[300], ipp, ipt, ifict_pixel__;
    static real diff;
    static integer mask, init, ippp, ispec, msave;
    static real zpmax;
    static integer iq_end__, iwshed;
    extern doublereal maxval_(real *, integer *);
    static integer iempty, ic_dist__;
    extern integer iminval_(integer *, integer *);

/* / */
/* / ------------------------------------------------------------------- / */
/* / local parameters */
/* / */
/* -------------------------------------------------------------------- / */
/* 0.  initializations */

    /* Parameter adjustments */
    --ind;
    --imi;
    --imo;
    --zp;
    neigh -= 10;

    /* Function Body */
    mask = -2;
    init = -1;
    iwshed = 0;
    ic_label__ = 0;
    ifict_pixel__ = -100;

    iq_start__ = 1;
    iq_end__ = 1;

    zpmax = maxval_(&zp[1], nspec);
    i__1 = *nspec;
    for (i__ = 1; i__ <= i__1; ++i__) {
	imo[i__] = init;
	imd[i__ - 1] = 0;
    }

/* -------------------------------------------------------------------- / */
/* 1.  loop over levels */

    m = 1;

    i__1 = *ihmax;
    for (ih = 1; ih <= i__1; ++ih) {
	msave = m;

/* 1.a pixels at level ih */

L10:
	ip = ind[m];
	if (imi[ip] != ih) {
	    goto L11;
	}

/*     flag the point, if it stays flagge, it is a separate minimum. */

	imo[ip] = mask;

/*     consider neighbors. if there is neighbor, set distance and add */
/*     to queue. */

	i__2 = neigh[ip * 9 + 9];
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ipp = neigh[i__ + ip * 9];
	    if (imo[ipp] > 0 || imo[ipp] == iwshed) {
		imd[ip - 1] = 1;
		fifo_add__(&ip, iq, &iq_end__, nspec);
		goto L11;
	    }
	}

	if (m + 1 > *nspec) {
	    goto L11;
	} else {
	    ++m;
	}

	goto L10;

L11:
/* 1.b process the queue */

	ic_dist__ = 1;
	fifo_add__(&ifict_pixel__, iq, &iq_end__, nspec);

L20:
	fifo_first__(&ip, iq, &iq_start__, nspec);

/*     check for end of processing */

	if (ip == ifict_pixel__) {
	    fifo_empty__(&iempty, &iq_start__, &iq_end__);
	    if (iempty == 1) {
		goto L21;
	    } else {
		fifo_add__(&ifict_pixel__, iq, &iq_end__, nspec);
		++ic_dist__;
		fifo_first__(&ip, iq, &iq_start__, nspec);
	    }
	}

/*     process queue */

	i__2 = neigh[ip * 9 + 9];
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ipp = neigh[i__ + ip * 9];

/*     check for labeled watersheds or basins */

	    if (imd[ipp - 1] < ic_dist__ && (imo[ipp] > 0 || imo[ipp] == 
		    iwshed)) {

		if (imo[ipp] > 0) {

		    if (imo[ip] == mask || imo[ip] == iwshed) {
			imo[ip] = imo[ipp];
		    } else if (imo[ip] != imo[ipp]) {
			imo[ip] = iwshed;
		    }

		} else if (imo[ip] == mask) {

		    imo[ip] = iwshed;

		}

	    } else if (imo[ipp] == mask && imd[ipp - 1] == 0) {

		imd[ipp - 1] = ic_dist__ + 1;
		fifo_add__(&ipp, iq, &iq_end__, nspec);

	    }

	}

	goto L20;
L21:

/* 1.c check for mask values in imo to identify new basins */

	m = msave;

L30:
	ip = ind[m];
	if (imi[ip] != ih) {
	    goto L31;
	}
	imd[ip - 1] = 0;

	if (imo[ip] == mask) {

/* ... new label for pixel */

	    ++ic_label__;
	    fifo_add__(&ip, iq, &iq_end__, nspec);
	    imo[ip] = ic_label__;

/* ... and all connected to it ... */

L40:
	    fifo_empty__(&iempty, &iq_start__, &iq_end__);
	    if (iempty == 1) {
		goto L41;
	    }
	    fifo_first__(&ipp, iq, &iq_start__, nspec);

	    i__2 = neigh[ipp * 9 + 9];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ippp = neigh[i__ + ipp * 9];
		if (imo[ippp] == mask) {
		    fifo_add__(&ippp, iq, &iq_end__, nspec);
		    imo[ippp] = ic_label__;
		}
	    }

	    goto L40;
L41:

	    ;
	}

	if (m + 1 > *nspec) {
	    goto L31;
	} else {
	    ++m;
	}

	goto L30;
L31:

	;
    }

/* -------------------------------------------------------------------- / */
/* 2.  find nearest neighbor of 0 watershed points and replace */
/*     use original input to check which group to affiliate with 0 */
/*     soring changes first in imd to assure symetry in adjustment. */

    for (j = 1; j <= 5; ++j) {
	i__1 = *nspec;
	for (ispec = 1; ispec <= i__1; ++ispec) {
	    imd[ispec - 1] = imo[ispec];
	}
	i__1 = *nspec;
	for (jl = 1; jl <= i__1; ++jl) {
	    ipt = -1;
	    if (imo[jl] == 0) {
		ep1 = zpmax;
		i__2 = neigh[jl * 9 + 9];
		for (jn = 1; jn <= i__2; ++jn) {
		    diff = (r__1 = zp[jl] - zp[neigh[jn + jl * 9]], dabs(r__1)
			    );
		    if (diff <= ep1 && imo[neigh[jn + jl * 9]] != 0) {
			ep1 = diff;
			ipt = jn;
		    }
		}
		if (ipt > 0) {
		    imd[jl - 1] = imo[neigh[ipt + jl * 9]];
		}
	    }
	}
	i__1 = *nspec;
	for (ispec = 1; ispec <= i__1; ++ispec) {
	    imo[ispec] = imd[ispec - 1];
	}
	if (iminval_(&imo[1], nspec) > 0) {
	    goto L60;
	}
    }
L60:

    *npart = ic_label__;

    return 0;

} /* pt_fld__ */

/* / ------------------------------------------------------------------- / */
/* Subroutine */ int fifo_add__(integer *iv, integer *iq, integer *iq_end__, 
	integer *nspec)
{

/*     add point to fifo queue. */


    /* Parameter adjustments */
    --iq;

    /* Function Body */
    iq[*iq_end__] = *iv;

    ++(*iq_end__);
    if (*iq_end__ > *nspec) {
	*iq_end__ = 1;
    }

    return 0;
} /* fifo_add__ */

/* / ------------------------------------------------------------------- / */
/* Subroutine */ int fifo_empty__(integer *iempty, integer *iq_start__, 
	integer *iq_end__)
{

/*     check if queue is empty. */


    if (*iq_start__ != *iq_end__) {
	*iempty = 0;
    } else {
	*iempty = 1;
    }

    return 0;
} /* fifo_empty__ */

/* / ------------------------------------------------------------------- / */
/* Subroutine */ int fifo_first__(integer *iv, integer *iq, integer *
	iq_start__, integer *nspec)
{

/*     get point out of queue. */


    /* Parameter adjustments */
    --iq;

    /* Function Body */
    *iv = iq[*iq_start__];

    ++(*iq_start__);
    if (*iq_start__ > *nspec) {
	*iq_start__ = 1;
    }

    return 0;
} /* fifo_first__ */

