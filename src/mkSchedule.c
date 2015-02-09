/* This program is a combination between gendaymtx and epw2wea.
 * Both of those programs have copyright, see below.
 * Also, gendaymtx borrowed a lot of code from Ian Ashdown's code... etc. etc.
 * Read below for more information and copyrights.
 * 
 * An improvement over the predecessors of this program is that this one can
 * interpolate the .epw file in order to create sub-hourly simulations.
 *
 * This program also allow to consider only the first "n" lines in the .epw file.
 *
 * Also, with this program it is possible to choose a different location (latitude and
 * longitude).
 *
 * It might be important to note that this code does not close the year; that is,
 * interpolate BETWEEN the first year and the final one, and looses the hour that 
 * connects them (the hour 8760)... I thought it was not a big deal... it can be 
 * modified I guess.
 *
 * Written by Germán Molina, Pontificia Universidad Católica de Chile.
 * gmolina1@uc.cl
 *

 *********************************************************************

 *  Copyright (c) 2003
 *  National Research Council Canada
 *  written by Christoph Reinhart
 *
 *
 * epw2wea: daylight analysis subprogram of DAYSIM
 * Program converts EnergyPlus weather format (*.ppw) into DAYSIM format (*.wea)

 *********************************************************************

 * #ifndef lint
 * static const char RCSid[] = "$Id: gendaymtx.c,v 2.11 2013/04/30 17:05:27 greg Exp $";
 * #endif
 *
 *  epw2daymtx.c
 *  
 *  Generate a daylight matrix based on Perez Sky Model.
 *
 *  Most of this code is borrowed (see copyright below) from Ian Ashdown's
 *  excellent re-implementation of Jean-Jacques Delaunay's gendaylit.c
 *
 *  Created by Greg Ward on 1/16/13.
 *

 *********************************************************************
 *
 *  H32_gendaylit.CPP - Perez Sky Model Calculation
 *
 *  Version:    1.00A
 *
 *  History:    09/10/01 - Created.
 *				11/10/08 - Modified for Unix compilation.
 *				11/10/12 - Fixed conditional __max directive.
 *				1/11/13 - Tweaks and optimizations (G.Ward)
 *
 *  Compilers:  Microsoft Visual C/C++ Professional V10.0    
 *
 *  Author:     Ian Ashdown, P.Eng.
 *              byHeart Consultants Limited
 *              620 Ballantree Road
 *              West Vancouver, B.C.
 *              Canada V7S 1W3
 *              e-mail: ian_ashdown@helios32.com
 *
 *  References:	Perez, R., P. Ineichen, R. Seals, J. Michalsky, and R.
 *				Stewart. 1990. ìModeling Daylight Availability and
 *				Irradiance Components from Direct and Global
 *				Irradiance,î Solar Energy 44(5):271-289.
 *
 *				Perez, R., R. Seals, and J. Michalsky. 1993.
 *				ìAll-Weather Model for Sky Luminance Distribution -
 *				Preliminary Configuration and Validation,î Solar Energy
 *				50(3):235-245.
 *
 *				Perez, R., R. Seals, and J. Michalsky. 1993. "ERRATUM to
 *				All-Weather Model for Sky Luminance Distribution -
 *				Preliminary Configuration and Validation,î Solar Energy
 *				51(5):423.
 *
 *  NOTE:		This program is a completely rewritten version of
 *				gendaylit.c written by Jean-Jacques Delaunay (1994).
 *
 *  Copyright 2009-2012 byHeart Consultants Limited. All rights
 *  reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted for personal and commercial purposes
 *  provided that redistribution of source code must retain the above
 *  copyright notice, this list of conditions and the following
 *  disclaimer:
 *
 *    THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED
 *    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 *    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *    DISCLAIMED. IN NO EVENT SHALL byHeart Consultants Limited OR
 *    ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *    USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *    ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *    POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************/


/* Include files ... from dctimestep and gendaymtx */
#define	_USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "./ray/src/common/standard.h"
#include "./ray/src/common/platform.h"
#include "./ray/src/common/paths.h"
#include "./ray/src/common/rtmath.h"
#include "./ray/src/common/color.h"
#include "./ray/src/common/resolu.h"
#include "./ray/src/common/bsdf.h"
#include "./ray/src/common/bsdf_m.h"

/* ============ Lua Requirements ============= */
/* =========================================== */

#include "/usr/include/lua5.2/lua.h" 
#include "/usr/include/lua5.2/luaconf.h"
#include "/usr/include/lua5.2/lualib.h" 
#include "/usr/include/lua5.2/lauxlib.h"


/* Sky information */
float	*mtx_data = NULL;	/* our matrix data */

/* the Lua interpreter */
lua_State* L;

/* ============ state variables  ============= */
/* =========================================== */

#define MAXsensor 5
#define MAXlum 5
#define MAXwin 5
#define MAXpos 10

float lum_state[MAXlum];
float sensor[MAXsensor][1];
int win_state[MAXwin];
int n_dirs[MAXwin];
int is_day;


/* ============= input defaults  ============= */
/* =========================================== */

int t_bool=0;
int l_bool=0;
int s_bool=0;
int dc_bool=0;
int lua_bool=0;
int w_bool=0;
int sen_bool=0;
int sch_bool=0;


int nwin=0;
int nsen=0;
int nlum=0;


char var_name[1024]; /*this is used to pull the variables from the Lua script */	

/* A color coefficient matrix -- vectors have ncols==1 */
typedef struct {
	int	nrows, ncols;
	COLORV	cmem[3];		/* extends struct */
} CMATRIX;

CMATRIX *Lum_array[MAXlum];
CMATRIX *DC3[MAXwin][MAXpos];

/* ============================================ */
/*          gendaymtx functions                 */
/* ============================================ */

char *progname;								/* Program name */
char errmsg[128];							/* Error message buffer */
const double DC_SolarConstantE = 1367.0;	/* Solar constant W/m^2 */
const double DC_SolarConstantL = 127.5;		/* Solar constant klux */

double altitude;			/* Solar altitude (radians) */
double azimuth;				/* Solar azimuth (radians) */
double apwc;				/* Atmospheric precipitable water content */
double dew_point = 11.0;		/* Surface dew point temperature (deg. C) */
double diff_illum;			/* Diffuse illuminance */
double diff_irrad;			/* Diffuse irradiance */
double dir_illum;			/* Direct illuminance */
double dir_irrad;			/* Direct irradiance */
int julian_date;			/* Julian date */
double perez_param[5];			/* Perez sky model parameters */
double sky_brightness;			/* Sky brightness */
double sky_clearness;			/* Sky clearness */
double solar_rad;			/* Solar radiance */
double sun_zenith;			/* Sun zenith angle (radians) */
int	input = 0;				/* Input type */
int	output = 0;				/* Output type */

/* other */
double ext_drybulb;
double	dir, dif;		/* direct and diffuse values */
double	hr;			/* hour (local standard time) */


extern double dmax( double, double );
extern double CalcAirMass();
extern double CalcDiffuseIllumRatio( int );
extern double CalcDiffuseIrradiance();
extern double CalcDirectIllumRatio( int );
extern double CalcDirectIrradiance();
extern double CalcEccentricity();
extern double CalcPrecipWater( double );
extern double CalcRelHorzIllum( float *parr );
extern double CalcRelLuminance( double, double );
extern double CalcSkyBrightness();
extern double CalcSkyClearness();
extern int CalcSkyParamFromIllum();
extern int GetCategoryIndex();
extern void CalcPerezParam( double, double, double, int );
extern void CalcSkyPatchLumin( float *parr );
extern void ComputeSky( float *parr );

/* Degrees into radians */
#define DegToRad(deg)	((deg)*(PI/180.))

/* Radiuans into degrees */
#define RadToDeg(rad)	((rad)*(180./PI))


/* Perez sky model coefficients */

/* Reference:	Perez, R., R. Seals, and J. Michalsky, 1993. "All- */
/*				Weather Model for Sky Luminance Distribution - */
/*				Preliminary Configuration and Validation," Solar */
/*				Energy 50(3):235-245, Table 1. */

static const double PerezCoeff[8][20] =
{
	/* Sky clearness (epsilon): 1.000 to 1.065 */
	{   1.3525,  -0.2576,  -0.2690,  -1.4366,   -0.7670,
	    0.0007,   1.2734,  -0.1233,   2.8000,    0.6004,
	    1.2375,   1.0000,   1.8734,   0.6297,    0.9738,
	    0.2809,   0.0356,  -0.1246,  -0.5718,    0.9938 },
    /* Sky clearness (epsilon): 1.065 to 1.230 */
	{  -1.2219,  -0.7730,   1.4148,   1.1016,   -0.2054,
	    0.0367,  -3.9128,   0.9156,   6.9750,    0.1774,
		6.4477,  -0.1239,  -1.5798,  -0.5081,   -1.7812,
		0.1080,   0.2624,   0.0672,  -0.2190,   -0.4285 },
    /* Sky clearness (epsilon): 1.230 to 1.500 */
	{  -1.1000,  -0.2515,   0.8952,   0.0156,    0.2782,
	   -0.1812, - 4.5000,   1.1766,  24.7219,  -13.0812,
	  -37.7000,  34.8438,  -5.0000,   1.5218,    3.9229,
	   -2.6204,  -0.0156,   0.1597,   0.4199,   -0.5562 },
    /* Sky clearness (epsilon): 1.500 to 1.950 */
	{  -0.5484,  -0.6654,  -0.2672,   0.7117,   0.7234,
	   -0.6219,  -5.6812,   2.6297,  33.3389, -18.3000,
	  -62.2500,  52.0781,  -3.5000,   0.0016,   1.1477,
	    0.1062,   0.4659,  -0.3296,  -0.0876,  -0.0329 },
    /* Sky clearness (epsilon): 1.950 to 2.800 */
	{  -0.6000,  -0.3566,  -2.5000,   2.3250,   0.2937,
	    0.0496,  -5.6812,   1.8415,  21.0000,  -4.7656 ,
	  -21.5906,   7.2492,  -3.5000,  -0.1554,   1.4062,
	    0.3988,   0.0032,   0.0766,  -0.0656,  -0.1294 },
    /* Sky clearness (epsilon): 2.800 to 4.500 */
	{  -1.0156,  -0.3670,   1.0078,   1.4051,   0.2875,
	   -0.5328,  -3.8500,   3.3750,  14.0000,  -0.9999,
	   -7.1406,   7.5469,  -3.4000,  -0.1078,  -1.0750,
	    1.5702,  -0.0672,   0.4016,   0.3017,  -0.4844 },
    /* Sky clearness (epsilon): 4.500 to 6.200 */
	{  -1.0000,   0.0211,   0.5025,  -0.5119,  -0.3000,
	    0.1922,   0.7023,  -1.6317,  19.0000,  -5.0000,
		1.2438,  -1.9094,  -4.0000,   0.0250,   0.3844,
		0.2656,   1.0468,  -0.3788,  -2.4517,   1.4656 },
    /* Sky clearness (epsilon): 6.200 to ... */
	{  -1.0500,   0.0289,   0.4260,   0.3590,  -0.3250,
	    0.1156,   0.7781,   0.0025,  31.0625, -14.5000,
	  -46.1148,  55.3750,  -7.2312,   0.4050,  13.3500,
	    0.6234,   1.5000,  -0.6426,   1.8564,   0.5636 }
};

/* Perez irradiance component model coefficients */

/* Reference:	Perez, R., P. Ineichen, R. Seals, J. Michalsky, and R. */
/*				Stewart. 1990. ìModeling Daylight Availability and */
/*				Irradiance Components from Direct and Global */
/*				Irradiance,î Solar Energy 44(5):271-289. */

typedef struct
{
	double lower;	/* Lower bound */
	double upper;	/* Upper bound */
} CategoryBounds;

/* Perez sky clearness (epsilon) categories (Table 1) */
static const CategoryBounds SkyClearCat[8] =
{
	{ 1.000, 1.065 },	/* Overcast */
	{ 1.065, 1.230 },
	{ 1.230, 1.500 },
	{ 1.500, 1.950 },
	{ 1.950, 2.800 },
	{ 2.800, 4.500 },
	{ 4.500, 6.200 },
	{ 6.200, 12.00 }	/* Clear */
};

/* Luminous efficacy model coefficients */
typedef struct
{
	double a;
	double b;
	double c;
	double d;
} ModelCoeff;

/* Diffuse luminous efficacy model coefficients (Table 4, Eqn. 7) */
static const ModelCoeff DiffuseLumEff[8] =
{
	{  97.24, -0.46,  12.00,  -8.91 },
	{ 107.22,  1.15,   0.59,  -3.95 },
	{ 104.97,  2.96,  -5.53,  -8.77 },
	{ 102.39,  5.59, -13.95, -13.90 },
	{ 100.71,  5.94, -22.75, -23.74 },
	{ 106.42,  3.83, -36.15, -28.83 },
	{ 141.88,  1.90, -53.24, -14.03 },
	{ 152.23,  0.35, -45.27,  -7.98 }
};

/* Direct luminous efficacy model coefficients (Table 4, Eqn. 8) */
static const ModelCoeff DirectLumEff[8] =
{
	{  57.20, -4.55, -2.98, 117.12 },
	{  98.99, -3.46, -1.21,  12.38 },
	{ 109.83, -4.90, -1.71,  -8.81 },
	{ 110.34, -5.84, -1.99,  -4.56 },
	{ 106.36, -3.97, -1.75,  -6.16 },
	{ 107.19, -1.25, -1.51, -26.73 },
	{ 105.75,  0.77, -1.26, -34.44 },
	{ 101.18,  1.58, -1.10,  -8.29 }
};

#ifndef NSUNPATCH
#define	NSUNPATCH	4		/* max. # patches to spread sun into */
#endif

extern int jdate(int month, int day);
extern double stadj(int  jd);
extern double sdec(int  jd);
extern double salt(double sd, double st);
extern double sazi(double sd, double st);
					/* sun calculation constants */
extern double  s_latitude;
extern double  s_longitude;
extern double  s_meridian;

int		nsuns = NSUNPATCH;	/* number of sun patches to use */
double		fixed_sun_sa = -1;	/* fixed solid angle per sun? */

int		verbose = 0;		/* progress reports to stderr? */

int		outfmt = 'a';		/* output format */

int		rhsubdiv = 1;		/* Reinhart sky subdivisions */

COLOR		skycolor = {.96, 1.004, 1.118};	/* sky coloration */
COLOR		suncolor = {1., 1., 1.};	/* sun color */
COLOR		grefl = {.2, .2, .2};		/* ground reflectance */

int		nskypatch;		/* number of Reinhart patches */
float		*rh_palt;		/* sky patch altitudes (radians) */
float		*rh_pazi;		/* sky patch azimuths (radians) */
float		*rh_dom;		/* sky patch solid angle (sr) */

#define		vector(v,alt,azi)	(	(v)[1] = tcos(alt), \
						(v)[0] = (v)[1]*tsin(azi), \
						(v)[1] *= tcos(azi), \
						(v)[2] = tsin(alt) )

#define		rh_vector(v,i)		vector(v,rh_palt[i],rh_pazi[i])

#define		rh_cos(i)		tsin(rh_palt[i])

extern int	rh_init(void);
extern float *	resize_dmatrix(float *mtx_data, int nsteps, int npatch);
extern void	AddDirect(float *parr);


/* ============================================== */
/*              dctimestep functions              */
/* ============================================== */

	
/* =========================================== */
/* =========================================== */

char	*progname;			/* global argv[0] */

/* Data types for file loading */
enum {DTfromHeader, DTascii, DTfloat, DTdouble, DTrgbe, DTxyze};



#define COLSPEC	(sizeof(COLORV)==sizeof(float) ? "%f %f %f" : "%lf %lf %lf")

#define cm_lval(cm,r,c)	((cm)->cmem + 3*((r)*(cm)->ncols + (c)))

#define cv_lval(cm,i)	((cm)->cmem + 3*(i))

/* Allocate a color coefficient matrix */
static CMATRIX *
cm_alloc(int nrows, int ncols)
{
	CMATRIX	*cm;

	if ((nrows <= 0) | (ncols <= 0))
		error(USER, "attempt to create empty matrix");
	cm = (CMATRIX *)malloc(sizeof(CMATRIX) +
				3*sizeof(COLORV)*(nrows*ncols - 1));
	if (cm == NULL)
		error(SYSTEM, "out of memory in cm_alloc()");
	cm->nrows = nrows;
	cm->ncols = ncols;
	return(cm);
}

#define cm_free(cm)	free(cm)

/* Resize color coefficient matrix */
static CMATRIX *
cm_resize(CMATRIX *cm, int nrows)
{
	if (nrows == cm->nrows)
		return(cm);
	if (nrows <= 0) {
		cm_free(cm);
		return(NULL);
	}
	cm = (CMATRIX *)realloc(cm, sizeof(CMATRIX) +
			3*sizeof(COLORV)*(nrows*cm->ncols - 1));
	if (cm == NULL)
		error(SYSTEM, "out of memory in cm_resize()");
	cm->nrows = nrows;
	return(cm);
}

/* Load header to obtain data type */
static int
getDT(char *s, void *p)
{
	char	fmt[32];
	
	if (formatval(fmt, s)) {
		if (!strcmp(fmt, "ascii"))
			*((int *)p) = DTascii;
		else if (!strcmp(fmt, "float"))
			*((int *)p) = DTfloat;
		else if (!strcmp(fmt, "double"))
			*((int *)p) = DTdouble;
		else if (!strcmp(fmt, COLRFMT))
			*((int *)p) = DTrgbe;
		else if (!strcmp(fmt, CIEFMT))
			*((int *)p) = DTxyze;
	}
	return(0);
}

static int
getDTfromHeader(FILE *fp)
{
	int	dt = DTfromHeader;
	
	if (getheader(fp, getDT, &dt) < 0)
		error(SYSTEM, "header read error");
	if (dt == DTfromHeader)
		error(USER, "missing data format in header");
	return(dt);
}

/* Allocate and load a matrix from the given file (or stdin if NULL) */
static CMATRIX *
cm_load(const char *fname, int nrows, int ncols, int dtype)
{
	FILE	*fp = stdin;
	CMATRIX	*cm;

	if (ncols <= 0)
		error(USER, "Non-positive number of columns");
	if (fname == NULL)
		fname = "<stdin>";
	else if ((fp = fopen(fname, "r")) == NULL) {
		sprintf(errmsg, "cannot open file '%s'", fname);
		error(SYSTEM, errmsg);
	}
#ifdef getc_unlocked
	flockfile(fp);
#endif
	if (dtype != DTascii)
		SET_FILE_BINARY(fp);		/* doesn't really work */
	if (dtype == DTfromHeader)
		dtype = getDTfromHeader(fp);
	switch (dtype) {
	case DTascii:
	case DTfloat:
	case DTdouble:
		break;
	default:
		error(USER, "unexpected data type in cm_load()");
	}
	if (nrows <= 0) {			/* don't know length? */
		int	guessrows = 147;	/* usually big enough */
		if ((dtype != DTascii) & (fp != stdin)) {
			long	startpos = ftell(fp);
			if (fseek(fp, 0L, SEEK_END) == 0) {
				long	endpos = ftell(fp);
				long	elemsiz = 3*(dtype==DTfloat ?
					    sizeof(float) : sizeof(double));

				if ((endpos - startpos) % (ncols*elemsiz)) {
					sprintf(errmsg,
					"improper length for binary file '%s'",
							fname);
					error(USER, errmsg);
				}
				guessrows = (endpos - startpos)/(ncols*elemsiz);
				if (fseek(fp, startpos, SEEK_SET) < 0) {
					sprintf(errmsg,
						"fseek() error on file '%s'",
							fname);
					error(SYSTEM, errmsg);
				}
				nrows = guessrows;	/* we're confident */
			}
		}
		cm = cm_alloc(guessrows, ncols);
	} else
		cm = cm_alloc(nrows, ncols);
	if (cm == NULL)					/* XXX never happens */
		return(NULL);
	if (dtype == DTascii) {				/* read text file */
		int	maxrow = (nrows > 0 ? nrows : 32000);
		int	r, c;
		for (r = 0; r < maxrow; r++) {
		    if (r >= cm->nrows)			/* need more space? */
			cm = cm_resize(cm, 2*cm->nrows);
		    for (c = 0; c < ncols; c++) {
		        COLORV	*cv = cm_lval(cm,r,c);
			if (fscanf(fp, COLSPEC, cv, cv+1, cv+2) != 3)
				if ((nrows <= 0) & (r > 0) & !c) {
					cm = cm_resize(cm, maxrow=r);
					break;
				} else
					goto EOFerror;
		    		
		    }
		}
		while ((c = getc(fp)) != EOF)
			if (!isspace(c)) {
				sprintf(errmsg,
				"unexpected data at end of ascii file %s",
						fname);
				error(WARNING, errmsg);
				break;
			}
	} else {					/* read binary file */
		if (sizeof(COLORV) == (dtype==DTfloat ? sizeof(float) :
							sizeof(double))) {
			int	nread = 0;
			do {				/* read all we can */
				nread += fread(cm->cmem + 3*nread,
						3*sizeof(COLORV),
						cm->nrows*cm->ncols - nread,
						fp);
				if (nrows <= 0) {	/* unknown length */
					if (nread == cm->nrows*cm->ncols)
							/* need more space? */
						cm = cm_resize(cm, 2*cm->nrows);
					else if (nread && !(nread % cm->ncols))
							/* seem to be  done */
						cm = cm_resize(cm, nread/cm->ncols);
					else		/* ended mid-row */
						goto EOFerror;
				} else if (nread < cm->nrows*cm->ncols)
					goto EOFerror;
			} while (nread < cm->nrows*cm->ncols);

		} else if (dtype == DTdouble) {
			double	dc[3];			/* load from double */
			COLORV	*cvp = cm->cmem;
			int	n = nrows*ncols;

			if (n <= 0)
				goto not_handled;
			while (n--) {
				if (fread(dc, sizeof(double), 3, fp) != 3)
					goto EOFerror;
				copycolor(cvp, dc);
				cvp += 3;
			}
		} else /* dtype == DTfloat */ {
			float	fc[3];			/* load from float */
			COLORV	*cvp = cm->cmem;
			int	n = nrows*ncols;

			if (n <= 0)
				goto not_handled;
			while (n--) {
				if (fread(fc, sizeof(float), 3, fp) != 3)
					goto EOFerror;
				copycolor(cvp, fc);
				cvp += 3;
			}
		}
		if (fgetc(fp) != EOF) {
				sprintf(errmsg,
				"unexpected data at end of binary file %s",
						fname);
				error(WARNING, errmsg);
		}
	}
	if (fp != stdin)
		fclose(fp);
#ifdef getc_unlocked
	else
		funlockfile(fp);
#endif
	return(cm);
EOFerror:
	sprintf(errmsg, "unexpected EOF reading %s", fname);
	error(USER, errmsg);
not_handled:
	error(INTERNAL, "unhandled data size or length in cm_load()");
	return(NULL);	/* gratis return */
}

/* Extract a column vector from a matrix */
static CMATRIX *
cm_column(const CMATRIX *cm, int c)
{
	CMATRIX	*cvr;
	int	dr;

	if ((c < 0) | (c >= cm->ncols))
		error(INTERNAL, "column requested outside matrix");
	cvr = cm_alloc(cm->nrows, 1);
	if (cvr == NULL)
		return(NULL);
	for (dr = 0; dr < cm->nrows; dr++) {
		const COLORV	*sp = cm_lval(cm,dr,c);
		COLORV		*dp = cv_lval(cvr,dr);
		dp[0] = sp[0];
		dp[1] = sp[1];
		dp[2] = sp[2];
	}
	return(cvr);
}

/* Scale a matrix by a single value */
static CMATRIX *
cm_scale(const CMATRIX *cm1, const COLOR sca)
{
	CMATRIX	*cmr;
	int	dr, dc;

	cmr = cm_alloc(cm1->nrows, cm1->ncols);
	if (cmr == NULL)
		return(NULL);
	for (dr = 0; dr < cmr->nrows; dr++)
	    for (dc = 0; dc < cmr->ncols; dc++) {
	        const COLORV	*sp = cm_lval(cm1,dr,dc);
		COLORV		*dp = cm_lval(cmr,dr,dc);
		dp[0] = sp[0] * sca[0];
		dp[1] = sp[1] * sca[1];
		dp[2] = sp[2] * sca[2];
	    }
	return(cmr);
}

/* Multiply two matrices (or a matrix and a vector) and allocate the result */
static CMATRIX *
cm_multiply(const CMATRIX *cm1, const CMATRIX *cm2)
{
	char	*rowcheck=NULL, *colcheck=NULL;
	CMATRIX	*cmr;
	int	dr, dc, i;

	if ((cm1->ncols <= 0) | (cm1->ncols != cm2->nrows))
		error(INTERNAL, "matrix dimension mismatch in cm_multiply()");
	cmr = cm_alloc(cm1->nrows, cm2->ncols);
	if (cmr == NULL)
		return(NULL);
				/* optimization: check for zero rows & cols */
	if (((cm1->nrows > 5) | (cm2->ncols > 5)) & (cm1->ncols > 5)) {
		static const COLOR	czero;
		rowcheck = (char *)calloc(cmr->nrows, 1);
		for (dr = cm1->nrows*(rowcheck != NULL); dr--; )
		    for (dc = cm1->ncols; dc--; )
			if (memcmp(cm_lval(cm1,dr,dc), czero, sizeof(COLOR))) {
				rowcheck[dr] = 1;
				break;
			}
		colcheck = (char *)calloc(cmr->ncols, 1);
		for (dc = cm2->ncols*(colcheck != NULL); dc--; )
		    for (dr = cm2->nrows; dr--; )
			if (memcmp(cm_lval(cm2,dr,dc), czero, sizeof(COLOR))) {
				colcheck[dc] = 1;
				break;
			}
	}
	for (dr = 0; dr < cmr->nrows; dr++)
	    for (dc = 0; dc < cmr->ncols; dc++) {
		COLORV	*dp = cm_lval(cmr,dr,dc);
		dp[0] = dp[1] = dp[2] = 0;
		if (rowcheck != NULL && !rowcheck[dr])
			continue;
		if (colcheck != NULL && !colcheck[dc])
			continue;
		for (i = 0; i < cm1->ncols; i++) {
		    const COLORV	*cp1 = cm_lval(cm1,dr,i);
		    const COLORV	*cp2 = cm_lval(cm2,i,dc);
		    dp[0] += cp1[0] * cp2[0];
		    dp[1] += cp1[1] * cp2[1];
		    dp[2] += cp1[2] * cp2[2];
		}
	    }
	if (rowcheck != NULL) free(rowcheck);
	if (colcheck != NULL) free(colcheck);
	return(cmr);
}

/* print out matrix as ASCII text -- no header */
static void
cm_print(const CMATRIX *cm, FILE *fp)
{
	int		r, c;
	const COLORV	*mp = cm->cmem;
	
	for (r = 0; r < cm->nrows; r++) {
		for (c = 0; c < cm->ncols; c++, mp += 3)
			fprintf(fp, "\t%.6e %.6e %.6e", mp[0], mp[1], mp[2]);
		fputc('\n', fp);
	}
}

/* Convert a BSDF to our matrix representation */
static CMATRIX *
cm_bsdf(const COLOR bsdfLamb, const COLOR specCol, const SDMat *bsdf)
{
	CMATRIX	*cm = cm_alloc(bsdf->nout, bsdf->ninc);
	int	nbadohm = 0;
	int	nneg = 0;
	int	r, c;
					/* loop over incident angles */
	for (c = 0; c < cm->ncols; c++) {
		const double	dom = mBSDF_incohm(bsdf,c);
					/* projected solid angle */
		nbadohm += (dom <= 0);

		for (r = 0; r < cm->nrows; r++) {
			float	f = mBSDF_value(bsdf,c,r);
			COLORV	*mp = cm_lval(cm,r,c);
					/* check BSDF value */
			if ((f <= 0) | (dom <= 0)) {
				nneg += (f < -FTINY);
				f = .0f;
			}
			copycolor(mp, specCol);
			scalecolor(mp, f);
			addcolor(mp, bsdfLamb);
			scalecolor(mp, dom);
		}
	}
	if (nneg | nbadohm) {
		sprintf(errmsg,
		    "BTDF has %d negatives and %d bad incoming solid angles",
				nneg, nbadohm);
		error(WARNING, errmsg);
	}
	return(cm);
}

/* Convert between input and output indices for reciprocity */
static int
recip_out_from_in(const SDMat *bsdf, int in_recip)
{
	FVECT	v;

	if (!mBSDF_incvec(v, bsdf, in_recip+.5))
		return(in_recip);		/* XXX should be error! */
	v[2] = -v[2];
	return(mBSDF_outndx(bsdf, v));
}

/* Convert between output and input indices for reciprocity */
static int
recip_in_from_out(const SDMat *bsdf, int out_recip)
{
	FVECT	v;

	if (!mBSDF_outvec(v, bsdf, out_recip+.5))
		return(out_recip);		/* XXX should be error! */
	v[2] = -v[2];
	return(mBSDF_incndx(bsdf, v));
}

/* Convert a BSDF to our matrix representation, applying reciprocity */
static CMATRIX *
cm_bsdf_recip(const COLOR bsdfLamb, const COLOR specCol, const SDMat *bsdf)
{
	CMATRIX	*cm = cm_alloc(bsdf->ninc, bsdf->nout);
	int	nbadohm = 0;
	int	nneg = 0;
	int	r, c;
					/* loop over incident angles */
	for (c = 0; c < cm->ncols; c++) {
		const int	ro = recip_out_from_in(bsdf,c);
		const double	dom = mBSDF_outohm(bsdf,ro);
					/* projected solid angle */
		nbadohm += (dom <= 0);

		for (r = 0; r < cm->nrows; r++) {
			const int	ri = recip_in_from_out(bsdf,r);
			float		f = mBSDF_value(bsdf,ri,ro);
			COLORV		*mp = cm_lval(cm,r,c);
					/* check BSDF value */
			if ((f <= 0) | (dom <= 0)) {
				nneg += (f < -FTINY);
				f = .0f;
			}
			copycolor(mp, specCol);
			scalecolor(mp, f);
			addcolor(mp, bsdfLamb);
			scalecolor(mp, dom);
		}
	}
	if (nneg | nbadohm) {
		sprintf(errmsg,
		    "BTDF has %d negatives and %d bad incoming solid angles",
				nneg, nbadohm);
		error(WARNING, errmsg);
	}
	return(cm);
}

/* Load and convert a matrix BSDF from the given XML file */
static CMATRIX *
cm_loadBSDF(char *fname, COLOR cLamb)
{
	CMATRIX		*Tmat;
	char		*fpath;
	int		recip;
	SDError		ec;
	SDData		myBSDF;
	SDSpectralDF	*tdf;
	COLOR		bsdfLamb, specCol;
					/* find path to BSDF file */
	fpath = getpath(fname, getrlibpath(), R_OK);
	if (fpath == NULL) {
		sprintf(errmsg, "cannot find BSDF file '%s'", fname);
		error(USER, errmsg);
	}
	SDclearBSDF(&myBSDF, fname);	/* load XML and check type */
	ec = SDloadFile(&myBSDF, fpath);
	if (ec)
		error(USER, transSDError(ec));
	ccy2rgb(&myBSDF.tLamb.spec, myBSDF.tLamb.cieY/PI, bsdfLamb);
	recip = (myBSDF.tb == NULL);
	tdf = recip ? myBSDF.tf : myBSDF.tb;
	if (tdf == NULL) {		/* no non-Lambertian transmission? */
		if (cLamb != NULL)
			copycolor(cLamb, bsdfLamb);
		SDfreeBSDF(&myBSDF);
		return(NULL);
	}
	if (tdf->ncomp != 1 || tdf->comp[0].func != &SDhandleMtx) {
		sprintf(errmsg, "unsupported BSDF '%s'", fpath);
		error(USER, errmsg);
	}
					/* convert BTDF to matrix */
	ccy2rgb(&tdf->comp[0].cspec[0], 1., specCol);
	Tmat = recip ? cm_bsdf_recip(bsdfLamb, specCol, (SDMat *)tdf->comp[0].dist)
			: cm_bsdf(bsdfLamb, specCol, (SDMat *)tdf->comp[0].dist);
	if (cLamb != NULL)		/* Lambertian is included */
		setcolor(cLamb, .0, .0, .0);
					/* free BSDF and return */
	SDfreeBSDF(&myBSDF);
	return(Tmat);
}

/* Sum together a set of images and write result to fout */
static int
sum_images(const char *fspec, const CMATRIX *cv, FILE *fout)
{
	int	myDT = DTfromHeader;
	COLOR	*scanline = NULL;
	CMATRIX	*pmat = NULL;
	int	myXR=0, myYR=0;
	int	i, y;

	if (cv->ncols != 1)
		error(INTERNAL, "expected vector in sum_images()");
	for (i = 0; i < cv->nrows; i++) {
		const COLORV	*scv = cv_lval(cv,i);
		char		fname[1024];
		FILE		*fp;
		int		dt, xr, yr;
		COLORV		*psp;
							/* check for zero */
		if ((scv[RED] == 0) & (scv[GRN] == 0) & (scv[BLU] == 0) &&
				(myDT != DTfromHeader) | (i < cv->nrows-1))
			continue;
							/* open next picture */
		sprintf(fname, fspec, i);
		if ((fp = fopen(fname, "r")) == NULL) {
			sprintf(errmsg, "cannot open picture '%s'", fname);
			error(SYSTEM, errmsg);
		}
		SET_FILE_BINARY(fp);
		dt = getDTfromHeader(fp);
		if ((dt != DTrgbe) & (dt != DTxyze) ||
				!fscnresolu(&xr, &yr, fp)) {
			sprintf(errmsg, "file '%s' not a picture", fname);
			error(USER, errmsg);
		}
		if (myDT == DTfromHeader) {		/* on first one */
			myDT = dt;
			myXR = xr; myYR = yr;
			scanline = (COLOR *)malloc(sizeof(COLOR)*myXR);
			if (scanline == NULL)
				error(SYSTEM, "out of memory in sum_images()");
			pmat = cm_alloc(myYR, myXR);
			memset(pmat->cmem, 0, sizeof(COLOR)*myXR*myYR);
							/* finish header */
			fputformat(myDT==DTrgbe ? COLRFMT : CIEFMT, fout);
			fputc('\n', fout);
			fprtresolu(myXR, myYR, fout);
			fflush(fout);
		} else if ((dt != myDT) | (xr != myXR) | (yr != myYR)) {
			sprintf(errmsg, "picture '%s' format/size mismatch",
					fname);
			error(USER, errmsg);
		}
		psp = pmat->cmem;
		for (y = 0; y < yr; y++) {		/* read it in */
			int	x;
			if (freadscan(scanline, xr, fp) < 0) {
				sprintf(errmsg, "error reading picture '%s'",
						fname);
				error(SYSTEM, errmsg);
			}
							/* sum in scanline */
			for (x = 0; x < xr; x++, psp += 3) {
				multcolor(scanline[x], scv);
				addcolor(psp, scanline[x]);
			}
		}
		fclose(fp);				/* done this picture */
	}
	free(scanline);
							/* write scanlines */
	for (y = 0; y < myYR; y++)
		if (fwritescan((COLOR *)cm_lval(pmat, y, 0), myXR, fout) < 0)
			return(0);
	cm_free(pmat);					/* all done */
	return(fflush(fout) == 0);
}

/* check to see if a string contains a %d or %o specification */
static int
hasNumberFormat(const char *s)
{
	if (s == NULL)
		return(0);

	while (*s) {
		while (*s != '%')
			if (!*s++)
				return(0);
		if (*++s == '%') {		/* ignore "%%" */
			++s;
			continue;
		}
		while (isdigit(*s))		/* field length */
			++s;
						/* field we'll use? */
		if ((*s == 'd') | (*s == 'i') | (*s == 'o') |
					(*s == 'x') | (*s == 'X'))
			return(1);
	}
	return(0);				/* didn't find one */
}


/* ================================== */
/*           New Functions            */
/* ================================== */

void
fill(int in, int out, int Input[in], int Output[out],int inter){
	/* is like Interpolate, but for days and months */
	int i,j;
	for(i=0;i<in; i++){
		int A=Input[i];

		for(j=0;j<inter;j++){
			Output[i*inter+j]=A;
		}
	}
}

void
interp_hour(int in, int out, float Input[in], float Output[out], int inter){
	/* Interpolate the hour. */
	int i,j;
	for(i=0;i<in; i++){
		float A=Input[i];
		float B=Input[i+1];
		float delta;
		if(B<A){ /* transition to next day */
			delta=24-A+B;
		}else{
			delta=B-A;
		}
				
		for(j=0;j<inter;j++){
			float O=A+j*delta/inter;
			if(O>=24){
				O=24-O;
			}	
			
			Output[i*inter+j]=O;	
		}
	}
}

void
interpolate(int in, int out, float Input[in], float Output[out], int inter){
	/* Interpolate the Input on "inter" subdivisions to return the Output */
	int i,j;
	for(i=0;i<in; i++){
		float A=Input[i];
		float delta=Input[i+1]-A;
		
		for(j=0;j<inter;j++){
			Output[i*inter+j]=A+j*delta/inter;	
		}
	}
}

/*float
date2n(int Month, int Day, int Hour_in, int Minute){
	float N;
	int b;
	
    if (Month==1){
        b=0;
    }else if(Month==2){
        b=31;
    }else if(Month==3){
        b=59;
    }else if(Month==4){
        b=90;
    }else if(Month==5){
        b=120;
    }else if(Month==6){
        b=151;
    }else if(Month==7){
        b=181;
    }else if(Month==8){
        b=212;
    }else if(Month==9){
        b=243;
    }else if(Month==10){
        b=273;
    }else if(Month==11){
        b=304;
    }else if(Month==12){
        b=334;
    }else{
        fprintf(stderr,"Incorrect month: input number not between 1 and 12");                                       
    }
    
    N=b+Day+(Hour_in*1.0-Minute*(0.5/60))/24;
	
	
	return N;
}

*/


/* add two matrices and allocate the result */
static CMATRIX *
cm_add(const CMATRIX *cm1, const CMATRIX *cm2)
{

	CMATRIX	*cmr;
	int	dr, dc;
	if ((cm1->ncols <= 0) | ((cm1->ncols != cm2->ncols) && (cm1->nrows != cm2->nrows)))
		error(INTERNAL, "matrix dimension mismatch in cm_add()");
	cmr = cm_alloc(cm2->nrows, cm2->ncols);
	if (cmr == NULL) {return(NULL);}

	for (dr = 0; dr < cmr->nrows; dr++) {

	    for (dc = 0; dc < cmr->ncols; dc++) {
		
			COLORV	*dp = cm_lval(cmr,dr,dc);
			dp[0] = dp[1] = dp[2] = 0;
		
		    const COLORV	*cp1 = cm_lval(cm1,dr,dc);
		    const COLORV	*cp2 = cm_lval(cm2,dr,dc);
		    dp[0] = cp1[0] + cp2[0];
		    dp[1] = cp1[1] + cp2[1];
		    dp[2] = cp1[2] + cp2[2]; 
		}
	}

	return(cmr);
}

/* scale a matrix by a certain scalar */
static CMATRIX *
cm_scale2(const CMATRIX *cm, float k)
{

	CMATRIX	*cmr;
	int	dr, dc;

	cmr = cm_alloc(cm->nrows, cm->ncols);
	if (cmr == NULL)
		return(NULL);
	
	for (dr = 0; dr < cmr->nrows; dr++) {
	    for (dc = 0; dc < cmr->ncols; dc++) {
		
			COLORV	*dp = cm_lval(cmr,dr,dc);
			dp[0] = dp[1] = dp[2] = 0;
		
		    const COLORV	*cp = cm_lval(cm,dr,dc);
		    dp[0] = k*cp[0];
		    dp[1] = k*cp[1];
		    dp[2] = k*cp[2];
		}
	}

	return(cmr);
}

/* Transform color matrix into illuminance values */
void cm_to_light(int cols, int rows, float light[rows][cols], const CMATRIX *cm)
{	
	
	int row;
	int col;
	
	for(row=0;row<cm->nrows;row++){
		for(col=0;col<cm->ncols;col++){
			const COLORV	*cp = cm_lval(cm,row,col);
			light[row][col]=179*(0.265*cp[0]+0.670*cp[1]+0.065*cp[2]);
		}
	}
}


void print_matrix(int cols, int rows, float mat[rows][cols]) {
	int i,j;
	
	for(i=0;i<rows;i++){
		for(j=0;j<cols;j++){
			printf("%f\t",mat[i][j]);
		}
		printf("\n");
	}
}


/* print out matrix as ASCII text -- no header */
void
print_cm(const CMATRIX *cm)
{
	int		r, c;
	const COLORV	*mp = cm->cmem;
	
	for (r = 0; r < cm->nrows; r++) {
		for (c = 0; c < cm->ncols; c++, mp += 3)
			printf("\t%.6e %.6e %.6e", mp[0], mp[1], mp[2]);
			printf("\n");
	}
}

/* Allocate a color coefficient matrix */
static CMATRIX *
zeroes(int nrows, int ncols)
{
	CMATRIX	*cmr;
	int	dr, dc;

	cmr = cm_alloc(nrows, ncols);
	
	for (dr = 0; dr < cmr->nrows; dr++) {
	    for (dc = 0; dc < cmr->ncols; dc++) {
		
			COLORV	*dp = cm_lval(cmr,dr,dc);
			dp[0] = dp[1] = dp[2] = 0.0;
		}
	}

	return(cmr);
}


/* ================== FUNCTIONS FOR LUA SCRIPTS =============== */
/* ============================================================ */

int 
now(lua_State *L) {
	/* Returns the actual time */
	lua_pushnumber(L,hr);
	return(1);
}

int
irradiance(lua_State *L){

	
	int args=lua_gettop(L);
	
	if(args==3 && lua_isnumber(L,-args)){	
		
		if(is_day==1){

			/*  sky_brightness, sky_clearness and sun_zenith are global*/
			float f11, f12, f13, f21, f22, f23, F1, F2, a, b, cosBeta, sinBeta, cosTheta, rho, R;
			
			
			float sun_direction[3];			
			vector(sun_direction,altitude,azimuth);
			
			float x=(float)lua_tonumber(L,-args); /* input direction */
			float y=(float)lua_tonumber(L,-args+1);
			float z=(float)lua_tonumber(L,-args+2);
			float norm=sqrt(x*x+y*y+z*z);
			
			x=x/norm;
			y=y/norm;
			z=z/norm;
			
			if(x*x+y*y < 0.001){
				cosBeta=1;
				sinBeta=0;
			}else{
				cosBeta=(x*x+y*y)/sqrt(x*x+y*y);
				sinBeta=sin(acos(cosBeta));
			}
			cosTheta=x*sun_direction[0]+y*sun_direction[1]+z*sun_direction[2]; /* dot product */

			if(sky_clearness <= 1.065){
				f11=-0.008;	f12=0.588; f13=-0.062; f21=-0.060; f22=0.072; f23=-0.022;
			}else if(sky_clearness > 1.065 && sky_clearness <= 1.230){
				f11=0.130;	f12=0.683; f13=-0.151; f21=-0.019; f22=0.066; f23=-0.029;	
			}else if(sky_clearness > 1.230 && sky_clearness <= 1.500){
				f11=0.330;	f12=0.487; f13=-0.221; f21=0.055; f22=-0.064; f23=-0.026;	
			}else if(sky_clearness > 1.500 && sky_clearness <= 1.950){
				f11=0.568;	f12=0.187; f13=-0.295; f21=0.109; f22=-0.152; f23=0.014;	
			}else if(sky_clearness > 1.950 && sky_clearness <= 2.800){
				f11=0.873;	f12=-0.392; f13=-0.362; f21=0.226; f22=-0.462; f23=0.001;	
			}else if(sky_clearness > 2.800 && sky_clearness <= 4.500){
				f11=1.132;	f12=-1.237; f13=-0.412; f21=0.288; f22=-0.823; f23=0.056;	
			}else if(sky_clearness > 4.500 && sky_clearness <= 6.200){
				f11=1.060;	f12=-1.6; f13=-0.359; f21=0.264; f22=-1.127; f23=0.131;	
			}else{
				f11=0.678;	f12=-0.327; f13=-0.250; f21=0.156; f22=-1.377; f23=0.251;	
			}
		
			a=dmax(0,cosTheta);
			b=dmax(.08715574274765836,cos(sun_zenith));
			F1=dmax(0,f11+f12*sky_brightness+f13*sun_zenith);
			F2=f21+f22*sky_brightness+f23*sun_zenith;
		
			rho=0.265*grefl[0]+0.670*grefl[1]+0.065*grefl[2];
			
			/* add diffuse and ground */
			R=dif*((1-F1)*0.5*(1+cosBeta)+F1*a/b+F2*sinBeta)+(dif+cos(sun_zenith)*dir)*rho*(1-cosBeta)*0.5;
			
			/* add direct if not altitude is over 0 */
			if (altitude>0){
				R=R+dir*cosTheta;
			}
			
			lua_pushnumber(L,R);
		}else{
			lua_pushnumber(L,0);
		}
		
		return(1);
	}else{
		fprintf(stderr,"warning: Bad input on 'irradiance' function on Lua script\n");
		lua_pushnumber(L,-1);
		return(1);
	}
		


}

int 
night(lua_State *L) {
	 /* Recieves one integer (the number of the sensor)
	and returns the value registered by that sensor during the night with all the 
	lights at full power */
	
	int args=lua_gettop(L);
	
	if(args==1 && lua_isnumber(L,-args) && (int)lua_tonumber(L,-args)>=0){		
		int s=(int)lua_tonumber(L,-args);
		
		CMATRIX *Tot;
		Tot = Lum_array[0];
		float night[MAXsensor][1];
		
		
		/* Calculate artificial lighting */
		if(l_bool==1){
			int j;
			for(j=1;j<nlum;j++){
				Tot=cm_add(Tot,Lum_array[j]);
			}
		}		

		/* Add everything and transform into Lux */
		cm_to_light(1, MAXsensor, night, Tot);
		lua_pushnumber(L,night[s-1][0]);

		
		return(1);
	}else{
		fprintf(stderr,"warning: Bad input on 'night' function on Lua script\n");
		lua_pushnumber(L,-1);
		return(1);
	}
}


int
update(lua_State *L){
	/* Updates the positions for new calculations within the script, allowing
	iteratively finding new states (i.e. add more and more luminaire power
	until we get 500 Lux) */
	
	/*First, we get and store the already defined values */
	
	int j;
	for(j=0;j<nwin;j++){ /* for each fenestration system, pull the position */
			/* Build the name of the variable */
			sprintf(var_name, "win%d", j+1);
			/* put results in the stack */
			lua_getglobal(L, var_name);
			if(lua_isnumber(L,1)){
				/* assign the value */
				win_state[j]=(int)lua_tonumber(L,1)-1;
			} 
			/*if the returned value is not a number (or does not exist), 
			use the previous value (do not modify) */
			lua_pop(L,1);
		} /* end loop collecting shading positions */
		
		if(l_bool){
			for(j=0;j<nlum;j++){ /* for each luminaire, pull the relative power */			
				/* Build the name of the variable */
				sprintf(var_name, "lum%d", j+1);
				/* put results in the stack */
				lua_getglobal(L, var_name);
				if(lua_isnumber(L,1)){
					/* assign the value */
					lum_state[j]=(float)lua_tonumber(L,1);
					/*if(lum_state[j]<0.01){
						lum_state[j]=0.0;
					}
					if(lum_state[j]>0.99){
						lum_state[j]=1.0;
					}*/
				} 
				/*if the returned value is not a number (or does not exist), 
				use the previous value (do not modify) */
				lua_pop(L,1);
			}	
		} /*end of collection luminaire power*/
	
	/* Then, we calculate the daylight, if needed */
	
	CMATRIX *Tot;
	Tot = zeroes(MAXsensor, 1);
			
	double		sda, sta;		
		mtx_data = resize_dmatrix(mtx_data, 1, nskypatch);
		if (dif > 1e-4) { /* if it is not night, make daylight calculations */
				
			dir_irrad = dir;
			diff_irrad = dif;
		
		
			/* =========================================================== */
			/*                COMPUTE THE 3PHASE SKY VECTOR                */
			/* =========================================================== */
	
			/* compute sky patch values */
			ComputeSky(mtx_data);
			AddDirect(mtx_data);
			
			/*******/
			CMATRIX *S3;
			S3 = cm_alloc(nskypatch, 1);
			
			/* Fill the Color sky vector */
			for (j = 0; j < nskypatch; j++) {
				int off=3*j;
				COLORV	*dp = cm_lval(S3,j,0);

				dp[0] = mtx_data[off];
				dp[1] = mtx_data[off+1]; 
				dp[2] = mtx_data[off+2];
				
			}
			
			/* =========================================================== */
			/*          CALCULATE SENSORS USING THE ACTUAL STATE           */
			/* =========================================================== */
	
			CMATRIX *fDC3;
			fDC3=DC3[0][win_state[0]];
			
			for(j=1;j<nwin;j++){
				fDC3=cm_add(fDC3,DC3[j][win_state[j]]);			
			}
			Tot=cm_multiply(fDC3, S3);
			
			
		}/*else {  if it is night, there is no daylight */
			/*for(j=0;j<nsen;j++){
				COLORV	*dp = cm_lval(Tot,j,0);
				dp[0] = dp[1] = dp[2] = 0;
			} Se inicia con Ceros... no es necesario.  
		} end if is day, else is night */ 	
			
			
		
		if(l_bool==1){
			/* Calculate artificial light contribution */
			for(j=0;j<nlum;j++){
				if(lum_state[j]>0){
					Tot=cm_add(Tot,cm_scale2(Lum_array[j],lum_state[j]));
				}
			}
		}		

		cm_to_light(1, MAXsensor, sensor, Tot);
	
	
	return(0);
}

int 
sensor_val(lua_State *L) {
	 /* Recieves one integer (the number of the sensor)
	and returns the value registered by that sensor using the previous positions.
	The idea is to know the new illuminance values when nothing is modified (position
	of shading systems and/or luminaire power). */
	
	int args=lua_gettop(L);
	
	if(args==1 && lua_isnumber(L,-args) && (int)lua_tonumber(L,-args)>=0){		
		int s=(int)lua_tonumber(L,-args);
		lua_pushnumber(L,sensor[s-1][0]);
		return(1);
	}else{
		fprintf(stderr,"warning: Bad input on 'sensor' function on Lua script\n");
		lua_pushnumber(L,-1);
		return(1);
	}
}

int 
ext_temp(lua_State *L) {
	 /* Returns the exterior dry-bulb temperature */
	/* No inputs required, no error when there are some */
	/*int args=lua_gettop(L);*/
	lua_pushnumber(L,ext_drybulb);
	return(1);
}

int 
power_fraction(lua_State *L) {
	 /* Recieves one integer (the number of the luminaire)
	and returns the actual power fraction (between 0 and 1). */
	
	int args=lua_gettop(L);
	
	if(args==1 && lua_isnumber(L,-args) && (int)lua_tonumber(L,-args)>=0){	/*Checks that only one argument, that is a number and greater or equal to 0 */	
		int s=(int)lua_tonumber(L,-args);
		
		lua_pushnumber(L,lum_state[s-1]);
		
		return(1);
	}else{
		fprintf(stderr,"warning: Bad input on 'lum' function on Lua script\n");
		lua_pushnumber(L,-1);
		return(1);
	}
}

int 
win_position(lua_State *L) {
	 /* Recieves one integer (the number of the window)
	and returns the actual position */
	
	int args=lua_gettop(L);
	
	if(args==1 && lua_isnumber(L,-args) && (int)lua_tonumber(L,-args)>=0){	/*Checks that only one argument, that is a number and greater or equal to 0 */	
		int s=(int)lua_tonumber(L,-args);
		
		lua_pushnumber(L,win_state[s-1]);
		
		return(1);
	}else{
		fprintf(stderr,"warning: Bad input on 'win' function on Lua script\n");
		lua_pushnumber(L,-1);
		return(1);
	}
}

int 
check_day(lua_State *L) {
	lua_pushnumber(L,is_day);
	return(1);
}

int 
sun_azimuth(lua_State *L) {
	lua_pushnumber(L,azimuth);
	return(1);
}

int 
sun_altitude(lua_State *L) {
	lua_pushnumber(L,altitude);
	return(1);
}

int 
f_sun_zenith(lua_State *L) {
	lua_pushnumber(L,1.5705-altitude);
	return(1);
}

/* ================================== */
/* ================================== */
/* ================================== */


int
main(int argc, char *argv[])
{
	char	buf[256];
	double	rotation = 0;		/* site rotation (degrees) */
	int	dir_is_horiz;		/* direct is meas. on horizontal? */

	int	ntsteps = 0;		/* number of rows in matrix */
	int	last_monthly = 0;	/* month of last report */
	int	year, mo, da, hour_in,minute;			/* month (1-12) and day (1-31) */
	int	mtx_offset;
	int	i, j;
	/* radiometric quantities... all the .wea files have a 1, according to
	 the code I based this program on. GML */
	input = 1;		
	dir_is_horiz = 0;
	
	float dir_norm_rad, dif_or_rad, dummy_float, dry_bulb;
    char keyword[2000]="";
	char epw_file[200]="";
	char wea_file[200] ="";
	char city[200] ="";
	char country[200] ="";
	char latitude[200] ="",longitude[200] ="",time_zone[200] ="",elevation[200] ="";
	char *epw_name;
	int nsteps=10000, inter=1, Nsteps; 
	int fbool=0, lbool=0, tbool=0, nTs=0;
	
	char *v_dir;
	char *t_dirs[MAXwin];
	char *d_dir;
	char *l_dir;
	char *s_fname;
	char *dc_fname;
	char *UsrScript;
	char *SchFile;
	
	char v_fname[1024];
	char t_fname[1024];
	char d_fname[1024];
	char l_fname[1024];


	/* In this program, in general, the lower case arrays are those before interpolating, 
	and those in upper case are those after the interpolation (longer ones) */
	
	float dir_norm[8784], dif_hor[8784], n[8784], tout[8784], hour[8784], in_lat, in_lon, in_t;
	int month[8784], day[8784];
	
	progname = argv[0];
					/* get options */

		
		
	for (i=1; i<argc ;i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'u':
					UsrScript=argv[i+1];
					lua_bool=1;
				break;
				
				case 'l':
					nlum=atoi(argv[i+1]);
				break;
				case 'L':
					if(hasNumberFormat(argv[i+1])) {
						l_dir=argv[i+1];
						l_bool=1;
					}
					else{
						fprintf(stderr, "%s input error: -V, -T, -D and -L need a numeric format. Ex: -V /V/Vmx%s.vmx\n",progname,"%d");
						return(1);
						}
				break;
				case 'x':
					nsen=atoi(argv[i+1]);
					sen_bool=1;
				break;
				case 'w':
					nwin=atoi(argv[i+1]);
					w_bool=1;
				break;
				case 'V': /*View matrix */
					if(hasNumberFormat(argv[i+1])) {
						v_dir=argv[i+1];
					}
					else{
						fprintf(stderr, "%s input error: -V, -T, -D and -L need a numeric format. Ex: -V /V/Vmx%s.vmx\n",progname,"%d");
						return(1);
						}
				break;

				case 'T': /* Transmission matrix */
					j=1;
					int pp=0;
					while(i+j<argc && argv[i+j][0] != '-'){
						
						if(hasNumberFormat(argv[i+j+1])) {
							nTs++;
							n_dirs[pp]=atoi(argv[i+j]);
							t_dirs[nTs-1]=argv[i+j+1];
							j++;
							j++;
							pp++;
							t_bool=1;
						}	
						else{
							fprintf(stderr, "%s input error: -V, -T, -D and -L need a numeric format. Ex: -V /V/Vmx%s.vmx\n",progname,"%d");
							return(1);
							}
					}	
				break;
						
				case 'D': /* Daylight matrices */
					if(hasNumberFormat(argv[i+1])) {
						d_dir=argv[i+1];
					}
					else{
						fprintf(stderr, "%s input error: -V, -T, -D and -L need a numeric format. Ex: -V /V/Vmx%s.vmx\n",progname,"%d");
						return(1);
						}
				break;
				
				case 'f':			/* Energy Plus Weather file*/
					fbool=1;
					epw_name=argv[i+1];
					break;
					
				case 'p':			/* Latitude, Longitude */
					in_lat = atof(argv[++i]);
					in_lon = atof(argv[++i]);
					lbool=1;
					break;	
				
				case 't':			/* timezone (GMT) */
					in_t=atoi(argv[i+1]);
					tbool=1;
					break;	
		
				case 'n':			/* nsteps of the weather file to simulate */
					nsteps=atoi(argv[i+1]);
					break;	
	
				case 'i':			/* interpolation for subhourly simulations */
					inter=atoi(argv[i+1]);
					break;
		
				case 'g':			/* ground reflectance */
					grefl[0] = atof(argv[++i]);
					grefl[1] = atof(argv[++i]);
					grefl[2] = atof(argv[++i]);
					break;
				case 'v':			/* verbose progress reports */
					verbose++;
					break;
					
				case 'o':			/* output schedule file */
					SchFile=argv[i+1];
					sch_bool=1;
					break;
				case 'm':			/* Reinhart subdivisions */
					rhsubdiv = atoi(argv[++i]);
					break;
				case 'c':			/* sky color */
					skycolor[0] = atof(argv[++i]);
					skycolor[1] = atof(argv[++i]);
					skycolor[2] = atof(argv[++i]);
					break;
				case 'd':			/* solar (direct) only */
					skycolor[0] = skycolor[1] = skycolor[2] = 0;
					if (suncolor[1] <= 1e-4)
						suncolor[0] = suncolor[1] = suncolor[2] = 1;
					break;
				case 's':			/* sky only (no direct) */
					suncolor[0] = suncolor[1] = suncolor[2] = 0;
					if (skycolor[1] <= 1e-4)
						skycolor[0] = skycolor[1] = skycolor[2] = 1;
					break;
				case 'r':			/* rotate distribution */
					if (argv[i][2] && argv[i][2] != 'z')
						fprintf(stderr, "Usage: %s -f File.epw [-n steps][-i int.][-l Lat Lon][-t time_zone][-v][-d|-s][-r deg][-m N][-g r g b][-c r g b][-o{f|d}]\n",
						progname);
						exit(1);
					rotation = atof(argv[++i]);
					break;
				case '5':			/* 5-phase calculation */
					nsuns = 1;
					fixed_sun_sa = 6.797e-05;
					break;
				/*default:
					goto userr;*/ /* for some reason this was causing troubles */
	
			}
		}
	}

	   char cwd[1024];
       if (getcwd(cwd, sizeof(cwd)) == NULL)
           fprintf(stdout, "Warning: Current working directory could not be found\n");
       

	
	if (fbool==0){
		fprintf(stderr, "Usage: %s -f File.epw -o Out.sch [-n steps][-i interpolation][-l Lat Lon][-t time_zone][-v][-d|-s][-r deg][-m N][-g r g b][-c r g b]\n",
				progname);
		exit(1);
	}
	if (sch_bool==-2){
		fprintf(stderr, "Usage: %s -f File.epw -o Out.sch [-n steps][-i interpolation][-l Lat Lon][-t time_zone][-v][-d|-s][-r deg][-m N][-g r g b][-c r g b]\n",
				progname);
		exit(1);
	}		
	
	/* If the number of systems assigned is different of the number of windows
	throw error. However, if it is just one, the same system will be applied for
	all windows */
	if (nTs!=nwin) {/* this was causing segmentation faults */
		/*if(nTs==1){
			int y;
			for(y=1;y<nwin;y++){
				printf("aaa\n");
				t_dirs[y]=t_dirs[0];
			}
		}
		else {*/
			fprintf(stderr, "%s input error: number of Ts do not match the number of windows\n",progname);
			exit(1);
		/*}*/
	}
	
	
	
		
	/* ====================================================== */
	/*                 initialize sky patches                 */
	/* ====================================================== */
	/* This does not make any calculations */

	rh_init();	

	if (verbose) {
		fprintf(stderr, "%s: location '%s'\n", progname, buf);
		fprintf(stderr, "%s: (lat,long)=(%.1f,%.1f) degrees north, west\n",
				progname, s_latitude, s_longitude);
		fprintf(stderr, "%s: %d sky patches per time step\n",
				progname, nskypatch);
		if (rotation != 0)
			fprintf(stderr, "%s: rotating output %.0f degrees\n",
					progname, rotation);
	}
					/* convert quantities to radians */
	
	
	/* ====================================================== */
	/*    LOAD LUMINAIRE, VIEW, BTDF AND DAYLIGHT MATRICES    */
	/* ====================================================== */
	/* Assemble and store the DC matrices for avoiding continuos matrix multiplication */
	
	

	COLOR	tLamb;
	for(i=0;i<nwin;i++){
		CMATRIX	*V3, *D3;
		
		int npos=n_dirs[i];
		sprintf(v_fname, v_dir, i+1);
		sprintf(d_fname, d_dir, i+1);

		V3 = cm_load(v_fname,nsen,145, DTfromHeader);
		D3 = cm_load(d_fname,145,nskypatch, DTfromHeader);
		for(j=0;j<npos;j++){

			CMATRIX *T, *TD;
			sprintf(t_fname, t_dirs[i], j+1);
			T = cm_loadBSDF(t_fname, tLamb);
			TD=cm_multiply(T,D3);
			cm_free(T);
			DC3[i][j] = cm_multiply(V3,TD);
			cm_free(TD);
		}
		cm_free(V3);cm_free(D3);
	}

	
	for(i=0;i<nlum;i++){
		sprintf(l_fname, l_dir, i+1);
		Lum_array[i]=cm_load(l_fname,nsen,1, DTfromHeader);
	}
	
	

	/* ====================================================== */
	/*                      READ EPW FILE                     */
	/* ====================================================== */


	FILE *EPW_FILE=fopen(epw_name,"r");
	
	fscanf(EPW_FILE,"%[^,]s",keyword);
	if( !strcmp(keyword,"LOCATION") ){
		fscanf(EPW_FILE,",%[^,]s",city);
		fscanf(EPW_FILE,",%[^,]s",country);
		fscanf(EPW_FILE,",%[^,]s",country);
		sprintf(keyword,"%s_%s",city,country);

		fscanf(EPW_FILE,",%[^,]s",country);
		fscanf(EPW_FILE,",%[^,]s",country);
		fscanf(EPW_FILE,",%[^,]s",latitude);

		fscanf(EPW_FILE,",%[^,]s",longitude);
		fscanf(EPW_FILE,",%[^,]s",time_zone);
	
		fscanf(EPW_FILE,",%s[^\n]",elevation);

		

		fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
		fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
		fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
		fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
		fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
		fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
		fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");
		fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");


		/* Scan for outside temperature, "N" and the radiation data */
		i=0;
		while( EOF != fscanf(EPW_FILE,"%d,%d,%d,%d,%d",&year,&mo,&da, &hour_in,&minute) && i<nsteps){
			
			fscanf(EPW_FILE,",%[^,]s",city);			
			fscanf(EPW_FILE,",%f",&dry_bulb);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);
			fscanf(EPW_FILE,",%f",&dummy_float);

			fscanf(EPW_FILE,",%f,%f",&dir_norm_rad, &dif_or_rad);

			
			/*n[i]=date2n(mo, da,hour_in, minute);*/
			dir_norm[i]=dir_norm_rad;
			dif_hor[i]=dif_or_rad;
			tout[i]=dry_bulb;
			hour[i]=hour_in*1.0-minute*(0.5/60);
			month[i]=mo;
			day[i]=da;
			
			
			
			fscanf(EPW_FILE,"%*[^\n]");fscanf(EPW_FILE,"%*[\n\r]");

			i++;
		}
		nsteps=i; /*the number of steps, just in case it was not defined. */
	}else{
			fprintf(stderr, "%s: weather file error\n", progname);
			return(1);
	}

	fclose(EPW_FILE);
	/* End of EPW scan */
	
	if(lbool==0){
		s_latitude = DegToRad(atof(latitude));
		s_longitude = DegToRad(-atof(longitude));
	}else{
		s_latitude = DegToRad(in_lat);
		s_longitude = DegToRad(-in_lon);
	}
	
	if(tbool==0){
		s_meridian = DegToRad(-15*atof(time_zone));
	}else{
		s_meridian = DegToRad(-15*in_t);
	}
	
	/* ====================================================== */
	/*                 PRINT HEADER OF RESULTS                */
	/* ====================================================== */
	printf("Command,"); 
	for (i=0; i<argc ;i++) {
		if(i<argc-1){
			printf("%s,",argv[i]);
		}else{
			printf("%s",argv[i]);	
		}
	}
	printf("\n");
	printf("Sensors,%d\n",nsen);
	printf("Lums,%d\n",nlum);
	printf("Wins,%d\n",nwin);
	printf("Timezone,%f\n",-s_meridian*3.819718); //This number sums up the -15 multiplication and the conversion from Radians to Degrees
	printf("Location,%f, %f\n",s_latitude*180.0/3.141592654,-s_longitude*180.0/3.141592654);
	printf("Script,%s/%s\n",cwd,UsrScript);		
	printf("WeatherFile,%s\n",epw_name);
	printf("Dmatrix,%s/%s\n",cwd,d_dir);
	printf("Tmatrix,");
		
	for(i=0;i<nwin;i++){
		if(i<nwin-1){
			printf("%d,%s/%s,",n_dirs[i],cwd,t_dirs[i]);
		}else{
			printf("%d,%s/%s",n_dirs[i],cwd,t_dirs[i]);
		}		
	}
	printf("\nMonth,Day,Hour,Dif,Dir,");
	for(i=0;i<nsen;i++){
		printf("Sen%d,",i+1);
	}
	for(i=0;i<nwin;i++){
		printf("Win%d,",i+1);
	}
	if(l_bool==1){
		for(i=0;i<nlum;i++){
			if(i<nlum-1){
				printf("Lum%d,",i+1);
			}else{
				printf("Lum%d",i+1);
			}
		}
	}	
	printf("\n");
		

	/* ====================================================== */
	/*                      INTERPOLATION                     */
	/* ====================================================== */

	
	Nsteps=(nsteps-1)*inter+1;
	
	float Dir_norm[Nsteps], Dif_hor[Nsteps], N[Nsteps], Tout[Nsteps], Hour[Nsteps];
	int Month[Nsteps], Day[Nsteps]; 
	
	interpolate(nsteps, Nsteps, tout, Tout,inter);
	interpolate(nsteps, Nsteps, n, N,inter);
	interpolate(nsteps, Nsteps, dir_norm, Dir_norm,inter);	
	interpolate(nsteps, Nsteps, dif_hor, Dif_hor,inter);
	
	interp_hour(nsteps, Nsteps, hour, Hour,inter);
	
	fill(nsteps, Nsteps, day, Day,inter);
	fill(nsteps, Nsteps, month, Month,inter);
	
	
	/* ====================================================== */
	/*                     DO SOME LUA STUFF                  */
	/* ====================================================== */
	

	/* initialize Lua  and load libraries */
	L=luaL_newstate();
	luaL_openlibs(L);

	/* Register functions */
	lua_register(L,"sensor",sensor_val); /*returns illuminance on a sensor*/
	lua_register(L,"update",update); /*updates illuminance within a script */
	lua_register(L,"ext_temp",ext_temp); /*returns the exterior temperature from the EPW*/
	lua_register(L,"lum",power_fraction); /*returns the previous fraction of the lighting power of the luminaire*/
	lua_register(L,"win",win_position); /*returns the previous shading position*/	
	lua_register(L,"is_day",check_day); /* check if it is day */
	lua_register(L,"night",night); /*  Value of a sensor at night, with all the lights on.*/
	lua_register(L,"altitude",sun_altitude);  /*  Solar altitude position */
	lua_register(L,"azimuth",sun_azimuth); /*  Solar azimuth position */
	lua_register(L,"zenith",f_sun_zenith); /*  Solar azimuth position */
	lua_register(L,"irradiance",irradiance);  /*  Solar irradiance over a sensor... assumed outside and unshaded */
	lua_register(L,"hour",now);  /*  Hour of day */
	
	/* ====================================================== */
	/*            SET VARIABLES FOR 1st ITERATION             */
	/* ====================================================== */

	for(i=0;i<nwin;i++){ win_state[i]=0; }
	for(i=0;i<nlum;i++){ lum_state[i]=0; }
	
		
	/* ====================================================== */
	/*                   PROCESS TIMESTEPS                    */
	/* ====================================================== */
	


	
	for(i=0;i<Nsteps;i++){
		
			/* Declare the matrix that contains the illuminance of the sensors */
		CMATRIX *Tot;
		Tot = zeroes(nsen, 1);
		
		int tot_bool=0;
		
		/* =========================================================== */
		/*                  GET THE RADIATION VALUES                   */
		/* =========================================================== */
		
		mo=Month[i];
		da=Day[i];
		hr=Hour[i];
		dir=Dir_norm[i];
		dif=Dif_hor[i];


		double		sda, sta;
					/* make space for next time step */
			
		mtx_offset = 0;/*3*nskypatch*ntsteps++;*/
		mtx_data = resize_dmatrix(mtx_data, 1, nskypatch);
		if (dif <= 1e-4) { /* is night */
			is_day=0;
			/*for(j=0;j<nsen;j++){
				COLORV	*dp = cm_lval(Tot,j,0);
				dp[0] = dp[1] = dp[2] = 0;
			}
			continue;*/
		} else {
			is_day=1;
			if (verbose && mo != last_monthly)
				fprintf(stderr, "%s: stepping through month %d...\n",
							progname, last_monthly=mo);
						/* compute solar position */
			julian_date = jdate(mo, da);
			sda = sdec(julian_date);
			sta = stadj(julian_date);
			altitude = salt(sda, hr+sta);
			azimuth = sazi(sda, hr+sta) + PI - DegToRad(rotation);
						/* convert measured values */


			/* Some stuff was deleted here since the "input" is always 1, and 	
			dir_is_horiz is always 0 */
		
			if (dir_is_horiz && altitude > 0.)
				dir /= sin(altitude);
				
			dir_irrad = dir;
			diff_irrad = dif;
		
						/* compute sky patch values */
			ComputeSky(mtx_data+mtx_offset);
			AddDirect(mtx_data+mtx_offset);	
			/*******/
			CMATRIX *S3;
			S3 = cm_alloc(nskypatch, 1);	
			/* Fill the Color sky vector */
			for (j = 0; j < nskypatch; j++) {
				int off=3*j;
				COLORV	*dp = cm_lval(S3,j,0);

				dp[0] = mtx_data[off];
				dp[1] = mtx_data[off+1]; 
				dp[2] = mtx_data[off+2];				
			}
		
			/* =========================================================== */
			/*                     CALCULATE DAYLIGHT                      */
			/* =========================================================== */

			tot_bool=1;

			CMATRIX *fDC3;
			fDC3=DC3[0][win_state[0]];

			for(j=1;j<nwin;j++){
				fDC3=cm_add(fDC3,DC3[j][win_state[j]]);			
			}

			Tot=cm_multiply(fDC3, S3);
			cm_free(S3);
		}		
		/* =========================================================== */
		/*                 CALCULATE ARTIFICIAL LIGHT                  */
		/* =========================================================== */

		if(l_bool){
			/* Calculate artificial light contribution */
			for(j=0;j<nlum;j++){			
				if(lum_state[j]>0.01){
					Tot=cm_add(Tot,cm_scale2(Lum_array[j],lum_state[j]));
					tot_bool=1;
				}
			}
		}		

		
		/* =========================================================== */
		/*                        RUN LUA SCRIPT                       */
		/* =========================================================== */

		/* run the script */
    	luaL_dofile(L, UsrScript);

    	/* =========================================================== */
		/*                 GET THE NEW STATE VARIABLES                 */
		/* =========================================================== */
		
		
		for(j=0;j<nwin;j++){ /* for each fenestration system, pull the position */
			/* Build the name of the variable */
			sprintf(var_name, "win%d", j+1);
			/* put results in the stack */
			lua_getglobal(L, var_name);
			if(lua_isnumber(L,1)){
				/* assign the value */
				win_state[j]=(int)lua_tonumber(L,1)-1;
			} 
			/*if the returned value is not a number (or does not exist), 
			use the previous value (do not modify) */
			lua_pop(L,1);
		} /* end loop collecting shading positions */
		
		if(l_bool==1){
			for(j=0;j<nlum;j++){ /* for each luminaire, pull the relative power */			
				/* Build the name of the variable */
				sprintf(var_name, "lum%d", j+1);
				/* put results in the stack */
				lua_getglobal(L, var_name);
				if(lua_isnumber(L,1)){
					/* assign the value */
					lum_state[j]=(float)lua_tonumber(L,1);
					if(lum_state[j]<0.01){
						lum_state[j]=0.0;
					}
					if(lum_state[j]>0.99){
						lum_state[j]=1.0;
					}
				} 
				/*if the returned value is not a number (or does not exist), 
				use the previous value (do not modify) */
				lua_pop(L,1);
			}	
		} /*end of collection luminaire power*/
		  
	   
	   
	   
	   
	   	/* =========================================================== */
		/* =========================================================== */
	   	/* =========================================================== */
	   
	   	Tot = zeroes(nsen, 1); // Clean Tot matrix	   
	 			
		mtx_offset = 0;/*3*nskypatch*ntsteps++;*/
		mtx_data = resize_dmatrix(mtx_data, 1, nskypatch);
		if (dif <= 1e-4) { /* is night */
			is_day=0;
			/*for(j=0;j<nsen;j++){
				COLORV	*dp = cm_lval(Tot,j,0);
				dp[0] = dp[1] = dp[2] = 0;
			}
			continue;*/
		} else {
			is_day=1;
			if (verbose && mo != last_monthly)
				fprintf(stderr, "%s: stepping through month %d...\n",
							progname, last_monthly=mo);
						/* compute solar position */
			julian_date = jdate(mo, da);
			sda = sdec(julian_date);
			sta = stadj(julian_date);
			altitude = salt(sda, hr+sta);
			azimuth = sazi(sda, hr+sta) + PI - DegToRad(rotation);
						/* convert measured values */


			/* Some stuff was deleted here since the "input" is always 1, and 	
			dir_is_horiz is always 0 */
		
			if (dir_is_horiz && altitude > 0.)
				dir /= sin(altitude);
				
			dir_irrad = dir;
			diff_irrad = dif;
		
						/* compute sky patch values */
			ComputeSky(mtx_data+mtx_offset);
			AddDirect(mtx_data+mtx_offset);	

			/*******/
			CMATRIX *S3;
			S3 = cm_alloc(nskypatch, 1);
		
			/* Fill the Color sky vector */
			for (j = 0; j < nskypatch; j++) {
				int off=3*j;
				COLORV	*dp = cm_lval(S3,j,0);

				dp[0] = mtx_data[off];
				dp[1] = mtx_data[off+1]; 
				dp[2] = mtx_data[off+2];				
			}
		
			/* =========================================================== */
			/*                     CALCULATE DAYLIGHT                      */
			/* =========================================================== */

			tot_bool=1;
			CMATRIX *fDC3;
			fDC3=DC3[0][win_state[0]];
		
			for(j=1;j<nwin;j++){
				fDC3=cm_add(fDC3,DC3[j][win_state[j]]);			
			}
			Tot=cm_multiply(fDC3, S3);
			cm_free(S3);
		}		
		/* =========================================================== */
		/*                 CALCULATE ARTIFICIAL LIGHT                  */
		/* =========================================================== */
		if(l_bool){
			/* Calculate artificial light contribution */
			for(j=0;j<nlum;j++){
				if(lum_state[j]>0.01){
					Tot=cm_add(Tot,cm_scale2(Lum_array[j],lum_state[j]));
					tot_bool=1;
				}
			}
		}	
	   
	   	/* =========================================================== */
		/* =========================================================== */
	   	/* =========================================================== */
	   
	   
	   
	   
	   
	   
		/* =========================================================== */
		/*                PRINT THE NEW STATE VARIABLES                */
		/* =========================================================== */
		printf("%d,%d,%f,%f,%f,",mo,da,hr,dif,dir);
		if(tot_bool){
			cm_to_light(1, MAXsensor, sensor, Tot);
			for(j=0;j<nsen;j++){
				printf("%f,",sensor[j][0]);
			}
		}else{
			for(j=0;j<nsen;j++){
				printf("%f,",0.0);
			}
		}	
		for(j=0;j<nwin;j++){
			printf("%d,",win_state[j]+1);
		}
		
		if(l_bool==1){
			for(j=0;j<nlum;j++){
				if(j<nlum-1){
					printf("%f,",lum_state[j]);
				}else{
					printf("%f",lum_state[j]);
				}
			}
		}	
		printf("\n");
		
		
		cm_free(Tot);
			   		
	}/* end loop timesteps */
	
				
	
	exit(0);
	
	writerr:
	{	fprintf(stderr, "%s: write error on output\n", progname);
		exit(1);
	}
	epwErr:
	{	fprintf(stderr,"FATAL ERROR - the epw file does not seem to be correct\n");
		exit(1);
	}
	
	NotNumericErr:
	{
		fprintf(stderr, "%s input error: -V, -T, -D and -L need a numeric format. Ex: -V /V/Vmx%s.vmx\n",progname,"%d");
		return(1);
	}
}


/* Return maximum of two doubles */
double dmax( double a, double b )
{ return (a > b) ? a : b; }


void
ComputeSky(float *parr)
{
	int index;			/* Category index */
	double norm_diff_illum;		/* Normalized diffuse illuimnance */
	int i;
	
	/* Calculate atmospheric precipitable water content */
	apwc = CalcPrecipWater(dew_point);

	/* Calculate sun zenith angle (don't let it dip below horizon) */
	/* Also limit minimum angle to keep circumsolar off zenith */
	if (altitude <= 0.0)
		sun_zenith = DegToRad(90.0);
	else if (altitude >= DegToRad(87.0))
		sun_zenith = DegToRad(3.0);
	else
		sun_zenith = DegToRad(90.0) - altitude;

	/* Compute the inputs for the calculation of the sky distribution */
	
	if (input == 0)					/* XXX never used */
	{
		/* Calculate irradiance */
		diff_irrad = CalcDiffuseIrradiance();
		dir_irrad = CalcDirectIrradiance();
		
		/* Calculate illuminance */
		index = GetCategoryIndex();
		diff_illum = diff_irrad * CalcDiffuseIllumRatio(index);
		dir_illum = dir_irrad * CalcDirectIllumRatio(index);
	}
	else if (input == 1)
	{
		sky_brightness = CalcSkyBrightness();
		sky_clearness =  CalcSkyClearness();

		/* Limit sky clearness */
		if (sky_clearness > 11.9)
			sky_clearness = 11.9;

		/* Limit sky brightness */
		if (sky_brightness < 0.01)
			sky_brightness = 0.01;

		/* Calculate illuminance */
		index = GetCategoryIndex();
		diff_illum = diff_irrad * CalcDiffuseIllumRatio(index);
		dir_illum = dir_irrad * CalcDirectIllumRatio(index);
	}
	else if (input == 2)
	{
		/* Calculate sky brightness and clearness from illuminance values */
		index = CalcSkyParamFromIllum();
	}

	if (output == 1) {			/* hack for solar radiance */
		diff_illum = diff_irrad * WHTEFFICACY;
		dir_illum = dir_irrad * WHTEFFICACY;
	}

	if (bright(skycolor) <= 1e-4) {			/* 0 sky component? */
		memset(parr, 0, sizeof(float)*3*nskypatch);
		return;
	}
	/* Compute ground radiance (include solar contribution if any) */
	parr[0] = diff_illum;
	if (altitude > 0)
		parr[0] += dir_illum * sin(altitude);
	parr[2] = parr[1] = parr[0] *= (1./PI/WHTEFFICACY);
	multcolor(parr, grefl);

	/* Calculate Perez sky model parameters */
	CalcPerezParam(sun_zenith, sky_clearness, sky_brightness, index);

	/* Calculate sky patch luminance values */
	CalcSkyPatchLumin(parr);

	/* Calculate relative horizontal illuminance */
	norm_diff_illum = CalcRelHorzIllum(parr);

	/* Normalization coefficient */
	norm_diff_illum = diff_illum / norm_diff_illum;

	/* Apply to sky patches to get absolute radiance values */
	for (i = 1; i < nskypatch; i++) {
		scalecolor(parr+3*i, norm_diff_illum*(1./WHTEFFICACY));
		multcolor(parr+3*i, skycolor);
	}
}

/* Add in solar direct to nearest sky patches (GW) */
void
AddDirect(float *parr)
{
	FVECT	svec;
	double	near_dprod[NSUNPATCH];
	int	near_patch[NSUNPATCH];
	double	wta[NSUNPATCH], wtot;
	int	i, j, p;

	if (dir_illum <= 1e-4 || bright(suncolor) <= 1e-4)
		return;
					/* identify nsuns closest patches */
	if (nsuns > NSUNPATCH)
		nsuns = NSUNPATCH;
	else if (nsuns <= 0)
		nsuns = 1;
	for (i = nsuns; i--; )
		near_dprod[i] = -1.;
	vector(svec, altitude, azimuth);
	for (p = 1; p < nskypatch; p++) {
		FVECT	pvec;
		double	dprod;
		rh_vector(pvec, p);
		dprod = DOT(pvec, svec);
		for (i = 0; i < nsuns; i++)
			if (dprod > near_dprod[i]) {
				for (j = nsuns; --j > i; ) {
					near_dprod[j] = near_dprod[j-1];
					near_patch[j] = near_patch[j-1];
				}
				near_dprod[i] = dprod;
				near_patch[i] = p;
				break;
			}
	}
	wtot = 0;			/* weight by proximity */
	for (i = nsuns; i--; )
		wtot += wta[i] = 1./(1.002 - near_dprod[i]);
					/* add to nearest patch radiances */
	for (i = nsuns; i--; ) {
		float	*pdest = parr + 3*near_patch[i];
		float	val_add = wta[i] * dir_illum / (WHTEFFICACY * wtot);

		val_add /= (fixed_sun_sa > 0)	? fixed_sun_sa 
						: rh_dom[near_patch[i]] ;
		*pdest++ += val_add*suncolor[0];
		*pdest++ += val_add*suncolor[1];
		*pdest++ += val_add*suncolor[2];
	}
}

/* Initialize Reinhart sky patch positions (GW) */
int
rh_init(void)
{
#define	NROW	7
	static const int	tnaz[NROW] = {30, 30, 24, 24, 18, 12, 6};
	const double		alpha = (PI/2.)/(NROW*rhsubdiv + .5);
	int			p, i, j;
					/* allocate patch angle arrays */
	nskypatch = 0;
	for (p = 0; p < NROW; p++)
		nskypatch += tnaz[p];
	nskypatch *= rhsubdiv*rhsubdiv;
	nskypatch += 2;
	rh_palt = (float *)malloc(sizeof(float)*nskypatch);
	rh_pazi = (float *)malloc(sizeof(float)*nskypatch);
	rh_dom = (float *)malloc(sizeof(float)*nskypatch);
	if ((rh_palt == NULL) | (rh_pazi == NULL) | (rh_dom == NULL)) {
		fprintf(stderr, "%s: out of memory in rh_init()\n", progname);
		exit(1);
	}
	rh_palt[0] = -PI/2.;		/* ground & zenith patches */
	rh_pazi[0] = 0.;
	rh_dom[0] = 2.*PI;
	rh_palt[nskypatch-1] = PI/2.;
	rh_pazi[nskypatch-1] = 0.;
	rh_dom[nskypatch-1] = 2.*PI*(1. - cos(alpha*.5));
	p = 1;				/* "normal" patches */
	for (i = 0; i < NROW*rhsubdiv; i++) {
		const float	ralt = alpha*(i + .5);
		const int	ninrow = tnaz[i/rhsubdiv]*rhsubdiv;
		const float	dom = 2.*PI*(sin(alpha*(i+1)) - sin(alpha*i)) /
						(double)ninrow;
		for (j = 0; j < ninrow; j++) {
			rh_palt[p] = ralt;
			rh_pazi[p] = 2.*PI * j / (double)ninrow;
			rh_dom[p++] = dom;
		}
	}
	return nskypatch;
#undef NROW
}



/* Resize daylight matrix (GW) */
float *
resize_dmatrix(float *mtx_data, int nsteps, int npatch)
{
	if (mtx_data == NULL)
		mtx_data = (float *)malloc(sizeof(float)*3*nsteps*npatch);
	else
		mtx_data = (float *)realloc(mtx_data,
					sizeof(float)*3*nsteps*npatch);
	if (mtx_data == NULL) {
		fprintf(stderr, "%s: out of memory in resize_dmatrix(%d,%d)\n",
				progname, nsteps, npatch);
		exit(1);
	}
	return(mtx_data);
}

/* Determine category index */
int GetCategoryIndex()
{
	int index;	/* Loop index */

	for (index = 0; index < 8; index++)
		if ((sky_clearness >= SkyClearCat[index].lower) &&
				(sky_clearness < SkyClearCat[index].upper))
			break;

	return index;
}

/* Calculate diffuse illuminance to diffuse irradiance ratio */

/* Reference:	Perez, R., P. Ineichen, R. Seals, J. Michalsky, and R. */
/*				Stewart. 1990. ìModeling Daylight Availability and */
/*				Irradiance Components from Direct and Global */
/*				Irradiance,î Solar Energy 44(5):271-289, Eqn. 7. */

double CalcDiffuseIllumRatio( int index )
{
	ModelCoeff const *pnle;	/* Category coefficient pointer */
	
	/* Get category coefficient pointer */
	pnle = &(DiffuseLumEff[index]);

	return pnle->a + pnle->b * apwc + pnle->c * cos(sun_zenith) +
			pnle->d * log(sky_brightness);
}

/* Calculate direct illuminance to direct irradiance ratio */

/* Reference:	Perez, R., P. Ineichen, R. Seals, J. Michalsky, and R. */
/*				Stewart. 1990. ìModeling Daylight Availability and */
/*				Irradiance Components from Direct and Global */
/*				Irradiance,î Solar Energy 44(5):271-289, Eqn. 8. */

double CalcDirectIllumRatio( int index )
{
	ModelCoeff const *pnle;	/* Category coefficient pointer */

	/* Get category coefficient pointer */
	pnle = &(DirectLumEff[index]);

	/* Calculate direct illuminance from direct irradiance */
	
	return dmax((pnle->a + pnle->b * apwc + pnle->c * exp(5.73 *
			sun_zenith - 5.0) + pnle->d * sky_brightness),
			0.0);
}

/* Calculate sky brightness */

/* Reference:	Perez, R., P. Ineichen, R. Seals, J. Michalsky, and R. */
/*				Stewart. 1990. ìModeling Daylight Availability and */
/*				Irradiance Components from Direct and Global */
/*				Irradiance,î Solar Energy 44(5):271-289, Eqn. 2. */

double CalcSkyBrightness()
{
	return diff_irrad * CalcAirMass() / (DC_SolarConstantE *
			CalcEccentricity());
}

/* Calculate sky clearness */

/* Reference:	Perez, R., P. Ineichen, R. Seals, J. Michalsky, and R. */
/*				Stewart. 1990. ìModeling Daylight Availability and */
/*				Irradiance Components from Direct and Global */
/*				Irradiance,î Solar Energy 44(5):271-289, Eqn. 1. */

double CalcSkyClearness()
{
	double sz_cubed;	/* Sun zenith angle cubed */

	/* Calculate sun zenith angle cubed */
	sz_cubed = sun_zenith*sun_zenith*sun_zenith;

	return ((diff_irrad + dir_irrad) / diff_irrad + 1.041 *
			sz_cubed) / (1.0 + 1.041 * sz_cubed);
}

/* Calculate diffuse horizontal irradiance from Perez sky brightness */

/* Reference:	Perez, R., P. Ineichen, R. Seals, J. Michalsky, and R. */
/*				Stewart. 1990. ìModeling Daylight Availability and */
/*				Irradiance Components from Direct and Global */
/*				Irradiance,î Solar Energy 44(5):271-289, Eqn. 2 */
/*				(inverse). */

double CalcDiffuseIrradiance()
{
	return sky_brightness * DC_SolarConstantE * CalcEccentricity() /
			CalcAirMass();
}

/* Calculate direct normal irradiance from Perez sky clearness */

/* Reference:	Perez, R., P. Ineichen, R. Seals, J. Michalsky, and R. */
/*				Stewart. 1990. ìModeling Daylight Availability and */
/*				Irradiance Components from Direct and Global */
/*				Irradiance,î Solar Energy 44(5):271-289, Eqn. 1 */
/*				(inverse). */

double CalcDirectIrradiance()
{
	return CalcDiffuseIrradiance() * ((sky_clearness - 1.0) * (1 + 1.041
			* sun_zenith*sun_zenith*sun_zenith));
}

/* Calculate sky brightness and clearness from illuminance values */
int CalcSkyParamFromIllum()
{
	double test1 = 0.1;
	double test2 = 0.1;
	int	counter = 0;
	int index = 0;			/* Category index */

	/* Convert illuminance to irradiance */
	diff_irrad = diff_illum * DC_SolarConstantE /
			(DC_SolarConstantL * 1000.0);
	dir_irrad = dir_illum * DC_SolarConstantE /
			(DC_SolarConstantL * 1000.0);

	/* Calculate sky brightness and clearness */
	sky_brightness = CalcSkyBrightness();
	sky_clearness =  CalcSkyClearness(); 

	/* Limit sky clearness */
	if (sky_clearness > 12.0)
		sky_clearness = 12.0;

	/* Limit sky brightness */
	if (sky_brightness < 0.01)
			sky_brightness = 0.01; 

	while (((fabs(diff_irrad - test1) > 10.0) ||
			(fabs(dir_irrad - test2) > 10.0)) && !(counter == 5))
	{
		test1 = diff_irrad;
		test2 = dir_irrad;	
		counter++;
	
		/* Convert illuminance to irradiance */
		index = GetCategoryIndex();
		diff_irrad = diff_illum / CalcDiffuseIllumRatio(index);
		dir_irrad = dir_illum / CalcDirectIllumRatio(index);
	
		/* Calculate sky brightness and clearness */
		sky_brightness = CalcSkyBrightness();
		sky_clearness =  CalcSkyClearness();

		/* Limit sky clearness */
		if (sky_clearness > 12.0)
			sky_clearness = 12.0;
	
		/* Limit sky brightness */
		if (sky_brightness < 0.01)
			sky_brightness = 0.01; 
	}

	return GetCategoryIndex();
}		

/* Calculate relative luminance */

/* Reference:	Perez, R., R. Seals, and J. Michalsky. 1993. */
/*				ìAll-Weather Model for Sky Luminance Distribution - */
/*				Preliminary Configuration and Validation,î Solar Energy */
/*				50(3):235-245, Eqn. 1. */

double CalcRelLuminance( double gamma, double zeta )
{
	return (1.0 + perez_param[0] * exp(perez_param[1] / cos(zeta))) *
		    (1.0 + perez_param[2] * exp(perez_param[3] * gamma) +
			perez_param[4] * cos(gamma) * cos(gamma));
}

/* Calculate Perez sky model parameters */

/* Reference:	Perez, R., R. Seals, and J. Michalsky. 1993. */
/*				ìAll-Weather Model for Sky Luminance Distribution - */
/*				Preliminary Configuration and Validation,î Solar Energy */
/*				50(3):235-245, Eqns. 6 - 8. */

void CalcPerezParam( double sz, double epsilon, double delta,
		int index )
{
	double x[5][4];		/* Coefficents a, b, c, d, e */
	int i, j;			/* Loop indices */

	/* Limit sky brightness */
	if (epsilon > 1.065 && epsilon < 2.8)
	{
		if (delta < 0.2)
			delta = 0.2;
	}

	/* Get Perez coefficients */
	for (i = 0; i < 5; i++)
		for (j = 0; j < 4; j++)
			x[i][j] = PerezCoeff[index][4 * i + j];

	if (index != 0)
	{
		/* Calculate parameter a, b, c, d and e (Eqn. 6) */
		for (i = 0; i < 5; i++)
			perez_param[i] = x[i][0] + x[i][1] * sz + delta * (x[i][2] +
					x[i][3] * sz);
	}
	else
	{
		/* Parameters a, b and e (Eqn. 6) */
		perez_param[0] = x[0][0] + x[0][1] * sz + delta * (x[0][2] +
				x[0][3] * sz);
		perez_param[1] = x[1][0] + x[1][1] * sz + delta * (x[1][2] +
				x[1][3] * sz);
		perez_param[4] = x[4][0] + x[4][1] * sz + delta * (x[4][2] +
				x[4][3] * sz);

		/* Parameter c (Eqn. 7) */
		perez_param[2] = exp(pow(delta * (x[2][0] + x[2][1] * sz),
				x[2][2])) - x[2][3];

		/* Parameter d (Eqn. 8) */
		perez_param[3] = -exp(delta * (x[3][0] + x[3][1] * sz)) + 
				x[3][2] + delta * x[3][3];
	}
}

/* Calculate relative horizontal illuminance (modified by GW) */

/* Reference:	Perez, R., R. Seals, and J. Michalsky. 1993. */
/*				ìAll-Weather Model for Sky Luminance Distribution - */
/*				Preliminary Configuration and Validation,î Solar Energy */
/*				50(3):235-245, Eqn. 3. */

double CalcRelHorzIllum( float *parr )
{
	int i;
	double rh_illum = 0.0;	/* Relative horizontal illuminance */

	for (i = 1; i < nskypatch; i++)
		rh_illum += parr[3*i+1] * rh_cos(i) * rh_dom[i];

	return rh_illum;
}

/* Calculate earth orbit eccentricity correction factor */

/* Reference:	Sen, Z. 2008. Solar Energy Fundamental and Modeling  */
/*				Techniques. Springer, p. 72. */

double CalcEccentricity()
{
	double day_angle;	/* Day angle (radians) */
	double E0;			/* Eccentricity */

	/* Calculate day angle */
	day_angle  = (julian_date - 1.0) * (2.0 * PI / 365.0);

	/* Calculate eccentricity */
	E0 = 1.00011 + 0.034221 * cos(day_angle) + 0.00128 * sin(day_angle)
			+ 0.000719 * cos(2.0 * day_angle) + 0.000077 * sin(2.0 *
			day_angle);

	return E0;
}

/* Calculate atmospheric precipitable water content */

/* Reference:	Perez, R., P. Ineichen, R. Seals, J. Michalsky, and R. */
/*				Stewart. 1990. ìModeling Daylight Availability and */
/*				Irradiance Components from Direct and Global */
/*				Irradiance,î Solar Energy 44(5):271-289, Eqn. 3. */

/* Note:	The default surface dew point temperature is 11 deg. C */
/*			(52 deg. F). Typical values are: */

/*			Celsius 	Fahrenheit	 	Human Perception */
/*			> 24 		> 75 			Extremely uncomfortable */
/*			21 - 24 	70 - 74 		Very humid */
/*			18 - 21		65 - 69		 	Somewhat uncomfortable */
/*			16 - 18 	60 - 64 		OK for most people */
/*			13 - 16 	55 - 59		 	Comfortable */
/*			10 - 12 	50 - 54		 	Very comfortable */
/*			< 10 		< 49		 	A bit dry for some */

double CalcPrecipWater( double dpt )
{ return exp(0.07 * dpt - 0.075); }

/* Calculate relative air mass */

/* Reference:	Kasten, F. 1966. "A New Table and Approximation Formula */
/*				for the Relative Optical Air Mass," Arch. Meteorol. */
/*				Geophys. Bioklimataol. Ser. B14, pp. 206-233. */

/* Note:		More sophisticated relative air mass models are */
/*				available, but they differ significantly only for */
/*				sun zenith angles greater than 80 degrees. */

double CalcAirMass()
{
	return (1.0 / (cos(sun_zenith) + 0.15 * pow(93.885 -
			RadToDeg(sun_zenith), -1.253)));
}

/* Calculate Perez All-Weather sky patch luminances (modified by GW) */

/* NOTE: The sky patches centers are determined in accordance with the */
/*       BRE-IDMP sky luminance measurement procedures. (See for example */
/*       Mardaljevic, J. 2001. "The BRE-IDMP Dataset: A New Benchmark */
/*       for the Validation of Illuminance Prediction Techniques," */
/*       Lighting Research & Technology 33(2):117-136.) */

void CalcSkyPatchLumin( float *parr )
{
	int i;
	double aas;				/* Sun-sky point azimuthal angle */
	double sspa;			/* Sun-sky point angle */
	double zsa;				/* Zenithal sun angle */

	for (i = 1; i < nskypatch; i++)
	{
		/* Calculate sun-sky point azimuthal angle */
		aas = fabs(rh_pazi[i] - azimuth);

		/* Calculate zenithal sun angle */
		zsa = PI * 0.5 - rh_palt[i];

		/* Calculate sun-sky point angle (Equation 8-20) */
		sspa = acos(cos(sun_zenith) * cos(zsa) + sin(sun_zenith) *
				sin(zsa) * cos(aas));

		/* Calculate patch luminance */
		parr[3*i] = CalcRelLuminance(sspa, zsa);
		if (parr[3*i] < 0) parr[3*i] = 0;
		parr[3*i+2] = parr[3*i+1] = parr[3*i];
	}
}
