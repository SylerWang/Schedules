#ifndef lint
static const char RCSid[] = "$Id: pabopto2bsdf.c,v 2.27 2014/10/03 21:57:06 greg Exp $";
#endif
/*
 * Load measured BSDF data in PAB-Opto format.
 *
 *	G. Ward
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "platform.h"
#include "bsdfrep.h"
#include "resolu.h"
				/* global argv[0] */
char			*progname;

typedef struct {
	const char	*fname;		/* input file path */
	double		theta, phi;	/* input angles (degrees) */
	double		up_phi;		/* azimuth for "up" direction */
	int		igp[2];		/* input grid position */
	int		isDSF;		/* data is DSF (rather than BSDF)? */
	long		dstart;		/* data start offset in file */
} PGINPUT;

PGINPUT		*inpfile;	/* input files sorted by incidence */
int		ninpfiles;	/* number of input files */

/* Compare incident angles */
static int
cmp_indir(const void *p1, const void *p2)
{
	const PGINPUT	*inp1 = (const PGINPUT *)p1;
	const PGINPUT	*inp2 = (const PGINPUT *)p2;
	int		ydif = inp1->igp[1] - inp2->igp[1];
	
	if (ydif)
		return(ydif);

	return(inp1->igp[0] - inp2->igp[0]);
}

/* Prepare a PAB-Opto input file by reading its header */
static int
init_pabopto_inp(const int i, const char *fname)
{
	FILE	*fp = fopen(fname, "r");
	FVECT	dv;
	char	buf[2048];
	int	c;
	
	if (fp == NULL) {
		fputs(fname, stderr);
		fputs(": cannot open\n", stderr);
		return(0);
	}
	inpfile[i].fname = fname;
	inpfile[i].isDSF = -1;
	inpfile[i].up_phi = 0;
	inpfile[i].theta = inpfile[i].phi = -10001.;
				/* read header information */
	while ((c = getc(fp)) == '#' || c == EOF) {
		char	typ[32];
		if (fgets(buf, sizeof(buf), fp) == NULL) {
			fputs(fname, stderr);
			fputs(": unexpected EOF\n", stderr);
			fclose(fp);
			return(0);
		}
		if (sscanf(buf, "sample_name \"%[^\"]\"", bsdf_name) == 1)
			continue;
		if (sscanf(buf, "format: theta phi %s", typ) == 1) {
			if (!strcasecmp(typ, "DSF")) {
				inpfile[i].isDSF = 1;
				continue;
			}
			if (!strcasecmp(typ, "BSDF") ||
					!strcasecmp(typ, "BRDF") ||
					!strcasecmp(typ, "BTDF")) {
				inpfile[i].isDSF = 0;
				continue;
			}
		}
		if (sscanf(buf, "upphi %lf", &inpfile[i].up_phi) == 1)
			continue;
		if (sscanf(buf, "intheta %lf", &inpfile[i].theta) == 1)
			continue;
		if (sscanf(buf, "inphi %lf", &inpfile[i].phi) == 1)
			continue;
		if (sscanf(buf, "incident_angle %lf %lf",
				&inpfile[i].theta, &inpfile[i].phi) == 2)
			continue;
	}
	inpfile[i].dstart = ftell(fp) - 1;
	fclose(fp);
	if (inpfile[i].isDSF < 0) {
		fputs(fname, stderr);
		fputs(": unknown format\n", stderr);
		return(0);
	}
	if ((inpfile[i].theta < -10000.) | (inpfile[i].phi < -10000.)) {
		fputs(fname, stderr);
		fputs(": unknown incident angle\n", stderr);
		return(0);
	}
				/* convert to Y-up orientation */
	inpfile[i].phi += 90.-inpfile[i].up_phi;
				/* convert angle to grid position */
	dv[2] = sin(M_PI/180.*inpfile[i].theta);
	dv[0] = cos(M_PI/180.*inpfile[i].phi)*dv[2];
	dv[1] = sin(M_PI/180.*inpfile[i].phi)*dv[2];
	dv[2] = sqrt(1. - dv[2]*dv[2]);
	pos_from_vec(inpfile[i].igp, dv);
	return(1);
}

/* Load a set of measurements corresponding to a particular incident angle */
static int
add_pabopto_inp(const int i)
{
	FILE		*fp = fopen(inpfile[i].fname, "r");
	double		theta_out, phi_out, val;
	int		n, c;
	
	if (fp == NULL || fseek(fp, inpfile[i].dstart, 0) == EOF) {
		fputs(inpfile[i].fname, stderr);
		fputs(": cannot open\n", stderr);
		return(0);
	}
					/* prepare input grid */
	if (!i || cmp_indir(&inpfile[i-1], &inpfile[i])) {
		if (i)			/* process previous incidence */
			make_rbfrep();
#ifdef DEBUG
		fprintf(stderr, "New incident (theta,phi)=(%.1f,%.1f)\n",
					inpfile[i].theta, inpfile[i].phi);
#endif
		new_bsdf_data(inpfile[i].theta, inpfile[i].phi);
	}
#ifdef DEBUG
	fprintf(stderr, "Loading measurements from '%s'...\n", inpfile[i].fname);
#endif
					/* read scattering data */
	while (fscanf(fp, "%lf %lf %lf\n", &theta_out, &phi_out, &val) == 3)
		add_bsdf_data(theta_out, phi_out+90.-inpfile[i].up_phi,
				val, inpfile[i].isDSF);
	n = 0;
	while ((c = getc(fp)) != EOF)
		n += !isspace(c);
	if (n) 
		fprintf(stderr,
			"%s: warning: %d unexpected characters past EOD\n",
				inpfile[i].fname, n);
	fclose(fp);
	return(1);
}

#ifndef TEST_MAIN
/* Read in PAB-Opto BSDF files and output RBF interpolant */
int
main(int argc, char *argv[])
{
	extern int	nprocs;
	int		i;
						/* start header */
	SET_FILE_BINARY(stdout);
	newheader("RADIANCE", stdout);
	printargs(argc, argv, stdout);
	fputnow(stdout);
	progname = argv[0];			/* get options */
	while (argc > 2 && argv[1][0] == '-') {
		switch (argv[1][1]) {
		case 'n':
			nprocs = atoi(argv[2]);
			break;
		default:
			goto userr;
		}
		argv += 2; argc -= 2;
	}
						/* initialize & sort inputs */
	ninpfiles = argc - 1;
	if (ninpfiles < 2)
		goto userr;
	inpfile = (PGINPUT *)malloc(sizeof(PGINPUT)*ninpfiles);
	if (inpfile == NULL)
		return(1);
	for (i = 0; i < ninpfiles; i++)
		if (!init_pabopto_inp(i, argv[i+1]))
			return(1);
	qsort(inpfile, ninpfiles, sizeof(PGINPUT), cmp_indir);
						/* compile measurements */
	for (i = 0; i < ninpfiles; i++)
		if (!add_pabopto_inp(i))
			return(1);
	make_rbfrep();				/* process last data set */
	build_mesh();				/* create interpolation */
	save_bsdf_rep(stdout);			/* write it out */
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-n nproc] meas1.dat meas2.dat .. > bsdf.sir\n",
					progname);
	return(1);
}
#else
/* Test main produces a Radiance model from the given input file */
int
main(int argc, char *argv[])
{
	PGINPUT	pginp;
	char	buf[128];
	FILE	*pfp;
	double	bsdf, min_log;
	FVECT	dir;
	int	i, j, n;

	progname = argv[0];
	if (argc != 2) {
		fprintf(stderr, "Usage: %s input.dat > output.rad\n", progname);
		return(1);
	}
	ninpfiles = 1;
	inpfile = &pginp;
	if (!init_pabopto_inp(0, argv[1]) || !add_pabopto_inp(0))
		return(1);
						/* reduce data set */
	if (make_rbfrep() == NULL) {
		fprintf(stderr, "%s: nothing to plot!\n", progname);
		exit(1);
	}
#ifdef DEBUG
	fprintf(stderr, "Minimum BSDF = %.4f\n", bsdf_min);
#endif
	min_log = log(bsdf_min*.5 + 1e-5);
#if 1						/* produce spheres at meas. */
	puts("void plastic yellow\n0\n0\n5 .6 .4 .01 .04 .08\n");
	n = 0;
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++)
		if (dsf_grid[i][j].sum.n > 0) {
			ovec_from_pos(dir, i, j);
			bsdf = dsf_grid[i][j].sum.v /
			   ((double)dsf_grid[i][j].sum.n*output_orient*dir[2]);
			if (bsdf <= bsdf_min*.6)
				continue;
			bsdf = log(bsdf + 1e-5) - min_log;
			ovec_from_pos(dir, i, j);
			printf("yellow sphere s%04d\n0\n0\n", ++n);
			printf("4 %.6g %.6g %.6g %.6g\n\n",
					dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf,
					.007*bsdf);
		}
#endif
#if 1						/* spheres at RBF peaks */
	puts("void plastic red\n0\n0\n5 .8 .01 .01 .04 .08\n");
	for (n = 0; n < dsf_list->nrbf; n++) {
		RBFVAL	*rbf = &dsf_list->rbfa[n];
		ovec_from_pos(dir, rbf->gx, rbf->gy);
		bsdf = eval_rbfrep(dsf_list, dir);
		bsdf = log(bsdf + 1e-5) - min_log;
		printf("red sphere p%04d\n0\n0\n", ++n);
		printf("4 %.6g %.6g %.6g %.6g\n\n",
				dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf,
				.011*bsdf);
	}
#endif
#if 1						/* output continuous surface */
	puts("void trans tgreen\n0\n0\n7 .7 1 .7 .04 .04 .9 1\n");
	fflush(stdout);
	sprintf(buf, "gensurf tgreen bsdf - - - %d %d", GRIDRES-1, GRIDRES-1);
	pfp = popen(buf, "w");
	if (pfp == NULL) {
		fprintf(stderr, "%s: cannot open '| %s'\n", progname, buf);
		return(1);
	}
	for (i = 0; i < GRIDRES; i++)
	    for (j = 0; j < GRIDRES; j++) {
		ovec_from_pos(dir, i, j);
		bsdf = eval_rbfrep(dsf_list, dir);
		bsdf = log(bsdf + 1e-5) - min_log;
		fprintf(pfp, "%.8e %.8e %.8e\n",
				dir[0]*bsdf, dir[1]*bsdf, dir[2]*bsdf);
	    }
	if (pclose(pfp) != 0)
		return(1);
#endif
	return(0);
}
#endif
