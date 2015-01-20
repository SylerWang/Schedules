#ifndef lint
static const char RCSid[] = "$Id: hemi.c,v 2.5 2013/11/26 17:33:55 greg Exp $";
#endif
/*
 *  Plot 3-D BSDF output based on scattering interpolant or XML representation
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "rtprocess.h"
#include "bsdfrep.h"

const float	colarr[6][3] = {
		.7, 1., .7,
		1., .7, .7,
		.7, .7, 1.,
		1., .5, 1.,
		1., 1., .5,
		.5, 1., 1.
	};

char	*progname;

/* Produce a Radiance model plotting the indicated incident direction(s) */
int
main(int argc, char *argv[])
{
	int	showPeaks = 0;
	int	doTrans = 1;
	int	inpXML = -1;
	RBFNODE	*rbf = NULL;
	FILE	*fp;
	char	buf[128];
	SDData	myBSDF;
	double	bsdf, min_log;
	FVECT	idir, odir;
	int	i, j, n, r, x, y, xp, yp;
	float rho;
						/* check arguments */
	progname = argv[0];
	if (argc > 1 && !strcmp(argv[1], "-p")) {
		++showPeaks;
		++argv; --argc;
	}
	if (argc > 1 && !strcmp(argv[1], "-r")) {
		r=atoi(argv[2]);
		++argv; ++argv; --argc; --argc;
	}
	if (argc >= 1 && (n = strlen(argv[1])-4) > 0) {
		if (!strcasecmp(argv[1]+n, ".xml"))
			inpXML = 1;
		else if (!strcasecmp(argv[1]+n, ".sir"))
			inpXML = 0;
	}
	if (inpXML < 0) {
		fprintf(stderr, "Usage: %s [-p] bsdf.sir theta1 phi1 .. > output.rad\n", progname);
		fprintf(stderr, "   Or: %s [-t] bsdf.xml theta1 phi1 .. > output.rad\n", progname);
		return(1);
	}
						/* load input */
	if (inpXML) {
		SDclearBSDF(&myBSDF, argv[1]);
		if (SDreportError(SDloadFile(&myBSDF, argv[1]), stderr))
			return(1);
		bsdf_min = 1./M_PI;
		if (myBSDF.rf != NULL && myBSDF.rLambFront.cieY < bsdf_min*M_PI)
			bsdf_min = myBSDF.rLambFront.cieY/M_PI;
		if (myBSDF.rb != NULL && myBSDF.rLambBack.cieY < bsdf_min*M_PI)
			bsdf_min = myBSDF.rLambBack.cieY/M_PI;
		if ((myBSDF.tf != NULL) | (myBSDF.tb != NULL) &&
				myBSDF.tLamb.cieY < bsdf_min*M_PI)
			bsdf_min = myBSDF.tLamb.cieY/M_PI;
		if (doTrans && (myBSDF.tf == NULL) & (myBSDF.tb == NULL)) {
			fprintf(stderr, "%s: no transmitted component in '%s'\n",
					progname, argv[1]);
			return(1);
		}
	} else {
		fp = fopen(argv[1], "rb");
		if (fp == NULL) {
			fprintf(stderr, "%s: cannot open BSDF interpolant '%s'\n",
					progname, argv[1]);
			return(1);
		}
		if (!load_bsdf_rep(fp))
			return(1);
		fclose(fp);
	}
#ifdef DEBUG
	//fprintf(stderr, "Minimum BSDF set to %.4f\n", bsdf_min);
#endif
	min_log = log(bsdf_min*.5);
						/* output BSDF rep. */
	printf("#?RADIANCE\n\n-Y %d +X %d\n",r+1,2*r+1);
	for (yp = r; yp>=0; yp--) {
		for(xp=-r; xp<=r; xp++){
		
			int  t, e;
			float d;

			x=xp;
			y=yp;
			rho=sqrt(x*x+y*y);
			

			input_orient = 1;
			output_orient = -input_orient;
		
			float theta, phi, azim;
			theta=90-rho*90/r;
			if(y){
				azim=atan(x/y)*180/3.141592654;
			}else{
				azim=0;
			}
		
			phi=90-azim;
		
			idir[2] = sin((M_PI/180.)*theta);
			idir[0] = idir[2] * cos((M_PI/180.)*phi);
			idir[1] = idir[2] * sin((M_PI/180.)*phi);
			idir[2] = input_orient * sqrt(1. - idir[2]*idir[2]);


			rbf = advect_rbf(idir, 15000);
			if(rho>r/1.){
				d=1;
			}else{
				if(abs(rho-r)<2){
					d=0;
				}else{
					d=SDdirectHemi(idir, SDsampSp|SDsampDf |	(output_orient > 0 ?  SDsampR : SDsampT), &myBSDF);
				}
			}
			float tau=d;
		
			d=d*9.0;

			if (d <= 1e-32) {
				t = 0;
				e = 0;
			}else{
				d = frexp(d, &e) * 255.9999 / d;

				if (tau > 0.0)
					t = tau * d;
				else
					t = 0;
				
				e = e + 128;
			}


						
			printf("%u %u %u %u\n", t,t,t,e);


		} //Loop X		
	} //Loop Y
	
	return(0);
}
