/****************************************************************************

 SDSS-spec: Read in and analyse SDSS spectra

 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "SDSS-spec.h"
#include <stats.h>
#include <const.h>
#include <file.h>
#include <fit.h>
#include <memory.h>
#include <error.h>

/* Global declarations */
char *progname = NULL;
plotenv plenv;
/* Default settings */
setblock set =
		{ .medcont = MEDCONT, .ncontsig = NCONSIG, .forest = 0,
				.keep = "\0", .pspec = 0, .outfile = "\0", .supress = 0,
				.nocomp = 0, .restframe = 0, .nforest = 1, .noext = 0,
				.pfilter = 0, .debug = 0 };

/****************************************************************************
 * Print the usage message
 ****************************************************************************/

void usage(void) {

	fprintf(stderr, "\n%s: Read in and analyse SDSS spectra\n", progname);

	fprintf(stderr, "By Michael Murphy (http://www.ast.cam.ac.uk/~mim)\n\
\n");
	fprintf(stderr,
			"With QSO absorption and D_a calculations additions\n\
By Jon Ouellet (4th May 2008)\n\nVersion: %4.3lf\n",
			VERSION);

	fprintf(stderr, "\nUsage: %s [OPTIONS] FILE\n", progname);

	fprintf(stderr,
			"\nOptions:\n\
 -p, --plot             =    : Plot each spectrum.\n\
 -n, --nforest [INT]    =    : Number of Lyman forests to compute.\n\
 -k, --keep [FILE]      =    : Keep a list of the QSOs used in this Composite in FILE\n\
 -o, --outfile [FILE]   =    : An output file [stdout]\n\
 -S, --supress          =    : Supress Output\n\
 -s, --silent           =    : Silence\n\
 -N, --null             =    : Do not compute a Composite Spectrum\n\
     --medcont [INT]    =    : Specify size of array for median fitting, [%i].\n\
     --ncontsig [FLOAT] =    : Rejection threshhold for cont. fit, [%3.2f].\n\
     --forest [INT]     =    : Specify which forest to analyze, [LYA].\n\
     --restspec [FILE,HDU]=  : Specify restframe spectra rather than powerlaw\n\
 -d, --debug            =    : Activate debugging output\n\
	 --metalcorr		=    : Activate metal corrections \n\
	 --noext            =	 : Deactivate dust corrections\n\
 -h, --help             =    : Print this message.\n\n",
			MEDCONT, NCONSIG);
	exit(3);
}

/****************************************************************************
 * The main program

 ****************************************************************************/

int main(int argc, char *argv[]) {

	/* DECLARATIONS and INITIALIZAIONS */
	/* Number of spectra read in from infile */
	int nspec = 0, speccount = 0;
	/* ID of pgplot device */
	int pg_id = 0;
	/* Counting Integers */
	int i = 0, j = 0, step = 0;
	// int miss = 0;
	/* Flags */
	int pgset = 0, flag = 0;
	// int moveback = 0;
	time_t t;
	/* Files for outputs */
	FILE *keepfile = NULL;
	FILE *files_used = NULL;
	FILE *files_not_used = NULL;
//  FILE     *status=fopen("/dev/null","w");
	files_used = fopen("files_used.txt", "w");
	files_not_used = fopen("files_not_used.txt", "w");
	/* Array of Spectra */
	spectrum *spec = NULL;
	int alpha_wrong_count = 0;


	/* Composite Spectrum Block */
	compspec *cspec = NULL;
	/* RestSpectrum If Needed */
	restspec *rspec = NULL;

	printf("Defined variables, now away we go....\n");
	/* Define the program name from the command line input */
	progname =
			((progname = strrchr(argv[0], '/')) == NULL ) ?
					argv[0] : progname + 1;
	/*Define start time */
	printf("Program name is %s\n", progname);
	t = time(NULL );
	/* Must be at least one argument */
	if (argc == 1)
		usage();

	/* Scan command line for options */
	while (++i < argc) {
		if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
			usage();
		else if (!strcmp(argv[i], "--ncontsig")) {
			i++;
			set.ncontsig = atof(argv[i]);
		} else if (!strcmp(argv[i], "--medcont")) {
			i++;
			set.medcont = ((j = atoi(argv[i])) % 2 == 1 ? j : j + 1);
		} else if (!strcmp(argv[i], "--noext"))
			set.noext = 1;
		else if (!strcmp(argv[i], "--metalcorr")){
			set.metal_corr = 1;
		}
		else if (!strcmp(argv[i], "--forest"))
			set.forest = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--outfile") || !strcmp(argv[i], "-o")) {
			if (!access(argv[++i], W_OK) || access(argv[i], F_OK)) {
				if (strlen(argv[i]) <= NAMELEN)
					strcpy(set.outfile, argv[i]);
				else
					errormsg("Output File Name too long: %s", argv[i]);
			} else
				errormsg("Outfile not write accessable");
		} else if (!strcmp(argv[i], "-N") || !strcmp(argv[i], "--null")) {
			set.nocomp = 1;
		} else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--plot")) {
			set.pspec = 1;
//      status=fopen("/dev/null","w"); 
		} else if (!strcmp(argv[i], "--pfilter"))
			set.pfilter = 1;
		else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
//      status=stderr;
		} else if (!strcmp(argv[i], "-S") || !strcmp(argv[i], "--supress")) {
			set.supress = 1;
		} else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--debug")) {
			set.debug = 1;
		} else if (!strcmp(argv[i], "-Y") || !strcmp(argv[i], "--synflux"))
			set.synflux = 1;
		else if (!strcmp(argv[i], "--restspec"))
			strcpy(set.restspec, argv[++i]);
		else if (!strncmp(argv[i], "--restspec=", 10))
			strcpy(set.restspec, strchr(argv[i], '=') + 1);
		else if (!strcmp(argv[i], "-k")) {
			if (!access(argv[++i], W_OK) || access(argv[i], F_OK)) {
				if (strlen(argv[i]) <= NAMELEN)
					strcpy(set.keep, argv[i]);
				else
					errormsg("Quasar Record Name too long: %s", argv[i]);
			}
		} else if (!strcmp(argv[i], "--nforest") || !strcmp(argv[i], "-n"))
			set.nforest = atoi(argv[++i]);
		else if (!access(argv[i], R_OK)) {
			if (strlen(argv[i]) <= VLNGSTRLEN)
				strcpy(set.infile, argv[i]);
			else
				errormsg("Input File Name Too Long: %s", argv[i]);
		} else
			usage();
	}
	printf("Input options processed\n");

	/* Make sure an input was specified */
	while (!strncmp(set.infile, "\0", 1))
		get_input("Please specify an input file", "%s", set.infile);
	/* Make sure an output file was specified */
	printf("Input file specified \n");
	while (!strncmp(set.outfile, "\0", 1) && !set.pspec)
		get_input("Please specify output file?", "%s", set.outfile);


	// put output folder in front of output file
	char temp[NAMELEN];
	strcpy(temp, "./OutFits/");
	strcat(temp, set.outfile);
	strcpy(set.outfile, temp);
	printf("%s \n\n", set.outfile);

	printf("Reading input file ...\n");
	/* Read input file */
	if (!SDSS_rinputfile(&spec, &nspec, set.infile))
		errormsg("Unknown error returned from SDSS_rinputfile()");
	/* Initialize Composite Arrays */
	if (!SDSS_compinit(&cspec, &set))
		errormsg("Could not initialize composite spectra arrays");

	/* If we are using Restframe spectra initialize them */
	if (strcmp(set.restspec, "\0")) {
		if (!SDSS_restinit(&rspec, &set))
			errormsg("Could not initialize rest spectra arrays");
	}
	//fprintf(status,"Operating Spectra..\nCompletion: 00%%");
	printf("Operating spectra...\n");

	/*===============*/
	/* START OF LOOP */
	/*===============*/

	/* Read in and operate in individual spectra, one at a time to save memory */
	printf("Reading in spectra...\n");
	for (i = 0; i < nspec; i++) {
		/* Just the completion counter */
		//if( rint(100.0*i/nspec) > rint(100.0*(i-1)/nspec))
		//printf("\b\b\b%2i%%",(int) 100*i/nspec);
		printf("Working on spectrum number %d...\n", i + 1);

		/* Open and Read The Fit File */
		if (spec[i].PlateSpec == 0) {
			if (!SDSS_rspSpec(&(spec[i]), &set)) 	// From an spSpec file
				errormsg("Unknown error returned from SDSS_rspSpec()");
		} else if (spec[i].PlateSpec == 1) {
			spec[i].fiberID = i+1;
			printf("Number nspec: %d\n", nspec);
			if (!SDSS_rplate(&(spec[i]), &set)) 	// From an spPlate file
				errormsg("Unknown error returned from SDSS_rplate()");
		} else if (spec[i].PlateSpec == 2) {
			if (!SDSS_rspec(&(spec[i]), &set))      // From an spec File
				errormsg("Unknown error returned from SDSS_rspec()");
		} else
			errormsg("PlateSpec error at line %i", i + j + 1);

		/* Enter Any specifications on the quasars */
		// 2.2 to 5.3
		if (spec[i].zem >= 2.2 && spec[i].zem <= 5.3) {

			if (!set.noext) {
				/* Correct For Galactic Extinction */
				printf("doing Gal_corr\n");
				if (!SDSS_Ebv(&(spec[i])))
					errormsg("Error returned from SDSS_Ebv()");

				if (!SDSS_Gal_extinct_correct(&(spec[i])))
					errormsg(
							"Unknown error returned from SDSS_Gal_extinct_correct()");
			}

			if (!strcmp(set.restspec, "\0")) {
				/* Build The Median Array */
				if (!(medianrun(spec[i].fl, spec[i].co, spec[i].st, spec[i].np,
						set.medcont)))
					errormsg(
							"SDSS-spec(): Error when constructing med. filt. cont of size %d\n\
                   \tfor %s",
							set.medcont, spec[i].abfile);

				/* Measure Spectroscopic Colors */
				if (!SDSS_colors(&(spec[i]), set.pfilter))
					errormsg(
							"SDSS-spec(): Error returned from spectroscopic colors");

				/* Fit a Power Law */
				if (!SDSS_red_powerlaw(&(spec[i])))
					nferrormsg(
							"Problem fitting red powerlaw to spectrum in file\n\t%s",
							spec[i].file);
			} else {
				if (!SDSS_restspec(&(spec[i]), rspec))
					errormsg("Problem fitting rest spectrum");
			}

			/**** Restrict Spectral Index Range ****/
			if ((spec[i].alpha < TOPCUT && spec[i].alpha > LOWCUT)
					&& !set.nocomp) {

				/* Build the Composite Spectrum */
				// Save variables of cspec locations
				//double *cspec_wl = NULL, *cspec_fl = NULL, *cspec_er = NULL, *cspec_sum = NULL, *cspec_sum2 = NULL, *cspec_wsum = NULL;
				//double *cspec_wsum2 = NULL;
				if ((flag = SDSS_compspec(cspec, &(spec[i]), &set)) == 0) // PROBLEM: incompatible wavelength range
					errormsg("Could not run comp spec at file %i", i);

				/* If the quasar contributes add it to the list */
				if (flag != 9) {
					speccount++;
					spec[i].usespec = 1;
					fprintf(files_used, "%s\n", spec[i].abfile);
				}
			}
			else{
				fprintf(files_not_used, "alpha\t%f\t%s\n", spec[i].alpha, spec[i].file);
			}
		}


		/* Initialize plotting if we are plotting */
		if (!pgset && set.pspec) {
			cpgqid(&pg_id);
			SDSS_pgenv_init();
			pg_open(&plenv, "?\0", progname, 1);
			pgset++;
		}
		if (set.pspec) {
			if ((step = SDSS_plot_spec(&(spec[i]))) == -9)
				errormsg("Problem plotting spectrum in file\n\t%s",
						spec[i].file);
			else if (step == -1 && i > 0 && i != 1)
				i -= 2;
			else if (step == -1 && i == 1){
				i -= 1;
			}
		}
		if (set.pspec && i == nspec - 1)
			//if (set.pspec > 400)						// employed for incompatibility problem - solved (for now)
			cpgclos();

		/** Free up major use of memory **/

		if(!set.pspec){
			free(spec[i].wl);
			free(spec[i].fl);
			free(spec[i].er);
			free(spec[i].co);
			free(spec[i].pc);
			free(spec[i].sn);
			free(spec[i].st);
			//free(spec[i].ps);
			//free(spec[i].eps);
		}
		if (step == 10)
			i = nspec - 1;

		printf("Done with spectrum number %d.\n", i + 1);
		if (spec[i].alpha == -999){
			alpha_wrong_count++;
		}
	}
	/*===============*/
	/*  END OF LOOP  */
	/*===============*/

	/* Do some number crunching */
	if (!SDSS_statistics(cspec, spec, nspec, &set))
		errormsg("Problem computing statistics");

	/* Print the data to a FITS file */
	if (strcmp(set.outfile, "\0")) {
//    fprintf(status,"\nBuilding FITS file..\n");
		printf("\nBuilding FITS file..\n");
		if (!SDSS_compfits(cspec, spec, nspec, speccount, &set))
			errormsg("Problem making fits file");
	}

	/* Print the used Spectra to a File */
	if (strcmp(set.keep, "\0")) {
//    fprintf(status,"Printing Used Spectra..\n");
		printf("Used Spectra..\n");
		if ((keepfile = fopen(set.keep, "w")) == NULL )
			errormsg("Could not open %s", set.keep);
		for (i = 0; i < nspec; i++) {
			if (spec[i].usespec)
				fprintf(keepfile, "%s\n", spec[i].file);
		}
		fclose(keepfile);
	}

	/* Clean up */
//  fprintf(status,"Clearing Memory Space..\n");
	printf("Clearing Memory Space..\n");
	free(spec);
	free(cspec->wl);
	free(cspec->sum);
	free(cspec->sum2);
	free(cspec->wsum);
	free(cspec->wsum2);
	free(cspec->nhist);
	free(cspec->fl);
	free(cspec->er);
	free(cspec);
	SDSS_colors(NULL, 0);
//  fprintf(status,"Used %i/%i spectra.\n",speccount,nspec);
	printf("Used %i/%i spectra.\n", speccount, nspec);
	printf("There are %i objects with alpha -999", alpha_wrong_count);
//  fprintf(status,"%s Completed in %.0f Seconds.\n",
	printf("%s Completed in %.0f Seconds.\n", progname,
			difftime(time(NULL ), t));
//  fclose(status);
	fclose(files_not_used);
	return 1;

}
