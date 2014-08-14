/****************************************************************************

 SDSS_rinputfile: Decide on & read type of input file.  Assign nspec & list.
 Allocate memory for nspec spectra.

 ****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "SDSS-spec.h"
#include "file.h"
#include "const.h"
#include "error.h"
#include "fitsio.h"


// Function which determines if file is of type spSpec, spPlate
int SDSS_filetype(char *infile, int *nspec){
	// return value:
	// 0 = spSpec
	// 1 = spPlate
	// 2 = spec
	int val = 0;
	int hdunum_exp = 0;

	fitsfile *infits;
	int status = 0, hdunum = 0;
	//printf("Name of File: %s\n", infile);
	fits_open_file(&infits, infile, READONLY, &status);
	fits_get_num_hdus(infits, &hdunum, &status);
	//printf("HDUNumber: %i\n\n", hdunum);


	// Number of HDUs not necessarily 9. In some cases might be different
	if( strstr(infile, "spSpec") != NULL ){
		val = 0;
		hdunum_exp = 7;
		//if(hdunum != hdunum_exp) warnmsg("Number of HDUs is not %d as expected from spSpec", hdunum_exp);
	}
	if( strstr(infile, "spPlate") != NULL ){
		val = 1;
		hdunum_exp = 9;
		//if(hdunum != hdunum_exp) warnmsg("Number of HDUs is not %d as expected from spPlate", hdunum_exp);
		fits_read_key(infits, TINT, "NAXIS2", nspec, NULL, &status);
	}
	if( strstr(infile, "spec-") != NULL ){
		val = 2;
		hdunum_exp = 14;
		//if(hdunum != hdunum_exp) warnmsg("Number of HDUs is not %d as expected from spec", hdunum_exp);
	}

	// report Error. TODO: Decide if want to use errormsg instead.
	fits_close_file(infits, &status);
	fits_report_error(stderr, status);

	return val;
}


int SDSS_rinputfile(spectrum **spec, int *nspec, char *infile) {

	int i = 0, j = 0, flag = 0;
	int dum1, dum2, dum3;
	char buffer[HUGESTRLEN] = "\0";
	char *cptr, *dir = NULL;
	int PlateSpec_tmp = 0;
	int nspec_tmp = 0;
	FILE *data_file = NULL;

	/* Get Directory Info */
	dir = getenv("SPECTRO_DATA");

	/* Open Data File */
	if ((data_file = fopen(infile, "r")) == NULL )
		errormsg("Could not open file %s\n", infile);

	/* Read First Line */
	if (fgets(buffer, HUGESTRLEN, data_file) == NULL ) {
		fclose(data_file);
		errormsg("SDSS_rinputfile: Problem reading file %s on line %d", infile,
				1);
	}
	if (!strncmp(buffer, "SIMPLE  =", 8)) {
		/* Single file specified that looks suspiciously like a FITS file */
		(*nspec) = 1;
		flag = 1;

		// Simply determine nspec by looking at number of rows in spPlate file, IF spPlate!
		// Additionally determine what kind of Fits file it is
		nspec_tmp = 0;
		// Determine the filetype and set nspec, if spPlate.
		PlateSpec_tmp = SDSS_filetype(infile, &nspec_tmp);
		if( PlateSpec_tmp == 1 ) *nspec = nspec_tmp;
	} else {
		printf("No FITS file\n");
		i = 1;
		while (fgets(buffer, HUGESTRLEN, data_file) != NULL ) {
			if (sscanf(buffer, "%i %i %i", &dum1, &dum2, &dum3) == 3)
				printf("%i %i %i\n", dum1, dum2, dum3);
			else if (sscanf(buffer, "%s", buffer) && !access(buffer, R_OK)) ;
			else
				errormsg("Could not open file/read mjd/plate/fiber on line %i",
						i);
			i++;
		}

		if (!feof(data_file)) {
			fclose(data_file);
			errormsg("SDSS_rinputfile: Problem reading file %s on line %d",
					infile, i + 1);
		} else
			(*nspec) = i;
	}
	rewind(data_file);

	/** Alocate memory for nspec spectra **/
	// TODO: Currently not implemented that if list file contains spPlate files, because in that case we would need N*nspec spectra,
	// where N is the number of spectra in one spSplate file. That means: spPlate currently only workable if single file.
	if (!((*spec) = (spectrum *) malloc((size_t) ((*nspec) * sizeof(spectrum)))))
		errormsg("SDSS_rinputfile: Could not allocate memory for header array\n\\t of size %d",	(*nspec));

	/** Initialisations **/
	for (i = 0; i < *nspec; i++) {
		(*spec)[i].fidx = i;
		(*spec)[i].tsf.primtarg = -1;
		(*spec)[i].beta = (*spec)[i].alpha = (*spec)[i].ealpha =
				(*spec)[i].delta = -999.0;
		(*spec)[i].usespec = 0;
		// smag array was not initialised. caused memory leak. Thinkable to put initialisation at point in code where used.
		for(j = 0; j < 5; j++){
			(*spec)[i].smag[j] = 0;
		}
	}

	if (flag) {
		for(i = 0; i < *nspec; i++){
			strcpy((*spec)[i].file, infile);
			cptr = ((cptr = strrchr((*spec)[i].file, '/')) == NULL ) ?
					(*spec)[i].file : cptr + 1;
			strcpy((*spec)[i].abfile, cptr);

			// Move down to end somewhere. After reading all files from list file still need to determine the kind of fits files within the list file
			// Set PlateSpec of all spectra accordingly
			(*spec)[i].PlateSpec = PlateSpec_tmp;
		}
	}

	else { // Comment for me: Checklist, index c
		for (i = 0; i < *nspec; i++) {
			fgets(buffer, HUGESTRLEN, data_file);
			strcpy((*spec)[i].file, buffer);
			if (sscanf(buffer, "%i\t%i\t%i", &dum1, &dum2, &dum3) == 3) {
				(*spec)[i].mjd = dum1;
				(*spec)[i].plate = dum2;
				(*spec)[i].fiberID = dum3;
				sprintf(buffer, "%s/%04i/spPlate-%04i-%05i.fits", dir, dum2,
						dum2, dum1);
				if (!access(buffer, R_OK)) {
					strcpy((*spec)[i].file, buffer);
					(*spec)[i].PlateSpec = 1;
				} else
					sprintf(buffer, "%s/%04i/spSpec-%05i-%04i-%03i.fit", dir,
							dum2, dum1, dum2, dum3);
				if (!access(buffer, R_OK))
					strcpy((*spec)[i].file, buffer);
				else
					errormsg("Can not find file to match: %05i %04i %03i", dum1,
							dum2, dum3);
			} else {
				if (sscanf(buffer, "%s %lf %lf %lf %lf %lf", (*spec)[i].file,
						&((*spec)[i].mag[0]), &((*spec)[i].mag[1]),
						&((*spec)[i].mag[2]), &((*spec)[i].mag[3]),
						&((*spec)[i].mag[4])) == 6)
					;
				else
					(*spec)[i].mag[0] = (*spec)[i].mag[1] = (*spec)[i].mag[2] =
							(*spec)[i].mag[3] = (*spec)[i].mag[4] = 0;
			}
			//strcpy((*spec)[i].file, infile);
			cptr = ((cptr = strrchr((*spec)[i].file, '/')) == NULL ) ?
					(*spec)[i].file : cptr + 1;
			strcpy((*spec)[i].abfile, cptr);

			// Put call to filetype() here and after that set s

			// Simply determine nspec by looking at number of rows in spPlate file, IF spPlate!
			// Additionally determine what kind of Fits file it is

			nspec_tmp = *nspec;
			// Determine the filetype and set nspec, if spPlate.
			PlateSpec_tmp = SDSS_filetype((*spec)[i].file, &nspec_tmp);
			//if( PlateSpec_tmp == 1 ) *nspec = nspec_tmp;

			// Move down to end somewhere. After reading all files from list file still need to determine the kind of fits files within the list file
			// Set PlateSpec of all spectra accordingly
			(*spec)[i].PlateSpec = PlateSpec_tmp;
		}
	}


	// Maybe here after all file names are determined and written to abfile, we can run over
	// all again and determine their individual properties? spSpec, spPlate or spec?


	//(*spec)[0].PlateSpec = 1;

	fclose(data_file);

	return 1;
}
