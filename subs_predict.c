#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "interface.h"
#include "subs_inoutput.h"
#include "subs_common_string.h"
#include "subs_fits.h"
#include "subs_asciifile.h"
#include "subs_memory.h"
#include "subs_lambert.h"
#include "subs_predict.h"

/******************************************************************************/
/* Fortran wrapper */
void DECLARE(fort_predict_thermal)
  (void  *  pNGal,
   float *  pGall,
   float *  pGalb,
   float *  pNu,
   char  *  pIPath,
   char  *  pResName,
   char  *  pUnitsName,
   void  *  pModelNum,
   void  *  pQInterp,
   void  *  pQNoloop,
   void  *  pQVerbose,
   float *  pOutput)
{
   int      iChar;
   int      modelNum;
   int      qInterp;
   int      qNoloop;
   int      qVerbose;
   long     iGal;
   long     nGal;
   float *  pTemp;

   /* Truncate the Fortran-passed strings with a null,
    * in case they are padded with spaces */
   for (iChar=0; iChar < 80; iChar++)
    if (pIPath[iChar] == ' ') pIPath[iChar] = '\0';
   for (iChar=0; iChar < 10; iChar++)
    if (pResName[iChar] == ' ') pResName[iChar] = '\0';
   for (iChar=0; iChar < 10; iChar++)
    if (pUnitsName[iChar] == ' ') pUnitsName[iChar] = '\0';

   /* Select the 4-byte words passed by a Fortran call */
   if (sizeof(short) == 4) {
      nGal = *((short *)pNGal);
      modelNum = *((short *)pModelNum);
      qInterp = *((short *)pQInterp);
      qNoloop = *((short *)pQNoloop);
      qVerbose = *((short *)pQVerbose);
   } else if (sizeof(int) == 4) {
      nGal = *((int *)pNGal);
      modelNum = *((int *)pModelNum);
      qInterp = *((int *)pQInterp);
      qNoloop = *((int *)pQNoloop);
      qVerbose = *((int *)pQVerbose);
   } else if (sizeof(long) == 4) {
      nGal = *((long *)pNGal);
      modelNum = *((long *)pModelNum);
      qInterp = *((long *)pQInterp);
      qNoloop = *((long *)pQNoloop);
      qVerbose = *((long *)pQVerbose);
   }
 
   pTemp = predict_thermal(nGal, pGall, pGalb, pNu, pIPath, pResName,
    pUnitsName, modelNum, qInterp, qNoloop, qVerbose);
 
   /* Copy results into Fortran-passed location for "pOutput",
    * assuming that memory has already been allocated */
   for (iGal=0; iGal < nGal; iGal++) pOutput[iGal] = pTemp[iGal];
}
 
/******************************************************************************/
/* Fortran wrapper */
void DECLARE(fort_predict_sync)
  (void  *  pNGal,
   float *  pGall,
   float *  pGalb,
   float *  pNu,
   char  *  pIPath,
   char  *  pUnitsName,
   void  *  pQInterp,
   void  *  pQNoloop,
   void  *  pQVerbose,
   float *  pOutput)
{
   int      iChar;
   int      qInterp;
   int      qNoloop;
   int      qVerbose;
   long     iGal;
   long     nGal;
   float *  pTemp;

   /* Truncate the Fortran-passed strings with a null,
    * in case they are padded with spaces */
   for (iChar=0; iChar < 80; iChar++)
    if (pIPath[iChar] == ' ') pIPath[iChar] = '\0';
   for (iChar=0; iChar < 10; iChar++)
    if (pUnitsName[iChar] == ' ') pUnitsName[iChar] = '\0';

   /* Select the 4-byte words passed by a Fortran call */
   if (sizeof(short) == 4) {
      nGal = *((short *)pNGal);
      qInterp = *((short *)pQInterp);
      qNoloop = *((short *)pQNoloop);
      qVerbose = *((short *)pQVerbose);
   } else if (sizeof(int) == 4) {
      nGal = *((int *)pNGal);
      qInterp = *((int *)pQInterp);
      qNoloop = *((int *)pQNoloop);
      qVerbose = *((int *)pQVerbose);
   } else if (sizeof(long) == 4) {
      nGal = *((long *)pNGal);
      qInterp = *((long *)pQInterp);
      qNoloop = *((long *)pQNoloop);
      qVerbose = *((long *)pQVerbose);
   }
 
   pTemp = predict_sync(nGal, pGall, pGalb, pNu, pIPath,
    pUnitsName, qInterp, qNoloop, qVerbose);
 
   /* Copy results into Fortran-passed location for "pOutput",
    * assuming that memory has already been allocated */
   for (iGal=0; iGal < nGal; iGal++) pOutput[iGal] = pTemp[iGal];
}

/******************************************************************************/
float * predict_thermal
  (long     nGal,
   float *  pGall,
   float *  pGalb,
   float *  pNu,
   char  *  pIPath,
   char  *  pResName,
   char  *  pUnitsName,
   int      modelNum,
   int      qInterp,
   int      qNoloop,
   int      qVerbose)
{
   int      ii;
   int      im;
   int      iz1 = -1; /* crash if Zindx not tabulated for alpha1 */
   int      iz2 = -1; /* crash if Zindx not tabulated for alpha2 */
   int      imap;
   int      iGal;
   float *  pI100;
   float *  pRmapval;
   float *  pInu;
   float    alpha1;
   float    alpha2;
   float    f1;
   float    q1q2;
   float    RfitA[6];
   float    lnR;
   float    lnRpow;
   float    T1;
   float    T2;
   float    lnT1;
   float    lnT2;
   float    tcoeff;
   const float nu100 = 2997.92458; /* Frequency in GHz for 100-microns */
   const float h_Pl = 6.6261e-27; /* cm^2 g s^-1 */
   const float k_B = 1.3806e-16;  /* erg K^-1 */

   /* Declarations for command-line keyword names */
   char pText_MJy[] = "MJy";
   char pText_microK[] = "microK";
   char pText_thermo[] = "thermo";

   /* Declarations for data file names */
   char     pFileN[MAX_FILE_NAME_LEN];
   char     pFileS[MAX_FILE_NAME_LEN];
   struct   mapParms {
      char *   pName;
      char *   pFile1;
      char *   pFile2;
   } ppMapAll[] = {
     { "D1024", "SFD_d100_1024_ngp.fits", "SFD_d100_1024_sgp.fits" },
     { "I1024", "SFD_i100_1024_ngp.fits", "SFD_i100_1024_sgp.fits" },
     { "I2048", "SFD_i100_2048_ngp.fits", "SFD_i100_2048_sgp.fits" },
     { "I4096", "SFD_i100_4096_ngp.fits", "SFD_i100_4096_sgp.fits" }
   };
   const int nmap = sizeof(ppMapAll) / sizeof(ppMapAll[0]);
   char * ppRmapFile[] =
     { "FINK_Rmap_ngp.fits" , "FINK_Rmap_sgp.fits" };

   /* Set model parameters */
   const float alpha1vec[] = {1.50, 1.70, 2.00, 2.20, 1.50, 2.00, 1.50, 1.67};
   const float alpha2vec[] = {0.00, 0.00, 0.00, 0.00, 2.60, 2.00, 2.60, 2.70};
   const float f1vec[]     = {1.00, 1.00, 1.00, 1.00, 0.25, 0.00261, 0.0309, 0.0363};
   const float q1q2vec[]   = {1.00, 1.00, 1.00, 1.00, 0.61, 2480.0, 11.2, 13.0};
   /* const int N_MODEL = sizeof(alpha1vec) / sizeof(alpha1vec[0]); */
 
   /* Rfita contains fit coefficients for T2_of_R */
   const float RfitAarr[][6] =
    {{2.9268E+00, 3.8419E-01, 5.0233E-02, 1.0852E-02, 3.0738E-03, 5.0595E-04},
     {2.8483E+00, 3.8044E-01, 4.6584E-02, 9.0938E-03, 2.7038E-03, 5.4664E-04},
     {2.7334E+00, 3.7537E-01, 4.1712E-02, 6.8839E-03, 2.0316E-03, 6.0311E-04},
     {2.6556E+00, 3.7377E-01, 3.9898E-02, 5.7662E-03, 1.4638E-03, 6.3723E-04},
     {2.9206E+00, 2.3254E-01, 2.3506E-02, 4.0781E-03, 1.0048E-03, 1.2004E-04},
     {2.9900E+00, 2.5041E-01, 2.9688E-02, 6.5641E-03, 1.5688E-03, 1.6542E-04},
     {2.8874E+00, 2.4172E-01, 2.9369E-02, 4.7867E-03, 9.7237E-04, 1.1410E-04},
     {2.8723E+00, 2.4071E-01, 2.9625E-02, 4.7196E-03, 9.3207E-04, 1.1099E-04} };
 
   /* Zeta integrals for alpha=[1.50, 1.67, 1.70, 2.00, 2.20, 2.60, 2.70]
      from equn (15) of Finkbeiner et al. */
   const float Zindx[] = {1.50, 1.67, 1.70, 2.00, 2.20, 2.60, 2.70};
   const float Zintegral[] = {5.3662E+01, 7.0562E+01, 7.4100E+01, 1.2208E+02,
    1.7194E+02, 3.4855E+02, 4.1770E+02};
   const int N_ZINDEX = sizeof(Zindx) / sizeof(Zindx[0]);

   /* Test that inputs are valid */
   if (nGal == 0 || pGall == NULL || pGalb == NULL || pNu == NULL) {
      printf("ERROR: Must specify coordinates and frequencies.\n");
      return NULL;
   }

   /* Select parameters for this model */
   alpha1 = alpha1vec[modelNum-1];
   alpha2 = alpha2vec[modelNum-1];
   f1 = f1vec[modelNum-1];
   q1q2 = q1q2vec[modelNum-1];
   for (ii=0; ii < 6; ii++) RfitA[ii] = RfitAarr[modelNum-1][ii];

   /* Determine the file names to use */
   for (imap=0; imap < nmap; imap++) {
      if (strcmp(pResName,ppMapAll[imap].pName) == 0) {
         sprintf(pFileN, "%s/%s", pIPath, ppMapAll[imap].pFile1);
         sprintf(pFileS, "%s/%s", pIPath, ppMapAll[imap].pFile2);
      }
   }

   /* Read the 100-micron map */
   pI100 = lambert_getval(pFileN, pFileS, nGal, pGall, pGalb,
    qInterp, qNoloop, qVerbose);

   /* Read the I100/240 ratio map */
   sprintf(pFileN, "%s/%s", pIPath, ppRmapFile[0]);
   sprintf(pFileS, "%s/%s", pIPath, ppRmapFile[1]);
   pRmapval = lambert_getval(pFileN, pFileS, nGal, pGall, pGalb,
    qInterp, qNoloop, qVerbose);

   /* Allocate memory for output array */
   pInu = ccvector_build_(nGal);

   if (modelNum <=4) {
      /* SINGLE-COMPONENT MODEL: Evaluate equn (1) from Finkbeiner et al */
      for (iGal=0; iGal < nGal; iGal++) {

         /* Compute ln(T1) from ln(Rmap) */
         lnR = log(pRmapval[iGal]);
         lnRpow = 1.0;
         lnT1 = RfitA[0];
         for (ii=1; ii < 6; ii++) {
            lnRpow *= lnR;
            lnT1 += RfitA[ii] * lnRpow;
         }
         T1 = exp(lnT1);

         pInu[iGal] =
          pI100[iGal] * pow(pNu[iGal]/nu100,alpha1) * planck(T1,pNu[iGal]) /
          ( planck(T1,nu100) * kfactor(alpha1,T1) );
      }
   } else {
      /* TWO-COMPONENT MODEL: Evaluate equn (6) from Finkbeiner et al */

      /* Find Zintegral index for the model values of "alpha" */
      for (im=0; im < N_ZINDEX; im++) {
         if (fabs(Zindx[im] - alpha1) < 1.e-4) iz1 = im;
         if (fabs(Zindx[im] - alpha2) < 1.e-4) iz2 = im;
      }
      tcoeff = pow( (Zintegral[iz2] / (q1q2*Zintegral[iz1]))
       * pow(h_Pl*nu100*1.e+9/k_B,alpha1-alpha2), 1./(4.+alpha1) );

      for (iGal=0; iGal < nGal; iGal++) {

         /* Compute ln(T2) from ln(Rmap) */
         lnR = log(pRmapval[iGal]);
         lnRpow = 1.0;
         lnT2 = RfitA[0];
         for (ii=1; ii < 6; ii++) {
            lnRpow *= lnR;
            lnT2 += RfitA[ii] * lnRpow;
         }
         T2 = exp(lnT2);

         /* Compute T1 as a function of T2; equn (13) of Finkbeiner et al. */
         T1 = tcoeff * pow( T2, ((4+alpha2)/(4+alpha1)) );

         pInu[iGal] = pI100[iGal] *
          ( f1 * q1q2 * pow(pNu[iGal]/nu100,alpha1) * planck(T1,pNu[iGal])
             + (1-f1) * pow(pNu[iGal]/nu100,alpha2) * planck(T2,pNu[iGal]) ) / 
          ( f1 * q1q2 * planck(T1,nu100) * kfactor(alpha1,T1)
          + (1-f1) * planck(T2,nu100) * kfactor(alpha2,T2) );
      }
   }

   /* Convert units */
   if (strcmp(pUnitsName,pText_MJy) == 0) {
      /* MJy/sr */
   } else if (strcmp(pUnitsName,pText_microK) == 0) {
      /* brightness temp micro-K */
      for (iGal=0; iGal < nGal; iGal++)
       pInu[iGal] *= fac_flux2temp(pNu[iGal]);
   } else if (strcmp(pUnitsName,pText_thermo) == 0) {
      /* thermodynamic micro-K */
      for (iGal=0; iGal < nGal; iGal++)
       pInu[iGal] *= fac_flux2temp(pNu[iGal]) * planckcorr(pNu[iGal]);
   } else {
      printf("ERROR: Invalid units name.\n");
      for (iGal=0; iGal < nGal; iGal++)
       pInu[iGal] = 0.0;
   }

   ccvector_free_(pI100);
   ccvector_free_(pRmapval);

   return pInu;
}

/******************************************************************************/
float * predict_sync
  (long     nGal,
   float *  pGall,
   float *  pGalb,
   float *  pNu,
   char  *  pIPath,
   char  *  pUnitsName,
   int      qInterp,
   int      qNoloop,
   int      qVerbose)
{
   int      iGal;
   float *  pAmap;
   float *  pBmap;
   float *  pInu;

   /* Declarations for command-line keyword names */
   char pText_MJy[] = "MJy";
   char pText_microK[] = "microK";
   char pText_thermo[] = "thermo";

   /* Declarations for data file names */
   char     pFileN[MAX_FILE_NAME_LEN];
   char     pFileS[MAX_FILE_NAME_LEN];
   char * ppHaslamFile[] =
     { "Haslam_clean_ngp.fits" , "Haslam_clean_sgp.fits" };
   char * ppBetaFile[] =
     { "Synch_Beta_ngp.fits" , "Synch_Beta_sgp.fits" };

   /* Test that inputs are valid */
   if (nGal == 0 || pGall == NULL || pGalb == NULL || pNu == NULL) {
      printf("ERROR: Must specify coordinates and frequencies.\n");
      return NULL;
   }

   /* Read the Haslam map */
   sprintf(pFileN, "%s/%s", pIPath, ppHaslamFile[0]);
   sprintf(pFileS, "%s/%s", pIPath, ppHaslamFile[1]);
   pAmap = lambert_getval(pFileN, pFileS, nGal, pGall, pGalb,
    qInterp, qNoloop, qVerbose);

   /* Read the Beta ratio map */
   sprintf(pFileN, "%s/%s", pIPath, ppBetaFile[0]);
   sprintf(pFileS, "%s/%s", pIPath, ppBetaFile[1]);
   pBmap = lambert_getval(pFileN, pFileS, nGal, pGall, pGalb,
    qInterp, qNoloop, qVerbose);

   /* Allocate memory for output array */
   pInu = ccvector_build_(nGal);

   /* microK brightness temp (Beta map is actually negative of spectral
    * index, and ranges from roughly 2.5 < beta < 3.0)
    */
   for (iGal=0; iGal < nGal; iGal++) {
      pInu[iGal] = 1.0e6 * pAmap[iGal] * pow(0.408/pNu[iGal], pBmap[iGal]);
   }

   /* Convert units */
   if (strcmp(pUnitsName,pText_MJy) == 0) {
      /* MJy/sr */
      for (iGal=0; iGal < nGal; iGal++)
       pInu[iGal] /= fac_flux2temp(pNu[iGal]);
   } else if (strcmp(pUnitsName,pText_microK) == 0) {
      /* brightness temp micro-K */
   } else if (strcmp(pUnitsName,pText_thermo) == 0) {
      /* thermodynamic micro-K */
      for (iGal=0; iGal < nGal; iGal++)
       pInu[iGal] *= planckcorr(pNu[iGal]);
   } else {
      printf("ERROR: Invalid units name.\n");
      for (iGal=0; iGal < nGal; iGal++)
       pInu[iGal] = 0.0;
   }

   return pInu;
}

/******************************************************************************/
/* Return DIRBE color-correction (K-factor) for 100-micron map.
 */
float kfactor
  (float    alpha,
   float    temp)
{
   int      ia = -1; /* crash if Kfit not tabulated for alpha */
   int      im;
   int      k;
   float    log10T;
   float    log10Tpow;
   float    sum1;
   float    sum2;
   float    Kvalue;

   /* We tabulate the K-factor for only the following emissivity profiles */
   const float Kfitindx[] = {1.50, 1.67, 1.70, 2.00, 2.20, 2.60, 2.70};
   const float KfitAarr[][4] =
    { {  1.00000,  2.08243, -4.72422,  2.29118 },
      {  1.00000,  2.15146, -4.84539,  2.35210 },
      {  1.00000,  2.14106, -4.83639,  2.35919 },
      {  1.00000,  2.18053, -4.89849,  2.38060 },
      {  1.00000,  2.55941, -5.41290,  2.57867 },
      {  1.00000,  3.16383, -6.23131,  2.86900 },
      {  1.00000,  3.31600, -6.43306,  2.93939 } };
   const float KfitBarr[][4] =
    { { -0.88339,  4.10104, -4.43324,  1.76240 },
      { -0.87985,  4.10909, -4.43404,  1.76591 },
      { -0.93625,  4.19278, -4.46069,  1.77103 },
      { -0.80409,  3.95436, -4.27972,  1.70919 },
      { -0.80318,  4.20361, -4.55598,  1.80207 },
      { -0.50356,  4.07226, -4.70080,  1.87416 },
      { -0.41568,  4.02002, -4.72432,  1.88865 } };
   const int N_AINDEX = sizeof(KfitAarr) / sizeof(KfitAarr[0]);

   /* Find Kfit index for the given value of "alpha" */
   for (im=0; im < N_AINDEX; im++)
    if (fabs(Kfitindx[im] - alpha) < 1.e-4) ia = im;

   log10T = log10(temp);
   log10Tpow = 1.0;
   sum1 = KfitAarr[ia][0];
   sum2 = KfitBarr[ia][0];
   for (k=1; k < 4; k++) {
      log10Tpow *= log10T;
      sum1 += KfitAarr[ia][k] * log10Tpow;
      sum2 += KfitBarr[ia][k] * log10Tpow;
   }
   Kvalue = sum1 / sum2 ;

   return Kvalue;
}

/******************************************************************************/
/* Return Planck function in MJy/sr.
 * This is based upon the IDL procedure PLANCK() written by Rich Isaacman
 * for the COBE analysis software.
 */
float planck
  (float    temp, /* temperature in Kelvin */
   float    nu) /* frequency in GHz */
{
   float    val;
   float    chi;
   float    fhz;

   const float cspeed = 299792.458e+0; /* speed of light in km/s */
   const float hck = 14387.69e+0;      /* h*c/k          */
   const float SBconst = 2.7794795e-13; /* Stephan-Boltzmann * 15/pi^5 */

   /* Variable fhz is a scale factor used to go from units of nu_I_nu units to
      MJy/sr if the keyword is set.  In that case, its value is
      frequency in Hz * w/cm^2/Hz ==> MJy */
   fhz = nu * 1.0e-15;

   /* Introduce dimensionless variable chi, used to check whether we are on
      Wien or Rayleigh Jeans tails */
   chi = hck * nu / (cspeed * temp);

   if (chi < 0.001) {
      /* Rayleigh-Jeans side */
      val = SBconst * pow(temp,4.0) * pow(chi,3.0) / fhz;
   } else if (chi > 50.) {
      /* Wein tail */
      val = SBconst * pow(temp*chi,4.0) * exp(-chi) / fhz;
   } else {
      /* Exact solution */
      val = SBconst * pow(temp*chi,4.0) / ( (exp(chi) - 1.0) * fhz);
   }

   return val;
}
/******************************************************************************/
/* Compute factor to convert from flux/sr to brightness temp.
   Return conversion factor (MJy/sr) / (micro-K)
*/
float fac_flux2temp
  (float    nu) /* frequency in GHz */
{
   float    fac;
   const float k_b = 1.3806e-16;  /* erg K^-1 */

   fac = 4.5e-9 / (k_b * nu * nu);
   return fac;
}
/******************************************************************************/
/* Compute factor to convert from brightness temp to thermodynamic temp
   in micro-K.
*/
float planckcorr
  (float    nu) /* frequency in GHz */
{
   float    x;
   float    result;
   const float k_b = 1.3806e-23; /* J/K */
   const float h_Pl = 6.6262e-34;   /* J*s */
   const float T_cmb = 2.73;     /* K   */
 
   x = h_Pl * nu * 1.e9 / (k_b * T_cmb);
   result = pow(exp(x)-1.,2.0) / (x*x * exp(x));
 
   return result;
}

