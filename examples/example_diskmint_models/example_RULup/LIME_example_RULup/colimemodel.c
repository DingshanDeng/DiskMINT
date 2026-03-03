/*
 *  model.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"
 /* NLIN is first line of limegrid.dat; NRS is 1st line of rgrid */
// lines below modified by UG and DD
#define NLIN 20000 // number of r cells X number of Z cells
#define NCOL 11    // number of variables
#define NRS  100   // number of radial grid cells

   long double modeldata[NLIN][NCOL];
   long double radii[NRS];
   int rlbeg[NRS],rlend[NRS];
   double gtdreadin;

/* can include other headers here */

/******************************************************************************/
/*   Functions to read data and compute nn here */
/******************************************************************************/

/*  Read in the model data file */

void readdata(){

   FILE *file = fopen("limegrid.dat","r");
   int nlines, ncols;
   fscanf(file,"%d %d",&nlines,&ncols);

   int i=0,j=0;
   for (i=0; i<nlines; i++){
      for (j=0; j<ncols; j++) {
        fscanf(file,"%Lf",&modeldata[i][j]); 
   }
   }

   fclose(file);

   file = fopen("rgrid.dat","r");
   fscanf(file,"%d",&nlines);

   i=0,j=0;
   for (i=0; i<nlines; i++){
        fscanf(file,"%Lf %d %d",&radii[i],&rlbeg[i],&rlend[i]); 
   }

   fclose(file);
    
   file = fopen("ratio_g2d.dat", "r");
   fscanf(file, "%lf", &gtdreadin);
   printf("gtd = %lf", gtdreadin);
}

/******************************************************************************/
/*  Find the nearest gridpoint for some (x,y,z)    */

void findnn(double x, double y, double z, int *icol) {

     double radius=sqrt(x*x+y*y)/AU;
     double ht=fabsl(z)/AU;

/*    Defaults  */
        int irow=0;
        *icol=-1;

/*    Out of grid bounds */
/*     if(radius < radii[0] ) {
        return;
     }     
     if(radius > radii[NRS-1] ) {
        return;
     } */     
     if(radius < radii[0] ) {
        return;
     }     
     if(radius > radii[NRS-1] ) {
        return;
     }     

     double dif,dummy;
     int i;

/*    Nearest R zone */     
     dif=fabsl(radius-radii[0]);
     irow=0;
     for(i=0;i<(NRS-1);i++){
       dummy=fabsl(radius-radii[i]);
       if(dummy <= dif) {
           dif=dummy;
           irow=i;
       }
     }

      int ib=rlbeg[irow]-1; // ib means at which line this r begins
      int ie=rlend[irow]-1; // ie means at which line this r ends
      dif=1.0e30;
      for(i=ib;i<=ie;i++){
       dummy=fabsl(ht-modeldata[i][1]);
       if(dummy <= dif) {
           dif=dummy;
           *icol=i; 
       }
     }
}

/******************************************************************************/

void
input(inputPars *par, image *img){
  /* Read the data file for this model setup */
  /* * Basic parameters. See cheat sheet for details.  */
  readdata();  // reads in the model parameter/structure file
  par->radius		    = 500*AU; // 500*AU; // dist from star to max z at max r
  par->minScale	   	    = 0.04*AU; // 0.04*AU; // smallest scale
  par->pIntensity    	= 1e5; // 1e5; // minimum # of lines in tgas.dat(# total points)
  par->sinkPoints    	= 5e3; // 2e4; // 5e3; // min 0.01-0.1 of the above (# of surface points)
  par->dust             = "dustkappa_currentlyusing.inp";
  par->moldatfile[0] 	= "molecule_c18o.inp";
  //par->antialias	= 4;
  par->sampling		    = 2; // 1 uniform, 2 log dist. for radius,directions distr. uniformly on sphere
//   par->outputfile 	   = "populations.pop";
//   par->binoutputfile 	= "restart.pop";
//   par->gridfile		   = "grid.vtk";
  par->nThreads         = 25; // number of cores
/*  par->pregrid          = "model01"; */
  par->nSolveIters      = 17;
  par->init_lte         = 1;
  par->lte_only         = 1; // 0, set 1 for lte when checking, but lose emission from outer disk
  // LTE or non-LTE are generally very close to each other, to save time, just run lte_only = 1
  // par->collPartIds[0]           = 1; // 1 for h2; 2 for p-h2; 3 for o-h2
  // par->nMolWeights[0]           = 1.0; // the 100*percentage of the collisional partner (fraction weight)
  // // par->collPartNames[0]         = "H2"; // hard link the partner to H2 in the LAMDA table
  // par->collPartMolWeights[0]    = 2.0159; // 2.0159;
  // par->dustWeights[0]   = 0; // 1.0e-22
  par->collPartIds[0]   = 2; // p-h2
  par->nMolWeights[0]   = 0.25;
  par->collPartMolWeights[0]    = 2.0159;
  // par->dustWeights[0]   = 0;
  par->collPartIds[1]   = 3; // o-h2
  par->nMolWeights[1]   = 0.75;
  par->collPartMolWeights[1]    = 2.0159;
  // par->dustWeights[1]   = 0; // 1.0e-22
  //par->collPartIds[2]   =4; // e-
  //par->nMolWeights[2]   = 0.0;
  //par->dustWeights[2]   = 0.0;

  /* * Definitions for image #0.  Add blocks for additional images.  */
  img[0].nchan		= 241;	    // Number of channels; 81, 241 for including the middle channel
  img[0].velres     = 83.3;     // Channel resolution in m/s 250, 83.3
  img[0].trans		= 1;        // 1 for (2-1), 2 for (3-2) // zero-indexed J quantum number
  img[0].pxls		= 151;      // Pixels per dimension
  img[0].imgres	    = 0.04;	    // Resolution in arc seconds in RULup Huang data 0.04
  img[0].distance	= 158.9*PC; // source distance in m
  img[0].source_vel = 0;        // source velocity in m/s 4500.0
  // the img.unit is changed between 1 and 4 to make different files
  img[0].unit		= 1;	    // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[0].filename	= "LIME_image_c18o_2-1.fits";	// Output filename
  img[0].azimuth    = 0.0;
  // the incl and PA from Huang+2020
  img[0].incl       = 18.8*(3.14/180.0); // incl angle
  img[0].posang     = -121*(3.14/180.0);
  
  //  /* * Definitions for image #1.  Add blocks for additional images.  */
  // img[1].freq       = 219.5603541e9; // in Hz for C18O
  // // img[1].nchan		= 81;	   // Number of channels; 80
  // // img[1].velres     = 250;    // Channel resolution in m/s 250; 83.3
  // // img[1].trans		= 1;     // zero-indexed J quantum number
  // img[1].pxls		= 151;   // Pixels per dimension 1001
  // img[1].imgres	    = 0.04;	// Resolution in arc seconds in RULup Huang data 0.04
  // img[1].distance	= 158.9*PC; // source distance in m
  // // img[1].source_vel = 0; // source velocity in m/s 4500.0
  // // the img.unit is changed between 1 and 4 to make different files
  // img[1].unit		= 1;	  // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  // img[1].filename	= "LIME_image_c18o_ct.fits";	// Output filename
  // img[1].azimuth    = 0.0;
  // // the incl and PA from Huang+2020
  // img[1].incl       = 18.8*(3.14/180.0); // incl angle
  // img[1].posang     = -121*(3.14/180.0);
}
/******************************************************************************/

void
density(double x, double y, double z, double *density){
    /*
    double r;             // Define variable for radial coordinate
    r=sqrt(x*x+y*y+z*z);  //  Calculate radial distance from origin 


    * Calculate a spherical power-law density profile
    * (Multiply with 1e6 to go to SI-units)
    density[0] = 1.5e6*pow(r/(300*AU),-1.5)*1e6;
    */
    double dummy;
    int i;
    findnn(x, y, z, &i);
    if(i<0) {
     density[0]=1.0e-30;
     // density[1]=1.0e-30;
     // density[2]=1.0e-30;
    }
    else {
    dummy=modeldata[i][2]*modeldata[i][5]; // n_H * X(H2) 
    density[0] = dummy; // n_H * X(H2) for all H2
    // density[0] = dummy*0.25; //*modeldata[i][9]; // n_H * X(H2)*0.25 for p
    // density[1] = dummy*0.75; // modeldata[i][10];  // n_H * X(H2)*0.75 for o
    //density[2] = dummy*modeldata[i][7];  // n_H * X(e-) */

    }
    /*   density[2] = modeldata[i][2]*modeldata[i][7];  // n_H * X(e-) */
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){

    int i;
    findnn(x, y, z, &i);
    if(i<0) {
     temperature[0]=2.7;
     temperature[1]=2.7;
    }
    else {
    temperature[0] = modeldata[i][3]; // gas temperature
    //if(temperature[0]<20) temperature[0]=30;
    temperature[1] = modeldata[i][4]; // dust temperature
    }

}

/******************************************************************************/

// void
// abundance(double x, double y, double z, double *abundance){
// // The abundance is the fractional abundance with respect to 
// // a weighted sum of the densities supplied for the collision partners. 
// // If the user does not supply the weights via the nMolWeights parameter, 
// // the code will try to guess them.

//     int i;
//     double radius=sqrt(x*x+y*y)/AU;
//     findnn(x, y, z, &i);
//     if (i<0) {
//     abundance[0]=0.0;
//     }
//     else {
//     abundance[0] = modeldata[i][8] * 2.0; // abundance of C18O [NOT density!!], typically 1.4e-7
//     // this should be the abundance related to the H2 density (collider)
//     }
// }

/******************************************************************************/

void
molNumDensity(double x, double y, double z, double *nmol){
// As an alternative to the abundance function, 
// the user is now able to supply a function which specifies directly 
// the number density of each of the radiating species.
// The densities are number densities, that is, 
// the number of molecules per unit volume (in cubic meters, not cubic centimeters).
    
    
    int i;
    double radius=sqrt(x*x+y*y)/AU;
    findnn(x, y, z, &i);
    if (i<0) {
    nmol[0] = 0.0;
    }
    else {
    nmol[0] = modeldata[i][2]*modeldata[i][8]; // n_H * X(C18O)
    }
    
    // nmol[1] = f1(x,y,z);
    // ...
    // nmol[n] = fn(x,y,z);
}


/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){

  // *doppler = 10.0;  // in m/s
    int i;
    double tg;
    findnn(x, y, z, &i);
    if(i<0) {
      tg=2.7;
    }
    else {
    tg = modeldata[i][3];
    }

    *doppler = (0.0 + 0.13*sqrt(tg/30))*1000.0; // in m/s
    // *doppler = 0.0; // in m/s, no broadening

}

/******************************************************************************/

void
gasIIdust(double x, double y, double z, double *gtd) {

    // *gtd = 100.0;
    *gtd = gtdreadin;

}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){

    double r, rToUse, ffSpeed, mstar;
    const double rMin = 0.01*AU; /* This cutoff should be chosen smaller than par->minScale but greater than zero (to avoid a singularity at the origin). */

    /*
    * Calculate radial distance from origin
    */
    r = sqrt(x*x+y*y);
    if(r>rMin)
      rToUse = r;
    else
      rToUse = rMin;
    
    mstar  = 1.0; // 0.7;
    vel[0] = sqrt(6.67e-11*mstar*1.99e30/rToUse)*sin(atan2(y,x));
    vel[1] = -sqrt(6.67e-11*mstar*1.99e30/rToUse)*cos(atan2(y,x));
    vel[2] = 0.0;

    /* Keplerian velocity */
    // vel[0] = sqrt(6.67e-11*1.2*1.99e30/sqrt(x*x+y*y))*sin(atan2(y,x));
    // vel[1] = -sqrt(6.67e-11*1.2*1.99e30/sqrt(x*x+y*y))*cos(atan2(y,x));
    // vel[2] = 0.0;
}

/******************************************************************************/












