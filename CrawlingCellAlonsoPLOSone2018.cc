// Program by Sergio Alonso
//Alonso, Sergio, Maike Stange, and Carsten Beta. 
//"Modeling random crawling, membrane deformation and intracellular polarity of motile amoeboid cells." 
//PloS one 13.8 (2018): e0201977.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
using namespace std;
   
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define pi 3.1415592653589793

				      
// Parameters for the integration
#define nx 300         // total pixels in the x direction.100 
#define ny 300         // total pixels in the y direction.100  
#define nsw 1         // Number of cells   For 200x200 nsw max is 81
float dt=0.0020;     // temporal step 0.0025;   
float dx=0.15;         // spatial step in um // old 0.1  
float stoptime=50.0; // total time of the integration.
int xof=0;
int yof=0;
float xcenter_mov=0;
float ycenter_mov=0;
float xmembrane=0;
float ymembrane=0;
int motion=0;    // 1 motion 0 no motion

// Parameter for Sawai Model


float D_b   = 0.5 ;   // in um2/s. 
float xdecay = 0.1;    // 
float kact  = 2.0;   // 5.0
float react = 1.0;  


// Parameter for the Noise
float tauxi  =  10.0;    // in   s
float sigma2 = 0.0150;       // 
float amplitude= sqrt(2.0*sigma2*dt)/(dx);
float tauxiinv= 1.0/tauxi;
float sigma= sqrt(sigma2);
int seed=3498705;
int MembraneNoise=1;
int GlobalNoise=0;
//long id=90348;


// Characteristics of the simulations
int    cercle=0;     /* Circular system of radius (nx/2-5) 
			with no-flux boundary conditions. 
			if 0 periodic boundary conditions in 
			a square system   */
int    pbc=1;     /* periodic boundary conditions    */
int    channel=0;     /* channel of size (20) 
			with no-flux boundary conditions. 
			if 0 periodic boundary conditions in 
			a square system   */
#define nv  500     // Time steps before saving data



// Phase field parameters
#define rodx  40
float   ro =  (rodx * dx  );
#define  threshold 0.001
float vol=(pi*rodx*rodx);


// Parameters of the activator conservation
float deltao2 = 0.5;
float delta2  = 0.5;
float vol2    = (0.25*(pi*rodx*rodx));  
float ro2     = 0.5*ro;

 
#define alphao 3.0     //          Active force
float alpha=alphao;
#define epsilon 0.750    //  1.50   
#define gamma 2.0     // 300.0    // Tension of the membrane
#define tau 2.0         // 3.0         Time scale of the whole equation for the membrane dynamics
#define ma  0.50       //             Volume conservation
#define mb  1.0       // 500   //    Repulsion among cells



// Concentrations and phase field matrix 

float exp_dt_tauxi = exp(-dt/tauxi);
float exp_2_dt_tauxi = exp(-2.*dt/tauxi);


float bct[nx][ny];      // Concentration variable
float bct2[nx][ny];     // Concentration variable for integration
float dx2bct[nx][ny];   // Laplacian of the Concentration 

float phis[nx][ny];          // Contribution of the boundaries of the system: 1 outside 
float phi[nx][ny];      // Phase field
float phi2[nx][ny];     // Phase field for integration
float gradphi[nx][ny];  // Gradient of the Phase field
float dphi2[nx][ny];    // Laplacian of the Phase field





// Matrix for the periodic boundary conditions
int ixyp1[nx];
int ixym1[nx];

//Averages 

float ctotal;
float ctotal2;
float xcenter,ycenter;
float xcenter_old,ycenter_old, velx, vely;
float anglex,angley;
int xframe,yframe;

// **************************************************************
// Subroutine generation of the index for boundary conditions 
//                 nx=ny!!!!   
// **************************************************************
void generation_index()
{
int ix,iz;
// for (ix=0;ix<nx-1;ix++){ixyp1[ix]=ix+1;} ixyp1[nx-1]=nx-1;
// for (ix=1;ix<nx;ix++){ixym1[ix]=ix-1;} ixym1[0]=0;
 for (ix=0;ix<nx-1;ix++){ixyp1[ix]=ix+1;} ixyp1[nx-1]=0;
 for (ix=1;ix<nx;ix++){ixym1[ix]=ix-1;} ixym1[0]=nx-1;
}
// **************************************************************
// **************************************************************
// **************************************************************


// **************************************************************
//		 Generation of the phase field: circular domain
// **************************************************************
void generation_phi(float (*phi)[ny], int ixo, int iyo)
{
int ix,iy,points;
float r,sum;

for (ix=0;ix<nx;ix++)
  for (iy=0;iy<ny;iy++)
  { 
  r= dx*sqrt (float ((ix-ixo)*(ix-ixo) +(iy-iyo)*(iy-iyo) )) ;
  phi[ix][iy] = 0.5+0.5*tanh((ro-r)/(dx*epsilon));
  }

}
// **************************************************************
// **************************************************************
// **************************************************************

// **************************************************************
// 	Generation of the initial concentration (random)
// **************************************************************
void generation_act_phi(int ixo, int iyo)
{
int ix,iy,points;
float r,sum;

for (ix=0;ix<nx;ix++)
  for (iy=0;iy<ny;iy++)
  {
	
    	bct[ix][iy] = 0.0  + sqrt(sigma)*( rand() / ( (float)RAND_MAX )  - 0.5);
  }
}
// **************************************************************
// **************************************************************
// **************************************************************

// **************************************************************
// 	Generation of the initial concentration (cercle)
// **************************************************************
void generation_act_phi_cercle(int ixo, int iyo)
{
int ix,iy,points;
float r,sum;

for (ix=0;ix<nx;ix++)
  for (iy=0;iy<ny;iy++)
  { 
  r= dx*sqrt (float ((ix-ixo)*(ix-ixo) +(iy-iyo)*(iy-iyo) )) ;
  bct[ix][iy] = 0.5+0.5*tanh((ro2-r)/(dx*epsilon));
  }


}
// **************************************************************
// **************************************************************
// **************************************************************


// **************************************************************
//	Generation of the phase field: circular closed domain
// **************************************************************
void generation_phis_Cercle(float (*phis)[ny])
{
int ix,iy,points;
float r,sum;

for (ix=0;ix<nx;ix++)
  for (iy=0;iy<ny;iy++)
  { 
  r= sqrt (float ((ix-nx/2)*(ix-nx/2) +(iy-ny/2)*(iy-ny/2) )) ;
  phis[ix][iy] = 0.5+0.5*tanh((nx/2-5-r)/(dx*epsilon));
  }
}
// **************************************************************
// **************************************************************
// **************************************************************



// **************************************************************
//	Generation of the phase field: channel domain
// **************************************************************
void generation_phis_channel(float (*phis)[ny])
{
int ix,iy,points;
float r,sum;

for (ix=0;ix<nx;ix++)
  for (iy=0;iy<ny;iy++)
  { 
  r= sqrt (float ((ix-nx/2)*(ix-nx/2) )) ;
  phis[ix][iy] = 0.5+0.5*tanh((20-r)/(dx*epsilon));
  }

}
// **************************************************************
// **************************************************************
// **************************************************************

// **************************************************************
// Derivate phase field model:   nabla*(phi nabla u)  
// **************************************************************
void Laplacian_b_phi_N (float (*c)[ny],float (*dx2c)[ny],float (*phi)[ny]){
int ix,iy,iz,isw;
float Dphiup,Dphidown,Dphiright,Dphileft;
float dudphi;
for (ix=0;ix<nx;ix++){
for (iy=0;iy<ny;iy++){

Dphiup    = (phi[ix][iy] + phi[ix][ixyp1[iy]]);
Dphidown  = (phi[ix][iy] + phi[ix][ixym1[iy]]);
Dphiright = (phi[ix][iy] + phi[ixyp1[ix]][iy]);
Dphileft  = (phi[ix][iy] + phi[ixym1[ix]][iy]);

dx2c[ix][iy]= 0.5 * 
( Dphiup    * (c[ix][ixyp1[iy]] - c[ix][iy]) - Dphidown * (c[ix][iy] - c[ix][ixym1[iy]])
+ Dphiright * (c[ixyp1[ix]][iy] - c[ix][iy]) - Dphileft * (c[ix][iy] - c[ixym1[ix]][iy]) 
) /(dx*dx);
}}

}
// **************************************************************
// **************************************************************
// **************************************************************







// **************************************************************
// Laplacian  
// **************************************************************
void Laplacian (float (*c)[ny],float (*dx2c)[ny]){
int ix,iy,iz;
float dudphi;
for (ix=0;ix<nx;ix++){
for (iy=0;iy<ny;iy++){
dx2c[ix][iy]= (c[ixyp1[ix]][iy]  + c[ixym1[ix]][iy]
	     + c[ix][ixyp1[iy]]  + c[ix][ixym1[iy]]
	     - 4.*c[ix][iy]  )/(dx*dx);
}}
}
// **************************************************************
// **************************************************************
// **************************************************************

// **************************************************************
// Gradient of phi 
// **************************************************************
void Gradient (float (*c)[ny],float (*gradu)[ny]){
int ix,iy,iz;
float dudphi;
for (ix=0;ix<nx;ix++){
for (iy=0;iy<ny;iy++){
gradu[ix][iy]= sqrt ( (c[ixyp1[ix]][iy] - c[ixym1[ix]][iy]) * (c[ixyp1[ix]][iy] - c[ixym1[ix]][iy]) 
		+     (c[ix][ixyp1[iy]] - c[ix][ixym1[iy]]) * (c[ix][ixyp1[iy]] - c[ix][ixym1[iy]]) 
		) / (2. * dx) ;
}}
}
// **************************************************************
// **************************************************************
// **************************************************************



// **************************************************************
// Gaussian noise 
// **************************************************************
float gasdev(long *idum)
{
float ran1(long *idum);
static int iset=0;
static float gset;
float fac,rsq,v1,v2;
if (iset == 0) {
 do { v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
 } while ( rsq >= 1.0 || rsq == 0.0);
 
 fac=sqrt(-2.0*log(rsq)/rsq);
 gset=v1*fac;
 iset=1;
 return v2*fac;
 } else {
 iset=0;
 return gset;
 }
}


float ran1(long *dseed){
int a=69069,c=1,xm;
float rm,rand_xx;
xm=int (pow( 2., 32 ));
rm=1./float (xm);
*dseed=*dseed*a+c % xm;
rand_xx=float (*dseed) * rm;
return(rand_xx);
}



float GN_generator()
{/* Generates additive white Gaussian Noise samples with zero mean and a standard deviation of 1. */
  float temp1;
  float temp2;
  float result;
  int p;
  p = 1;
  while( p > 0 )
  {
    temp2 = ( rand() / ( (float)RAND_MAX ) ); 
    if ( temp2 == 0 ) {p = 1;}
    else {p = -1;}
  }
  temp1 = cos( ( 2.0 * pi ) * rand() / ( (float)RAND_MAX ) );
  result = sqrt( -2.0 * log( temp2 ) ) * temp1;
  return result; 
}
 
 
 
float ou(float xi_xy)
{
float result;
//result = xi_xy * exp(-dt/tauxi) + sqrt( (sigma/tauxi)*(1.-exp(-2.*dt/tauxi)) ) * GN_generator();
result = xi_xy * exp_dt_tauxi + sqrt( (sigma/tauxi)*(1.- exp_2_dt_tauxi) ) * GN_generator();
return result;
}


// **************************************************************
// **************************************************************
// **************************************************************


// AIXO ES LA -G'(phi)
float f_phi(float phi)
{
float fu;
fu= 72.0 * ( phi* (1.- phi)*(phi - 0.5));
return fu;
}


// **************************************************************
// ***** Save Data
// **************************************************************
void Save_Parameters ()
{
ofstream filenoise ("Valuedates.dat");
   filenoise << "Diffusion             :  " << D_b    << "\n";
   filenoise << "Tau_xi                :  " << tauxi  << "\n";
   filenoise << "White noise amplitude :  " << sigma2  << "\n";
   filenoise << "White noise variance  :  " << sqrt(sigma2)  << "\n";
   filenoise << "Sigma2                :  " << sigma2  << "\n";
   filenoise << "Sigma2/tauxi          :  " << sigma2/tauxi  << "\n";
   filenoise << "Decay                 :  " << xdecay  << "\n";
   filenoise << "Nx and Ny             :  " << nx  << "\n";
   filenoise << "dt                    :  " << dt  << "\n";
   filenoise << "dx                    :  " << dx  << "\n";
   filenoise << "Total time            :  " << stoptime  << "\n";
   filenoise << "Active force          :  " << alpha  << "\n";
   filenoise << "length                :  " << epsilon  << "\n";
   filenoise << "Membrane Tension      :  " << gamma  << "\n";
   filenoise << "Time scale membrane   :  " << tau  << "\n";
   filenoise << "Volume conservation   :  " << ma  << "\n";
   filenoise << "Repulsion among cells :  " << mb  << "\n"; 
 filenoise.close();   
}
// **************************************************************
// **************************************************************
// **************************************************************





// **************************************************************
// ***** Save Picture
// **************************************************************
void save_picture_Con (float (*w)[ny], float  (*vv)[ny], int step)
{
int ix,iy,isw,numero=1000000;
double ww[nx][ny];
char title[20],titlec[20],titleact[20],titlegnu[20];


for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) { ww[ix][iy]=w[ix][iy]*phi[ix][iy];
if (ww[ix][iy]>=1.0) ww[ix][iy]=1.0;
if (ww[ix][iy]<=0.0) ww[ix][iy]=0.0;
}
numero=numero+step;
sprintf(titlec,"Color2_%d.ppm",numero);
ofstream filec (titlec);
filec  << "P3" << "\n";  filec << ny << "\n";  filec << nx << "\n";  filec << 999 << "\n";
for (iy=0;iy<ny;iy++){ 
for (ix=0;ix<nx;ix++){
  filec << int (0) << " " <<  int ( (950*(phi[ix][iy]))  - 950*ww[ix][iy]*phi[ix][iy])<< " " <<  int ( 950* ww[ix][iy])   <<"\n";  
}
}
filec.close();
}
// **************************************************************
// **************************************************************
// **************************************************************


// **************************************************************
//Generation of the phase field: no domain periodic boundary conditions
// **************************************************************
void generation_phis(float (*phi)[ny])
{
int ix,iy,points;
float r,sum;
  
for (ix=0;ix<nx;ix++)
  for (iy=0;iy<ny;iy++)
  { 
  phi[ix][iy] = 1.0;
  }
 
}
// **************************************************************
// **************************************************************
// **************************************************************




// **************************************************************
//		 Generation of Active swimmers
// **************************************************************
void  generation_Active_swimmers (float (*phi_bc)[ny],float xcenter,float ycenter)
{
int ix,iy,isw,isw2,xsw,ysw,points,xran,yran;
int bandera;
int xx,yy;
float r,rad_o,dist;

  xsw= nx/2; 
  ysw= ny/2;
  generation_phi(phi,xsw,ysw); 
  xx=xsw; yy=ysw;  xcenter=xsw; ycenter=ysw;
  xran= - 4 + int  (8.*(rand()/(RAND_MAX+1.0)));
  yran= - 4 + int  (8.*(rand()/(RAND_MAX+1.0)));  
    generation_act_phi_cercle(xsw+xran,ysw+yran);  
}
// **************************************************************
// **************************************************************
// **************************************************************







void Movement_x_m(float (*noise)[ny])
{
int ix,iy;
float acces[nx][ny];
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=phi[ixyp1[ix]][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { phi[ix][iy]=acces[ix][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=bct[ixyp1[ix]][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { bct[ix][iy]=acces[ix][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=noise[ixyp1[ix]][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { noise[ix][iy]=acces[ix][iy];}
}

void Movement_x_p(float (*noise)[ny])
{
int ix,iy;
float acces[nx][ny];
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=phi[ixym1[ix]][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { phi[ix][iy]=acces[ix][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=bct[ixym1[ix]][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { bct[ix][iy]=acces[ix][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=noise[ixym1[ix]][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { noise[ix][iy]=acces[ix][iy];}
}

void Movement_y_m(float (*noise)[ny])
{
int ix,iy;
float acces[nx][ny];
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=phi[ix][ixyp1[iy]];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { phi[ix][iy]=acces[ix][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=bct[ix][ixyp1[iy]];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { bct[ix][iy]=acces[ix][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=noise[ix][ixyp1[iy]];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { noise[ix][iy]=acces[ix][iy];}
}

void Movement_y_p(float (*noise)[ny])
{
int ix,iy;
float acces[nx][ny];
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=phi[ix][ixym1[iy]];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { phi[ix][iy]=acces[ix][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=bct[ix][ixym1[iy]];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { bct[ix][iy]=acces[ix][iy];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { acces[ix][iy]=noise[ix][ixym1[iy]];}
	for (ix=0;ix<nx;ix++)  for (iy=1;iy<ny-1;iy++) { noise[ix][iy]=acces[ix][iy];}
}



float Perimeter_particles(float (*gradph)[ny]) 
{
int ix,iy;
float xperimeter=0.0; 
for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) { xperimeter+=gradph[ix][iy];}
return xperimeter;
}




float Size_particles(float (*ph)[ny])
{
int ix,iy,isw;
float phtot=0.0; 
for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) { phtot+=ph[ix][iy];}
return phtot;
}


void reaction_particles (float (*ac)[ny], float (*ph)[ny], float (phit))
{
int ix,iy,isw;    
float tot,tot2;
     tot=0.0; tot2=0.0;
     for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
     	tot+=ac[ix][iy]*ph[ix][iy];
     	tot2+=ac[ix][iy]*ac[ix][iy]*ph[ix][iy];
     }
     ctotal=tot/phit;
     ctotal2=tot2/phit;
     delta2 = 0.5 + 0.001*(tot-vol2);
	if (delta2 > 1.0) delta2=1.0;
}



void intial_values(float (*phi_b)[ny],float(*phi_t)[ny],float(*phis)[ny],float(*phi)[ny],float(*xi)[ny])
{ 
int ix,iy,isw;    
    for (ix=0;ix<nx;ix++)
     for (iy=0;iy<ny;iy++){
      phi_b[ix][iy] = phis[ix][iy] ; 
      phi_t[ix][iy] = 1.-phis[ix][iy] ;
      xi[ix][iy]=0.0; 
    }

}



void Center_particles_movement (float (phtot), float (*ph)[ny], int stp)
{ 
int ix,iy,isw;
float tot_a; 
	for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) { 
			xcenter_mov+=ix*ph[ix][iy];
			ycenter_mov+=iy*ph[ix][iy];
			}
	xcenter_mov=xcenter_mov/phtot;
	ycenter_mov=ycenter_mov/phtot;
}


void Center_vel_particles (float (phtot), float (*ph)[ny], float (*ac)[ny], int stp)
{ 
int ix,iy,isw;
float tot_a; 
		xcenter=0.0; anglex=0.0; tot_a=0.0;
		if (xcenter_old<=nx/2) {  
		for (ix=0;ix<nx-2*rodx;ix++) for (iy=0;iy<ny;iy++) { 
			xcenter+=ix*ph[ix][iy];
			anglex+= ix*ph[ix][iy]*ac[ix][iy];
			tot_a+=ph[ix][iy]*ac[ix][iy];
			}
		for (ix=nx-2*rodx;ix<nx;ix++) for (iy=0;iy<ny;iy++) { 
			xcenter+=(ix-nx+1)*ph[ix][iy]; 
			anglex+=(ix-nx+1)*ph[ix][iy]*ac[ix][iy]; 
			tot_a+=ph[ix][iy]*ac[ix][iy];
			}
		}
		
		else {
		for (ix=0;ix<2*rodx;ix++) for (iy=0;iy<ny;iy++) { 
			xcenter+=(nx-1+ix)*ph[ix][iy]; 
			anglex+=(nx-1+ix)*ph[ix][iy]*ac[ix][iy]; 
			tot_a+=ph[ix][iy]*ac[ix][iy];
			}
		for (ix=2*rodx;ix<nx;ix++) for (iy=0;iy<ny;iy++) { 
			xcenter+=ix*ph[ix][iy]; 
			anglex+=ix*ph[ix][iy]*ac[ix][iy]; 
			tot_a+=ph[ix][iy]*ac[ix][iy];
			}
		}
		xcenter=xcenter/phtot;
		anglex=anglex/tot_a;
			
		ycenter=0.0; angley=0.0;
		if (ycenter_old<=ny/2) { 
		for (ix=0;ix<nx;ix++) for (iy=0;iy<ny-2*rodx;iy++) {  
			ycenter+=iy*ph[ix][iy];
			angley+=iy*ph[ix][iy]*ac[ix][iy];
		}
		for (ix=0;ix<nx;ix++) for (iy=ny-2*rodx;iy<ny;iy++) { 
			ycenter+=(iy-ny+1)*ph[ix][iy];
			angley+=(iy-ny+1)*ph[ix][iy]*ac[ix][iy];
		}
		} 
		
		else {
		for (ix=0;ix<nx;ix++) for (iy=0;iy<2*rodx;iy++) { 
			ycenter+=(ny-1+iy)*ph[ix][iy];
			angley+=(ny-1+iy)*ph[ix][iy]*ac[ix][iy];
			}
		for (ix=0;ix<nx;ix++) for (iy=2*rodx;iy<ny;iy++) { 
			ycenter+=iy*ph[ix][iy];
			angley+=iy*ph[ix][iy]*ac[ix][iy];
			}
		}
		ycenter=ycenter/phtot;
		angley=angley/tot_a;
				
		velx= dx*(xcenter + xof- xcenter_old)/(nv*dt);
		vely= dx*(ycenter + yof- ycenter_old)/(nv*dt);
		
		anglex=dx*(anglex-xcenter);
		angley=dx*(angley-ycenter);
		
		if (xcenter<=0)    { xcenter=xcenter+nx-1; xframe--;}
		if (xcenter>=nx-1) { xcenter=xcenter-nx+1; xframe++;}
		if (ycenter<=0)    { ycenter=ycenter+ny-1; yframe--;}
		if (ycenter>=ny-1) { ycenter=ycenter-ny+1; yframe++;}
			
		xcenter_old=xcenter + xof; 
		ycenter_old=ycenter + yof;
		if(stp <= 100 ) {velx= 0.0; vely= 0.0;}
}



// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************
//    MAIN PROGRAM  MAIN PROGRAM  MAIN PROGRAM  MAIN PROGRAM  
// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************


int main()
{
  float time=0.,phistot=0.0;
  float phitot;    // Area of each particle
  float perimeter; // Perimeter of each cell   
  float perimeter_o; // Initial Perimeter of each cell  
  float phi_b[nx][ny];  // 1 for free area
  float phi_t[nx][ny];  // 1 for ocupied area (swimmer or Bound. cond.)
  float xi[nx][ny];          
  float xcen_fr,ycen_fr;
  float dphidt,fabsawai;
  float alphab=0.0;
  int step=0,nstep=0,nstep2=0;
  int ix,iy,iz,isw,isw2,nn,ixcen,iycen;
  int xsw,ysw;	


//******************************************************************************************************************   
//**** Files opening ***********************************************************************************************  

  ofstream file1 ("Motion_1.dat");
  ofstream file2 ("Speed.dat");
  ofstream file3 ("Positions.dat");
  ofstream file4 ("Speed_xy.dat");
  ofstream file6 ("Angle_xy.dat");
  ofstream filecon ("Concentration.dat");
  ofstream fileexc ("Excentricity.dat");
  ofstream filephixy ("Phi_XY_center.dat");
  ofstream filephixy2 ("Phi_XY_Membrane.dat"); 
  Save_Parameters();
 
  
//******************************************************************************************************************   
//**** Initial processes *******************************************************************************************  

srand(seed);     
generation_index();
if (cercle==1)  {generation_phis_Cercle(phis);}                         // generates the domain
if (pbc==1)     {generation_phis(phis);} 				// generates the domain
if (channel==1) {generation_phis_channel(phis);}			// generates the domain

generation_Active_swimmers(phis,xcenter_old,ycenter_old);  		// generates the cell
intial_values(phi_b,phi_t,phis,phi,xi);  


phitot=vol; velx=0.0; vely=0.0; xframe=0; yframe=0; xcen_fr=xcenter_old; ycen_fr=ycenter_old;

Gradient(phi,gradphi);
perimeter_o=Perimeter_particles(gradphi); // Perimeter initial 
cout << "Perimeter   : " << (2.*pi*rodx) << "   Calculated   : " << perimeter_o << " Relation Constant   "  << (2.*pi*rodx)/perimeter_o << "\n";

//******************************************************************************************************************   
//**** Process *****************************************************************************************************   

  for(step=0;time<=stoptime;step++)
    {
    time=time+dt;
    if (step==1000) alphab= 0.01;
	Laplacian(phi,dphi2);Gradient(phi,gradphi);
	Laplacian_b_phi_N(bct,dx2bct,phi);  // version 1

    for (ix=0;ix<nx;ix++)
    for (iy=0;iy<ny;iy++)   {
	dphidt=  (gamma /tau) * (f_phi(phi[ix][iy])/(epsilon*epsilon)  + dphi2[ix][iy]) 
	-  (mb/tau) *  (phi[ix][iy]*phi_t[ix][iy]) * gradphi[ix][iy] 
	-  (ma/tau) *  (phitot-vol) * gradphi[ix][iy]	
	+ (alpha/tau) * (bct[ix][iy]*phi[ix][iy]) * gradphi[ix][iy]
	;
	
	phi2[ix][iy] = phi[ix][iy] + dt * dphidt ;
        xi[ix][iy]= xi[ix][iy] * (1.0 - dt * tauxiinv) + amplitude*GN_generator();   // version 1

	if (phi[ix][iy] > threshold) {
	bct2[ix][iy]= bct[ix][iy] 
	+ dt * D_b * dx2bct[ix][iy]/phi[ix][iy] 
	- dt * (react) * xdecay * bct[ix][iy] 
	- dt*(1.0-react)*phi[ix][iy]*xi[ix][iy] 
	+ dt*  (react) * kact* bct[ix][iy]*(1.-bct[ix][iy])*(bct[ix][iy]-delta2) 
	+ dt* MembraneNoise * (10.0 * react) *phi[ix][iy]*xi[ix][iy]*(1.-phi[ix][iy]) 
	+ dt * GlobalNoise * phi[ix][iy]*xi[ix][iy]; 


				}
	else {	bct2[ix][iy] = bct[ix][iy] - dt * bct[ix][iy] * xdecay; }
}
    
//******************************************************************************************************************   
//**** Next step and control of the data ***************************************************************************
    for (ix=0;ix<nx;ix++)
    for (iy=0;iy<ny;iy++) {
	if (phi2[ix][iy] < 0.00)  {  cout << "!:  Value of phi: " << phi2[ix][iy] << ",  Value of b: " << bct2[ix][iy] << ",  Value of xi: " << xi[ix][iy] << ",  Value of grad Phi: " << gradphi[ix][iy] << ", ix: " << ix << ", iy: " << iy <<   "\n" ;return 0;}
	if (phi2[ix][iy] > 1.00)   { phi2[ix][iy]= 1.0; }	
	
	bct[ix][iy] = bct2[ix][iy];
	if (bct2[ix][iy] > 1.00)   { bct[ix][iy] = 1.0;}
	if (bct2[ix][iy] < -1.00)  { bct[ix][iy] =-1.0;}
	
	phi_b[ix][iy] = (1.-phi2[ix][iy]);
	phi_t[ix][iy] = (1.-phis[ix][iy]);
	phi[ix][iy]=phi2[ix][iy];
    }
     
     
//******************************************************************************************************************   
//**** Calculation and saving data  ********************************************************************************

// Calculation of the total size of each particle
     phitot=Size_particles(phi);
     perimeter=Perimeter_particles(gradphi);

// Save picture		
	if(step % 500 == 0) {nstep2=nstep2+1; 
				save_picture_Con (bct,phi_b,nstep2); 
				}	

//Moviment of Frame
	if (motion==1){	
		if(step % 5 == 0) {
			Center_particles_movement (phitot, phi, step); 
			if (xcenter_mov>=(2*nx/3)) {Movement_x_m(xi);xof++;}
			if (xcenter_mov<=(nx/3))   {Movement_x_p(xi);xof--;}
			if (ycenter_mov>=(2*ny/3)) {Movement_y_m(xi);yof++;}
			if (ycenter_mov<=(ny/3))   {Movement_y_p(xi);yof--;}
		}
	}


//calculation of center of each particle and the velocity
	if(step % nv == 0) { 
	Center_vel_particles (phitot, phi, bct, step);
        file1 << time << " " << xcenter + xof + float(1.0*(xframe*nx)) << " " << ycenter + yof + float (1.0*(yframe*ny)) << " " << velx << " " << vely<< " "  
	<<   sqrt((velx*velx)+(vely*vely)) << " " << phitot << "\n";
	filephixy << time << " " << phi[int (xcenter)][int (ycenter)]<< "\n";
	filephixy2 << time << " " << phi[int (xmembrane)][int (ymembrane)]<< "\n";
			}



// Save velocities and positions in files			
	if(step % nv == 0) {   
		file2 << time  <<  " " <<  sqrt((velx*velx)+(vely*vely)) << "\n";
		file3 << xcenter + xof<< " " << ycenter + yof <<  " " << time << "\n"; 
		file4 << velx<< " " << vely <<  " " << time << "\n";
		file6 << anglex<< " " << angley <<  " "<< time << "\n";
	}			

 
// Calculation of the total concentration on each cell
	if(step % 100 == 0) {
                reaction_particles (bct,phi,phitot); 
}
// Save data  
 	if(step % nv == 0) {
		ixcen=xcenter; iycen=ycenter; 
		cout << time << ", Perimeter : " <<  perimeter << ", Exc : " <<  perimeter/perimeter_o << ", Phi : " <<  phitot<< ", Conc : " <<  ctotal<< ", Conc2 :  " << ctotal2<< ", xcenter : " << ixcen  << ", ycenter : " << iycen << ", Conccenter : " << bct[ixcen][iycen] << ", Frame : " << xframe << ", " << xcenter + xof + float(1.0*(xframe*nx)) << ", Frame " << yframe  <<  " , Delta   : " << delta2 << "\n";
		filecon << time << " " <<  ctotal << " " <<  ctotal2<< " " << xi[nx/2][ny/2] << "\n";
		fileexc << time << " " <<  perimeter/perimeter_o << " " << perimeter << "\n";
}


}


//******************************************************************************************************************   
//**** Closing files and end ***************************************************************************************

file1.close(); file2.close(); file3.close(); file4.close(); file6.close(); filecon.close(); fileexc.close();
return 0;
}      
