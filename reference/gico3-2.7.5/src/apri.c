#include <stdio.h>
#include <math.h>
#include "sofam.h"
#include "header.h"
#include "common.h"
double UTC2MJD(time_t epoch) {double mjd0,mjd;struct tm tm;gmtime_r(&epoch,&tm);iauCal2jd(tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,&mjd0,&mjd);return mjd+(epoch%86400)/86400.0;}
void  scale(double v[3][3],double r){int ix,iy;for(iy=0;iy<3;iy++)for(ix=0;ix<3;ix++) v[iy][ix]=v[iy][ix]*r;}


// 黄経 mean ecliptic longitude
// by kaijou-hoantyo- 2006 PP439
// rerun mean obliquity of sun by Simon et al 1994 kaijo-hoantyoP439
double MEL(double jd){
  // JY = julianyear
  double jy = jd/36525.0;
  double mel= 280.4664573+0.98564735800624*jd+0.0003032*jy*jy;
 return (mel*M_PI/180.0); // [rad]

}
// rerun mean obliquity of sun by Lieske et al 1977 kaijo-hoantyoP439
double Epsilon(double jd){
  double jy = jd/36525.0;
  double ep = 23.439291+jy*(-0.01300417+jy*(-0.000000164+0.0000005036*jy));
  return ep/180.0*M_PI;
 } 
// rerun Kinjituten-oukei of sun by Simon et al 1994 kaijo-hoantyoP439
double Omega(double jd){
  double jy = jd/36525.0;
  double omega = 357.529109+ jd*0.9856002818 - jy*jy*0.0001537;
  return omega/180.0*M_PI;
}
double Eccentricity(double jd){
  double jy = jd/36525.0;
  double e =  0.01670863-jy*(0.000042037+0.000000127*jy);
  return e/180.0*M_PI;
}
void Aberration(double A[3],double *B, double jd){
  // aberration accourding to tentaino ichi keisan by nagasawa takumi
  //jd = Juli an day;
  double kappa=9.936491e-5;
  double ek = Eccentricity(jd)*kappa;
  double L =   A[0], M=A[1], N = A[2];
  double mel = MEL(jd), eps = Epsilon(jd), omg = Omega(jd);
  double t[6];
   t[0]= sin(mel), t[1]=cos(mel)*cos(eps), t[2]=cos(mel)*sin(eps),t[3]=sin(omg),t[4]=cos(omg)*cos(eps),t[5]=cos(omg)*sin(eps);
  double w[3][3];
  w[0][0] = (1-L*L), w[0][1] = -L*M , w[0][2] = -L*N ;
  w[1][0] = -L*M   , w[1][1] = 1-M*M, w[1][2] = -M*N ;
  w[2][0] = -L*N   , w[2][1] = -M*N , w[2][2] = 1-N*N;

  B[0]=L+kappa*(w[0][0]*t[0] -w[0][1]*t[1] -w[0][2]*t[2])+ek*(-w[0][0]*t[3] +w[0][1]*t[4] + w[0][2]*t[5]);
  B[1]=M+kappa*(w[1][0]*t[0] -w[1][1]*t[1] -w[1][2]*t[2])+ek*(-w[1][0]*t[3] +w[1][1]*t[4] + w[1][2]*t[5]);
  B[2]=N+kappa*(w[2][0]*t[0] -w[2][1]*t[1] -w[2][2]*t[2])+ek*(-w[2][0]*t[3] +w[2][1]*t[4] + w[2][2]*t[5]);
  double norm=sqrt(pow(B[0],2)+pow(B[1],2)+pow(B[2],2));
  int i;
for ( i=0; i<3; i++) { B[i]= B[i]/norm;  } // Normalized    
}

// adopted by kondo-san lib_apri_a.c
double Atm(double Height, double El /*[rad]*/)
{
  // Xatm[sec] caused by atmospehric delay
/* !--------------------------------------------------------------
   ! **** ATM model is CHAO                                      !
   !      Zenith excess path is caluculated standard model       !
   !           Hight=0   P=1013.25 mb     15\DFC                   !
   !           P=P0-0.115*Height                                  !
   !           WV PRESS = 8.5                                    !
   !           ZP=0.002277*(P+37.5)=2.393-2.62E-4*Height (m)      !
   !             =7.98E-9-8.74E-13*Height (sec)                   !
   !-------------------------------------------------------------*/
    double Zp;
    Zp=7.98e-9 - 8.74e-13*Height;
    return (Zp/(sin(El)+0.00143/(tan(El)+0.0445)));
}

// adopted by kondo-san lib_apri_a.c
void Xyz_to_llh(Station xyz, double *Longi, 
                double *Latit, double *Height)
{
      double Xlamd, Epsp, Epsh, d1, d2, d3;
      double Rearth, Eflat, Esquar;
      double P,H, Xn, X, Y, Z, Dx, Dy, Dz, Xm, Dp, Dh;

      //printf("%f %f %f \n",Xyz[0],Xyz[1],Xyz[2]);
      d1=xyz.pos_x;
      d2=xyz.pos_y;
      d3=xyz.pos_z;

      Xlamd=atan2(d2,d1);
      Epsp=1.0E-6;
      Epsh=1.0E-4;
      /*
   !----------------------------
   !   VLBI coordinate system  !
   !----------------------------
   !Definition of Earth ellipsoid : GRS80(ITRF) and WGS84(GPS)
      */
      Rearth=6378137.0;
      //Eflat=1.0/298.257;
      Eflat=1.0/298.257222101;   // GRS80
      //Eflat=1.0/298.257223563;   // WGS84
      Esquar=(2.0-Eflat)*Eflat;
      /*
   !----------------------------------------
   !      INITIAL PARAMETER CALCULATED     !
   !----------------------------------------
      */

      P=atan(d3/sqrt(d1*d1+d2*d2)*(1.0-Esquar));
      H=0.0;
      /*     Iteration      */

 Loop_llh:             //  Loop iteration
         d1=sin(P);
         Xn=Rearth/sqrt(1.0-Esquar*pow(d1,2));
         X=(Xn+H)*cos(P)*cos(Xlamd);
         Y=(Xn+H)*cos(P)*sin(Xlamd);
         Z=(Xn*(1.0-Esquar)+H)*d1;
         Dx=xyz.pos_x-X;
         Dy=xyz.pos_y-Y;
         Dz=xyz.pos_z-Z;
         d2=1.0-Esquar*pow(d1,2);
         d3=pow(d2,1.5);
         Xm=Rearth*(1.0-Esquar)/d3;
         Dp=(-(Dx*cos(Xlamd)+Dy*sin(Xlamd))*sin(P)+Dz*cos(P))/(Xm+H);
         Dh=(Dx*cos(Xlamd)+Dy*sin(Xlamd)*cos(P)+Dz*sin(P));
         if (fabs(Dp) < Epsp && fabs(Dh) < Epsh) {
             goto L3626;
         }
         P=P+Dp;
         H=H+Dh;
         goto Loop_llh;

L3626:
      *Latit=P;
      *Longi=Xlamd;
      *Height=H;
}


void xyz2azel(double longitude,double latitude /*[rad]*/, Source source , time_t epoch, double msec,double* az , double* el ){
  double phs,mjd,omega,NPB[3][3],v[5][3][3],star[3],x[3],z[3]; double ut1utc=0;
  
  omega=2*M_PI/86400.0*1.002737909350795;
  mjd=UTC2MJD(epoch)+(ut1utc+msec)/86400.0; //sign of UT1UTC is no check !!!!
  // (2400000.5 is MJD method) iaupnm06a returnd bias-precession-nutation matrix NPB
  iauPnm06a(2400000.5,mjd,NPB);
  // Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.
  // Returned     double        Greenwich apparent sidereal time (radians)
  phs=iauGst06(2400000.5,mjd,2400000.5,mjd,NPB);
  //  cout<<"GAST: "<<phs<<"[rad] "<<phs*180/M_PI<<"[deg]"<<endl;
  star[0]=cos(source.right_ascension)*cos(source.declination);
  star[1]=sin(source.right_ascension)*cos(source.declination);
  star[2]=                            sin(source.declination);
  
  // iauIr( initialize 3dim vector ),iauRz(rotate matirix about Z-axis) rotates 3dim vector with Z-axis
  iauIr(v[0]);iauRz(-phs+0*M_PI/2,v[0]);v[0][2][2]=1;
  iauRxp(NPB,star,x);

  phs += longitude;
  iauIr(v[0]);iauRz(phs,v[0]); // rotate + direction not minus direction
  iauRxp(v[0],x,z);    
  iauIr(v[0]);iauRy(M_PI*0.5-latitude,v[0]); // rotate + direction not minus direction
  iauRxp(v[0],z,x);    

  *el = asin(x[2]);

  double tmp = atan2(x[1],x[0]); 
  
  if(tmp<0){ *az = -1.0*tmp ;}else{ *az = 2*M_PI- tmp;}

}
double atm_delay(double height, double longitude,double latitude,Source source,time_t epoch ,double msec ){

  double az, el;
  xyz2azel(longitude,latitude,source,epoch,msec ,&az,&el);

  return   Atm(height,el);
}


Clock apri_old(Station station,Source source,time_t epoch,double msec)
{
  char time_code[32];struct tm tm;double julian,phs,mjd,omega,NPB[3][3],v[5][3][3],pos[3],star[3],x[3],y[3],z[3];Clock gdata; double ut1utc=0.0; //-0.0720;
  omega=2*M_PI/86400.0*1.002737909350795;

  mjd=UTC2MJD(epoch)+(ut1utc+msec)/86400.0; //sign of UT1UTC is no check !!!!
  julian = (mjd+2400000.5-2451545.0); // change MJD to julian day

  // (2400000.5 is MJD method) iaupnm06a returnd bias-precession-nutation matrix NPB
  iauPnm06a(2400000.5,mjd,NPB);

  // Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.
  // Returned     double        Greenwich apparent sidereal time (radians)
  // phs=iauGst06a(2400000.5,mjd,2400000.5,mjd);
  phs=iauGst06(2400000.5,mjd,2400000.5,mjd,NPB);

  pos[0]=station.pos_x;star[0]=cos(source.right_ascension)*cos(source.declination);
  pos[1]=station.pos_y;star[1]=sin(source.right_ascension)*cos(source.declination);
  pos[2]=station.pos_z;star[2]=                            sin(source.declination);
  // iauIr initialization , Z axis rotation, bibun keisu initialization and scale 
  iauIr(v[0]);iauRz(phs-0*M_PI/2,v[0]); v[0][2][2]=1; scale(v[0],pow(-omega,0));
  iauIr(v[1]);iauRz(phs-1*M_PI/2,v[1]); v[1][2][2]=0; scale(v[1],pow(-omega,1));
  iauIr(v[2]);iauRz(phs-2*M_PI/2,v[2]); v[2][2][2]=0; scale(v[2],pow(-omega,2));
  iauIr(v[3]);iauRz(phs-3*M_PI/2,v[3]); v[3][2][2]=0; scale(v[3],pow(-omega,3));
  iauIr(v[4]);iauRz(phs-4*M_PI/2,v[4]); v[4][2][2]=0; scale(v[4],pow(-omega,4));
  // x = NPB*star, star is 3dimension vector for geocenter
  iauRxp(NPB,star,x);  
  // transform aberration to the star vector
  Aberration(x,y,julian); 

  //z = V*Y vecotr
  iauRxp(v[0],y,z);gdata.delay = -(z[0]*pos[0]+z[1]*pos[1]+z[2]*pos[2])/299792458.0;
  iauRxp(v[1],y,z);gdata.rate  = -(z[0]*pos[0]+z[1]*pos[1]+z[2]*pos[2])/299792458.0;
  iauRxp(v[2],y,z);gdata.acel  = -(z[0]*pos[0]+z[1]*pos[1]+z[2]*pos[2])/299792458.0;
  iauRxp(v[3],y,z);gdata.jerk  = -(z[0]*pos[0]+z[1]*pos[1]+z[2]*pos[2])/299792458.0;
  iauRxp(v[4],y,z);gdata.snap  = -(z[0]*pos[0]+z[1]*pos[1]+z[2]*pos[2])/299792458.0;

  gmtime_r(&epoch,&tm);strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",&tm);

  return gdata; // Clock;

}

Clock apri(Station station,Source source,time_t epoch,double msec){

  double xlon,xlat,height,xdelay,atm; Clock tmp,gdata;
  Xyz_to_llh(station,&xlon,&xlat,&height);

  double az,el;
  //  xyz2azel(xlon,xlat,source,epoch,&az,&el);
  //  xdelay = ATM(height,el);
  //  printf("az el h, %f,%f,%f\n",az*180/M_PI, el*180/M_PI,xheight);
  double t0,t1,t2,t3,t4;
  tmp = apri_old( station, source, epoch-2, msec);   xyz2azel(xlon,xlat,source,epoch-2,msec,&az,&el);xdelay = Atm( height,el );  t0 = tmp.delay+xdelay;
  tmp = apri_old( station, source, epoch-1, msec);   xyz2azel(xlon,xlat,source,epoch-1,msec,&az,&el);xdelay = Atm( height,el );  t1 = tmp.delay+xdelay;
  tmp = apri_old( station, source, epoch+0, msec);   xyz2azel(xlon,xlat,source,epoch+0,msec,&az,&el);xdelay = Atm( height,el );  t2 = tmp.delay+xdelay;  atm = xdelay;
  tmp = apri_old( station, source, epoch+1, msec);   xyz2azel(xlon,xlat,source,epoch+1,msec,&az,&el);xdelay = Atm( height,el );  t3 = tmp.delay+xdelay;
  tmp = apri_old( station, source, epoch+2, msec);   xyz2azel(xlon,xlat,source,epoch+2,msec,&az,&el);xdelay = Atm( height,el );  t4 = tmp.delay+xdelay; 
 //if(message>=2)  printf("atm delay, %0.15f, %0.15f\n",xdelay,atm);
  gdata.delay =    t2+(12*(t1+t3-2*t2)-3*(t0+t4-2*t2))/35 ; 
  gdata.rate  =   (t0-t4-8*(t1-t3))/12;
  gdata.acel  =   (2*(t0-t2+t4)-t1-t3)/7;
  gdata.jerk  =   (2*(t1-t3)+t4-t0)/2;
  gdata.snap  =   gdata.jerk*gdata.rate;
  gdata.sec=epoch; gdata.nsec=msec*1e+9; //printf("msec: %f",msec);
if(message>=2) fprintf(stdout,"GICO3 : [level-2] apri[%8s@ %+.7e,%+.7e,%+.7e,%+.7e,%+.7e]\n",station.name,gdata.delay,gdata.rate,gdata.acel,gdata.jerk,gdata.snap);
// if(message>=2) printf(stdout,"GICO3 : [level-2] apri[%8s@ %+.7e,%+.7e,%+.7e]\n",station.name,gdata.delay,gdata.rate,atm);

  return gdata;
}

Clock apri_org(Station station,Source source,time_t epoch)
{
  char time_code[32];struct tm tm;double phs,mjd,omega,NPB[3][3],v[5][3][3],pos[3],star[3],x[3],y[3];Clock gdata;
  omega=2*M_PI/86400.0*1.002737909350795;mjd=UTC2MJD(epoch)+ut1utc/86400.0; //sign of UT1UTC is no check !!!!
  iauPnm06a(2400000.5,mjd,NPB);phs=iauGst06(2400000.5,mjd,2400000.5,mjd,NPB);
  pos[0]=station.pos_x;star[0]=cos(source.right_ascension)*cos(source.declination);
  pos[1]=station.pos_y;star[1]=sin(source.right_ascension)*cos(source.declination);
  pos[2]=station.pos_z;star[2]=                            sin(source.declination);
  iauIr(v[0]);iauRz(-phs+0*M_PI/2,v[0]);v[0][2][2]=1;scale(v[0],pow(-omega,0));
  iauIr(v[1]);iauRz(-phs+1*M_PI/2,v[1]);v[1][2][2]=0;scale(v[1],pow(-omega,1));
  iauIr(v[2]);iauRz(-phs+2*M_PI/2,v[2]);v[2][2][2]=0;scale(v[2],pow(-omega,2));
  iauIr(v[3]);iauRz(-phs+3*M_PI/2,v[3]);v[3][2][2]=0;scale(v[3],pow(-omega,3));
  iauIr(v[4]);iauRz(-phs+4*M_PI/2,v[4]);v[4][2][2]=0;scale(v[4],pow(-omega,4));
  iauRxp(NPB,star,x);gdata.sec=epoch;gdata.nsec=0;
  iauRxp(v[0],pos,y);gdata.delay=-(x[0]*y[0]+x[1]*y[1]+x[2]*y[2])/299792458.0;
  iauRxp(v[1],pos,y);gdata.rate =-(x[0]*y[0]+x[1]*y[1]+x[2]*y[2])/299792458.0;
  iauRxp(v[2],pos,y);gdata.acel =-(x[0]*y[0]+x[1]*y[1]+x[2]*y[2])/299792458.0;
  iauRxp(v[3],pos,y);gdata.jerk =-(x[0]*y[0]+x[1]*y[1]+x[2]*y[2])/299792458.0;
  iauRxp(v[4],pos,y);gdata.snap =-(x[0]*y[0]+x[1]*y[1]+x[2]*y[2])/299792458.0;
  
  gmtime_r(&epoch,&tm);strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",&tm);
  if(message>=2) fprintf(stdout,"GICO3 : [level-2] apri[%8s@%s %+.7e,%+.7e,%+.7e,%+.7e,%+.7e]\n",station.name,time_code,gdata.delay,gdata.rate,gdata.acel,gdata.jerk,gdata.snap);
  return gdata;
}
