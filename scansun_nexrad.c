/******************************************************************************/
/*This program read NEXRAD Level II files and reads the radials of the volume */
/*scan from the binary file. For each radial the noise level in the horizontal*/
/* and vertical receiver channels are extracted. The noise levels represent   */
/*the along radial, constant noise level and can be used to detect solar      */
/*interferences. The algorithm used to determine the noise levels is described*/
/* in Ivic et al, JTECH (2013), 30, p2737.                                    */
/*The format of the Level II files is described in the Interface Control      */
/*Documents (ICD), 2620002M.pdf and 2620010E.pdf.                             */
/*For each radial also the date, time, elevation, and azimuth are printed. In */
/*addition the solar elevation and azimuth are provided. Using available      */
/*metadate in the Level II file, the horizontal and vertical noise powers are */
/*converted to received solar flux at the antenna feed.                       */
/********************************************************************************/

/*Program: scansun_nexrad.c*/
/*Author: Iwan Holleman*/
/*Created: 17 February 2017*/

/*27 Feb 2017, receiver losses changed to receiver gain.*/
/*28 Feb 2017, dynamic reading of antenna area and beam width.*/
/*2 Mar 2017, bug with ICAO in volume header fixed.*/
/*10 Mar 2017, dynamic reading and use of receiver system gains.*/
/*10 Mar 2017, dynamic reading of TX pulse width / matched RX width.*/
/*18 Aug 2017, reformatted output fitting OPERA format.*/
/*18 Jul 2018, reading of Zdr Bias (calibration) is added.*/
/*27 Apr 2019, Julian-Day dependency in Solar_elev_azim function removed.*/
/*28 Apr 2019, reading of VCP data from Message 5 or 7 added.*/
/*11 May 2019, static tables added with PRFs for surveillance and Doppler modes.*/
/*19 May 2019, waveform type dependent processing and corerctions added.*/
/*24 May 2019, missing minus-sign added in correction formula for azimuth.*/
/*29 May 2019, correction of '-1' error in reading of VCP data.*/
/*8 Jun 2019, introduced noise level threshold for sun flux data.*/
/*8 Jun 2019, introduced noise level correction for sun flux data.*/
/*18 Jul 2019, corrected bug with printing longitutde and latitude.*/
/*15 Aug 2019, correction of azimuth width in batch mode (WFT4) added.*/
/*25 Aug 2019, arrays for decompression enlarged and error handling adjusted.*/
/*5 Sep 2019, check for corrupted volume header record included.*/
/*5 Sep 2019, bug fixed with non-31-type message size in data records.*/
/*6 Sep 2019, removed fixed record size, messsges are decompressed one-by-one.*/
/*15 Sep 2019, incorporated effect of Von Hann windowing on super-resolution width.*/
/*22 Sep 2019, erroneous decompression of message (header) data will cause exit.*/
/*4 Oct 2019, beam width not read from file anymore but calculated.*/
/*31 Jan 2020, big/litte endian conversion handled by explicit formulas.*/
/*6 Feb 2020, fix of BZ2-lib initializing and partial restructuring.*/
/*13 Jun 2020, reading of matched filter losses added.*/
/*20 Sep 2020, no super resolution for (WFT1) noise level data.*/
/*20 Sep 2020, fixed matched filter bandwidths for all radars.*/
/*26 Sep 2020, corrected calculation of RXgain from I0 and MF Loss.*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <bzlib.h>

/******************************************************************************/
/*Definition of standard parameters.                                          */
/******************************************************************************/

#define DEG2RAD     (0.017453293) /*Degrees to radians*/
#define RAD2DEG     (57.29578)    /*Radians to degrees*/
#define LSTR        (128)         /*Length of all strings used*/

/******************************************************************************/
/*NEXRAD level II dataformat parameters                                       */
/******************************************************************************/

#define LVOLHEAD    (24)          /*Length of Archive II Volume Header Record*/
#define NMETA       (134)         /*Number of the metadata messages*/
#define LMETA       (2432)        /*Length of the metadata messages*/
#define NDATA       (120)         /*Number of the Type-31 messages*/
#define LDATA       (32768)       /*Maximum length Type-31 messages*/
#define LTCMHEAD    (12)          /*Length of TCM Message Header*/
#define LMSGHEAD    (16)          /*Length of Message Header Data*/

/******************************************************************************/
/*Physical parameters of NEXRAD radars.                                       */
/*www.roc.noaa.gov/wsr88d/Engineering/NEXRADTechInfo.aspx                     */
/******************************************************************************/

#define CLIGHT      (299792458.0) /*Speed of light in vacuum in m/s*/
#define SUNDIAM     (0.5708)      /*Radio diameter of sun in deg*/
#define ANTDIAM     (8.5)         /*Aperture diameter of NEXRAD antennas in m*/
#define BANDW_SP    (0.714)       /*Short pulse matched filter bandwidth in MHz*/
#define BANDW_LP    (0.252)       /*Long pulse matched filter bandwidth in MHz*/
                                  /*Both numbers from Alan Free, mail 16-Sep-2020*/

/******************************************************************************/
/*Structure for Volume Coverage Pattern (VCP) data (Table XI, see Manual)     */
/******************************************************************************/

#define VCPMAX      (128) /*Maximum number of elevations in VCP.*/
struct vcp {
    float elev;           /*Elevation angle of scan in deg.*/
    unsigned char WF;     /*Waveform Type of scan, flag 0...5 (see Table XI).*/
    unsigned char SNUM;   /*Surveillance PRF Number, flag 0...8.*/
    float SPRF;           /*Surveillance PRF frequency in Hz.*/
    int SCNT;             /*Surveillance PRF pulse count per radial.*/
    unsigned char DNUM;   /*Doppler PRF Number, flag 0...8.*/
    float DPRF;           /*Doppler PRF frequency in Hz.*/
    int DCNT;             /*Doppler PRF pulse count per radial.*/
    int Nmsg;             /*Counted number of Type-31 messages at elevation.*/
};
typedef struct vcp VCP;

/******************************************************************************/
/*Pulse Repetition Frequencies from VCP definition files (Brandon Taylor).    */
/******************************************************************************/

static float SurvPRFvalues[9]={0,321.89,350.26,389.36,446.43,512.80,602.41,719.39,862.07};
static float DoppPRFvalues[9]={0,0,803.20,857.14,928.79,1013.51,1094.89,1181.10,1282.05};

/******************************************************************************/
/*Parameters for output selection.                                            */
/******************************************************************************/

#define NTHRS       (2.0)      /*Noise level threshold for solar hits in dB.*/
#define EDIFX       (5.0)      /*Maximum elevation devation from sun in deg.*/
#define ADIFX       (10.0)     /*Maximum azimuth devation from sun in deg.*/

/******************************************************************************/
/*Prototypes of local functions:                                              */
/******************************************************************************/

float ScanningLosses(float wbeam,float wazim);
float erfcc(float x);
float code2elev(int code);
void solar_elev_azim(float lon,float lat,int yyyymmdd,int hhmmss,float *elev,float *azim);
int datecalc(int yyyymmdd,int nday);

/******************************************************************************/
/*Main program:                                                               */
/******************************************************************************/
int main(int argc, char *argv[])
{
char *nexfile,*name,first=1,prthead=1,prtdata=0;
unsigned char volhead[LVOLHEAD],rechead[4],msgmeta[NMETA*LMETA],msgdata[NDATA*LDATA];
unsigned char *record,*message;
int julday,date,sec,time,recsize,datsize,bz_error,ipulse;
int msgsize,msgseg,msgtype,vcp,nelev,freq,blck,ielev,ngate,ival,n,m;
float noisethr=NTHRS,elevdif=EDIFX,azimdif=ADIFX,latrad,lonrad,RXhgain,RXvgain,factor;
float azim,elev,azimsun,elevsun,Gain,AntArea,wband,wbeam,wazim,dazim,hflux,vflux;
float hnoise,vnoise,hnoise_s,hnoise_l,hnoise0,hI0,vnoise_s,vnoise_l,vnoise0,vI0,ZdrBias;
float loss_mf_lp,loss_mf_sp,loss_mf;
VCP VCPdata[VCPMAX];
FILE *fp;
bz_stream strm;

/*Obtaining parameters from command line.*/

if (argc<2) {
   printf("Usage: %s <Nexrad Level II file> [-nthrs<x>/-delev<x>/-dazim<x>]!\n",argv[0]);
   exit(1);
}
nexfile=argv[1];
for (n=2 ; n<argc ; n++) {
   sscanf(argv[n],"-nthrs%f",&noisethr);
   sscanf(argv[n],"-delev%f",&elevdif);
   sscanf(argv[n],"-dazim%f",&azimdif);
}

/*Opening of NEXRAD Level II file.*/

printf("#Reading NEXRAD Level II file: %s\n",nexfile);
fp=fopen(nexfile,"rb");
if (fp==NULL) {
   printf("##Error: NEXRAD file could not be opened!\n");
   exit(2);
}

/*Reading of volume header message.*/

fread(volhead,1,LVOLHEAD,fp);

/*Decoding of volume header record (24 bytes).*/

julday=((volhead[12]*256+volhead[13])*256+volhead[14])*256+volhead[15];
date=datecalc(19700101,julday-1);
sec=(((volhead[16]*256+volhead[17])*256+volhead[18])*256+volhead[19])/1000;
time=10000*(sec/3600)+100*((sec/60)%60)+sec%60;

/*Printing header information.*/

printf("#Nexrad Archive II filename  : %12.12s\n",volhead);
printf("#Date of NEXRAD volume file  : %08d\n",date);
printf("#Time of NEXRAD volume file  : %06d\n",time);
printf("#ICAO name of NEXRAD radar   : %4.4s\n",volhead+20);
 
/*Reading of the compressed METADATA record.*/

fread(rechead,1,4,fp);
recsize=abs(((rechead[0]*256+rechead[1])*256+rechead[2])*256+rechead[3]);
record=(unsigned char *)malloc(recsize);
if (!record) {
   printf("##Error: record (size: %d) could not be allocated!\n",recsize);
   exit(4);
}
if (!fread(record,1,recsize,fp)) {
   printf("##Error: metadata record could not be read!\n");
   exit(5);
}

/*Decompressing all messages in the METADATA record.*/

strm.bzalloc=NULL;
strm.bzfree=NULL;
strm.opaque=NULL;
BZ2_bzDecompressInit(&strm,0,0);
strm.avail_in=recsize;
strm.next_in=(char *)record;
strm.avail_out=NMETA*LMETA;
strm.next_out=(char *)msgmeta;
if (BZ2_bzDecompress(&strm)!=BZ_STREAM_END) {
   printf("##Error: metadata record could not be decompressed!\n");
   exit(6);
}
BZ2_bzDecompressEnd(&strm);
free(record);
 
/*Reading of METADATA messages.*/

for (n=0 ; n<NMETA ; n++) {
   message=msgmeta+n*LMETA;
   message+=LTCMHEAD;
   msgtype=(int)message[3];
   msgseg=(int)(message[14]*256+message[15]);
   message+=LMSGHEAD;
   
/*Reading of calibration data: noise levels at RX input and ouput, and Zdr Bias.*/
   
   if (msgtype==3&&msgseg==1) {
      ival=((message[700]*256+message[701])*256+message[702])*256+message[703];
      hnoise_s=*(float *)&ival; 
      ival=((message[704]*256+message[705])*256+message[706])*256+message[707];
      hnoise_l=*(float *)&ival;
      ival=((message[712]*256+message[713])*256+message[714])*256+message[715];
      vnoise_s=*(float *)&ival;
      ival=((message[716]*256+message[717])*256+message[718])*256+message[719];
      vnoise_l=*(float *)&ival;
      ival=((message[764]*256+message[765])*256+message[766])*256+message[767];
      hI0=*(float *)&ival;
      ival=((message[768]*256+message[769])*256+message[770])*256+message[771];
      vI0=*(float *)&ival;
      ival=((message[800]*256+message[801])*256+message[802])*256+message[803];
      ZdrBias=*(float *)&ival;
   }
   
/*Reading characteristics of VCP (incl. short or long TX pulse flag).*/
   
   if ((msgtype==5||msgtype==7)&&msgseg==1) {
      vcp=message[4]*256+message[5];
      nelev=message[6]*256+message[7];
      if (nelev>VCPMAX) {
         printf("##Error: Number of elevations (%d) is higher than %d!\n",nelev,VCPMAX);
         exit(7);
      }
      ipulse=message[11];
      for (m=0 ; m<nelev ; m++) {
         VCPdata[m].elev=code2elev(message[2*(11+m*23)]*256+message[2*(11+m*23)+1]);
         VCPdata[m].WF=message[2*(11+m*23)+3];
         VCPdata[m].SNUM=message[2*(11+m*23)+5];
         VCPdata[m].SPRF=SurvPRFvalues[VCPdata[m].SNUM];
         VCPdata[m].SCNT=message[2*(11+m*23)+6]*256+message[2*(11+m*23)+7];
         VCPdata[m].DNUM=message[2*(11+m*23)+25];
         VCPdata[m].DPRF=DoppPRFvalues[VCPdata[m].DNUM];
         VCPdata[m].DCNT=message[2*(11+m*23)+26]*256+message[2*(11+m*23)+27];
         VCPdata[m].Nmsg=0;
      }
   }
   
/*Reading of radar frequency, matched filter losses, antenna beamwidth, and antenna gain.*/
   
   if (msgtype==18&&msgseg==1) {
      freq=((message[1092]*256+message[1093])*256+message[1094])*256+message[1095];
      ival=((message[1120]*256+message[1121])*256+message[1122])*256+message[1123];
      loss_mf_lp=*(float *)&ival;
      ival=((message[1124]*256+message[1125])*256+message[1126])*256+message[1127];
      loss_mf_sp=*(float *)&ival;
      ival=((message[1136]*256+message[1137])*256+message[1138])*256+message[1139];
      Gain=*(float *)&ival;
   }
}

/*Selection calibration noise levels (short (2) or long (4) pulse). Calculation of receiver */
/*system gains and matched RX bandwidth.*/

if (ipulse==2) {
   hnoise0=hnoise_s;
   vnoise0=vnoise_s;
   wband=BANDW_SP;
   loss_mf=loss_mf_sp;
}
else {
   hnoise0=hnoise_l;
   vnoise0=vnoise_l;
   wband=BANDW_LP;
   loss_mf=loss_mf_lp;
}
        
/*Calculation of RX gain, using noise level, I0 and Matched Filter Losses.*/

RXhgain=hnoise0-(hI0+loss_mf);
RXvgain=vnoise0-(vI0+loss_mf);

/*Calculation of effective antenna area (m2) and antenna beam width (deg).*/

AntArea=CLIGHT/(freq*1e6);
AntArea=AntArea*AntArea;
AntArea*=exp(0.1*log(10)*Gain)/(4*M_PI);
wbeam=1.27*RAD2DEG*CLIGHT/(freq*1e6);
wbeam/=ANTDIAM;

/*Printing of all relevant metadata to stdout.*/

printf("#Volume Coverage Pattern     : %d\n",vcp);
printf("#Number of elevations in VCP : %d\n",nelev);
printf("#Short (2) or long (4) pulse : %d\n",ipulse);
printf("#Matched RX band width (MHz) : %5.3f\n",wband);
printf("#Matched RX filter loss (dB) : %6.3f\n",loss_mf);
printf("#Horizontal noise level (dBm): %6.3f\n",hnoise0);
printf("#Vertical noise level (dBm)  : %6.3f\n",vnoise0);
printf("#Horizontal RX gain (dB)     : %6.3f\n",RXhgain);
printf("#Vertical RX gain (dB)       : %6.3f\n",RXvgain);
printf("#Beam width of antenna (deg) : %5.3f\n",wbeam);
printf("#Effective antenna area (m2) : %5.2f\n",AntArea);
printf("#Calibration of Zdr bias (dB): %6.3f\n",ZdrBias);

/*Reading of all Type-31 messages in the compressed DATA records.*/

while (fread(rechead,1,4,fp)) {
   recsize=abs(((rechead[0]*256+rechead[1])*256+rechead[2])*256+rechead[3]);
   record=(unsigned char *)malloc(recsize);
   if (!record) {
      printf("##Error: record (size: %d) could not be allocated!\n",recsize);
      exit(8);
   }
   if (!fread(record,1,recsize,fp)) {
      printf("##Error: data record could not be read!\n");
      exit(9);
   }

/*Decompressing all messages in the DATA record.*/

   strm.bzalloc=NULL;
   strm.bzfree=NULL;
   strm.opaque=NULL;
   BZ2_bzDecompressInit(&strm,0,0);
   strm.avail_in=recsize;
   strm.next_in=(char *)record;
   strm.avail_out=NDATA*LDATA;
   strm.next_out=(char *)msgdata;
   bz_error=BZ2_bzDecompress(&strm);
   BZ2_bzDecompressEnd(&strm);
   free(record);
   if (bz_error==BZ_OK) printf("##Warning: record could not be decompressed fully!\n");
   if (bz_error!=BZ_STREAM_END) {
      printf("##Warning: record could not be decompressed!\n");
      continue;
   }
   datsize=NDATA*LDATA-(int)strm.avail_out;
   
/*Reading all Digital Radar Data messages (type 31) in the record.*/

   message=msgdata;
   while (message-msgdata<datsize) {
      msgtype=message[LTCMHEAD+3];
      if (msgtype!=31) {
         message+=LMETA;
         continue;
      }
      message+=LTCMHEAD;
      msgsize=2*(message[0]*256+message[1]);
      message+=LMSGHEAD;
   
/*Reading data message header: date and time*/

      sec=(((message[4]*256+message[5])*256+message[6])*256+message[7])/1000;
      time=10000*(sec/3600)+100*((sec/60)%60)+sec%60;
      julday=message[8]*256+message[9];
      date=datecalc(19700101,julday-1);

/*Reading message header: azimuth, azimuthal averaging, elevation, and elevation number*/

      ival=((message[12]*256+message[13])*256+message[14])*256+message[15];
      azim=*(float *)&ival;
      ielev=(int)message[22];       /*Elevation number within volume scan, 1...25.*/
      ielev--;                      /*Change elevation number range, 0...24.*/
      VCPdata[ielev].Nmsg++;
      ival=((message[24]*256+message[25])*256+message[26])*256+message[27];
      elev=*(float *)&ival;
 
/*Reading Volume Data Constant Type: longitude, latitude*/

      blck=((message[32]*256+message[33])*256+message[34])*256+message[35];
      name=(char *)(message+blck);
      if (strncmp(name,"RVOL",4)==0) {
         ival=((message[blck+8]*256+message[blck+9])*256+message[blck+10])*256+message[blck+11];
         latrad=*(float *)&ival;
         ival=((message[blck+12]*256+message[blck+13])*256+message[blck+14])*256+message[blck+15];
         lonrad=*(float *)&ival;
      }
      else latrad=lonrad=0;

/*Printing geo-information of radar (only once).*/

      if (first) {
         printf("#Latitude of NEXRAD radar    : %8.3f\n",latrad);
         printf("#Longitude of NEXRAD radar   : %8.3f\n",lonrad);
         first=0;
      }
      
/*Reading Radial Data Constant Type: h/v noise levels (!)*/

      blck=((message[40]*256+message[41])*256+message[42])*256+message[43];
      name=(char *)(message+blck);
      if (strncmp(name,"RRAD",4)==0) {
         ival=((message[blck+8]*256+message[blck+9])*256+message[blck+10])*256+message[blck+11];	
         hnoise=*(float *)&ival;
	      ival=((message[blck+12]*256+message[blck+13])*256+message[blck+14])*256+message[blck+15];
         vnoise=*(float *)&ival;
      }
      else {
         hnoise=hnoise0;
         vnoise=vnoise0;
      }

/*Reading Data Block (reflectivity): number of range gates*/

      blck=((message[44]*256+message[45])*256+message[46])*256+message[47];
      name=(char *)(message+blck);
      if (strncmp(name,"DREF",4)==0) {
         ngate=message[blck+8]*256+message[blck+9];
      }
      else ngate=0;

/*Calculating solar elevation and azimuth from radar location and date/time.*/

      solar_elev_azim(lonrad,latrad,date,time,&elevsun,&azimsun);

/*Selecting sun hits via strength of (noise) signal and angular proximity to calculated sun position.*/

      if (hnoise > hnoise0+noisethr && vnoise > vnoise0+noisethr &&
                        fabs(elev-elevsun)<elevdif && fabs(azim-azimsun)<azimdif) {

/*Depending on waveform type flags are set and corrections are made.*/
/*WFT1 is mostly recorded in super resolution (0.5 instead of 1.0 deg).*/
/*However, the Von Hann window is not applied to the noise level data and thus */
/*these data have a normal resolution (1 deg).*/
          
         switch (VCPdata[ielev].WF) {
            case 1: case 3:
               dazim=0.0;
               wazim=1.0;
               prtdata=1;
               break;
            case 4:
               factor=VCPdata[ielev].DCNT*VCPdata[ielev].SPRF+VCPdata[ielev].SCNT*VCPdata[ielev].DPRF;
               if (factor>0.0) factor=(VCPdata[ielev].SCNT*VCPdata[ielev].DPRF)/factor;
               else factor=1.0;
               dazim=1.0*(factor-1.0)/2.0;
               wazim=1.0*factor;
               prtdata=1;
               break;
            case 2: case 5: default:
               dazim=0.0;
               wazim=1.0;
               prtdata=0;
               break;
         }

/*Performing corrections (noise level, bandwidth, etc) and printing data to stdout.*/

         if (prthead) {
            printf("#   Date   Time Elevatn Azimuth ElevSun AzimSun Ngate H_dBsfu H_Noise H_sd V_dBsfu V_sd Refl ZDR WnPnRx.xx\n");
            prthead=0;
         }
         if (prtdata) {
            if (hnoise>hnoise0) hflux=10.0*log10(exp(0.1*log(10)*hnoise)-exp(0.1*log(10)*hnoise0));
            if (vnoise>vnoise0) vflux=10.0*log10(exp(0.1*log(10)*vnoise)-exp(0.1*log(10)*vnoise0));
            hflux+=130-10*log10(wband)-RXhgain-10*log10(AntArea)+ScanningLosses(wbeam,wazim)+10*log10(2);
            vflux+=130-10*log10(wband)-RXvgain-10*log10(AntArea)+ScanningLosses(wbeam,wazim)+10*log10(2);
            printf("%08d %06d %7.3f %7.3f %7.3f %7.3f %5d %7.3f %7.2f  0.0 %7.3f  0.0   TH  TV W%1dP%1dR%4.2f\n",
                            date,time,elev,azim+dazim,elevsun,azimsun,ngate,hflux,hnoise,vflux,VCPdata[ielev].WF,
                            ipulse,wazim);
         }
      }
      message+=msgsize-LMSGHEAD;
   }
}

/*Printing details of Volume Coveage Patterrn and number of Type 31 messages per elevation.*/

printf("#NVCP Elev WaveFT SRVNUM SRVPRF SRVCNT DOPNUM DOPPRF DOPCNT Nmsg\n");
for (n=0 ; n<nelev ; n++) {
   printf("#%-4d %4.1f %6d %6d %6.1f %6d %6d %6.1f %6d %4d\n",n+1,VCPdata[n].elev,VCPdata[n].WF,
          VCPdata[n].SNUM,VCPdata[n].SPRF,VCPdata[n].SCNT,VCPdata[n].DNUM,VCPdata[n].DPRF,VCPdata[n].DCNT,VCPdata[n].Nmsg);
}

/*Closing of NEXRAD file.*/

fclose(fp);

/*End of main.*/

exit(0);
}

/******************************************************************************/
/*LOCAL FUNCTIONS                                                             */
/******************************************************************************/

/******************************************************************************/
/*This function calculates the losses due to the imperfect overlap between    */
/*the sun and the antenna beam pattern due to different angular sizes and the */
/*losses due to the fact that the antenna is scanning in azimuthal direction. */
/*The formulas are derived in Holleman et al, JTECH 27 (2010), p159.          */
/******************************************************************************/
float ScanningLosses(float wbeam,float wazim)
{
float xval,L0,L1;

/*Calculation of L0, the overlap factor between sun and antenna (Eq. 13).*/

xval=SUNDIAM/wbeam;
xval=log(2)*xval*xval;
L0=(1-exp(-xval))/xval;

/*Integration over azimuth due to scanning (Eq. 14).*/

L1=0.5*L0*sqrt(M_PI/log(2))*wbeam/wazim;
L1*=1-erfcc(sqrt(log(2))*wazim/wbeam);

/*Returning attenuation due to antenna-sun overlap and azimuthal scanning.*/

return -10*log10(L1);
}

/******************************************************************************/
/*This function returns the complementary error function erfc(x) with         */
/*fractional error everywhere less than 1.2e-7.                               */
/*For details see Numerical Recipes, section 6.2                              */
/******************************************************************************/
float erfcc(float x)
{
float t,z,ans;
z=fabs(x);
t=1.0/(1.0+0.5*z);
ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
	t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
	t*(-0.82215223+t*0.17087277)))))))));
return (x>=0.0)?(ans):(2.0-ans);
}

/******************************************************************************/
/*Conversion of code-value for elevation to actual float-value. The format is */
/*defined in Table III-A of the Interface Control Dicument for the RDA/RPG.   */
/******************************************************************************/
float code2elev(int code)
{
float elev=0.0;

/*Conversion according to Table III-A*/

if (code&0x8000) elev+=180.0;
if (code&0x4000) elev+=90.0;
if (code&0x2000) elev+=45.0;
if (code&0x1000) elev+=22.5;
if (code&0x0800) elev+=11.25;
if (code&0x0400) elev+=5.625;
if (code&0x0200) elev+=2.8125;
if (code&0x0100) elev+=1.40625;
if (code&0x0080) elev+=0.70313;
if (code&0x0040) elev+=0.35156;
if (code&0x0020) elev+=0.17578;
if (code&0x0010) elev+=0.08789;
if (code&0x0008) elev+=0.043945;

/*Returing the result.*/

return elev;
}

/******************************************************************************/
/*This function calculates the solar elevation and azimuth using the          */
/*geographical position, date, and time. The equations and constants are taken*/
/*from the WMO guide on Meteorological Instruments and Methods of Observations*/
/*(CIMO, WMO no. 8), annex 7.D. The equations have been slightly modified and */
/*extended to include the calculation of both the sine and cosine of the      */
/*azimuth.                                                                    */
/******************************************************************************/
void solar_elev_azim(float lon,float lat,int yyyymmdd,int hhmmss,float *elev,float *azim)
{
float MeanLon,MeanAnom,EclipLon,Obliquity,RightAsc,Declinat;
float GMST,angleH;
double days,hour;
struct tm timestruct;
time_t time1,time0;

/*Conversion of lon,lat.*/

lon*=DEG2RAD;
lat*=DEG2RAD;

/*Initializing of time structure with actual date and time.*/

timestruct.tm_year=(yyyymmdd/10000)-1900;
timestruct.tm_mon=((yyyymmdd/100)%100)-1;
timestruct.tm_mday=(yyyymmdd%100);
timestruct.tm_hour=(int)(hhmmss/10000);
timestruct.tm_min=(int)((hhmmss/100)%100);
timestruct.tm_sec=(int)(hhmmss%100);
timestruct.tm_isdst=0;
time1=mktime(&timestruct);

/*Initializing of time structure with reference (noon 1 Jan 2000).*/

timestruct.tm_year=100;
timestruct.tm_mon=0;
timestruct.tm_mday=1;
timestruct.tm_hour=12;
timestruct.tm_min=0;
timestruct.tm_sec=0;
timestruct.tm_isdst=0;
time0=mktime(&timestruct);

/*Calculation of fractional days.*/

days=difftime(time1,time0)/(24*3600);

/*Calculation of eclips coordinates.*/

MeanLon=280.460+0.9856474*days;
MeanAnom=357.528+0.9856003*days;
EclipLon=MeanLon+1.915*sin(MeanAnom*DEG2RAD)+0.020*sin(2*MeanAnom*DEG2RAD);
EclipLon*=DEG2RAD;
Obliquity=23.439-0.0000004*days;
Obliquity*=DEG2RAD;

/*Calculation of the celestial coordinates of the sun.*/

RightAsc=atan2(cos(Obliquity)*sin(EclipLon),cos(EclipLon));
Declinat=asin(sin(Obliquity)*sin(EclipLon));

/*Calculation of current, local hour angle.*/

hour=(double)(hhmmss/10000)+((hhmmss/100)%100)/60.0+(hhmmss%100)/3600.0;
GMST=hour+6.697375+0.0657098242*days;
angleH=GMST*15*DEG2RAD+lon-RightAsc;

/*Calculation of elevation and azimuth.*/

*elev=asin(sin(Declinat)*sin(lat)+cos(Declinat)*cos(lat)*cos(angleH));
*azim=atan2(-sin(angleH),tan(Declinat)*cos(lat)-cos(angleH)*sin(lat));

/*Scaling and shifting of values.*/

(*elev)*=RAD2DEG;
(*azim)*=RAD2DEG;
if ((*azim)<0) (*azim)+=360;
}

/******************************************************************************/
/*Calculation of new date which is offset 'nday' from input date 'yyyymmdd'.  */
/******************************************************************************/
int datecalc(int yyyymmdd,int nday)
{
struct tm timestruct;
int newdate;

/*Initializing time structure.*/

timestruct.tm_year=(yyyymmdd/10000)-1900;
timestruct.tm_mon=((yyyymmdd/100)%100)-1;
timestruct.tm_mday=(yyyymmdd%100)+nday;
timestruct.tm_hour=0;
timestruct.tm_min=0;
timestruct.tm_sec=0;

/*Calculation of new date.*/

timegm(&timestruct);
newdate=(timestruct.tm_year+1900)*10000;
newdate+=(timestruct.tm_mon+1)*100;
newdate+=timestruct.tm_mday;

/*Return new date.*/

return newdate;
}
