#include "atommqdtfinal.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
extern "C" {void dgeev_(char*,char*,int*,double*,int*,double*,double*,double*,int*,double*,int*,double*,int*,int*);}
int jindex(double j, int l, double s){
  if (j == l+ s && s==1){
    return(3);
  }
  else if (j== l && s==1){
    return(2);
  }
  else if (j== l-s && s==1){
    return(1);
  }
  else if (j== l && s==0){
    return(0);
  }
  else return(-1);
         };

double omega12(double m11,double m12,double m21,double m22,double j1primed,double j2primed){
  //originally for theta=0.0;
  double result=0.0;
  double k1,k2,j1,j2,totalk,doublek1,doubletotalk;
  k1=k2=1.0;
  j1=j2=2.0;
  totalk=k1+k2;
  doublek1=2.0*k1;
  doubletotalk=2*totalk;
  double cg1,cg2,cg3,cg4,cg5,fact1,fact2;
  cg1=cg2=cg3=cg4=cg5=fact1=fact2=0.0;
  factorial_(&doublek1,&fact1);
  factorial_(&doubletotalk,&fact2);
  double angular0=sqrt(4.0*pi*fact2/(2*k1+ 2*k2 + 1))/fact1;
  for (double m1primed=(-(j1primed)); m1primed<(j1primed+1); m1primed++){
    for (double m2primed=(-(j2primed)); m2primed<(j2primed+1); m2primed++){
      if (fabs(m1primed+m2primed) != fabs(m11 + m12)) continue;
      double angular3=0.0;
      for (double p1=-1.0;p1<2.0;p1++){
	for (double p2=-1.0;p2<2.0;p2++){
	  // for (double p= -2;p< 3;p++){
	  double p=p1+p2;
	  clebschgordan_(&k1,&p1,&k2,&p2,&totalk,&p,&cg1);
	  clebschgordan_(&j1,&m11,&k1,&p1,&j1primed,&m1primed,&cg2);
	  // printf("j1=%g m1=%g k1=%g p1=%g j1p=%g m1p=%g, cg1=%g\n",j1,m11,k1,p1,j1primed,m1primed,cg2);
	  clebschgordan_(&j2,&m12,&k2,&p2,&j2primed,&m2primed,&cg3);
	  // printf("j2=%g m2=%g k2=%g p2=%g j2p=%g m2p=%g, cg2=%g\n",j2,m12,k2,p2,j2primed,m2primed,cg3);
	  clebschgordan_(&j1,&m21,&k1,&p1,&j1primed,&m1primed,&cg4);
	  clebschgordan_(&j2,&m22,&k2,&p2,&j2primed,&m2primed,&cg5);
	  double costheta= 1.0;
	  double spherical= gsl_sf_legendre_sphPlm(k1+k2,(int)fabs(p), costheta);
	  // if (p < 0){
	  //   spherical= spherical*pow(-1,p);
	  // }
	  double entry=pow(spherical,2)*(4.0*pi*fact2/(2*k1+ 2*k2 + 1))*pow(cg1,2)*cg2*cg3*cg4*cg5/pow(fact1,2);//
	  // printf("m1primed=%g  m2primed=%g  entry=%g\n",m1primed,m2primed,entry);
	  result+= entry;
	  // printf("p1=%g p2=%g p=%g,entry=%g\n",p1,p2,p,entry);//,cg1=%g cg2=%g cg3=%g,cg1,cg4,cg5
	  // angular3+=cg1*cg2*cg3;//spherical*
	  // }
	}
      }
      // printf("ang3=%g m1=%g m2=%g m1p=%g m2p=%g\n",angular3,*m11,*m12,m1primed,m2primed);
    }
  }
  return(result);
}

int main(){
  const double Ry=109736.627;
  const double Is=45932.1982;
  const double Id = 60628.26;
  const double Ip = 70048.11;
  int nin,lin,na,nb;
  int sin;
  double energyin,jin, throwaway;
  ifstream defectsfile ("srenergies.csv");
  time_t tic,toc;
  double ***energyinitial;
  energyinitial= new double **[101];
  for (int i=0;i<101;i++){
    energyinitial[i]= new double *[4];
    for (int j=0;j<4;j++){
      energyinitial[i][j]= new double [4];
    }
  }
  //--------------------------------------------------------------------------------
  //Read in energy levels.
  printf("Loading quantum defects.\n");
  if (defectsfile.is_open())
    {
      while (defectsfile.good())
	{
	  nin=sin=lin=0;
	  energyin=jin=0.0;
	  defectsfile >> nin >> lin >> sin >> jin >> energyin;
	  int j= jindex(jin,lin,sin);
	  energyinitial[nin][lin][j]= (energyin-45932.1982)/219475.0 ;
	}
    }

  //--------------------------------------------------------------------------------
  //Read in MQDT state fractions.
  FILE * statefractionfile= fopen("singletstatefractions6channel.csv","r");//
  // FILE * statefractionfile= fopen("tripletstatefractions6channel.csv","r");//
  double *channel1fractions= new double[76];
  double *channel2fractions= new double[76];
  double *channel3fractions= new double[76];
  double *channel4fractions= new double[76];
  double *channel5fractions= new double[76];
  double *channel6fractions= new double[76];

  double *channelenergy= new double[76];
  int scanresult;
  for (int i=0;i<76;i++){
    double E;
    int n;
    scanresult=fscanf(statefractionfile,"%lf   %lf   %lf   %lf   %lf   %lf   %lf",&channelenergy[i],&channel1fractions[i],&channel2fractions[i],&channel3fractions[i],&channel4fractions[i],&channel5fractions[i],&channel6fractions[i]);
  }
  fclose(statefractionfile);

  FILE * pstatefractionfile= fopen("pstatefractionsextended.csv","r");//
  double *pchannel1fractions= new double[76];
  double *pchannel2fractions= new double[76];
  double *pchannelenergy= new double[76];
  for (int i=0;i<76;i++){
    scanresult=fscanf(pstatefractionfile,"%lf   %lf   %lf",&pchannelenergy[i],&pchannel1fractions[i],&pchannel2fractions[i]);
  }
  fclose(pstatefractionfile);

  FILE * fstatefractionfile= fopen("fstatefractionsextended.csv","r");//
  double *fchannel1fractions= new double[77];
  double *fchannel2fractions= new double[77];
  double *fchannelenergy= new double[77];
  for (int i=0;i<77;i++){
    scanresult=fscanf(fstatefractionfile,"%lf   %lf   %lf",&fchannelenergy[i],&fchannel1fractions[i],&fchannel2fractions[i]);
  }
  fclose(fstatefractionfile);

  FILE * tripletf3file= fopen("tripletf3statefractions.csv","r");
  double *tripletf3channel1fractions= new double[80];
  double *tripletf3channel2fractions= new double[80];
  double *tripletf3channel3fractions= new double[80];
  double *tripletf3channelenergy= new double[80];
  for (int i=0;i<80;i++){
    scanresult=fscanf(tripletf3file,"%lf   %lf   %lf   %lf",&tripletf3channelenergy[i],&tripletf3channel1fractions[i],&tripletf3channel2fractions[i],&tripletf3channel3fractions[i]);
  }
  fclose(tripletf3file);
  FILE * tripletf2file= fopen("tripletf2statefractions.csv","r");
  double *tripletf2channel1fractions= new double[80];
  double *tripletf2channel2fractions= new double[80];
  double *tripletf2channel3fractions= new double[80];
  double *tripletf2channelenergy= new double[80];
  for (int i=0;i<80;i++){
    scanresult=fscanf(tripletf2file,"%lf   %lf   %lf   %lf",&tripletf2channelenergy[i],&tripletf2channel1fractions[i],&tripletf2channel2fractions[i],&tripletf2channel3fractions[i]);
  }
  fclose(tripletf2file);

  FILE * tripletp2file= fopen("tripletp2statefractions.csv","r");
  double *tripletp2channel1fractions= new double[80];
  double *tripletp2channel2fractions= new double[80];
  double *tripletp2channel3fractions= new double[80];
  double *tripletp2channelenergy= new double[80];
  for (int i=0;i<80;i++){
    scanresult=fscanf(tripletp2file,"%lf   %lf   %lf   %lf",&tripletp2channelenergy[i],&tripletp2channel1fractions[i],&tripletp2channel2fractions[i],&tripletp2channel3fractions[i]);
  }
  fclose(tripletp2file);

  FILE * tripletp1file= fopen("tripletp1statefractions.csv","r");
  double *tripletp1channel1fractions= new double[80];
  double *tripletp1channel2fractions= new double[80];
  double *tripletp1channel3fractions= new double[80];
  double *tripletp1channelenergy= new double[80];
  for (int i=0;i<80;i++){
    scanresult=fscanf(tripletp1file,"%lf   %lf   %lf   %lf",&tripletp1channelenergy[i],&tripletp1channel1fractions[i],&tripletp1channel2fractions[i],&tripletp1channel3fractions[i]);
  }
  fclose(tripletp1file);

  int l=2;
  double j,m,s;
  j=m=s=0.0;

  double omegaq[6][4];
  for (int i=0;i<6;i++){
    for (int j=0;j<4;j++){
      omegaq[i][j]=0.0;
    }
  }
  omegaq[0][0]=-sqrt(2.0/3.0);
  omegaq[0][1]=sqrt(1.0/3.0);
  omegaq[0][2]= sqrt(1.0/90.0);
  omegaq[1][0]=sqrt(3.0/7.0);
  omegaq[1][1]=sqrt(1.0/3.0);
  omegaq[1][2]= sqrt(2.0/5.0);
  omegaq[2][0]=1.0/sqrt(10.0);
  omegaq[2][1]=1.0/sqrt(20.0);
  omegaq[2][2]=sqrt(9.0/40.0);
  omegaq[2][3]= -sqrt(3.0/20.0);
  omegaq[3][0]=-1.0/sqrt(2.0);
  omegaq[3][1]=0.5;
  omegaq[3][2]=-0.5;
  omegaq[3][3]=-1.0/sqrt(12.0);
  omegaq[4][0]= 1.0/sqrt(15.0);
  omegaq[4][1]= sqrt(7.0/45);
  omegaq[5][0]=sqrt(8.0/21.0);
  omegaq[5][1]=sqrt(8.0/27.0);

  // FILE *c6outputfile=fopen("c6dstates.csv","w");
  // FILE *singlechannelfile=fopen("c6dstatessinglechannel.csv","w");
  // int narray[18]={7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,56,75};
  // for (int i=0;i<18;i++){
  //   int n=narray[i];
    FILE *contributionsfile=fopen("c6multichannelcontributions12.csv","w");

  // for (int n=7;n<50;n++){
  int n=12;
  double interval=0.01;
    // if (n <= 20){
    //   interval=0.01;
    // }
    // else if (n <= 40 && n> 20){
    //   interval=0.05;
    // }
    // else if (n >40){
    //   interval=0.1;
    // }

  double targetenergy= energyinitial[n][2][0];//(channelenergy[n-4]-45932.1982)/219475.0;
    double targetenergyd= targetenergy + (Is - Id)/219475.0;
    double targetenergyp= targetenergy + (Is - Ip)/219475.0;
    printf("target energy is %g\n",targetenergy);
    // printf("energy=%g   target1=%g   targetd=%g, targetp=%g\n",channelenergy[n-4],targetenergy,targetenergyd,targetenergyp);
    double nustarget= sqrt(Ry/(45932.1982-targetenergy));
    double nud1target= sqrt(Ry/(60768.43-targetenergy));
    double nud2target= sqrt(Ry/(60488.09-targetenergy));
    //printf("energy=%g, nu1=%g, nu2=%g\n",energy,nu1,nu2);

    atomstate targetchannel1,targetchannel2,targetchannel3,targetchannel4;
    atomstate targetchannel5,targetchannel6;
    targetchannel1.loadstate(n,2,2.0,0.0,0.0,targetenergy);//singlets
    targetchannel2.loadstate(n,2,2.0,1.0,0.0,targetenergy);//triplet

    double dipole4d6s=6.48399986;
    double dipole5p2=-3.10055045;
    double dipole5p2ground=-32.24665918;
    double dipole4d2=-12.68701545;
    double dipole4d2f=28.26852701;

    double totaldipole[20][6];
    for (int i=0;i<20;i++){
      for (int j=0;j<6;j++) totaldipole[i][j]=0.0;
    }
    double energydifference[20][6];
    for (int i=0;i<20;i++){
      for (int j=0;j<6;j++) energydifference[i][j]=0.0;
    }
    double singledipole[20][2];
    for (int i=0;i<20;i++){
      for (int j=0;j<2;j++) singledipole[i][j]=0.0;
    }
    double dipole[20][2];
    for (int i=0;i<20;i++){
      for (int j=0;j<6;j++) dipole[i][j]=0.0;
    }

    //--------------------------------------------------------------------
    //------------------- Transition 1 -----------------------------------
    //---------------------1D2 to 1P1 ------------------------------------

    for (int nfinal=n+10;nfinal>n-11;nfinal--){
      if (nfinal <5) continue;
      int nindex= nfinal-(n-10);
      double energy= targetenergy-energyinitial[nfinal][1][0];
      energydifference[nindex][0]= energy;
      atomstate finalstate;
      finalstate.loadstate(nfinal,1,1.0,0.0,0.0,energyinitial[nfinal][1][0]);//
      dipmatels pstates(targetchannel1,finalstate);
      double dip1;
      dip1= pstates.dipmatel(interval);

      totaldipole[nindex][0]= (omegaq[0][0]*dip1*pchannel1fractions[nfinal-5]*channel1fractions[n-4]);
      totaldipole[nindex][0]+=(omegaq[0][1]*dipole4d6s*pchannel2fractions[nfinal-5]*channel3fractions[n-4]);
      totaldipole[nindex][0]+=(omegaq[0][2]*dipole5p2*pchannel2fractions[nfinal-5]*channel5fractions[n-4]);

      singledipole[nindex][0]=(omegaq[0][0]*dip1);
      dipole[nindex][0]=dip1;
    }
    //--------------------------------------------------------------------
    //------------------- Transition 2 -----------------------------------
    //---------------------1D2 to 1F3 ------------------------------------

    for (int nfinal=n+10;nfinal>n-11;nfinal--){
      if (nfinal <4) continue;
      int nindex= nfinal-(n-10);
      double energy= targetenergy-energyinitial[nfinal][3][0];
      energydifference[nindex][1]= energy;
      atomstate finalstate;
      finalstate.loadstate(nfinal,3,3.0,0.0,0.0,energyinitial[nfinal][3][0]);//
      dipmatels pstates(targetchannel1,finalstate);
      double dip1;
      dip1= pstates.dipmatel(interval);

      totaldipole[nindex][1]= (omegaq[1][0]*dip1*fchannel1fractions[nfinal-5]*channel1fractions[n-4]);
      totaldipole[nindex][1]+=(omegaq[1][1]*dipole4d6s*fchannel2fractions[nfinal-5]*channel3fractions[n-4]);
      totaldipole[nindex][1]+=(omegaq[1][2]*dipole5p2*fchannel2fractions[nfinal-5]*channel5fractions[n-4]);

      singledipole[nindex][1]=(omegaq[1][0]*dip1);
      dipole[nindex][1]=dip1;
  }
    //--------------------------------------------------------------------
    //------------------- Transition 3 -----------------------------------
    //---------------------3D2 to 3P2 ------------------------------------

    for (int nfinal=n+10;nfinal>n-11;nfinal--){
      if (nfinal <5) continue;
      int nindex= nfinal-(n-10);
      double energy= targetenergy-energyinitial[nfinal][1][3];
      energydifference[nindex][2]= energy;
      atomstate finalstate;
      finalstate.loadstate(nfinal,1,2.0,1.0,0.0,energyinitial[nfinal][1][3]);//
      dipmatels pstates(targetchannel2,finalstate);
      double dip1;
      dip1= pstates.dipmatel(interval);

      totaldipole[nindex][2]= (omegaq[2][0]*dip1*tripletp2channel1fractions[nfinal-5]*channel2fractions[n-4]);
      totaldipole[nindex][2]+=(omegaq[2][1]*dipole4d6s*tripletp2channel2fractions[nfinal-5]*channel4fractions[n-4]);
      totaldipole[nindex][2]+=(omegaq[2][2]*dipole4d2*tripletp2channel2fractions[nfinal-5]*channel6fractions[n-4]);
      totaldipole[nindex][2]+=(omegaq[2][3]*dipole4d2f*tripletp2channel3fractions[nfinal-5]*channel6fractions[n-4]);
  }
    //--------------------------------------------------------------------
    //------------------- Transition 4 -----------------------------------
    //---------------------3D2 to 3P1 ------------------------------------

    for (int nfinal=n+10;nfinal>n-11;nfinal--){
      if (nfinal <5) continue;
      int nindex= nfinal-(n-10);
      double energy= targetenergy-energyinitial[nfinal][1][2];
      energydifference[nindex][3]= energy;
      atomstate finalstate;
      finalstate.loadstate(nfinal,1,1.0,1.0,0.0,energyinitial[nfinal][1][2]);//
      dipmatels pstates(targetchannel2,finalstate);
      double dip1;
      dip1= pstates.dipmatel(interval);

      totaldipole[nindex][3]= (omegaq[3][0]*dip1*tripletp1channel1fractions[nfinal-5]*channel2fractions[n-4]);
      totaldipole[nindex][3]+=(omegaq[3][1]*dipole4d6s*tripletp1channel2fractions[nfinal-5]*channel4fractions[n-4]);
      totaldipole[nindex][3]+=(omegaq[3][2]*dipole4d2*tripletp1channel2fractions[nfinal-5]*channel6fractions[n-4]);
      totaldipole[nindex][3]+=(omegaq[3][3]*dipole4d2f*tripletp1channel3fractions[nfinal-5]*channel6fractions[n-4]);
  }
    //--------------------------------------------------------------------
    //------------------- Transition 5 -----------------------------------
    //---------------------3D2 to 3F2 ------------------------------------

    for (int nfinal=n+10;nfinal>n-11;nfinal--){
      if (nfinal <4) continue;
      int nindex= nfinal-(n-10);
      double energy= targetenergy-energyinitial[nfinal][3][1];
      energydifference[nindex][4]= energy;
      atomstate finalstate;
      finalstate.loadstate(nfinal,3,2.0,1.0,0.0,energyinitial[nfinal][3][1]);//
      dipmatels pstates(targetchannel2,finalstate);
      double dip1;
      dip1= pstates.dipmatel(interval);

      totaldipole[nindex][4]= (omegaq[4][0]*dip1*tripletf2channel1fractions[nfinal-5]*channel2fractions[n-4]);
      totaldipole[nindex][4]+=(omegaq[4][1]*dipole4d6s*tripletf2channel2fractions[nfinal-5]*channel4fractions[n-4]);
  }
    //--------------------------------------------------------------------
    //------------------- Transition 6 -----------------------------------
    //---------------------3D2 to 3F3 ------------------------------------

    for (int nfinal=n+10;nfinal>n-11;nfinal--){
      if (nfinal <4) continue;
      int nindex= nfinal-(n-10);
      double energy= targetenergy-energyinitial[nfinal][3][2];
      energydifference[nindex][5]= energy;
      atomstate finalstate;
      finalstate.loadstate(nfinal,3,3.0,1.0,0.0,energyinitial[nfinal][3][2]);//
      dipmatels pstates(targetchannel2,finalstate);
      double dip1;
      dip1= pstates.dipmatel(interval);

      totaldipole[nindex][5]= (omegaq[5][0]*dip1*tripletf3channel1fractions[nfinal-5]*channel2fractions[n-4]);
      totaldipole[nindex][5]+=(omegaq[5][1]*dipole4d6s*tripletf3channel2fractions[nfinal-5]*channel4fractions[n-4]);
    }

    //--------------------------------------------------------------------
    //-------------------All together now!--------------------------------
    //------------------------------- ------------------------------------
    double C6=0.0;
    double jprimed[6];
    jprimed[0]=1;
    jprimed[1]=3;
    jprimed[2]=2;
    jprimed[3]=1;
    jprimed[4]=2;
    jprimed[5]=3;
    double singletsinglet,singlettriplet,tripletsinglet,triplettriplet;
    singletsinglet=singlettriplet=tripletsinglet=triplettriplet=0.0;
    // (2,2)(2,2)
    double m1=2;
    double m2=2;
    for (int trans1=0;trans1<6;trans1++){
      for (int trans2=0;trans2<6;trans2++){
	for (int n1=0;n1<20;n1++){
	  for (int n2=0;n2<20;n2++){
	    double entry=omega12(m1,m2,m1,m2,jprimed[trans1],jprimed[trans2])*pow(totaldipole[n1][trans1]*totaldipole[n2][trans2],2);
	    double omega1,omega2,omegaboth, angular;
	    omega1= omegaq[trans1][0];
	    omega2= omegaq[trans2][0];
	    omegaboth=omega12(m1,m2,m1,m2,jprimed[trans1],jprimed[trans2]);
	    angular= omega1*omega2*omegaboth;
	    // printf("omega1=%g, omega2=%g, omega12=%g, product=%g\n",omega1,omega2,omegaboth,angular);
	    double energy=(energydifference[n1][trans1] + energydifference[n2][trans2]);
	    if (fabs(entry)< 1e-15 ||fabs(entry)> 1e25  || fabs(energy)< 1e-15 || fabs(energy)> 1.0) continue;
	    double c6entry=(entry)/energy;
	    C6+= c6entry;
	    int s1,s2;
	    if (trans1<2 && trans2<2) singletsinglet+=c6entry;
	    else if (trans1<2 && trans2>=2) singlettriplet+=c6entry;
	    else if (trans2<2 && trans1>=2)tripletsinglet+=c6entry;
	    else if (trans1>=2 && trans2>=2) triplettriplet+=c6entry;
	    fprintf(contributionsfile,"%d,%d,%d,%d,%g,%g,%g,%g,%g\n",trans1,trans2,n1,n2,c6entry,energy,omegaboth,totaldipole[n1][trans1],totaldipole[n2][trans2]);
	    if (C6 != C6) {
	      printf("entry=%g energy=%g trans1=%d  trans2=%d n1=%d n2=%d\n",entry,energy,trans1,trans2,n1,n2);
	      printf("angular=%g  radial1=%g radial2=%g\n",angular,totaldipole[n1][trans1],totaldipole[n2][trans2]);
	      exit(0);
	    }
	  }
	}
      }
    }
    printf("n=%d, C6=%g\n",n,C6);
    // fprintf(contributionsfile,"%d   %g   %g   %g   %g\n",n,singletsinglet,singlettriplet,tripletsinglet,triplettriplet);

     //------------------------------------------------------------------
     //-------------------single channel---------------------------------
     FILE *contributionsoutput= fopen("c6contributionssinglechannel12.csv","w");
    double C6single[4];
    for (int i=0;i<4;i++) C6single[i]=0.0;
    // (2,2)(2,2)
    m1=2;
    m2=2;
    for (int trans1=0;trans1<2;trans1++){
      for (int trans2=0;trans2<2;trans2++){
    	for (int n1=0;n1<20;n1++){
    	  for (int n2=0;n2<20;n2++){
    	    double entry=omega12(m1,m2,m1,m2,jprimed[trans1],jprimed[trans2])*pow(singledipole[n1][trans1]*singledipole[n2][trans2],2);
    	      double omegaboth=omega12(m1,m2,m1,m2,jprimed[trans1],jprimed[trans2]);
    	      // printf("omega12=%g, j1primed=%g, j2primed=%g\n",omegaboth,jprimed[trans1],jprimed[trans2]);
    	    double energy=(energydifference[n1][trans1] + energydifference[n2][trans2]);
    	    if (fabs(entry)< 1e-10 ||fabs(entry)> 1e25  || fabs(energy)< 1e-15 || fabs(energy)> 1.0) continue;
    	    C6single[0]+=entry/energy;
    	    fprintf(contributionsoutput,"%d,%d,%d,%d,%g,%g,%g,%g,%g\n",trans1,trans2,n1,n2,(entry/energy),energy,omegaboth,singledipole[n1][trans1],singledipole[n2][trans2]);
    	    if (C6single[0] != C6single[0]) exit(0);
    	  }
    	}
      }
    }
    fclose(contributionsoutput);
    printf("C6single[0]=%g\n",C6single[0]);


  // }
     fclose(contributionsfile);
}
