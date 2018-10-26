#include "atommqdtfinal.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

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


int main(){
  const double Ry=109736.627;
  const double Is=45932.1982;
  const double Id= 0.5*(60768.43+60488.09);
  const double Ip= 70048.11;
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
  FILE * statefractionfile= fopen("sstatefractionsnew.csv","r");//
  double *channel1fractions= new double[67];
  double *channel2fractions= new double[67];
  double *channel3fractions= new double[67];
  double *channelenergy= new double[67];
  for (int i=0;i<67;i++){
    double E;
    int n;
    fscanf(statefractionfile,"%lf   %lf   %lf   %lf",&channelenergy[i],&channel1fractions[i],&channel2fractions[i],&channel3fractions[i]);
  }

  FILE * pstatefractionfile= fopen("pstatefractionsextended.csv","r");//
  double *pchannel1fractions= new double[76];
  double *pchannel2fractions= new double[76];
  double *pchannelenergy= new double[76];
  for (int i=0;i<76;i++){
    fscanf(pstatefractionfile,"%lf   %lf   %lf",&pchannelenergy[i],&pchannel1fractions[i],&pchannel2fractions[i]);
  }

  FILE * tripletpstatefractionfile= fopen("tripletp1statefractions.csv","r");//
  double *tripletpchannel1fractions= new double[76];
  double *tripletpchannel2fractions= new double[76];
  double *tripletpchannel3fractions= new double[76];
  double *tripletpchannelenergy= new double[76];
  for (int i=0;i<76;i++){
    double throwaway;
    fscanf(tripletpstatefractionfile,"%lf   %lf   %lf   %lf",&tripletpchannelenergy[i],&tripletpchannel1fractions[i],&tripletpchannel2fractions[i],&tripletpchannel3fractions[i]);
  }
  fclose(tripletpstatefractionfile);
  
  int l=0;
  double j,m,s;
  j=m=s=0.0;
  // printf("Enter n:\n");
  // scanf("%d",&n);
  FILE * lifetimesoutputfile= fopen("lifetimessstates300K.csv","w");
  FILE * widthsoutputfile= fopen("partialwidthssstates.csv","w");
  FILE *lifetimescoulombfile=fopen("lifetimescoulomb.csv","w");
  double T=300.0;
  for (int n=40;n<75;n++){
  // int n=75;
    double targetenergy= (channelenergy[n-5]-Is)/219475.0;//energyinitial[n][l][0];
    double targetenergy2= (channelenergy[n-5] - Ip)/219475.0;
    double targetenergy3= (channelenergy[n-5]- Id)/219475.0;
    // printf("target energy1 %g   energy2 %g   energy3 %g\n",targetenergy,targetenergy2,targetenergy3);

    double nustarget= sqrt(Ry/(Is-targetenergy));
    double nud1target= sqrt(Ry/(Ip-targetenergy));
    double nud2target= sqrt(Ry/(Id-targetenergy));
    //printf("energy=%g, nu1=%g, nu2=%g\n",energy,nu1,nu2);

    atomstate targetchannel1,targetchannel2,targetchannel3,targetchannel4;
    targetchannel1.loadstate(n,0,0.0,0.0,0.0,targetenergy);//singlets
    double dipole5p2=-3.10055045;
    double dipole5p2ground=-32.24665918;
    double dipole4d2=-12.68701545;
    double dipole4d4f=28.26852701;

    double blackbody=0.0;
    //-------------------------------------------------------------------
    //------------------- Channel 1 -----------------------------------
    double dipolechannel;
    dipolechannel=0.0;
    double width=0.0;
    double width0K=0.0;
    double coulomb=0.0;
    double channel1,channel2,channel3;
    channel1=channel2=channel3=0.0;
    //---------------------------------------------------
    //Singlet P states
    for (int nfinal=n+5;nfinal>4;nfinal--){
      atomstate finalstate;
      finalstate.loadstate(nfinal,1,1.0,0.0,0.0,energyinitial[nfinal][1][0]);//
      dipmatels pstates(targetchannel1,finalstate);
      double dip;
      dip= pstates.dipmatel();
      // if (dip!=dip){
      // }
      double energy= targetenergy-energyinitial[nfinal][1][0];
      dipolechannel=-(pchannel1fractions[nfinal-5]*dip*channel1fractions[n-5]);//
      dipolechannel+=(pchannel2fractions[nfinal-5]*dipole4d2*channel2fractions[n-5]*sqrt(2.0/5.0));
      if (energy<0.0){
	blackbody= 1.0/(exp(fabs(energy)/(T*3.166812e-6))-1);
      }
      else {
	blackbody= 1.0 + 1.0/(exp(fabs(energy)/(T*3.166812e-6))-1);
	width0K+=pow(dipolechannel,2)*pow(fabs(energy),3);
	channel1+=pow(pchannel1fractions[nfinal-5]*dip*channel1fractions[n-5],2)*pow(fabs(energy),3);//
	channel2+= pow(pchannel2fractions[nfinal-5]*dipole4d2*channel2fractions[n-5],2)*pow(fabs(energy),3)*2.0/5.0;      }
      width+=pow(dipolechannel,2)*blackbody*pow(fabs(energy),3);
      coulomb+=pow(dip,2)*blackbody*pow(fabs(energy),3);
    }
    // printf("width1=%g\n",width);
    //--------------------------------------------------------
    //4d5p 1P1 state
    //41172.054   0.825708307641   0.174291692359
    
    
    double energyfinal=(41172.054-Is)/219475.0;//Is
    atomstate perturberfinalstate;
    perturberfinalstate.loadstate(6,1,1.0,0.0,0.0,energyfinal);//
    dipmatels pstates(targetchannel1,perturberfinalstate);
    double dip;
    dip= pstates.dipmatel();
    // if (dip!=dip){
    // }
    double energy= targetenergy-(41172.054-Is)/219475.0;//energyfinal
    dipolechannel=-(sqrt(0.825708307641)*dip*channel1fractions[n-5]);
    dipolechannel+=(sqrt(0.174291692359*2.0/5.0)*dipole4d2*channel2fractions[n-5]);

    if (energy<0.0){
      blackbody= 1.0/(exp(fabs(energy)/(T*3.166812e-6))-1);
    }
    else {
      blackbody= 1.0 + 1.0/(exp(fabs(energy)/(T*3.166812e-6))-1);
      channel1+=pow(0.825708307641*dip*channel1fractions[n-5],2)*pow(fabs(energy),3);
      channel2+=(2.0/5.0)*pow(0.174291692359*dipole4d2*channel2fractions[n-5],2)*pow(fabs(energy),3);
      width0K+=pow(dipolechannel,2)*pow(fabs(energy),3);
    }
    width+=blackbody*pow(dipolechannel,2)*pow(fabs(energy),3);
    coulomb+=blackbody*pow(dip,2)*pow(fabs(energy),3)*0.825708307641;
    // printf("width2=%g\n",dipolechannel);
    //-------------------------------------------------------------------
    //------------------- Channel 3 -----------------------------------
    // //--------------------------------------------------------
    //---------------------------------------------------
    // 3P1 states, 4d5p component (channel2)
    for (int nfinal=n+5;nfinal>4;nfinal--){
      double energy= targetenergy-(tripletpchannelenergy[nfinal-5]-Is)/219475.0;
      double dipole=-sqrt(3.0/10.0)*channel3fractions[n-5]*dipole4d2*tripletpchannel2fractions[nfinal-5];
      dipole+=-tripletpchannel3fractions[nfinal-5]*channel3fractions[n-5]*dipole4d4f/sqrt(5.0);
      if (energy<0.0){
	blackbody= 1.0/(exp(fabs(energy)/(T*3.166812e-6))-1);
      }
      else {
	channel3+=(3.0/10.0)*pow(channel3fractions[n-5]*dipole4d2*tripletpchannel2fractions[nfinal-5],2)*pow(fabs(energy),3);
	channel3+=pow(tripletpchannel3fractions[nfinal-5]*channel3fractions[n-5]*dipole4d4f,2)*pow(fabs(energy),3)/(5.0);
	blackbody= 1.0 + 1.0/(exp(fabs(energy)/(T*3.166812e-6))-1);
	width0K+=pow(dipole,2)*pow(fabs(energy),3);
      }
      double contr=-sqrt(3.0/10.0)*channel3fractions[n-5]*dipole4d2*tripletpchannel2fractions[nfinal-5];
      // printf("%d   %g   %g\n",n,channel3fractions[n-5],contr);
      width+=pow(dipole,2)*blackbody*pow(fabs(energy),3);//
    }
    // printf("width3=%g\n",width);
    //--------------------------------------------------------
    //4d5p state 3P1
    //37302.731   0.833117115407   0.552844880536   -0.0166856249825
    
    energy= targetenergy-(37302.731-Is)/219475.0;//finalenergy
    double finalenergy=(37302.731-Is)/219475.0;
    // printf("dip=%g, omega=%g\n",pow(dip,2),pow(fabs(energy),3));
    if (energy>0.0){
      double blackbody= 1.0 + 1.0/(exp(energy/(T*3.166812e-6))-1);
      double dipole=-sqrt(0.305637461934*3.0/10.0)*channel3fractions[n-5]*dipole4d2;
      dipole+=-sqrt(0.305637461934/5.0)*channel3fractions[n-5]*dipole4d4f;
      channel3+=pow(dipole,2)*pow(fabs(energy),3);
      width+=pow(dipole,2)*blackbody*pow(fabs(energy),3);//
      width0K+=pow(dipole,2)*pow(fabs(energy),3);
    }
    // printf("width4=%g\n",width);
    //--------------------------------------------------------
    //4d5p state 3D1
    //36264.151

    energy= targetenergy-(36264.151-Is)/219475.0;//finalenergy
    if (energy<0.0){
      blackbody= 1.0/(exp(fabs(energy)/(T*3.166812e-6))-1);
    }
    else {
      blackbody= 1.0 + 1.0/(exp(fabs(energy)/(T*3.166812e-6))-1);
      channel3+=pow(channel3fractions[n-5]*dipole4d2,2)*pow(fabs(energy),3)/10.0;
      width0K+=pow(channel3fractions[n-5]*dipole4d2,2)*pow(fabs(energy),3);
    }
    width+=blackbody*pow(channel3fractions[n-5]*dipole4d2,2)*pow(fabs(energy),3)/10.0;
    // printf("width5=%g\n",width);
    //----------------------------------------------------------
    //----------------------------------------------------------
    //----------------------------------------------------------
    //Final sum.
    coulomb*=4.0/(3.0*pow(137.036,3));
    width*=4.0/(3.0*pow(137.036,3));
    width0K*=4.0/(3.0*pow(137.036,3));
    channel1*=4.0/(3.0*pow(137.036,3));
    channel2*=4.0/(3.0*pow(137.036,3));
    channel3*=4.0/(3.0*pow(137.036,3));
    double lifetime;
    lifetime= 1.0/(2.0*pi*width*6.57969e9);
    double lifetimecoulomb;
    lifetimecoulomb= 1.0/(2.0*pi*coulomb*6.57969e9);

    printf("%d   %g   %g\n",n,lifetime,width);
    // printf("%d\n",n);
    fprintf(lifetimesoutputfile,"%d   %g\n",n,lifetime);
    fprintf(lifetimescoulombfile,"%d   %g\n",n,lifetimecoulomb);
    fprintf(widthsoutputfile,"%d   %g   %g   %g   %g\n",n,channel1,channel2,channel3,width0K);
    // printf("%g   %g   %g   %g   %g\n",width1,width2,width3,width,lifetime);
    // printf("n=%d   width=%g   lifetime=%g, energy=%g\n",n,width,lifetime,channelenergy[n-5]);
    // printf("------------------------------------------------\n");
    // fprintf(widthsoutputfile,"%d   %g\n",n,width1);
  }
  fclose(widthsoutputfile);
  fclose(lifetimesoutputfile);
  fclose(lifetimescoulombfile);
}

