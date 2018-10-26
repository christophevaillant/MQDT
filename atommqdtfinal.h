/* This header file includes classes that describe atomic states
   and their interactions for two electron atoms. Includes dipole and
   quadrupole interactions, as well as angular dependences.
   Author: Christophe Vaillant. */

#include <math.h>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_sf_gamma.h"
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <time.h>

#define pi 3.141592653589793238

using namespace std;

extern "C" {void clebschgordan_(double*,double*,double*,double*,double*,double*,double*);}
extern "C" {void sixjsymbol_(double*,double*,double*,double*,double*,double*,double*);}
extern "C" {void threejsymbol_(double*,double*,double*,double*,double*,double*,double*);}
extern "C" {void factorial_(double*, double*);}
extern "C" {void binomial_(int*,int*,double*);}
extern "C" {void __coulombmod_MOD_num(double*,int*, double*, double*, int*, double*, double*, double*, double*, int*);}
extern "C" {void __coulombmod_MOD_fgh(double*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);}

class atomstate{
  /* Contains the quantum numbers of the input states, and their
     energies that should be read from a quantum defects file.*/
  double max(double a, double b){
    double answer;
    if (a>b)
      answer=a;
    else
      answer=b;
    return(answer);
  };
 public:

  double squareintegral(double * integrand, int n, double h,double a, double b){
  double result=0.0;
  int lowerindex= (int) (a/h);
  int upperindex= (int) (b/h);
  /* printf("a=%g   lower=%d   b=%g   upper=%d\n",a,lowerindex,b,upperindex); */
  result= pow(integrand[lowerindex],2)+pow(integrand[upperindex-1],2);
  for (int i=lowerindex+2;i<upperindex;i+=2){
    result+= 2.0*pow(integrand[i],2) + 4.0*pow(integrand[i-1],2);
  }
  result*= h/3.0;
  if (result <0.0) {printf("%g\n",result);}
  return(result);
}
 int wavefunction(double * twf, int elements,double interval){
    /* printf("nu=%g\n",nu); */
   double *rwf=(double*) calloc(elements,sizeof(double));//regular wavefunction
   double tp,tp2;
   double energy= 2.0*E;//-1.0/pow(nu,2);
   double rcore=1.73484090163;
   tp= (-1 + sqrt(1.0 + l*(l+1)*energy))/energy;
   tp2= (-1 - sqrt(1.0 + l*(l+1)*energy))/energy;
   int nelem2= (int) (tp/interval)+1;
   /* printf("nelem2 should be %g, but is %d\n",tp/interval, nelem2); */
    if (tp<interval) tp=0.0;
    int nelem,rmax;
    double r0,g0,f0,ff,wfp,gf;
    r0=g0=f0=ff=wfp=0.0;
    int der=0;
    rmax= (int) 5*pow(n,2);
    nelem= elements-nelem2;
    double *rwf1= (double*) calloc(elements,sizeof(double));
    r0=(double) rmax;//2.0*((nelem2-30)*interval);
    /* printf("tp1=%g, tp2=%g, nelem=%d elements=%d, rmax=%d\n",tp,tp2,nelem,elements,rmax); */
    double thetamax= theta(nu,l,r0);
    double etamax= eta(nu,l,r0);
    double gammain= (double) l + 1.0 + nu;
    double gammaterm1;
    gammaterm1= gsl_sf_gamma(gammain);
    double gammaterm2= gsl_sf_gamma(nu-l);
    bool scale=false;
    double interval2=-interval;
    //regular solution:
    double A=1.0;
    for (int i=0;i<l+1;i++){
      A*=1 + pow(i,2)*energy;
    }

    f0= -cos(pi*nu)*thetamax/gammaterm1;
    f0+= sin(pi*nu)*gammaterm2*etamax/pi;
    f0*=pow(-1,n-l-1)*pow(nu,l+1);//
    double thetar= theta(nu,l,r0+interval2);
    double etar= eta(nu,l,r0+interval2);
    ff= -cos(pi*nu)*thetar/gammaterm1;
    ff+= sin(pi*nu)*gammaterm2*etar/pi;
    ff*=pow(-1,n-l-1)*pow(nu,l+1); //
    rwf1[0]=f0;
    rwf1[1]=ff;
    __coulombmod_MOD_num(&energy,&l,&r0,&interval2,&elements,&f0,&ff,rwf1,&wfp,&der);
    if (tp < interval){
      for (int i=1;i<elements;i++){
	//printf("i=%d, rwf1=%g\n",i,rwf1[nelem-i]);
	rwf[i]=rwf1[elements-i];
	if (rwf[i] != rwf[i]){
	  /* printf("i=%d, NaN for rwf (94)\n",i); */
	  rwf[i]=0.0;
	}

      };
      double r0= 1.0;
      int rindex= (int) (r0/interval);
      double f,fp,g,gp,h,hp,fbar,gbar,hbar,h0,hf;
      double ac=1e-10;
      __coulombmod_MOD_fgh(&energy,&l,&r0,&ac,&f0,&fp,&g0,&gp,&h0,&hp,&fbar,&gbar,&hbar);
      double twf0=sqrt(A/2.0)*f0*cos(pi*nu) + g0*sin(pi*nu)/sqrt(A*2.0);//f*cos(pi*nu) + h*sin(pi*nu);//
      double scalefactor= rwf[rindex]/twf0;
      for (int i=1;i<elements;i++){
	rwf[i]= rwf[i]/scalefactor;
      }
    }
    else {
      for (int i=0;i<elements;i++){
	rwf[elements-i-1]=rwf1[i];
	//printf("i=%d,nelem=%d\n",i,nelem);
	if (rwf[i] != rwf[i]){
	  rwf[i]=0.0;
	  /* printf("i=%d, NaN for rwf (110)\n",i);	 */
	}
	if (rwf[i]>1.0e154){
	  scale=true;
	}
      }
    }
    free(rwf1);
    double scalefactor=1.0;
    if (tp > interval){
      /* printf("neleme2=%d  tp=%g   interval=%g\n",nelem2,tp,interval); */
      double *rwf2=(double*) calloc(nelem2+1,sizeof(double));
      r0=tp;
      double ac=1e-10;
      double f,fp,g,gp,h,hp,fbar,gbar,hbar,h0,hf;
      double mu0= n-nu;
      /* double rchar=0.5*tp; */
      double rchar;
      if (tp > rcore) rchar= tp;
      else rchar=rcore;
      for (int i=1;i<nelem2;i++){
      	double r=i*interval;
	double mu=(1.0-exp(-pow(r/rchar,10)))*mu0;
	double nu=n- mu;//
	/* printf("r=%g   nu=%g, mu=%g\n",r,nu,mu); */
	__coulombmod_MOD_fgh(&energy,&l,&r,&ac,&f,&fp,&g,&gp,&h,&hp,&fbar,&gbar,&hbar);
      	rwf2[i]=pow(-1,l)*(sqrt(A/2.0)*f*cos(pi*nu) + g*sin(pi*nu)/sqrt(A*2.0));//f*cos(pi*nu) + g*sin(pi*nu);//
      }
      double r=interval*nelem2;
	__coulombmod_MOD_fgh(&energy,&l,&r,&ac,&f,&fp,&g,&gp,&h,&hp,&fbar,&gbar,&hbar);
      	rwf2[nelem2]=pow(-1,l)*(sqrt(A/2.0)*f*cos(pi*nu) + g*sin(pi*nu)/sqrt(A*2.0));//f*cos(pi*nu) + g*sin(pi*nu);//
      /* printf("nu=%g\n",nu); */
      scalefactor=rwf[nelem2]/rwf2[nelem2];//
    
      /* printf("1st scalefactor= %g, rwf=%g rwf2=%g\n",scalefactor,rwf[nelem2],rwf2[nelem2]); */
      /* printf("f0(%g)=%g   ff(%g)=%g\n",r0,f0,r,ff); */
      for (int i=0;i<nelem2;i++){
      	rwf2[i]*=scalefactor;
      	rwf[i]=rwf2[i];
      	if (rwf[i]>1.0e154){
      	  scale=true;
      	}
      }
      //printf("nelem+nelem2=%d, elements=%d\n",nelem+nelem2,elements);
    }
    if (scale){
      for (int i=0;i<elements;i++){
	rwf[i]/=1.0e200;
      }
    }
    for (int i=0;i<elements;i++){
      if (rwf[i]>1.0e154){// || iwf[i]>1.0e154
	scale=true;
      }
    }
    rwf[0]=0.0;
    for (int i=0;i<elements;i++){
      twf[i]= sqrt(2.0/pow(nu,3))*rwf[i]/scalefactor;//*fabs(interval)*i;//*cos(pi*nu)+gwf[i]*sin(pi*nu);
      if (twf[i] != twf[i]) {
	printf("NaN at twf[%d] wavefunction 1\n",i);
	printpair();
	exit(0);}
      if (scale){
	twf[i]/=1e50;
      }
    }

    double normalisation=0.0;
    normalisation= squareintegral(twf,elements,fabs(interval),0.0,rmax);
    /* printf("normalisation=%g\n",normalisation); */
    for (int i=0;i<elements;i++){
      twf[i]*=1.0/sqrt(normalisation);
    }
    /* if (l==1){ */
    /*   FILE * wavefunctionoutput; */
    /*   wavefunctionoutput= fopen("testwavefunctionmqdt.csv","w"); */
    /*   for (int i=0;i<elements;i++){ */
    /*   	fprintf(wavefunctionoutput,"%g   %g   %g   %g\n",i*interval,rwf[i],gwf[i],twf[i]); */
    /*   } */
    /*   fclose(wavefunctionoutput); */
    /*   //exit(0); */

    /* } */
    free(rwf);
    return(0);
  }

  double cursiveg(double energy,int l){
    double gsum=0.0;
    gsum= 1.0 + energy/10.0 + pow(energy,2)/21 + pow(energy,3)/20;
    gsum/=12.0;
    for (int i=0;i<l+1;i++){
      gsum+= i/(1.0+ pow(i,2)*energy);
    }
    gsum*= energy/pi;
    return(gsum);
  }


  double theta(double nu,int l, double r){
    double x= 2.0*r/nu;
    double answer= exp(-r/nu)*pow(x,nu);
    return(answer);
  }

  double eta(double nu,int l, double r){
    double x= 2.0*r/nu;
    double answer= exp(r/nu)*pow(x,-nu);
    return(answer);
  }

  atomstate (int,int,double,double,double,double);
  atomstate ();
  int n,l;
  double m,j,s;
  double nu;
  double E;
  int loadstate (int na, int la,double ja, double sa, double ma, double energya){
    n=na;
    l=la;
    m=ma;
    nu= sqrt(1.0/(-2.0*energya));
    E=energya;
    s=sa;
    j=ja;
    return(0);
  };

  int printpair(){
    printf("n=%d, l=%d, j=%g, s=%g, m=%g, nu=%g E=%g\n",n,l,j,s,m,nu,E);
  }

};

//Constructor:
atomstate::atomstate (int na, int la,double ja, double sa, double ma, double energya) {
  int n = na;
  int l=la;
  double m=ma;
  double nu=sqrt(1.0/(-2.0*energya));
  double s=sa;
  double j=ja;
  double E=energya;
};

atomstate::atomstate(){
 int n = 0;
 int l=0;
 double m=0.0;
 double nu=0.0;
 double s=0.0;
 double j=0.0;
 double E=0.0;

}; //alternate constructor.


double integral(double * integrand, int n, double h,double a, double b){
  //integrate between a and b
  int lowerindex= (int) (a/h);
  int upperindex= (int) (b/h);
  /* printf("a=%g   b=%g   lower=%d   upper=%d, n=%d\n",a,b,lowerindex,upperindex,n); */
  double result=0.0;
  result= integrand[lowerindex]+integrand[upperindex-1];
  for (int i=lowerindex+2;i<upperindex;i+=2){
    result+= 2.0*integrand[i] + 4.0*integrand[i-1];
  }
  result*= h/3.0;
  return(result);
}

double integral(double * integrand, int n, double h){
  double result=0.0;
  result= integrand[0]+integrand[n-1];
  for (int i=2;i<n;i+=2){
    result+= 2.0*integrand[i] + 4.0*integrand[i-1];
  }
  result*= h/3.0;
  return(result);
}

double squareintegral(double * integrand, int n, double h){
  double result=0.0;
  result= pow(integrand[0],2)+pow(integrand[n-1],2);
  for (int i=2;i<n;i+=2){
    result+= 2.0*pow(integrand[i],2) + 4.0*pow(integrand[i-1],2);
  }
  result*= h/3.0;
  return(result);
}


class dipmatels{
  /*Contains functions to calculate matrix elements.
    Dipole matrix elements calculated using Edmonds et al
  (1979)*/
  double max(double a, double b){
    double answer;
    if (a>b)
      answer=a;
    else
      answer=b;
    return(answer);
  };


 public: 
  dipmatels(atomstate, atomstate);

  double quadmatel(){
    double interval=0.01;
    double matrixelement;
    int elements,elementsinitial,elementsfinal;
    elementsinitial= (int) (5*pow(initial.n,2)/interval);
    elementsfinal= (int) (5*pow(final.n,2)/interval);
    elements=min(elementsinitial,elementsfinal);
    double *waveinitial,*integrand,*wavefinal;
    waveinitial=(double*) calloc(elementsinitial,sizeof(double));
    wavefinal=(double*) calloc(elementsfinal,sizeof(double));
    integrand=(double*) calloc(elements,sizeof(double));    
    initial.wavefunction(waveinitial,elementsinitial,interval);
    final.wavefunction(wavefinal,elementsfinal,interval);
    for (int i=0;i<elements;i++){
      integrand[i]= waveinitial[i]*wavefinal[i]*pow((i+1)*interval,2);
    }
    matrixelement= integral(integrand,elements,interval);
    free(integrand);
    free(waveinitial);
    free(wavefinal);
    return(matrixelement);
  }

  double quadmatel(double interval){
    double matrixelement;
    int elements,elementsinitial,elementsfinal;
    elementsinitial= (int) (5*pow(initial.n,2)/interval);
    elementsfinal= (int) (5*pow(final.n,2)/interval);
    elements=min(elementsinitial,elementsfinal);
    double *waveinitial,*integrand,*wavefinal;
    waveinitial=(double*) calloc(elementsinitial,sizeof(double));
    wavefinal=(double*) calloc(elementsfinal,sizeof(double));
    integrand=(double*) calloc(elements,sizeof(double));    
    initial.wavefunction(waveinitial,elementsinitial,interval);
    final.wavefunction(wavefinal,elementsfinal,interval);
    for (int i=0;i<elements;i++){
      integrand[i]= waveinitial[i]*wavefinal[i]*pow((i+1)*interval,2);
    }
    matrixelement= integral(integrand,elements,interval);
    free(integrand);
    free(waveinitial);
    free(wavefinal);
    return(matrixelement);
  }

  double dipmatel(){
    double interval=0.01;
    double matrixelement;
    int elements,elementsinitial,elementsfinal;
    elementsinitial= (int) (5*pow(initial.n,2)/interval);
    elementsfinal= (int) (5*pow(final.n,2)/interval);
    elements=min(elementsinitial,elementsfinal);
    double *waveinitial,*integrand,*wavefinal;
    waveinitial=(double*) calloc(elementsinitial,sizeof(double));
    wavefinal=(double*) calloc(elementsfinal,sizeof(double));
    integrand=(double*) calloc(elements,sizeof(double));
    initial.wavefunction(waveinitial,elementsinitial,interval);
    final.wavefunction(wavefinal,elementsfinal,interval);
    for (int i=0;i<elements;i++){
      integrand[i]= waveinitial[i]*(i+1)*interval*wavefinal[i];
    }
    matrixelement= integral(integrand,elements,interval);
    free(integrand);
    free(waveinitial);
    free(wavefinal);
    return(matrixelement);
  }

  double dipmatel(double interval){
    double matrixelement;
    int elements,elementsinitial,elementsfinal;
    elementsinitial= (int) (5*pow(initial.n,2)/interval);
    elementsfinal= (int) (5*pow(final.n,2)/interval);
    elements=min(elementsinitial,elementsfinal);
    double *waveinitial,*integrand,*wavefinal;
    waveinitial=(double*) calloc(elementsinitial,sizeof(double));
    wavefinal=(double*) calloc(elementsfinal,sizeof(double));
    integrand=(double*) calloc(elements,sizeof(double));
    initial.wavefunction(waveinitial,elementsinitial,interval);
    final.wavefunction(wavefinal,elementsfinal,interval);
    for (int i=0;i<elements;i++){
      integrand[i]= waveinitial[i]*(i+1)*interval*wavefinal[i];
    }
    matrixelement= integral(integrand,elements,interval);
    free(integrand);
    free(waveinitial);
    free(wavefinal);
    return(matrixelement);
  }

  double angular(){
    double angular0= max(initial.l,final.l)*(2.0*final.j + 1)*(2.0*initial.s + 1);
    double a,b,c,d,e,f,wigner;
    a=initial.j;
    b=1.0;
    c=final.j;
    d=final.l;
    e=initial.s;
    f=initial.l;
    sixjsymbol_(&a,&b,&c,&d,&e,&f,&wigner);
    double result=pow(wigner,2)*angular0;
    return(result);
  }

  atomstate initial,final;
};

//Constructor:
dipmatels::dipmatels(atomstate a, atomstate b) : initial(a), final(b){};
