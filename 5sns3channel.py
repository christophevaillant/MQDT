import scipy as sp
import numpy as np
from scipy.optimize import fminbound,fmin,fmin_bfgs
import pylab as pl

Ry=109736.627
Is=45932.1982
Ip=60488.09#d3/2
Id2= 60768.43#d5/2

def F(nup):
    nus= 1.0/np.sqrt((Is - Ip)/Ry + (1.0/nup**2))
    return(nus)

def H(nup):
    nud2= 1.0/np.sqrt((Id2 - Ip)/Ry + (1.0/nup**2))
    return(nud2)

def I(nus):
    nup= 1.0/np.sqrt((Ip - Is)/Ry + (1.0/nus**2))
    return(nup)

def detcondition1(nup,mu1initial,mu2initial,mu3initial,R12,R13,dmude1,dmude2,dmude3):#
    nud2= H(nup)
    E= Ip - Ry/np.power(nup,2)
    mu1final=mu1initial+dmude1*(45932.1982 - E)/45932.1982
    mu2final=mu2initial+dmude2*(45932.1982 - E)/45932.1982
    mu3final=mu3initial+dmude3*(45932.1982 - E)/45932.1982
    epsilon2= np.tan(np.pi*(nup))+mu2final
    epsilon3= np.tan(np.pi*(nud2))+mu3final
    num=epsilon3*np.power(R12,2) + epsilon2*np.power(R13,2)
    denom= epsilon2*epsilon3
    epsilon1= num/denom
    result=  np.arctan(epsilon1 - mu1final)/np.pi
    return(result%1)

def diff(nup,mu1,mu2,mu3,R12,R13,dmude1,dmude2,dmude3):
    ryd= F(nup)
    detcon= detcondition1(nup,mu1,mu2,mu3,R12,R13,dmude1,dmude2,dmude3)
    result= ryd%1-detcon%1
    return(result**2)

def chisquared(params):#,expE,error
    R12,R13,mu1,mu2,mu3,dmude1,dmude2,dmude3=params#,dmude2
    chi=0.0
    for i in range(len(expE)):
        nupexp= np.sqrt(Ry/(Ip-expE[i]))
        nupcalc= fminbound(diff,nupexp-0.0005,nupexp+0.0005,args=(mu1,mu2,mu3,R12,R13,dmude1,dmude2,dmude3),xtol=1e-12)#,dmude2
        Ecalc= Ip - Ry/nupcalc**2
        chi= chi + ((expE[i]-Ecalc)/error[i])**2 #
    print chi/(len(expE) - len(params)),params#
    return(chi)

R12=-0.02301876
R13=0.36992972
mu1=1.05410457
mu2=2.90825856
mu3=-0.66272033
dmude1=0.82792176
dmude2=-13.88066116
dmude3=0.50817284

initialparams= [R12,R13,mu1,mu2,mu3,dmude1,dmude2,dmude3]#
finalparams=initialparams
global expE
global error
expE= []
error=[]
n=[]

inputfile= open("expenergies.csv","r")
for line in inputfile:
    nin,a,err= line.split()
    n.append(int(nin))
    expE.append(float(a)+float(err))#
    error.append(float(err))

# Hessfun= nd.Hessian(chisquared)
# Hessfun.romberg_terms=5
# arguments=(initialparams,expE,error)
# hessian= Hessfun(initialparams)
# print hessian

# curvature= 0.5*hessian

# error= linalg.inv(curvature)
# print error

# for i in range(len(initialparams)):
#     print initialparams[i], np.sqrt(error[i,i])


#######Minimization
difference=1e8
oldchi=chisquared(finalparams)
while (difference > 0.0001):
    try:
        finalparams= fmin(chisquared,initialparams,maxfun=1e8,maxiter=1e8)
        chi=chisquared(finalparams)
        # parambounds= [(0.0,1.0),(0.0,1.0),(0.0,1.0),(0.0,1.0),(0.0,1.0),(0.0,1.0),(0.0,1.0),(-2.0,2.0),(-2.0,2.0)]
        # minoutput=fmin_l_bfgs_b(chisquared,finalparams,args=(expE,error),approx_grad=1,bounds=parambounds)
        # finalparams,chi,throw= minoutput
        difference= abs(oldchi-chi)
        oldchi=chi
        print difference
        print '\a'
        initialparams=finalparams
    except KeyboardInterrupt:
        break

initialparams= [R12,R13,mu1,mu2,mu3,dmude1,dmude2,dmude3]#,dmude2
for i in range(len(initialparams)):
    print initialparams[i], abs(initialparams[i]-finalparams[i])


for i in range(len(expE)):
    nupexp= np.sqrt(Ry/(Ip-expE[i]))
    nupcalc= fminbound(diff,nupexp-0.0005,nupexp+0.0005,args=(mu1,mu2,mu3,R12,R13,dmude1,dmude2,dmude3),xtol=1e-12)#,dmude2,dmude3
    Ecalc= Ip - Ry/nupcalc**2
    print Ecalc, expE[i], Ecalc-expE[i]

