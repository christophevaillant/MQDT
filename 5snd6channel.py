import numpy as np
import pylab as pl
from scipy.optimize import fmin, fminbound,brent
import numpy.linalg as linalg

Ry=109736.627
Is= 45932.1982
Id1= 60768.43
Id2= 60488.09
Ip= Id2
Id3= Id1

def F(nud1):
    nus= 1.0/np.sqrt((Is - Id1)/Ry + (1.0/np.power(nud1,2)))
    return(nus)

def G(nud1):
  nud2= 1.0/np.sqrt((Id2-Id1)/Ry + 1.0/(nud1**2))
  return(nud2)

def H(nud1):
  nup= 1.0/np.sqrt((Ip-Id1)/Ry + 1.0/(nud1**2))#
  return(nup)

def I(nud1):
  nud3= 1.0/np.sqrt((Id3-Id1)/Ry + 1.0/(nud1**2))#
  return(nud3)

def J(nus):
  nud1= 1.0/np.sqrt(-(45932.1982-Id1)/Ry + 1.0/(nus**2))
  return(nud1)

def detcondition1(nu3,params):
    mu1initial=params[10]#mu1#
    b=params[0]#R12#
    c=params[1]#R13#
    d=params[2]#R14#
    e=params[3]#R15#
    f=params[4]#R16#
    mu2initial=params[11]#mu2#
    h=params[5]#R23#
    i=params[6]#R24#
    j=params[7]#R25#
    k=params[8]#R26#
    mu3initial=params[12]#mu3#
    m=params[9]#R34#
    mu4initial=params[13]#mu4#
    mu5initial=params[14]#mu5#
    mu6initial=params[15]#mu6#
    dmude1=params[16]
    dmude2=params[17]
    dmude5=params[18]

    E= Id1 - Ry/nu3**2
    nu2=F(nu3)
    nu4= G(nu3)
    nu5=H(nu3)
    nu6=I(nu3)
    mu1= mu1initial + dmude1*(45932.1982 - E)/45932.1982
    mu2= mu2initial + dmude2*(45932.1982 - E)/45932.1982
    mu3= mu3initial# + dmude3*(45932.1982 - E)/45932.1982
    mu4= mu4initial# + dmude4*(45932.1982 - E)/45932.1982
    mu5= mu5initial + dmude5*(45932.1982 - E)/45932.1982
    mu6= mu6initial#+ dmude6*(45932.1982 - E)/45932.1982

    g= np.tan(np.pi*(nu2)) + mu2
    l= np.tan(np.pi*(nu3)) + mu3
    p= np.tan(np.pi*(nu4)) + mu4
    s= np.tan(np.pi*(nu5)) + mu5
    u= np.tan(np.pi*(nu6)) + mu6
    #Channel 5 coupled to channel 1:

    num= l*p*s*u*b**2 -s*u*(b*m)**2 -2*b*c*h*p*s*u + 2*b*c*i*m*s*u + 2*b*d*h*m*s*u
    num+= -2*b*d*i*l*s*u - 2*b*e*j*l*p*u + 2*b*e*j*u*m**2 - 2*b*f*k*l*p*s
    num+= 2*b*f*k*s*m**2 + g*p*s*u*c**2 - s*u*(c*i)**2 - p*u*(c*j)**2 - p*s*(c*k)**2
    num+= -2*c*d*g*m*s*u + 2*c*d*h*i*s*u + 2*c*d*m*u*j**2 + 2*c*d*m*s*k**2
    num+= 2*c*e*h*j*p*u - 2*c*e*i*j*m*u + 2*c*f*h*k*p*s - 2*c*f*i*k*m*s
    num+= g*l*s*u*d**2 - s*u*(d*h)**2 - l*u*(d*j)**2 - l*s*(d*k)**2 - 2*d*e*h*j*m*u
    num+= 2*d*e*i*j*l*u - 2*d*f*h*k*m*s + 2*d*f*i*k*l*s + g*l*p*u*e**2 - g*u*(e*m)**2
    num+= - p*u*(e*h)**2 + 2*h*i*m*u*e**2 - l*u*(e*i)**2 - l*p*(e*k)**2 + (e*k*m)**2
    num+= 2*e*f*j*k*l*p - 2*e*f*j*k*m**2 + g*l*p*s*f**2 - g*s*(f*m)**2 - p*s*(f*h)**2
    num+= 2*h*i*m*s*f**2 - l*s*(f*i)**2 - l*p*(f*j)**2 + (f*j*m)**2
    denom=g*l*p*s*u - g*s*u*m**2 - p*s*u*h**2 +2*h*i*m*s*u - l*s*u*i**2
    denom+= - l*p*u*j**2 + u*(j*m)**2 - l*p*s*k**2 + s*(k*m)**2
    epsilon1= num/denom
    nu1calc= (np.arctan(epsilon1- mu1))/np.pi+1.0
    return(nu1calc)

def diff(nud1,params):
    ryd= F(nud1)
    detcon= detcondition1(nud1,params)
    result= ryd%1-detcon%1
    return(np.power(result,2))

def chisquared(params,expE,error):
    chi=0.0
    for i in range(len(expE)):
        nud1exp= np.sqrt(Ry/(Id1-expE[i]))
        # if (i<10):
        #     difference=1e-3#7e-4
        # else:
        #     difference=3e-4
        difference=1e-4
        nud1calc= fminbound(diff,nud1exp-difference,nud1exp+difference,args=[params],xtol=1e-12)#
        Ecalc= Id1 - Ry/nud1calc**2
        #print expE[i], Ecalc
        chi= chi + ((expE[i]-Ecalc)/error[i])**2 #(nud1exp-nud1calc)**2#
    print chi/(len(expE)- len(params)),params# 
    return(chi)

def coeffs(nu3,params):
    mu1initial=params[10]#mu1#
    R12=params[0]#R12
    R13=params[1]#R13
    R14=params[2]
    R15=params[3]#R15
    R16=params[4]
    mu2initial=params[11]#mu2
    R23=params[5]#R23
    R24=params[6]#R24
    R25=params[7]#R25
    R26=params[8]#R26
    mu3initial=params[12]#mu3
    R34=params[9]#R34
    mu4initial=params[13]#mu4
    mu5initial=params[14]#mu5
    mu6initial=params[15]#mu6
    dmude1=params[16]
    dmude2=params[17]
    # dmude3=params[18]
    # dmude4=params[19]
    dmude5=params[18]
    # dmude6=params[19]
    E= Id1 - Ry/nu3**2
    nu2=F(nu3)
    nu4= G(nu3)
    nu5=H(nu3)
    nu6=I(nu3)
    mu1= mu1initial + dmude1*(45932.1982 - E)/45932.1982
    mu2= mu2initial + dmude2*(45932.1982 - E)/45932.1982
    mu3= mu3initial
    mu4= mu4initial
    mu5= mu5initial + dmude5*(45932.1982 - E)/45932.1982
    mu6= mu6initial
    epsilon1= np.tan(np.pi*(nu2)) + mu1
    epsilon2= np.tan(np.pi*(nu2)) + mu2
    epsilon3= np.tan(np.pi*(nu3)) + mu3
    epsilon4= np.tan(np.pi*(nu4)) + mu4
    epsilon5= np.tan(np.pi*(nu5)) + mu5
    epsilon6= np.tan(np.pi*(nu6)) + mu6
    
    Bnum=epsilon1*(R34**2 - epsilon3*epsilon4) - R13*R14*R34
    Bnum+= epsilon4*R13**2 + R14**2*(epsilon3*epsilon4-R34**2)/epsilon4
    Bnum+= (R14*R34)**2/epsilon4 - R13*R14*R34
    Bdenom= R12*(epsilon3*epsilon4-R34**2) + R13*(R24*R34-epsilon4*R23)
    Bdenom+= (R14*R24*R34**2 - R14*R24*R34**2)/epsilon4 - R14*R24*epsilon3
    B= Bnum/Bdenom

    C=R14*R34 - epsilon4*R13 + B*(R24*R34- R23*epsilon4)
    C/=(epsilon3*epsilon4-R34**2)

    D= -(R14 + B*R24 + C*R34)/epsilon4
    
    E= -(R15+ R25*B)/epsilon5

    A= -(R16 + R26*B)/epsilon6

    a1= 1.0/np.sqrt(np.power(nu1,3)*np.power(1.0/np.cos(np.pi*nu1),2) + np.power(nu2,3)*np.power(B/np.cos(np.pi*nu2),2) + np.power(nu3,3)*np.power(C/np.cos(np.pi*nu3),2) + np.power(nu4,3)*np.power(D/np.cos(np.pi*nu4),2) + np.power(nu5,3)*np.power(E/np.cos(np.pi*nu5),2) + np.power(nu6,3)*np.power(A/np.cos(np.pi*nu6),2))
    a2=np.sqrt(np.power(nu2,3))*B*a1/np.cos(np.pi*nu2)
    a3= np.sqrt(np.power(nu3,3))*C*a1/np.cos(np.pi*nu3)
    a4= np.sqrt(np.power(nu4,3))*D*a1/np.cos(np.pi*nu4)
    a5= np.sqrt(np.power(nu5,3))*E*a1/np.cos(np.pi*nu5)
    a6= np.sqrt(np.power(nu6,3))*A*a1/np.cos(np.pi*nu6)
    a1=np.sqrt(np.power(nu1,3))*a1/np.cos(np.pi*nu1)
    return([a1,a2,a3,a4,a5,a6])

    
params=[]


#Final paramters
R12=-1.17766788e-01
R13=-6.65233833e-01
R14=2.06597487e-04
R15=-4.29172343e-01
R16=-2.60642176e-02
R23=3.27673824e-01
R24=-6.55741156e-01
R25=1.49013823e-01
R26=9.24143883e-02
R34=2.13643266e-01
mu1initial=-7.69912252e-01
mu2initial=-5.99619507e-01
mu3initial=1.08190813e+00
mu4initial=1.17740030e+00
mu5initial=4.21094332e-01
mu6initial=2.62710223e+00
dmude1=3.44520643e+00
dmude2=-1.07064669e-02
dmude5=-3.95553066e-01

initialparams=[R12,R13,R14,R15,R16,R23,R24,R25,R26,R34,mu1initial,mu2initial,mu3initial,mu4initial,mu5initial,mu6initial,dmude1,dmude2,dmude5]
expE= []
error=[]
n=[]

inputfile= open("beigangextended.csv","r")
for line in inputfile:
    nin,E,err= line.split()
    n.append(int(nin))
    expE.append(float(E))
    error.append(float(err))

finalparams=initialparams
#######Minimization
difference=1e8
oldchi=chisquared(finalparams,expE,error)
while (difference > 1e-5):
    try:
        finalparams= fmin(chisquared,initialparams,args=(expE,error),maxfun=1e8,maxiter=1e8,ftol=1e-3)
        chi=chisquared(finalparams,expE,error)
        difference= abs(oldchi-chi)
        oldchi=chi
        print difference
        initialparams=finalparams
    except KeyboardInterrupt:
        break

difference=2e-4
print finalparams
R12,R13,R14,R15,R16,R23,R24,R25,R26,R34,mu1initial,mu2initial,mu3initial,mu4initial,mu5initial,mu6initial,dmude1,dmude2,dmude5=finalparams#
oldchi=chisquared(finalparams,expE,error)
print oldchi/(len(expE) - len(finalparams))
Ecalc=np.zeros_like(expE)
for i in range(len(expE)):
    nud1exp= np.sqrt(Ry/(Id1-expE[i]))
    nud1calc= fminbound(diff,nud1exp-difference,nud1exp+difference,args=(finalparams,),xtol=1e-12)
    Ecalc[i]= Id1 - Ry/nud1calc**2
    print n[i],nud1calc,nud1exp,nud1calc-nud1exp,difference

for i in range(len(expE)):
    nud1exp= np.sqrt(Ry/(Id1-expE[i]))
    nud1calc= fminbound(diff,nud1exp-difference,nud1exp+difference,args=(finalparams,),xtol=1e-12)
    Ecalc[i]= Id1 - Ry/nud1calc**2
    print n[i],Ecalc[i],expE[i],(Ecalc[i]-expE[i]),(Ecalc[i]-expE[i])/error[i]




singletfile=open("singletextended.csv","r")
nutermsinglet=[]
for line in singletfile:
    nin,E,err= line.split()
    nud1exp= np.sqrt(Ry/(Id1-float(E)))    
    nud1calc= fminbound(diff,nud1exp-difference,nud1exp+difference,args=(finalparams,),xtol=1e-10)
    nutermsinglet.append(nud1calc)
    nus= F(nud1exp)%1

tripletfile=open("tripletextended.csv","r")
nutermtriplet=[]
for line in tripletfile:
    nin,E,err= line.split()
    nud1exp= np.sqrt(Ry/(Id1-float(E)))    
    nud1calc= fminbound(diff,nud1exp-difference,nud1exp+difference,args=(finalparams,),xtol=1e-12)
    nutermtriplet.append(nud1calc)
    print nud1exp
    nus=F(nud1exp)%1
nutermtriplet=np.asarray(nutermtriplet)

Rmatrix= np.matrix([[mu1initial,R12,R13,R14,R15,R16],[R12,mu2initial,R23,R24,R25,R26],[R13,R23,mu3initial,R34,0,0],[R14,R24,R34,mu4initial,0,0],[R15,R25,0,0,mu5initial,0],[R16,R26,0,0,0,mu6initial]])
eigR= linalg.eigvals(Rmatrix)
mu= (1.0 + np.arctan(eigR)/np.pi)%1
print "mu="
print mu

nutermsinglet=np.asarray(nutermsinglet)
pl.subplot(311)
for nud1exp in nutermsinglet:
    pl.axvline(x=nud1exp,color='b')


V1=np.sqrt(0.6)
V2=np.sqrt(0.4)
a1=np.zeros_like(nutermsinglet)
a2=np.zeros_like(nutermsinglet)
a3=np.zeros_like(nutermsinglet)
a4=np.zeros_like(nutermsinglet)
a5=np.zeros_like(nutermsinglet)
a6=np.zeros_like(nutermsinglet)
for i in range(len(nutermsinglet)):
    E= Id1 - Ry/nutermsinglet[i]**2
    nud1calc= fminbound(diff,nutermsinglet[i]-difference,nutermsinglet[i]+difference,args=(finalparams,),xtol=1e-12)
    Ecalc= Id1 - Ry/nud1calc**2
    print E, Ecalc, E-Ecalc
    mu1= mu1initial + dmude1*(45932.1982 - E)/45932.1982
    mu2= mu2initial + dmude2*(45932.1982 - E)/45932.1982
    mu3= mu3initial
    mu4= mu4initial
    mu5= mu5initial + dmude5*(45932.1982 - E)/45932.1982
    mu6= mu6initial
    nu1= F(nud1calc)
    nu2=F(nud1calc)
    nu4= G(nud1calc)
    nu5=H(nud1calc)
    nu6=I(nud1calc)
    epsilon1= np.tan(np.pi*(nu1)) + mu1
    epsilon2= np.tan(np.pi*(nu2)) + mu2
    epsilon3= np.tan(np.pi*(nud1calc)) + mu3
    epsilon4= np.tan(np.pi*(nu4)) + mu4
    epsilon5= np.tan(np.pi*(nu5)) + mu5
    epsilon6= np.tan(np.pi*(nu6)) + mu6
    a= coeffs(nud1calc,finalparams)
    V= np.matrix([[-V2,V1,0,0,0,0],[V1,V2,0,0,0,0],[0,0,-V2,V1,0,0],[0,0,V1,V2,0,0],[0,0,0,0,-V2,V1],[0,0,0,0,V1,V2]])
    a= np.dot(V,a)
    a1[i]= a[0,0]**2
    a2[i]= a[0,1]**2
    a3[i]= a[0,2]**2
    a4[i]= a[0,3]**2
    a5[i]= a[0,4]**2
    a6[i]= a[0,5]**2

pl.plot(nutermsinglet,a1,'b-')
pl.plot(nutermsinglet,a2,'r-')
pl.plot(nutermsinglet,a3,'c-')
pl.plot(nutermsinglet,a4,'y-')
pl.plot(nutermsinglet,a5,'k-')
pl.plot(nutermsinglet,a6,'g-')
# pl.subplot(413)
pl.subplot(312)

for nud1exp in nutermtriplet:
    pl.axvline(x=nud1exp,color='g')


a1=np.zeros_like(nutermtriplet)
a2=np.zeros_like(nutermtriplet)
a3=np.zeros_like(nutermtriplet)
a4=np.zeros_like(nutermtriplet)
a5=np.zeros_like(nutermtriplet)
a6=np.zeros_like(nutermtriplet)
for i in range(len(nutermtriplet)):
    E= Id1 - Ry/nutermtriplet[i]**2

    nud1calc= fminbound(diff,nutermtriplet[i]-difference,nutermtriplet[i]+difference,args=(finalparams,),xtol=1e-12)
    Ecalc= Id1 - Ry/nud1calc**2
    print E, Ecalc, E-Ecalc
    mu1= mu1initial + dmude1*(45932.1982 - Ecalc)/45932.1982
    mu2= mu2initial + dmude2*(45932.1982 - Ecalc)/45932.1982
    mu3= mu3initial
    mu4= mu4initial
    mu5= mu5initial + dmude5*(45932.1982 - E)/45932.1982
    mu6= mu6initial
    nu1=F(nud1calc)
    nu2=F(nud1calc)
    nu4= G(nud1calc)
    nu5=H(nud1calc)
    nu6=I(nud1calc)
    epsilon1= np.tan(np.pi*(nu1)) + mu1
    epsilon2= np.tan(np.pi*(nu2)) + mu2
    epsilon3= np.tan(np.pi*(nud1calc)) + mu3
    epsilon4= np.tan(np.pi*(nu4)) + mu4
    epsilon5= np.tan(np.pi*(nu5)) + mu5
    epsilon6= np.tan(np.pi*(nu6)) + mu6
    a= coeffs(nud1calc,finalparams)
    V= np.matrix([[-V2,V1,0,0,0,0],[V1,V2,0,0,0,0],[0,0,-V2,V1,0,0],[0,0,V1,V2,0,0],[0,0,0,0,-V2,V1],[0,0,0,0,V1,V2]])
    a= np.dot(V,a)
    a1[i]= a[0,0]**2
    a2[i]= a[0,1]**2
    a3[i]= a[0,2]**2
    a4[i]= a[0,3]**2
    a5[i]= a[0,4]**2
    a6[i]= a[0,5]**2

pl.plot(nutermtriplet,a1,'b-')
pl.plot(nutermtriplet,a2,'r-')
pl.plot(nutermtriplet,a3,'c-')
pl.plot(nutermtriplet,a4,'y-')
pl.plot(nutermtriplet,a5,'k-')
pl.plot(nutermtriplet,a6,'g-')

##################################
#g factors!


pl.subplot(313)
gfactorfile= open("gfactorsenergies.csv","r")
gfactorsexpsinglet=[]
gfactorsexptriplet=[]
for line in gfactorfile:
    a,b,c,d,e,f,g = line.split()
    nusinglet=np.sqrt(Ry/(Id1-float(b)))
    nutriplet=np.sqrt(Ry/(Id1-float(e)))    
    gfactorsexpsinglet.append((nusinglet,float(c),float(d)))
    gfactorsexptriplet.append((nutriplet,float(f),float(g)))

gfactorssinglet=np.zeros_like(nutermsinglet)
gfactorstriplet=np.zeros_like(nutermtriplet)

for i in range(len(gfactorsexpsinglet)):
    E= Id1 - Ry/(gfactorsexpsinglet[i][0]**2)
    nud1calc= fminbound(diff,gfactorsexpsinglet[i][0]-difference,gfactorsexpsinglet[i][0]+difference,args=(finalparams,),xtol=1e-12)
    Ecalc= Id1 - Ry/nud1calc**2
    mu1= mu1initial + dmude1*(45932.1982 - E)/45932.1982
    mu2= mu2initial + dmude2*(45932.1982 - E)/45932.1982
    mu3= mu3initial
    mu4= mu4initial
    mu5= mu5initial + dmude5*(45932.1982 - E)/45932.1982
    mu6= mu6initial
    nu1= F(nud1calc)
    nu2=F(nud1calc)
    nu4= G(nud1calc)
    nu5=H(nud1calc)
    nu6=I(nud1calc)
    epsilon1= np.tan(np.pi*(nu1)) + mu1
    epsilon2= np.tan(np.pi*(nu2)) + mu2
    epsilon3= np.tan(np.pi*(nud1calc)) + mu3
    epsilon4= np.tan(np.pi*(nu4)) + mu4
    epsilon5= np.tan(np.pi*(nu5)) + mu5
    epsilon6= np.tan(np.pi*(nu6)) + mu6
    a= coeffs(nud1calc,finalparams)
    V= np.matrix([[-V2,V1,0,0,0,0],[V1,V2,0,0,0,0],[0,0,-V2,V1,0,0],[0,0,V1,V2,0,0],[0,0,0,0,-V2,V1],[0,0,0,0,V1,V2]])
    a= np.dot(V,a)
    a1= a[0,0]**2
    a2= a[0,1]**2
    a3= a[0,2]**2
    a4= a[0,3]**2
    a5= a[0,4]**2
    a6= a[0,5]**2
    gfactorssinglet[i]= 1.0*(a1+a3+a5) + 7.0*(a2+a4)/6.0  +1.5*a6

for i in range(len(gfactorsexptriplet)):
    E= Id1 - Ry/(gfactorsexptriplet[i][0]**2)
    nud1calc= fminbound(diff,gfactorsexptriplet[i][0]-difference,gfactorsexptriplet[i][0]+difference,args=(finalparams,),xtol=1e-12)
    Ecalc= Id1 - Ry/nud1calc**2
    mu1= mu1initial + dmude1*(45932.1982 - E)/45932.1982
    mu2= mu2initial + dmude2*(45932.1982 - E)/45932.1982
    mu3= mu3initial
    mu4= mu4initial
    mu5= mu5initial + dmude5*(45932.1982 - E)/45932.1982
    mu6= mu6initial
    nu1= F(nud1calc)
    nu2=F(nud1calc)
    nu4= G(nud1calc)
    nu5=H(nud1calc)
    nu6=I(nud1calc)
    epsilon1= np.tan(np.pi*(nu1)) + mu1
    epsilon2= np.tan(np.pi*(nu2)) + mu2
    epsilon3= np.tan(np.pi*(nud1calc)) + mu3
    epsilon4= np.tan(np.pi*(nu4)) + mu4
    epsilon5= np.tan(np.pi*(nu5)) + mu5
    epsilon6= np.tan(np.pi*(nu6)) + mu6
    Rmatrix= np.matrix([[epsilon1,R12,R13,R14,R15,R16],[R12,epsilon2,R23,R24,R25,R26],[R13,R23,epsilon3,R34,0,0],[R14,R24,R34,epsilon4,0,0],[R15,R25,0,0,epsilon5,0],[R16,R26,0,0,0,epsilon6]])
    guess=np.array([1,1,1,1,1,1,0])/np.sqrt(6.0)
    a= coeffs(nud1calc,finalparams)
    V= np.matrix([[-V2,V1,0,0,0,0],[V1,V2,0,0,0,0],[0,0,-V2,V1,0,0],[0,0,V1,V2,0,0],[0,0,0,0,-V2,V1],[0,0,0,0,V1,V2]])
    a= np.dot(V,a)
    a1= a[0,0]**2
    a2= a[0,1]**2
    a3= a[0,2]**2
    a4= a[0,3]**2
    a5= a[0,4]**2
    a6= a[0,5]**2

    gfactorstriplet[i]= (a1+a3+a5) + 7.0*(a2+a4)/6.0 +1.5*a6

chisquared=0.0
expsinglet=np.zeros(7)
exptriplet=np.zeros(7)
calcsinglet=np.zeros(7)
calctriplet=np.zeros(7)
for i in range(7):
    pl.errorbar(gfactorsexpsinglet[i][0],gfactorsexpsinglet[i][1],gfactorsexpsinglet[i][2],marker='o',mfc='b')#,fmt='bo'
    pl.errorbar(gfactorsexptriplet[i][0],gfactorsexptriplet[i][1],gfactorsexptriplet[i][2],marker='o',mfc='r')#'ro',
    chisinglet=((gfactorssinglet[i]-gfactorsexpsinglet[i][1])/gfactorsexpsinglet[i][2])
    chitriplet=((gfactorstriplet[i]-gfactorsexptriplet[i][1])/gfactorsexptriplet[i][2])
    # print chisinglet,chitriplet
    chisquared+= chisinglet**2 + chitriplet**2
    expsinglet[i]=gfactorsexpsinglet[i][0]
    exptriplet[i]=gfactorsexptriplet[i][0]
    calcsinglet[i]=gfactorssinglet[i]
    calctriplet[i]=gfactorstriplet[i]
    print calcsinglet[i], gfactorsexpsinglet[i][1],calctriplet[i], gfactorsexptriplet[i][1]

pl.plot(expsinglet,calcsinglet,'b-')
pl.plot(exptriplet,calctriplet,'r-')

chisquared/=16
print chisquared

pl.show()
