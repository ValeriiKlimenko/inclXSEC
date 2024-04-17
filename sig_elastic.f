      double precision function sig_elastic(Eb,theta,M,Zp)
      implicit none
      include 'constants8.inc'
      double precision Eb,theta,M,gel,gmag,cn2,Z
      double precision sn2,tn2,Epr,q2,gep2,gmp2,rmott,tau
      integer Zp
      Z=dble(Zp)
      sn2=0.5d0*(1.d0-cos(theta))
      cn2=1.d0-sn2
      tn2=sn2/cn2
c Energy of the scattered electron from recoil of the target: 
      Epr=Eb/(1.d0+2.d0*Eb*sn2/M)
c Q squared, neglecting electron mass:
      q2=4.d0*Eb*Epr*sn2
c Form factors
      gep2=(gel(q2))**2
      gmp2=(gmag(q2))**2
c Mott
      rmott=igev2mub*Z**2*cn2/((2.d0*Eb*sn2*ialpha)**2) ! microbarns/sr
c Recoil correction 
      rmott=rmott*(Epr/Eb)
c Form factor
      tau=q2/(4.d0*M**2)
      sig_elastic=rmott*((gep2+tau*gmp2)/(1.d0+tau)+2.d0*tau*(gmp2)*tn2)
      return
      end
c Proton Electric Form-factor
      double precision function gel(q2)
      implicit none
      double precision q2,q1
      q1=sqrt(q2)
c Dipole fit
c      gel=1.0/(1.0+q2/0.71)**2
c Bosted parameterization of form factors
c      gel=1.0/(1.0+0.62*q1+0.68*q2+2.8*q1**3+0.83*q1**4)
c Anghinolfi
c      if(q2.gt.(0.4)) then
c      gel=(1.0-0.25*(q2-0.4)/1.6)/(1.0+q2/0.71)**2
c      else
c      gel=1.0/(1.0+q2/0.71)**2
c      endif
c Jlab modification to Bosted parameterization
      if(q2.gt.(0.4)) then
      gel=(1.0-0.2*(q2-0.4)/1.6)/
     &(1.0+0.62*q1+0.68*q2+2.8*q1**3+0.83*q1**4)
      else
      gel=1.0/(1.0+0.62*q1+0.68*q2+2.8*q1**3+0.83*q1**4)
      endif
      return
      end
c Proton Magnetic Form-factor
      double precision function gmag(q2)
      implicit none
      double precision q2,q1,gel
      q1=sqrt(q2)
c Dipole fit
c      gmag=2.793*gel(q2)
c Bosted parameterization of form factors
      gmag=2.7928/(1.0+0.35*q1+2.44*q2+0.5*q1**3+1.04*q1**4+
     &0.34*q1**5)
c Jlab modification to dipol fit
c      gmag=2.793/(1.0+q2/0.71)**2
      return
      end
c Peak cross section
      double precision function peak(x,width)
      implicit none
      include 'constants8.inc'
      double precision x,norm,width,pow
      peak=0.0d+0
      norm=sqrt(2.d0*pi)*width
      pow=x**2/(2.d0*width**2)
      if(pow.lt.5.d1.and.norm.gt.1.d-10) then
      peak=exp(-pow)/norm
      endif
      return
      end
