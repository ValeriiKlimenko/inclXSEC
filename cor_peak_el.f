      subroutine smrad_el(Es,theta,smel)
      implicit none
      include 'constants8.inc'
      double precision Es,theta,M,smel
      double precision elasrad,elcor,sig_elastic
      double precision t,tiw,tfw,delta_E,elra,Ep3
      double precision dWp,sn2,beam_res,mom_res
      integer Z,A,Ziw,Zfw
      data dWp/0.08d0/,beam_res/1.d-4/,mom_res/1.0d-2/
c Straggling before and after scattering
c      data t/0.0056d0/,tiw/0.00023d0/,tfw/0.00012d0/
c Straggling only before scattering
      data t/0.0056d0/,tiw/0.00023d0/,tfw/0.d0/
c No straggling
c      data t/0.d0/,tiw/0.d0/,tfw/0.d0/
      data Z/1/,A/1/,Ziw/13/,Zfw/13/
c Init
      elra=0.d0
      smel=0.d0
c input values
      M=mp*dble(A)
      Ep3=Es/(1.d0+Es/M*(1.d0-cos(theta)))
      sn2=0.5d0*(1.d0-cos(theta))
c elastic peak width==beam + momentum resolution (cut two sigma)
c      delta_E=(dWp+dWp**2/(2.d0*M))/(1.d0+2.d0*Es/M*sn2)
C Valerii Change:
c      delta_E=2.d0*sqrt((beam_res*Es)**2+(mom_res*Ep3)**2)
      delta_E=0.00088
c cross section
      elcor=elasrad(Es,theta,t,delta_E,Z,A,Ziw,Zfw,tiw,tfw)
      elra=sig_elastic(Es,theta,M,Z)*elcor
      smel=elra
      return
      end
