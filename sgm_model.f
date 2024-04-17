      real function sgm_model(Ebeam,theta,Epr)
      implicit none
      include 'constants.inc'
      double precision Es,th,Ep,smel,smin
      real Ebeam,theta,Epr
      real nu,q2,x,s,w2
      sgm_model=0.e0
c Kinematic tests
      nu=Ebeam-Epr
      if(nu.lt.0.e0) return
      q2=2.e0*Ebeam*Epr*(1.e0-cos(theta))
      x=q2/(2.e0*mp*nu)
c      if(x.gt.1.0) return
      s=mp**2+2.e0*Ebeam*mp
      w2=mp**2+q2*(1.e0/x-1.e0)
c      print*,'w2= ',s,w2,nu,x,q2
      if(w2.lt.0.e0.or.w2.gt.s) return
      Es=dble(Ebeam)
      th=dble(theta)
      Ep=dble(Epr)
      call smrad(Es,th,Ep,smel,smin)
      sgm_model=real(smin)
      return
      end
