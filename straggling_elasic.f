      double precision function straggling_elasic(Zp,M,t,Es,Ep,theta,Ziw,Zfw,tiw,tfw)
      implicit none
      double precision Z,M,t,tiw,tfw,Es,Ep,theta,Eb,Epr,peakpos
      double precision eta1,eta2,Ib,sig_elastic
      integer Zp,Ziw,Zfw
      straggling_elasic=0.0d+0
      peakpos=Es/(1.d0+Es/M*(1.d0-cos(theta)))
c      if(Ep.lt.(M/(1.d0-cos(theta)))) then
      if(Ep.lt.peakpos) then
      Z=dble(Zp)
      eta1=1.d0/(1.d0-Ep/M*(1.d0-cos(theta)))
      eta2=1.d0+Es/M*(1.d0-cos(theta))
      Eb=Ep*eta1
      Epr=Es/eta2
c      if(Eb.gt.Ep.and.Ep.gt.(Eb/(1.d0+2.d0*Eb/M))) then
      if(Eb.gt.Ep) then
c Straggling before and after scattering
c      straggling_elasic=
c     &(Ib(Es,Eb,tiw,dble(Ziw))+Ib(Es,Eb,(t/2.d0),Z))*(eta1**2)*
c     &sig_elastic(Eb,theta,M,Zp)+
c     &(Ib(Epr,Ep,(t/2.d0),Z)+Ib(Epr,Ep,tfw,dble(Zfw)))*
c     &sig_elastic(Es,theta,M,Zp)
c Straggling only before scattering
      straggling_elasic=
     &(Ib(Es,Eb,tiw,dble(Ziw))+Ib(Es,Eb,(t/2.d0),Z))*(eta1**2)*
     &sig_elastic(Eb,theta,M,Zp)
      endif
      endif
c      endif
      return
      end
c Ib function
      double precision function Ib(Eb,Ep,t,Z)
      implicit none
      double precision Eb,Ep,t,Z,bfunc,bt
      bt=bfunc(Z)*t
      Ib=bt/(Eb-Ep)*(Ep/Eb+3.d0/4.d0*((Eb-Ep)/Eb)**2)*(log(Eb/Ep))**bt
      return
      end
