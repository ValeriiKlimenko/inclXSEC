      double precision function rad_tail_elastic(Es,Ep,M,theta)
      implicit none
      include 'constants8.inc'
      double precision s(0:3),p(0:3),pin(0:3),u(0:3)
      double precision Es,Ep,theta,sp,mods,modp,modu
      double precision costs,costp,sints,sintp,u2,phi
      double precision M,haha,prec,dgauss
      double precision M1,Es1,Ep1,theta1,up,down
      double precision bound(7),integ(8),Gp,Gs
      integer i1
      common/params/Es1,M1,Ep1,theta1,sp,mods,modp,modu,
     &u,costs,costp,sints,sintp,u2
      data phi/0.0d+0/,prec/1.0d-8/,down/-1.0d+0/,up/1.0d+0/
      external haha
      rad_tail_elastic=0.d0
      M1=M
      Es1=Es
      Ep1=Ep
      theta1=theta
      s(0)=Es
      s(1)=0.0d+0
      s(2)=0.0d+0
      s(3)=sqrt(Es**2-me**2)
      p(0)=Ep
      p(1)=sqrt(Ep**2-me**2)*sin(theta)*cos(phi)
      p(2)=sqrt(Ep**2-me**2)*sin(theta)*sin(phi)
      p(3)=sqrt(Ep**2-me**2)*cos(theta)
      pin(0)=M
      pin(1)=0.0d+0
      pin(2)=0.0d+0
      pin(3)=0.0d+0
      do i1=0,3
      u(i1)=s(i1)+pin(i1)-p(i1)
      enddo
      sp=s(0)*p(0)-s(1)*p(1)-s(2)*p(2)-s(3)*p(3)
      mods=sqrt(s(1)**2+s(2)**2+s(3)**2)
      modp=sqrt(p(1)**2+p(2)**2+p(3)**2)
      modu=sqrt(u(1)**2+u(2)**2+u(3)**2)
      if(mods.gt.0.d0.and.modp.gt.0.d0.and.modu.gt.0.d0) then
      costs=(s(1)*u(1)+s(2)*u(2)+s(3)*u(3))/mods/modu
      costp=(p(1)*u(1)+p(2)*u(2)+p(3)*u(3))/modp/modu
      sints=sqrt(max(0.d0,(1.d0-costs**2)))
      sintp=sqrt(max(0.d0,(1.d0-costp**2)))
c      sints=sin(acos(costs))
c      sintp=sin(acos(costp))
      u2=u(0)**2-u(1)**2-u(2)**2-u(3)**2
      if(Ep.gt.0.d0.and.Es.gt.0.d0) then
      Gp=sqrt(me/Ep)
      Gs=sqrt(me/Es)
c cos(theta+-gamma) should have theta_p>gamma_p=sqrt(me/Ep)
      if(acos(costp).gt.Gp) then
      bound(1)=costp*cos(Gp)-sintp*sin(Gp)
      bound(2)=costp
      bound(3)=costp*cos(Gp)+sintp*sin(Gp)
      else
      bound(1)=costp-(costp-down)/2.d0
      bound(2)=costp
      bound(3)=costp+(costs-costp)/2.d0
      endif
c cos(theta+-gamma) should have theta_s>gamma_s=sqrt(me/Es)
      if(acos(costs).gt.Gs) then
      bound(4)=costs*cos(Gs)-sints*sin(Gs)
      bound(5)=costs
      bound(6)=costs*cos(Gs)+sints*sin(Gs)
      bound(7)=up-(up-bound(4))/2.d0
      else
      bound(4)=costs-(costs-bound(3))/2.d0
      bound(5)=costs
      bound(6)=costs+(up-costs)/2.d0
      bound(7)=up-(up-bound(6))/2.d0
      endif
      
c      print*,'bound: ',bound(1),bound(2),bound(3),bound(4),bound(5),bound(6),bound(7)
c      print*,'intervals: ',(bound(1)-down),(bound(2)-bound(1)),(bound(3)-bound(2)),
c     &(bound(4)-bound(3)),(bound(5)-bound(4)),(bound(6)-bound(5)),(bound(7)-bound(6)),
c     &(up-bound(7))
c      print*,costs,costs*cos(Gs),sints*sin(Gs),acos(costs),Gs
      
      integ(1)=dgauss(haha,down,bound(1),prec)
c      print*,'i1: ',integ(1)
      integ(2)=dgauss(haha,bound(1),bound(2),prec)
c      print*,'i2: ',integ(2)
      integ(3)=dgauss(haha,bound(2),bound(3),prec)
c      print*,'i3: ',integ(3)
      integ(4)=dgauss(haha,bound(3),bound(4),prec)
c      print*,'i4: ',integ(4)
      integ(5)=dgauss(haha,bound(4),bound(5),prec)
c      print*,'i5: ',integ(5)
      integ(6)=dgauss(haha,bound(5),bound(6),prec)
c      print*,'i6: ',integ(6)
      integ(7)=dgauss(haha,bound(6),bound(7),prec)
c      print*,'i7: ',integ(7)
      integ(8)=dgauss(haha,bound(7),up,prec)
c      print*,'i8: ',integ(8)
c      stop
      
      rad_tail_elastic=igev2mub*(integ(1)+integ(2)+integ(3)+
     &integ(4)+integ(5)+integ(6)+integ(7)+integ(8))
      endif
      endif
      return
      end
c######################################################
      double precision function haha(costk)
      implicit none
      include 'constants8.inc'
      double precision u(0:3)
      double precision q2,Es,Ep,theta,sp,mods,modp,modu
      double precision costs,costp,sints,sintp,u2
      double precision M,costk,zero,vud
      double precision sintk,omega,a,b,apr,bpr,nu,nupr,norm,differ
      double precision md(5),nd(3),mdt,ndt
      double precision Fformfact,Gformfact
      common/params/Es,M,Ep,theta,sp,mods,modp,modu,
     &u,costs,costp,sints,sintp,u2
      zero=(Es*costp/mods-Ep*sints/modp)/sin(theta)
      haha=0.0d+0
c Zero divide by zero ambiguity point (Mo & Tsai, Rev.Mod.Phys.41, 1969, Appendix D)
      if(costk.ne.zero) then
      sintk=sqrt(max(0.d0,(1.d0-costk**2)))
      omega=(u2-M**2)/2.d0/(u(0)-modu*costk)
      q2=-2.d0*(me**2-sp-omega*(Es-Ep)+omega*modu*costk)
      a=omega*(Ep-modp*costp*costk)
      b=-omega*modp*sintp*sintk
      apr=omega*(Es-mods*costs*costk)
      bpr=-omega*mods*sints*sintk
      vud=(Ep*mods*sints-Es*modp*sintp+mods*modp*sin(theta)*costk)
c      if((apr**2-bpr**2).gt.0.d0) then
      if((apr**2-bpr**2).gt.0.d0.and.(a**2-b**2).gt.0.d0
     &.and.abs(vud).gt.0.d0.and.abs(omega).gt.0.d0) then
      
c      if(vud.gt.0.d0) print*,'vud<0: ',vud
c      if(omega.lt.0.d0) then
c      print*,'omega<0: ',omega
c      stop
c      endif
      
c      nu=-modp*sintp/omega/(Ep*mods*sints-
c     &Es*modp*sintp+mods*modp*sin(theta)*costk)
      nu=-modp*sintp/omega/vud
c      nupr=-mods*sints/omega/(Ep*mods*sints-
c     &Es*modp*sintp+mods*modp*sin(theta)*costk)
      nupr=-mods*sints/omega/vud
      norm=(alpha**3)*Ep/Es/((2.d0*pi)**2)/M
      differ=omega/2.d0/(q2**2)/(u(0)-modu*costk)
      md(1)=-2.d0*pi*a*(me**2)*(2.d0*Es*(Ep+omega)
     &-q2/2.d0)/((a**2-b**2)**1.5d0)
      md(2)=-2.d0*pi*apr*(me**2)*(2.d0*Ep*(Es-omega)
     &-q2/2.d0)/((apr**2-bpr**2)**1.5d0)
      md(3)=4.d0*pi*(nu/sqrt(a**2-b**2)-
     &nupr/sqrt(apr**2-bpr**2))*(me**2*(sp-omega**2)+sp*
     &(2.d0*Es*Ep-sp+omega*(Es-Ep)))
      md(4)=2.d0*pi*(2.d0*(Es*Ep+Es*omega+Ep**2)-q2/2.d0-
     &sp-me**2)/sqrt(a**2-b**2)
      md(5)=-2.d0*pi*(2.d0*(Es*Ep-Ep*omega+Es**2)-q2/2.d0-
     &sp-me**2)/sqrt(apr**2-bpr**2)
      mdt=M**2*Fformfact(q2)*(md(1)+md(2)-
     &4.d0*pi+md(3)+md(4)+md(5))
      nd(1)=(2.d0*pi*a/((a**2-b**2)**1.5d0)+2.d0*pi*apr/
     &((apr**2-bpr**2)**1.5d0))*(me**2)*(2.d0*me**2-q2)
      nd(2)=8.d0*pi*(nu/sqrt(a**2-b**2)-
     &nupr/sqrt(apr**2-bpr**2))*sp*(sp-2.d0*me**2)
      nd(3)=2.d0*pi*(1.d0/sqrt(a**2-b**2)-1.d0/
     &sqrt(apr**2-bpr**2))*(2.d0*sp+2.d0*me**2+q2)
      ndt=Gformfact(q2)*(nd(1)+8.d0*pi+nd(2)+nd(3))
      haha=norm*differ*(mdt+ndt)
c      print*,haha,q2,podchlen1c,(a**2-b**2),(apr**2-bpr**2),nu,nupr,
c     &(Ep*mods*sints-Es*modp*sintp+mods*modp*sin(theta)*costk),
c     &(Ep*mods*sints-Es*modp*sintp+mods*modp*sin(theta)*costk)
c      if(haha.ne.haha) print*,'err: ',vud,costk,zero
      
      endif
      endif
      return
      end
c F Form factor
      double precision function Fformfact(q2)
      implicit none
      include 'constants8.inc'
      double precision q2,gel,gmag,tau
      tau=q2/(4.d0*mp**2)
      Fformfact=4.d0*((gel(q2))**2+tau*(gmag(q2))**2)/(1.d0+tau)
      return
      end
c G Form-factor
      double precision function Gformfact(q2)
      implicit none
      double precision q2,gmag
      Gformfact=q2*(gmag(q2))**2
      return
      end
