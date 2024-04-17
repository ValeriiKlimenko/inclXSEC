      subroutine smrad(Es,theta,Ep,smel,smin)
      implicit none
      include 'constants8.inc'
      double precision Es,theta,Ep,M,eltail,smin,smel
      double precision peak,elasrad,elcor,innn,sig_elastic,varr
      double precision straggling_elasic,rad_tail_elastic
      double precision cor_cont,t,tiw,tfw,delta_E,delta,elra,Ep3,bfunc
      double precision dWp,sn2,weizscker,beam_res,mom_res,strg
      integer Z,A,Ziw,Zfw
      data dWp/0.01d0/,beam_res/1.d-4/,mom_res/0.5d-2/
c Straggling before and after scattering
c      data t/0.0056d0/,tiw/0.00023d0/,tfw/0.00012d0/
c      data t/0.0056d0/,tiw/0.00023d0/,tfw/0.00023d0/
c Straggling only before scattering
      data t/0.0056d0/,tiw/0.00023d0/,tfw/0.d0/
c No straggling
c      data t/0.d0/,tiw/0.d0/,tfw/0.d0/
      data Z/1/,A/1/,Ziw/13/,Zfw/13/
      common /results/elra,eltail,innn,strg
c Init
      elra=0.d0
      eltail=0.d0
      innn=0.d0
      strg=0.d0
      smel=0.d0
      smin=0.d0
c input values
c      M=mp*dble(A)
      M=weizscker(Z,A)
      Ep3=Es/(1.d0+Es/M*(1.d0-cos(theta)))
      sn2=0.5d0*(1.d0-cos(theta))
c Elastic peak width==beam + momentum resolution (cut two sigma)
c      delta_E=(dWp+dWp**2/(2.d0*M))/(1.d0+2.d0*Es/M*sn2)
c was before Valerii change:
c      delta_E=2.d0*sqrt((beam_res*Es)**2+(mom_res*Ep3)**2)

	  delta_E = 0.00088
c cross section
c      if(Ep.gt.(Ep3-delta_E)) then
      if(Ep.lt.(Ep3+2.d0*delta_E)) then
      varr=Ep3-Ep
      elcor=elasrad(Es,theta,t,delta_E,Z,A,Ziw,Zfw,tiw,tfw)
      elra=sig_elastic(Es,theta,M,Z)*elcor*peak(varr,(delta_E/2.d0))
      endif
c      else
c Small energy interval (Es,Es+delta), (Ep,Ep+delta)
c Should be large enough to safely neglect region IV of integration
      delta=Es*exp(-0.03d0/(bfunc(dble(z))*t/2.d0))
c      delta=Ep*exp(-0.03d0/(bfunc(dble(z))*t/2.d0))
c      print*,'d1: ',delta
c but on which the cross section is constant
      delta=max(delta,(0.015d0*((Ep3-delta_E)-Ep)))
c      print*,'d2: ',delta,(0.015d0*((Ep3-delta_E)-Ep))
      innn=cor_cont(Es,Ep,M,theta,t,delta,Z,Ziw,Zfw,tiw,tfw)
c Tails do not need Delta parameters
      if(Ep.lt.Ep3) then
      strg=straggling_elasic(Z,M,t,Es,Ep,theta,Ziw,Zfw,tiw,tfw)
      eltail=rad_tail_elastic(Es,Ep,M,theta)
     &+strg
      endif
c      endif
      smel=elra+eltail
      smin=smel+innn
      return
      end

c##################################################################
c Cole Smith program for radiative correction in the elastic peak			     
      double precision function elasrad(es,theta,t,delta_E,
     &Zvvv,Atnum,Ziw,Zfw,tiw,tfw)
      implicit none
      include 'constants8.inc'
      double precision spence,arg,weizscker
      double precision m,z,t,es,theta,delta_E
      double precision cst1,eta,eel,qs
      double precision znuc,deltac(27),del_mo,delta_t
      double precision arg11,arg15,arg19,arg23
      double precision epr,e1,e3,e4,beta4
      double precision radcor
      double precision snth,tiw,tfw
      double precision stragg_peak
      integer Zvvv,Atnum,Ziw,Zfw,idel
      data z/1.0d+0/
c      m=mp*dble(Atnum)
      m=weizscker(Zvvv,Atnum)
      znuc      = dble(Zvvv)
      snth	= sin(theta)
      cst1	= 1.d0-cos(theta)
      eel	= es/(1.d0+es*cst1/m)
      epr	= es+m-eel
      e1	= es
      e3	= eel
      e4	= epr
      beta4	= sqrt(e4**2-m**2)/e4
      eta	= es/eel
      qs	= 2.d0*es*eel*cst1
      deltac(1)=28.d0/9.d0-13.d0/6.d0*log(qs/(me**2))
      deltac(2)=(log(qs/(me**2))-1.d0+2.d0*znuc*log(eta))*
     &(2.d0*log(e1/delta_E)-3.d0*log(eta))
      arg=(e3-e1)/e3
      deltac(3)=-spence(arg)
      deltac(4)=znuc**2*log(e4/m)
      deltac(5)=znuc**2*log(m/eta/delta_E)*
     &(log((1.d0+beta4)/(1.d0-beta4))/beta4-2.d0)
      deltac(6)=znuc**2/beta4*(log((1.0+beta4)/(1.d0-beta4))*
     &log((e4+m)/2.d0/m)/2.d0)
      arg=sqrt((e4-m)*(1.d0+beta4)/(e4+m)/(1.d0-beta4))
      deltac(7)=-znuc**2*spence(-arg)/beta4
      arg=(e3-m)/e1
      deltac(8)=znuc*spence(arg)
      arg=(m-e3)*m/(2.d0*e3*e4-m*e1)
      deltac(9)=-znuc*spence(arg)
      arg=2.d0*e3*(m-e3)/(2.d0*e3*e4-m*e1)
      deltac(10)=znuc*spence(arg)
      arg11=(2.d0*e3*e4-m*e1)/e1/(m-2.d0*e3)
      arg11=abs(arg11)
      deltac(11)=znuc*log(arg11)*log(m/2.d0/e3)
      arg=(e3-e4)/e3
      deltac(12)=-znuc*spence(arg)
      arg=(e4-e3)*m/(2.d0*e1*e4-m*e3)
      deltac(13)=znuc*spence(arg)
      arg=2.d0*e1*(e4-e3)/(2.d0*e1*e4-m*e3)
      deltac(14)=-znuc*spence(arg)
      arg15=(2.d0*e1*e4-m*e3)/e3/(m-2.d0*e1)
      arg15=abs(arg15)
      deltac(15)=-znuc*log(arg11)*log(m/2.d0/e1)
      arg=(e1-m)/e1
      deltac(16)=-znuc*spence(arg)
      arg=(m-e1)/e1
      deltac(17)=znuc*spence(arg)
      arg=2.d0*(m-e1)/m
      deltac(18)=-znuc*spence(arg)
      arg19=abs(m/(2.d0*e1-m))
      deltac(19)=-znuc*log(arg19)*log(m/2.d0/e1)
      arg=(e3-m)/e3
      deltac(20)=znuc*spence(arg)
      arg=(m-e3)/e3
      deltac(21)=-znuc*spence(arg)
      arg=2.d0*(m-e3)/m
      deltac(22)=znuc*spence(arg)
      arg23=abs(m/(2.d0*e3-m))
      deltac(23)=znuc*log(arg23)*log(m/2.d0/e3)
      arg=(e1-e3)/e1
      deltac(24)=-spence(arg)
      arg=(e4-m)*(1.d0-beta4)/(e4+m)/(1.d0+beta4)
      arg=sqrt(arg)
      deltac(25)=znuc**2*spence(arg)/beta4
      arg=(e4-m)/(e4+m)
      arg=sqrt(arg)
      deltac(26)=-znuc**2*spence(arg)/beta4
      deltac(27)=znuc**2*spence(-arg)/beta4
      del_mo=0.0d+0
      do idel=1,27
      del_mo=del_mo+deltac(idel)
      enddo
      del_mo=-alpha*del_mo/pi
c Straggling correction delta_t
      delta_t=stragg_peak(znuc,t,es,eel,eta,delta_E,Ziw,Zfw,tiw,tfw)
c      delta_t=0.0d+0
      radcor=del_mo+delta_t
      elasrad=exp(radcor)
      return
      end
c######################################################
      double precision function stragg_peak(z,t,es,eel,eta,delta_E,ziw,zfw,tiw,tfw)
      implicit none
      double precision z,t,es,eel,eta,delta_E,bt,bfunc
      double precision tiw,tfw,btiw,btfw
      integer ziw,zfw
      bt=bfunc(z)*t
      btiw=bfunc(dble(ziw))*tiw
      btfw=bfunc(dble(zfw))*tfw
c Straggling before and after scattering
c      stragg_peak=-((btiw+bt/2.d0)*log(es/(eta**2*delta_E))+
c     &(btfw+bt/2.d0)*log(eel/delta_E))
c Straggling only before scattering
      stragg_peak=-((btiw+bt/2.d0)*log(es/(eta**2*delta_E)))
      return
      end
c B-function
      double precision function bfunc(z)
      implicit none
      double precision z,xi
      xi=log(1440.d0)-2.d0/3.d0*log(z)
      xi=xi/(log(183.d0)-log(z)/3.d0)
      bfunc=4.d0/3.d0*(1.d0+((z+1.d0)/(z+xi))/(log(183.d0)-
     &log(z)/3.d0)/9.d0)
      return
      end
c Spence function
      double precision function spence(x)
      double precision x,EPS,dgauss,VAL,null
      external VAL
      data EPS/1.0d-8/,null/0.0d+0/
      spence=dgauss(VAL,null,x,EPS)
      return
      end
c Integrand of Spence function
      double precision function VAL(y)
      double precision y
      VAL=-log(abs(1.d0-y))/y
      return
      end
