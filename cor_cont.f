      FUNCTION DGAUSS(F,A,B,EPS)
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DGAUSS')
      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
    3 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
    4 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
       H=0
       write(*,*)' DGAUSS: D103.1, TOO HIGH ACCURACY REQUIRED'
       stop
c       GO TO 99
      END IF
   99 DGAUSS=H
      RETURN
      END
*


      double precision function cor_cont(Es,Ep,M,theta,t,delta,Zp,Ziw,Zfw,tiw,tfw)
      implicit none
      include 'constants8.inc'
      double precision bfunc,z,b,bt,delta_t,Es,Ep,M,theta
      double precision norm,cn(3),spence
      double precision delta_r,tr,fs,fp,Esmin,Epmax
      double precision integ(2),sss,ppp,dgauss,sig_cont
      double precision t,delta,q2,int1es,Es1,Ep1,theta1
      double precision btiw,btfw,tiw,tfw,int2ep,prec,lgme
      integer Zp,Ziw,Zfw
      common/int/Es1,Ep1,theta1,norm,q2,lgme,btiw,btfw,bt,fs,fp
      external int1es,int2ep
      data prec/1.d-8/
      Es1=Es
      Ep1=Ep
      theta1=theta
      z=dble(Zp)
      b=bfunc(z)
      bt=b*t
      btiw=bfunc(dble(Ziw))*tiw
      btfw=bfunc(dble(Zfw))*tfw
      q2=2.d0*Es*Ep*(1.d0-cos(theta))
      norm=alpha/pi
      lgme=log(q2/(me**2))
c Straggling before and after scattering
c      delta_t=-((btiw+bt/2.d0)*log(Es/delta)+
c     &(btfw+bt/2.d0)*log(Ep/delta))
c Straggling only before scattering
      delta_t=-((btiw+bt/2.d0)*log(Es/delta))
c No straggling
c      delta_t=0.d0
c Schwinger correction
      cn(1)=28.d0/9.d0-13.d0/6.d0*lgme
      cn(2)=(log(Es/delta)+log(Ep/delta))*(lgme-1.d0)
      cn(3)=-spence(-(Es-Ep)/Ep)-spence((Es-Ep)/Es)
      delta_r=-norm*(cn(1)+cn(2)+cn(3))
c Bremsstahlung
      tr=norm/b*(lgme-1.d0)
      Esmin=(mpi**2+2.d0*M*mpi+2.d0*M*Ep)/2.d0/
     &(M-Ep*(1.d0-cos(theta)))
      Epmax=(2.d0*M*Es-2.d0*M*mpi-mpi**2)/2.d0/
     &(M+Es*(1.d0-cos(theta)))
      fs=b*tr+btiw+bt/2.d0
      fp=b*tr+btfw+bt/2.d0
c electron radiate before scattering
      if((Es-delta).gt.Esmin) then
      integ(1)=dgauss(int1es,Esmin,(Es-delta),prec)
      else
      integ(1)=0.d0
      endif
c electron radiate after scattering
      if((Ep+delta).lt.Epmax) then
      integ(2)=dgauss(int2ep,(Ep+delta),Epmax,prec)
      else
      integ(2)=0.d0
      endif
c Cross section
      sss=fs/2.d0
      ppp=fp/2.d0
      cor_cont=sig_cont(Es,Ep,theta)*exp(delta_t+delta_r)+
     &integ(1)*(delta/Ep)**ppp+integ(2)*(delta/Es)**sss
c      print*,cor_cont,sig_cont(Es,Ep,theta),delta_t,delta_r,integ(1),integ(2)
      return
      end
c Electron Radiates before scattering
      double precision function int1es(Espr)
      implicit none
      double precision Es,Ep,theta,norm,q2,lgme,btiw,btfw,bt,fs,fp
      double precision xs,ts,Espr,sig_cont
      common/int/Es,Ep,theta,norm,q2,lgme,btiw,btfw,bt,fs,fp
      xs=Espr/Es
      ts=norm*(0.5d0*(1.d0+xs**2)*lgme-xs)
      int1es=(ts+(btiw+bt/2.d0)*(xs+3.d0/4.d0*(1.d0-xs)**2))*
     &(log(1.d0/xs))**fs*sig_cont(Espr,Ep,theta)/(Es-Espr)
      return
      end
c Electron Radiates after scattering
      double precision function int2ep(Eppr)
      implicit none
      double precision Es,Ep,theta,norm,q2,lgme,btiw,btfw,bt,fs,fp
      double precision xp,tp,Eppr,sig_cont
      common/int/Es,Ep,theta,norm,q2,lgme,btiw,btfw,bt,fs,fp
      xp=Ep/Eppr
      tp=norm*(0.5d0*(1.d0+xp**2)*lgme-xp)
c Straggling before and after scattering
c      int2ep=(tp+(btfw+bt/2.d0)*(xp+3.d0/4.d0*(1.d0-xp)**2))*
c     &(log(1.d0/xp))**fp*sig_cont(Es,Eppr,theta)/(Eppr-Ep)
c Straggling only before scattering
      int2ep=(tp)*
     &(log(1.d0/xp))**fp*sig_cont(Es,Eppr,theta)/(Eppr-Ep)
      return
      end
