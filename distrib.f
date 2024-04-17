c SimplesT Event Generator
      program distrib
      implicit none
      include 'constants.inc'
      double precision elra,eltail,pureel,pureel_noWidth,innn,strg
      real degrad,raddeg
      real Ebeam
      real w_min,w_max,q2_gev
      real theta,p_el
      real w,nu,cs_th
      real sgm_model,sig_rc
      double precision sig_cont,sig
      integer i,np
      character*100 fileout
c      data np/200/
      common /results/elra,eltail,pureel,pureel_noWidth,innn,strg
      degrad=pi/180.e0
      raddeg=180.e0/pi
c Input
      print*,'Enter beam energy [GeV]:'
      read*,Ebeam
      print*,'Enter Q2 [GeV2]:'
      read*,q2_gev
      print*,'Enter W_min and W_max [GeV]:'
      read*,w_min,w_max
      print*,'Enter number of W points:'
      read*,np
c------------
      print*,'Enter output filename:'
      read(*,'(A100)') fileout
      open(33,file=fileout(1:index(fileout,' ')-1),status='unknown')
c Loop
      do i=0,np
       w=w_min+real(i)*(w_max-w_min)/real(np)
       nu=(w**2-mp**2+q2_gev)/(2.e0*mp)
       p_el=Ebeam-nu
       cs_th=1.e0-q2_gev/(2.e0*Ebeam*p_el)
       theta=0.e0
       if(abs(cs_th).le.1.e0) theta=acos(cs_th)
       
       sig=sig_cont(dble(Ebeam),dble(p_el),dble(theta))
       sig_rc=sgm_model(Ebeam,theta,p_el)
      
c       write(33,20) p_el,sig,elra,eltail,innn,strg
c       write(33,20) w,(theta*raddeg),p_el,((sig_rc-elra)/sig),sig,(sig_rc-elra),elra,eltail,innn,strg
c My update to match my RC calculations
       write(33,20) w,(theta*raddeg),p_el,sig,pureel,elra,eltail,(sig_rc-elra-eltail),sig_rc,pureel_noWidth

       
c       write(33,20) w,(theta*raddeg),p_el,(sig_rc/sig),sig,(sig_rc-elra),elra,eltail,innn,strg
      enddo
      close(33)
10    format(f5.1,1pe12.4,1pe12.4,1pe12.4,1pe12.4,1pe12.4)
20    format(1pe12.4,1pe12.4,1pe12.4,1pe12.4,1pe12.4,1pe12.4,1pe12.4,1pe12.4,1pe12.4,1pe12.4)
      END
