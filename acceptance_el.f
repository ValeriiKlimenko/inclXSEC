      subroutine acceptance_el(p,theta,t_current,delta_phi)
      implicit none
      real p,theta,t_current,accvb1,accvb1_new,pm,thm,delta_phi
      pm=p
      thm=theta
c      delta_phi=accvb1(pm,thm,t_current)
      delta_phi=accvb1_new(pm,thm,t_current)
      if(delta_phi.gt.30.e0) delta_phi=30.e0
      if(delta_phi.lt.0.e0) delta_phi=0.e0
      return
      end

      real function accvb1(p,theta,t_current)
      IMPLICIT NONE
      REAL p,theta,t_current,acc
      REAL t_max
      REAL phi0_el
      REAL theta0_el
      REAL thetas_el
      REAL p_shift,cel_ex,pel_ex
      REAL theta_min,theta_max,delta_phi,zexp
      REAL d2r
c Enlarged cuts
      data phi0_el/30.e0/
      data theta0_el/10.6e0/
      data thetas_el/15.e0/
      data theta_max/52.e0/
      data p_shift/0.25e0/
      data pel_ex/0.333e0/
      data cel_ex/0.25e0/
c Standard cuts
c      data phi0_el/30./
c      data theta0_el/12.2/
c      data thetas_el/21.5/
c      data theta_max/50./
c      data p_shift/0.15/
c      data pel_ex/0.416667/
c      data cel_ex/0.25/
c Init
      data t_max/3375.e0/
      data d2r/0.0174532925e0/
      
      acc=0.0e+0
      
        theta_min = theta0_el+thetas_el/(p*t_max/t_current+p_shift)
        if(theta.gt.theta_min.and.theta.lt.theta_max)then
          zexp = cel_ex*(p*t_max/t_current)**pel_ex
          delta_phi = phi0_el*sin((theta-theta_min)*d2r)**zexp
          Acc=delta_phi
        endif
c	print*,theta_min,delta_phi
      
      accvb1 = acc
      return
      end

      real function accvb1_new(p,theta,t_current)
      IMPLICIT NONE
      REAL p,theta,t_current,acc
      REAL t_max
      REAL phi0_el
      REAL theta0_el
      REAL thetas_el
      REAL p_shift,cel_ex,pel_ex,alpha_el
      REAL theta_min,theta_max,delta_phi,zexp
      REAL d2r
c Init
      data t_max/3375.e0/
      data d2r/0.0174532925e0/
c Enlarged cuts
      data phi0_el/23.e0/
      data theta0_el/6.5e0/
      data theta_max/52.e0/
      data thetas_el/36.5e0/
      data p_shift/0.8e0/
      data pel_ex/0.55e0/
      data cel_ex/0.40317e0/
      data alpha_el/5.2171e0/
c Standard cuts
c      data phi0_el/20.281e0/
c      data theta0_el/8.1908e0/
c      data theta_max/50.e0/
c      data thetas_el/36.601e0/
c      data p_shift/0.72590e0/
c      data pel_ex/0.63976e0/
c      data cel_ex/0.40317e0/
c      data alpha_el/5.2171e0/
      acc=0.0e+0
        theta_min=theta0_el+thetas_el/(p*t_max/t_current+p_shift)
        if(theta.gt.theta_min.and.theta.lt.theta_max)then
         zexp=cel_ex*(p*t_max/t_current)**pel_ex
	 if(alpha_el*(theta-theta_min).lt.90.e0) then
         delta_phi=phi0_el*(sin(alpha_el*(theta-theta_min)*d2r))**zexp
	 else
	 delta_phi=phi0_el
	 endif
         acc=delta_phi
        endif
      accvb1_new = acc
      return
      end

      subroutine newphi(phi,phinew,sec)
      IMPLICIT NONE
      real phinew,phi
      integer sec
        phinew = phi
        sec=1
      if (phi.gt.330.e0) then
        phinew = phi-360.e0
        sec=1
      elseif (phi.gt.30.e0.and.phi.le.90.e0) then
        phinew = phi-60.e0
        sec=2
      elseif (phi.gt.90.e0.and.phi.le.150.e0) then
        phinew = phi-120.e0
        sec=3
      elseif (phi.gt.150.e0.and.phi.le.210.e0) then
        phinew = phi-180.e0
        sec=4
      elseif (phi.gt.210.e0.and.phi.le.270.e0) then
        phinew = phi-240.e0
        sec=5
      elseif (phi.gt.270.e0.and.phi.le.330.e0) then
        phinew = phi-300.e0
        sec=6
      endif
      return
      end
