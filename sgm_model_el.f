      real function sgm_model_el(Ebeam,theta)
      implicit none
      include 'constants.inc'
      double precision Es,th,smel
      real Ebeam,theta
      sgm_model_el=0.e0
c Kinematic tests
      if(Ebeam.lt.0.1e0) return
      if(theta.lt.0.e0.or.theta.gt.pi) return
      Es=dble(Ebeam)
      th=dble(theta)
      call smrad_el(Es,th,smel)
      sgm_model_el=real(smel)
      return
      end
