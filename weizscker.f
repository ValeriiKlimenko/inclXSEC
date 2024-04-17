      double precision function weizscker(Z,A)
      implicit none
      include 'constants8.inc'
      double precision Ad,Zd,Nd,Ebnd
      double precision a_V,a_S,a_C,a_A,delta_wzr
      integer Z,A
c data from Rohlf: Modern Physics from A to Z, James William Rohlf, Wiley, 1994
c      data a_V/15.75d0/,a_S/17.8/,a_C/0.711d0/,a_A/23.7d0/
c fit
      data a_V/15.8d0/,a_S/18.3/,a_C/0.714d0/,a_A/23.2d0/
      Ad=dble(A)
      Zd=dble(Z)
      Nd=dble(A-Z)
      Ebnd=0.d0
      if(A.gt.1) then
c      if(A.gt.12) then
      if(A.gt.4) then
      Ebnd=(a_V*Ad-a_S*Ad**(2.d0/3.d0)-a_C*Zd*(Zd-1.d0)/(Ad**(1.d0/3.d0))
     &-a_A*(Ad-2.d0*Zd)**2/Ad
     &+delta_wzr(Z,A))*1.d-3
      else
      if(Z.eq.1.and.A.eq.2) Ebnd=0.d0
      if(Z.eq.2.and.A.eq.4) Ebnd=0.d0
      if(Z.eq.3.and.A.eq.7) Ebnd=0.d0
      if(Z.eq.4.and.A.eq.9) Ebnd=0.d0
      if(Z.eq.5.and.A.eq.11) Ebnd=0.d0
      if(Z.eq.6.and.A.eq.12) Ebnd=0.d0
      endif
      endif
      weizscker=Zd*mp+Nd*mn-Ebnd
      return
      end

      double precision function delta_wzr(Z,A)
      implicit none
      double precision delta_0,a_P
      integer Z,A
c data from Rohlf: Modern Physics from A to Z, James William Rohlf, Wiley, 1994
c      data a_P/11.18d0/
c Fit
      data a_P/12.d0/
      delta_wzr=0.d0
      if(mod(A,2).ne.0) return
c      delta_0=a_P
      delta_0=a_P/(dble(A)**(1.d0/3.d0))
      if(mod(Z,2).eq.0) then
      delta_wzr=delta_0
      else
      delta_wzr=-delta_0
      endif
      return
      end
