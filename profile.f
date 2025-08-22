C2345678901234567890123456789012345
     *    789012345678901234567890123456789012
C
      program profile
C
C  This program expand a profile from the zonal coefficients into the
C    anomaly, for a specified latitude of closest approach.
C    The input and output are in mgal.
C
C
      implicit none
      real*8 C1lm(5150),C2lm(5150),pi,th,norm,GM,a,cth,sth,fi
      real*8 N(5150),P(5150),R,j1,j2,j3,lat,latcr
      integer*4 l,lmin,lmax,m,i,dl,j,k
      character*20 inname,outname
C
      pi = dacos(-1.d0)
C	GM = (6.67241d-11)*(4.799844d22)*(1.d5)
C	a = 1560.8d3
      a = 1450.d3
C
C Comparing Gael's position vector to his C/A distance, we see a
C   variation between 1560 and 1565, with an average of 1562.8 km.
C
      write(*,*)'Enter the coefficient file name:'
      read(*,*)inname
      print *," "
      open(unit=10,file=inname,status='unknown',blank='null')
      lmax = 100
      print *,'Input file is degree-and-order ',lmax
      print *," "
      write(*,*)'Enter latitude of closest approach:'
      read(*,*)j1
      latcr = dble(j1)
      print *," "
C
C  Read in the coefficients from the file:
      do l = 1,lmax
        do m = 0,l
          i = l*(l+1)/2 + m
          read(10,*)j,k,C1lm(i),C2lm(i)
        enddo
      enddo
      close(10)
C
C	write(*,*)'Enter the maximum spherical harmonic degree:'
C	read(*,*)lmax
C	print *," "
C
C	write(*,*)'Enter the minimum spherical harmonic degree:'
C	read(*,*)lmin
      lmin = 6
C	print *," "
C
      write(*,*)'Enter the output file name:'
      read(*,*)outname
      print *," "
      open(unit=11,file=outname,status='unknown',blank='null')
C	write(*,*)'Enter the width of cosine taper:'
C	read(*,*)dl
      dl = 0
C	print *," "
C
C  Compute the normalization factors for the spherical harmonics:
      do l = 1,lmax
        do m = 0,l
          if (l.ge.lmin) then
            N(l*(l+1)/2+m) = norm(l,m)
          else
            N(l*(l+1)/2+m) = 0.d0
          endif
C	    if (l.ge.(lmax-dl)) then
C	      j1 = (dcos(pi*dble(l-lmax+dl)/dble(dl)) + 1.d0)/2.d0
C	      N(l*(l+1)/2+m) = N(l*(l+1)/2+m)*j1
C	    endif
        enddo
      enddo
C
C  The guts of this section are "barrowed" from Banerdt's code.
      do i = 1,361
        lat = (-45.d0 + dble(i-1)*0.25d0) + latcr
        th = (90.d0 - lat)*pi/180.d0
        cth = dcos(th)
        sth = dsin(th)
        R = 1590.d3/dcos((lat - latcr)*pi/180.d0)
        if (lmax.ge.2) then
C********* calculate degree 2 polynomials
          P(3)=(3.d0*cth*cth-1.d0)/2.d0
          P(4)=3.d0*sth*cth
          P(5)=3.d0*sth*sth
        endif
        if (lmax.ge.3) then
C********* calculate degree 3 polynomials
          P(6)=(5.d0*cth*P(3)-2.d0*cth)/3.d0
          P(7)=(5.d0*cth*P(4)-3.d0*sth)/2.d0
          P(8)=5.d0*sth*P(4)
          P(9)=5.d0*sth*P(5)
        endif
C********* calculate degree l>3 polynomials using recursion
        do l = 4,lmax
C********* calculate order 0 and 1 polynomials for each degree
          P(l*(l+1)/2+0)=(dble(2*l-1)*cth*P((l-1)*l/2+0)
     *                     -dble(l-1)*P((l-2)*(l-1)/2+0))/dble(l)
          P(l*(l+1)/2+1)=(dble(2*l-1)*cth*P((l-1)*l/2+1)
     *                     -dble(l)*P((l-2)*(l-1)/2+1))/dble(l-1)
C********* calculate order (l-1) and l polynomials for each degree
          P(l*(l+1)/2+(l-1))=dble(2*l-1)*sth*P((l-1)*l/2+(l-2))
          P(l*(l+1)/2+l)=dble(2*l-1)*sth*P((l-1)*l/2+(l-1))
C********* calculate polynomials of 1<order<(l-1) for each degree
          do m = 2,l-2
            P(l*(l+1)/2+m)=(dble(2*l-1)*cth*P((l-1)*l/2+m)
     *                     -dble(l+m-1)*P((l-2)*(l-1)/2+m))/dble(l-m)
          enddo
        enddo
        fi = 0.d0
        j2 = ((a/R)**dble(1+2))
        j1 = (N(1)*cth*C1lm(1) + 
     *        N(2)*sth*(C1lm(2)*dcos(fi) + C2lm(2)*dsin(fi)))*j2
        if (lmax.ge.2) then
          j2 = ((a/R)**dble(2+2))
          j1 = j1 + (N(3)*P(3)*C1lm(3) + N(4)*P(4)*(C1lm(4)*dcos(fi) +
     *          C2lm(4)*dsin(fi)) + N(5)*P(5)*(C1lm(5)*dcos(2.d0*fi) +
     *          C2lm(5)*dsin(2.d0*fi)))*j2
        endif
        if (lmax.ge.3) then
          j2 = ((a/R)**dble(3+2))
          j1 = j1 + (N(6)*P(6)*C1lm(6) + N(7)*P(7)*(C1lm(7)* 
     *          dcos(fi) + C2lm(7)*dsin(fi)) + N(8)*P(8)*(C1lm(8)*
     *          dcos(2.d0*fi) + C2lm(8)*dsin(2.d0*fi)) + N(9)*P(9)*
     *          (C1lm(9)*dcos(3.d0*fi) + C2lm(9)*dsin(3.d0*fi)))*j2
        endif
        do l = 4,lmax
          j2 = ((a/R)**dble(l+2))
          j1 = j1 + (N(l*(l+1)/2+0)*P(l*(l+1)/2+0)*C1lm(l*(l+1)/2+0) + 
     *          N(l*(l+1)/2+1)*P(l*(l+1)/2+1)*(C1lm(l*(l+1)/2+1)*
     *          dcos(fi) + C2lm(l*(l+1)/2+1)*dsin(fi)) +
     *          N(l*(l+1)/2+(l-1))*P(l*(l+1)/2+(l-1))*(
     *          C1lm(l*(l+1)/2+(l-1))*dcos(dble(l-1)*fi) + 
     *          C2lm(l*(l+1)/2+(l-1))*dsin(dble(l-1)*fi)) + 
     *          N(l*(l+1)/2+l)*P(l*(l+1)/2+l)*(C1lm(l*(l+1)/2+l)*
     *          dcos(dble(l)*fi) + C2lm(l*(l+1)/2+l)*
     *          dsin(dble(l)*fi)))*j2
          do m = 2,l-2
            j1 = j1 + (N(l*(l+1)/2+m)*P(l*(l+1)/2+m)*(
     *            C1lm(l*(l+1)/2+m)*dcos(dble(m)*fi) + 
     *            C2lm(l*(l+1)/2+m)*dsin(dble(m)*fi)))*j2
          enddo
        enddo
C	  if (lat.gt.90.d0) lat = 180.d0 - lat
        write(11,666)lat,(R-1560.d3)/1.d3,j1
      enddo
C
      close(unit=11)
666   format(2x,f7.3,2x,f10.3,2x,E18.12)
      end
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C  The guts of this section are "barrowed" from Banerdt's code.
C
      function norm(l,m)
C
      implicit none
      real*8 norm,nl0,fl
      integer*4 l,m
C
      nl0 = dsqrt(dble(2*l + 1))
      if (m.eq.0) then
        norm = nl0
      else
        if (m.eq.l) then
          norm = dsqrt(2.d0)*nl0*dexp(-fl(2*l)/2.d0)
        else
          norm = dsqrt(2.d0)*nl0*dexp((fl(l - m) - fl(l + m))/2.d0)
        endif
      endif
C
      return
      end    
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C  The guts of this section are "barrowed" from Banerdt's code.
C
      function fl(arg)
C
      implicit none
      real*8 fl
      integer*4 i,arg
C
      fl = 0.d0
      do i = 1,arg
        fl = fl + dlog(dble(i))
      enddo	    
C
      return
      end    
C
C23456789012345678901234567890123456789012345678901234567890123456789012
