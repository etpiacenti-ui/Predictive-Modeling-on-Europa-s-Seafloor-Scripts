C23456789012345678901234567890123456789012345678901234567890123456789012
C
	program ZonalCoeff
C
C  This program makes zonal spherical harmonics coefficients for an
C    input axisymmetric profile.
C
	implicit none
	real*8 colat(721),value(721),P(7380),N(7380)
	real*8 pi,th,dth,norm
	real*8 r,cth,sth
	integer*4 l,ll,lmin,lmax,m,i,j
	character*20 inname,outname
C
	pi = dacos(-1.d0)
C
	write(*,*)'Enter the input file name:'
	read(*,*)inname
	print *," "
	open(unit=10,file=inname,status='unknown',blank='null')
C
C  Read in the profile values from the file:
C    Assumes quarter degree resolution from pole to pole.
	do i = 1,721
	  read(10,*)colat(i),value(i)
	enddo
	close(10)
C
C	write(*,*)'Enter the maximum spherical harmonic degree:'
C	read(*,*)lmax
C	print *," "
	lmax = 120
C
C	write(*,*)'Enter the minimum spherical harmonic degree:'
C	read(*,*)lmin
C	print *," "
	lmin = 1
C
	write(*,*)'Enter the output file name:'
	read(*,*)outname
	print *," "
	open(unit=11,file=outname,status='unknown',blank='null')
C
	dth = 0.25d0
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
	do ll = lmin,lmax
	  r = 0.d0
	  i = ll*(ll+1)/2
	  do j = 1,721
	    th = dble(j-1)*dth*pi/180.d0
	    cth = dcos(th)
	    sth = dsin(th)
	    P(1) = cth
	    P(2) = sth
	    if (lmax.ge.2) then
********** calculate degree 2 polynomials
	      P(3)=(3.d0*cth*cth-1.d0)/2.d0
	      P(4)=3.d0*sth*cth
	      P(5)=3.d0*sth*sth
	    endif
	    if (lmax.ge.3) then
********** calculate degree 3 polynomials
	      P(6)=(5.d0*cth*P(3)-2.d0*cth)/3.d0
	      P(7)=(5.d0*cth*P(4)-3.d0*sth)/2.d0
	      P(8)=5.d0*sth*P(4)
	      P(9)=5.d0*sth*P(5)
	    endif
********** calculate degree l>3 polynomials using recursion
	    do l = 4,lmax
********** calculate order 0 and 1 polynomials for each degree
	      P(l*(l+1)/2+0)=(dble(2*l-1)*cth*P((l-1)*l/2+0)
     *                     -dble(l-1)*P((l-2)*(l-1)/2+0))/dble(l)
	      P(l*(l+1)/2+1)=(dble(2*l-1)*cth*P((l-1)*l/2+1)
     *                     -dble(l)*P((l-2)*(l-1)/2+1))/dble(l-1)
********** calculate order (l-1) and l polynomials for each degree
	      P(l*(l+1)/2+(l-1))=dble(2*l-1)*sth*P((l-1)*l/2+(l-2))
	      P(l*(l+1)/2+l)=dble(2*l-1)*sth*P((l-1)*l/2+(l-1))
********** calculate polynomials of 1<order<(l-1) for each degree
	      do m = 2,l-2
	        P(l*(l+1)/2+m)=(dble(2*l-1)*cth*P((l-1)*l/2+m)
     *                     -dble(l+m-1)*P((l-2)*(l-1)/2+m))/dble(l-m)
	      enddo
	    enddo
	    r = r + N(i)*P(i)*value(j)*dsin(th)*dth*pi/180.d0
	  enddo
	  r = r/2.d0
	  write(11,666)ll,0,r,0.0
	  do m = 1,ll
	    write(11,666)ll,m,0.0,0.0
	  enddo
	enddo
C
	close(unit=11)
666	format(2x,i3,2x,i3,2x,E18.12,2x,E18.12)
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
