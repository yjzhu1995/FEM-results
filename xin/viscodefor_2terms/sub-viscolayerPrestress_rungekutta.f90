module para
	implicit none
	integer :: np   ! number of postseismic snapshots
	real*8,parameter :: Re=6371.d3
	real*8 :: pi, s2y, g0
	real*8,parameter :: UGC=6.673d-11   ! universal gravitation constant
	real*8,parameter :: km2m=1.d3, Gpa2pa=1.d10
	real*8,parameter :: step=5.d3
	integer,parameter :: nstep = 5

contains
	    subroutine lufact( A )
        !... P * transpose(A) = L * U
        !... http://www.netlib.org/lapack/explore-3.1.1-html/dgetrf.f.html
        !... (note that the definition of P is opposite to that of the above page)

        complex*16, intent(inout) :: A(:,:)
        complex*16, allocatable, dimension(:,:) :: L, U, P
        integer, allocatable  :: ipiv(:)
        complex*16, allocatable :: row(:)
        integer :: i, n, m, info
        external ZGETRF

        n = size( A, 1 )
        m = size( A, 2 )
        allocate( L( n, m ), U( m, m ), P( n, n ), ipiv( m ), row( n ) )

        L = A
        call ZGETRF( n, m, L, n, ipiv, info )
        if ( info /= 0 ) stop "lufact: info /= 0"

        U = 0.0d0
        P = 0.0d0
        do i = 1, n
            U( i, i:m ) = L( i, i:m )
            L( i, i:m ) = 0.0d0
            if(i<=m) L( i, i ) = 1.0d0
            P( i, i ) = 1.0d0
        enddo

        !... Assuming that P = P[ipiv(n),n] * ... * P[ipiv(1),1]
        !... where P[i,j] is a permutation matrix for i- and j-th rows.
        do i = 1, m
            row = P( i, : )
            P( i, : ) = P( ipiv(i), : )
            P( ipiv(i), : ) = row
        enddo

        A = matmul(transpose(P),L)

    end subroutine lufact

end module para



	subroutine Spheroidalsol(lambda,mu,rho,N,r1,Yr)
!
!      Spheroidal fundamental matrix, 4x4 matrix consisting of ODE's solutions is obtained analytically.
!      Matrix contains singularity solutions.
!      In the matrix Yr, solutions in the first two column are sigularity at the center.
!      (y1,y3,y2,y4,y5,y6)^T
!      y1 -- vertical displacement
!      y3 -- horizontal displacement
!      y2 -- vertical component of vertical stress
!      y4 -- horizontal component of vertical stress
!      y5 -- gravity potential
!      y6 -- gravity flux density
!
		use para
		implicit none
		complex*16 :: lambda,mu
		real*8 :: r1,rho
		integer :: N
		real*8 :: NN
		integer,parameter :: nl=6,nc=3
		complex*16 :: Yr(nl,nc)
		complex*16 :: X,Z,beta,xi,vp2
		complex*16 :: f13,f23,f33,f43
		complex*16 :: M1, M2
		COMPLEX*16 :: CONHYP

		NN = dble(N)
		beta = lambda+2.d0*mu
		vp2 = beta/rho
		xi = 4.d0*pi*UGC*rho/3.d0
		X = r1*r1*xi/vp2/2.d0


		M1 = CONHYP(dcmplx(NN/2.d0),dcmplx(NN+0.5d0),X,0,10)
		M2 = CONHYP(dcmplx(NN/2.d0),dcmplx(NN+1.5d0),X,0,10)

		f13 = -NN+(2.d0*NN+1.d0)*M1-(NN+1.d0)*M2
		f23 = -1.d0+M2
		f33 = 2.d0*NN*(NN-1)*mu*vp2-(2.d0*NN+1.d0)*M1*(beta*xi*r1*r1-4.d0*mu*vp2) &
		     +(NN+1.d0)*M2*(beta*xi*r1*r1-2.d0*(NN+2.d0)*mu*vp2)
		f43 = -2.d0*(NN-1.d0)+2.d0*(2*NN+1.d0)*M1-2.d0*(NN+2.d0)*M2

		Yr = dcmplx(0.d0,0.d0)

		Yr(2,1) = dcmplx(-1.d0/NN/(NN+1.d0))
		Yr(3,1) = lambda/r1
		Yr(4,1) = -mu/(NN+1.d0)/r1
		Yr(5,1) = 3.d0*xi/2.d0/(2.d0*NN+3.d0)*r1
		Yr(6,1) = 3.d0*xi/2.d0

		Yr(1,2) = vp2*f13
		Yr(2,2) = vp2*f23
		Yr(3,2) = -f33/r1
		Yr(4,2) = vp2*mu*f43/r1
		Yr(5,2) = 3.d0*xi*vp2*f23*r1
		Yr(6,2) = 3.d0*xi*(NN+1.d0)*vp2*f23

		Yr(5,3) = dcmplx(1.d0)*r1
		Yr(6,3) = dcmplx(2.d0*NN+1.d0)
	end subroutine Spheroidalsol

! --------------------------------------------------------

	subroutine fundamatrixA(lambda,mu,Cg,xi1,xi2,N,rr,nm,Yr)
!
!   Spheroidal fundamental matrix, a 6x6 matrix consisting of ODE's equation.
!   dy/dr = Yr(lambda,mu,rr) y
!   (A1,A3,A2,A4,A5,A6)^T
!   xi2 = 4.d0*pi*UGC*rho*Re/g0
!
	implicit none
	complex*16 :: lambda,mu,Cg,xi1,xi2
	integer :: i,nm,N
	real*8 :: rr,NN
	complex*16 :: Yr(nm,nm)
	complex*16 :: beta,kappa

	NN = dble(N)

	beta = lambda+2.d0*mu
	kappa = lambda+2.d0*mu/3.d0

	Yr = dcmplx(0.d0,0.d0)
	
	Yr(1,1) = -2.d0*lambda/beta/rr
	Yr(1,2) =  lambda/beta/rr*NN*(NN+1.d0)
	Yr(1,3) =  1.d0/beta

	Yr(2,1) =  dcmplx(-1.d0/rr)
	Yr(2,2) =  dcmplx(1.d0/rr)
	Yr(2,4) =  1.d0/mu

	Yr(3,1) =  12.d0*kappa*mu/beta/rr/rr-2.d0*Cg*(1.d0+lambda/beta)/rr+xi1
	Yr(3,2) =  NN*(NN+1)*(Cg*lambda-6.d0*kappa*mu/rr)/beta/rr
	Yr(3,3) = -4.d0*mu/beta/rr+Cg/beta
	Yr(3,4) =  dcmplx(NN*(NN+1.d0)/rr)

	Yr(4,1) =  (Cg-6.d0*mu*kappa/beta/rr)/rr
	Yr(4,2) =  2.d0*mu/beta*((2.d0*NN*(NN+1.d0)-1.d0)*lambda+2.d0*(NN*(NN+1.d0)-1.d0)*mu)/rr/rr
	Yr(4,3) = -lambda/beta/rr
	Yr(4,4) = dcmplx(-3.d0/rr)

	Yr(5,1) = xi2
	Yr(5,5) = -(NN+1.d0)/rr
	Yr(5,6) = dcmplx(1.d0)

	Yr(6,1) = xi2*(NN+1.d0)/rr
	Yr(6,2) = -NN*Yr(6,1)
	Yr(6,6) = (NN-1.d0)/rr
	end subroutine fundamatrixA

! --------------------------------------------------------


! --------------------------------------------------------

	subroutine fundamatrixT(mu,NN,rr,rn,nm,Yr)
!
!      Toroidal fundamental matrix, 2x2 matrix consisting of ODE's solutions is obtained analytically.
!      Matrix contains singularity solutions.
!      In the matrix Yr, solutions in the last column are sigularity at the center.
!
	implicit none
	complex*16 :: mu
	integer :: nm,NN
	real*8 :: rr,rn
	complex*16 :: Yr(nm,nm)
	real*8 :: rd          ! rd=rr/Re guarantees r^n can be calculated

	rd = rr/rn

	Yr(1,1) = dcmplx(rd**NN)
	Yr(1,2) = dcmplx(rd**(-NN-1))

	Yr(2,1) =  mu*dble(NN-1)*rd**NN/rr
	Yr(2,2) = -mu*dble(NN+2)*rd**(-NN-1)/rr

	end subroutine fundamatrixT

! --------------------------------------------------------

	subroutine invfundamatrixT(mu,N,rr,rn,nm,Yr)
!
!      Invserse of Spheroidal fundamental matrix, 4x4 matrix consisting of ODE's solutions is obtained analytically.
!      Matrix contains singularity solutions.
!
	implicit none
	complex*16 :: mu
	integer :: i,nm,N
	real*8 :: rr,rn
	complex*16 :: Yr(nm,nm)
	real*8 :: rd          ! rd=rr/Re guarantees r^n can be calculated
	real*8 :: fl1

	rd = rr/rn

	fl1 = dble(2*N+1)
	Yr(1,1) = dcmplx(rd**(-N)*dble(N+2)/fl1)
	Yr(1,2) = rd**(-N)/fl1/mu*rr
	Yr(2,1) = dcmplx(rd**(N+1)*dble(N-1)/fl1)
	Yr(2,2) = -rd**(N+1)/fl1/mu*rr

	end subroutine invfundamatrixT
! --------------------------------------------------------

	SUBROUTINE sourceS(lambda,mu,N,RAD,S)
!
!   spheroidal source function
!   S(1)=y1(rs+)-y1(rs-)      ! vertical displacement
!   S(2)=y3(rs+)-y3(rs-)      ! horizontal displacement
!   S(3)=y2(rs+)-y2(rs-)      ! vertical stress
!   S(4)=y4(rs+)-y4(rs-)      ! horizontal stress
!
	use para
	implicit none
	complex*16 :: lambda, mu
	integer :: N
	real*8 :: RAD
	complex*16 :: S(6,4)
	
	S = dcmplx(0.d0,0.d0)

! ************************** Strike-slip **********************************

	S(4,1) = -(2.D0*N+1.D0)/8.D0/PI/N/(N+1.D0)/RAD**3*mu

! ***************************** Dip-slip **********************************

	S(2,2) = dcmplx((2.D0*N+1.D0)/8.D0/PI/N/(N+1.D0)/RAD**2)

! **************************** Horizontal tensile**************************

    S(1,3) = (2.D0*N+1.D0)/4.D0/PI/RAD/RAD*lambda/(lambda+2.D0*mu)
	S(3,3) = -(2.D0*N+1.D0)/2.D0/PI/RAD**3 &
     	*mu*(3.D0*lambda+2.D0*mu)/(lambda+2.D0*mu)
    S(4,3) = (2.D0*N+1.D0)/4.D0/PI/RAD**3 &
     	*mu*(3.D0*lambda+2.D0*mu)/(lambda+2.D0*mu)

! **************************** Vertical tensile***************************

    S(1,4) = dcmplx((2.D0*N+1.D0)/4.D0/PI/RAD/RAD)

! *************************************************************************
	end SUBROUTINE sourceS

! --------------------------------------------------------

	SUBROUTINE sourceT(mu,N,dept,ST)
!
!   toroidal source function
!
		use para
		implicit none
		complex*16 :: mu
		integer :: N
		real*8 :: dept,RAD
		complex*16 :: ST(2,2)

	RAD = Re-dept
    ST = dcmplx(0.d0,0.d0)
! ************************** Strike-slip **********************************

	ST(2,1) = -(2.D0*N+1.D0)*mu/8.D0/PI/N/(N+1.D0)/RAD**3

! ***************************** Dip-slip **********************************

	ST(1,2) = dcmplx((2.D0*N+1.D0)/8.D0/PI/N/(N+1.D0)/RAD**2)

! *************************************************************************
    end subroutine sourceT

! --------------------------------------------------------

	subroutine Burgers(s,nlayer,lambda,mu,tau,lambdas,mus)
		integer :: nlayer, k
		real*8 :: lambda(nlayer), mu(nlayer,2), tau(nlayer,2), bulk(nlayer)
		complex*16 :: lambdas(nlayer), mus(nlayer), s

		do k=1,nlayer
			bulk(k) = lambda(k) + 2.d0*mu(k,2)/3.d0
		if(abs(tau(k,1))<1.d-8 .and. dabs(tau(k,2))<1.d-8) then
			! elastic
			mus(k) = dcmplx(mu(k,2),0.d0)
		elseif(abs(mu(k,1))<1.d-8 .or. abs(tau(k,1))<1.d-8) then
			! Maxwell body
			mus(k) = s*mu(k,2)/(s+1.d0/tau(k,2))
		else
			! Burgers body
			mus(k) = s*mu(k,2)/(s+1.d0/tau(k,2)+s/(mu(k,1)/mu(k,2)+s*tau(k,1)))
		endif
		enddo
		lambdas = bulk - 2.d0*mus/3.d0
	end subroutine Burgers

! --------------------------------------------------------

	subroutine LOVESLAP(nlayer,rlayer,lambda,mu,rho,g,tau,t,NN,dept,HN,LN,KN)
! 
!   THIS SUBROUTINE TO CALCULATE spheroidal Love number -- HN,LN,KN
!   by inverse Laplace transform
! 
	use para
	integer :: nlayer
	real*8 :: rlayer(nlayer,2),lambda(nlayer), mu(nlayer,2), rho(nlayer), g(nlayer), tau(nlayer,2), t
	complex*16 :: lambdas(nlayer), mus(nlayer)
	INTEGER :: k,NN,M
	integer, PARAMETER :: ns=4
	REAL*8 :: dept,HN(ns),LN(ns),KN(ns)
	real*8 :: theta, sigma, ctan
	real*8 :: part1(3,ns),part2(3,ns)
	complex*16 :: HNs(ns),LNs(ns),KNs(ns),r,st

	M = precision(1.d0)
	r = 2.d0*M/5.d0/t
	call Burgers(r,nlayer,lambda,mu,tau,lambdas,mus)
	call LOVES(nlayer,rlayer,lambdas,mus,rho,g,NN,dept,HNs,LNs,KNs)
	part1(1,1:ns) = 0.5d0*real(HNs(1:ns)/r*zexp(r*t))
	part1(2,1:ns) = 0.5d0*real(LNs(1:ns)/r*zexp(r*t))
	part1(3,1:ns) = 0.5d0*real(KNs(1:ns)/r*zexp(r*t))

	part2 = 0.d0
	do k=1,M-1
		theta = k*pi/M
		ctan = 1.d0/dtan(theta)
		st = r*theta*(ctan+(0.d0,1.d0))
		sigma = theta+(theta*ctan-1.d0)*ctan
		call Burgers(st,nlayer,lambda,mu,tau,lambdas,mus)
		call LOVES(nlayer,rlayer,lambdas,mus,rho,g,NN,dept,HNs,LNs,KNs)
		part2(1,1:ns) = part2(1,1:ns)+real(zexp(t*st)*HNs(1:ns)/st*dcmplx(1.d0,sigma))
		part2(2,1:ns) = part2(2,1:ns)+real(zexp(t*st)*LNs(1:ns)/st*dcmplx(1.d0,sigma))
		part2(3,1:ns) = part2(3,1:ns)+real(zexp(t*st)*KNs(1:ns)/st*dcmplx(1.d0,sigma))
	end do

	HN(1:ns) = real(r)/M*(part1(1,1:ns)+part2(1,1:ns))
	LN(1:ns) = real(r)/M*(part1(2,1:ns)+part2(2,1:ns))
	KN(1:ns) = real(r)/M*(part1(3,1:ns)+part2(3,1:ns))

	end subroutine LOVESLAP

! --------------------------------------------------------

    SUBROUTINE LOVES(nlayer,rlayer,lambda,mu,rho,g,NN,dept,HN,LN,KN)
! 
!   THIS SUBROUTINE TO CALCULATE spheroidal Love number -- HN,LN,KN
!   in the Laplace domain
! 
	use para
	IMPLICIT NONE
	integer :: nlayer
	real*8 :: rlayer(nlayer,2), rho(nlayer), g(nlayer)
	complex*16 :: lambda(nlayer), mu(nlayer)
	INTEGER :: NN
	integer, PARAMETER :: nm=4,nc=3,ns=4
	REAL*8 :: dept
	complex*16 :: HN(ns),LN(ns),KN(ns)
	complex*16 :: d(nc,ns)

	CALL COEFFS(nlayer,rlayer,lambda,mu,rho,g,NN,dept,d)

	HN(1:ns) = d(1,1:ns)
	LN(1:ns) = d(2,1:ns)
	KN(1:ns) = d(3,1:ns)*g0
! 	HN(1:ns) = d(1,1:ns)*Re*Re**2
! 	LN(1:ns) = d(2,1:ns)*Re*Re**2
! 	KN(1:ns) = d(3,1:ns)*Re*g0*Re**2

!  	write(*,*) NN
!  	write(*,*) HN(1)
!  	write(*,*) LN(1)
!  	write(*,*) KN(1)
!  	read(*,*)

	end SUBROUTINE LOVES

! --------------------------------------------------------

! --------------------------------------------------------

	subroutine LOVETLAP(nlayer,rlayer,lambda,mu,tau,t,NN,dept,LNT)
! 
!   THIS SUBROUTINE TO CALCULATE spheroidal Love number -- HN,LN
!   by inverse Laplace transform
! 
	use para
	integer :: nlayer
	real*8 :: rlayer(nlayer,2),lambda(nlayer), mu(nlayer,2), tau(nlayer,2), t
	complex*16 :: lambdas(nlayer), mus(nlayer)
	INTEGER :: k,NN,M
	integer, PARAMETER :: ns=2
	REAL*8 :: dept,LNT(ns)
	real*8 :: theta, sigma, ctan
	real*8 :: part1(ns),part2(ns)
	complex*16 :: LNTs(ns),r,st

	M = precision(1.d0)
	r = 2.d0*M/5.d0/t
	call Burgers(r,nlayer,lambda,mu,tau,lambdas,mus)
	call LOVET(nlayer,rlayer,mus,NN,dept,LNTs)
	part1(1:ns) = 0.5d0*real(LNTs(1:ns)/r*zexp(r*t))

	part2 = 0.d0
	do k=1,M-1
		theta = k*pi/M
		ctan = 1.d0/dtan(theta)
		st = r*theta*(ctan+(0.d0,1.d0))
		sigma = theta+(theta*ctan-1.d0)*ctan
		call Burgers(st,nlayer,lambda,mu,tau,lambdas,mus)
		call LOVET(nlayer,rlayer,mus,NN,dept,LNTs)
		part2(1:ns) = part2(1:ns)+real(zexp(t*st)*LNTs(1:ns)/st*dcmplx(1.d0,sigma))
	end do
	
	LNT(1:ns) = real(r)/M*(part1(1:ns)+part2(1:ns))

	end subroutine LOVETLAP

! --------------------------------------------------------

    SUBROUTINE LOVET(nlayer,rlayer,mu,NN,dept,LNT)
! 
!   THIS SUBROUTINE TO CALCULATE Toroidal Love number -- HN,LN
! 
	use para
    implicit none
	integer :: nlayer
	real*8 :: rlayer(nlayer,2)
	complex*16 :: mu(nlayer)
	INTEGER :: i,j,NN
	integer,parameter :: nm=2,N=3,nc=1,ns=2
	real*8 :: dept
	complex*16 :: LNT(ns),d(nc,ns),C(nm,ns)
	complex*16 :: X(N,ns)
	complex*16 :: YT(nm,nm),Y(nc,nm)

	!  at the Earth's surface
	call fundamatrixT(mu(1),NN,Re,Re,nm,YT)
	Y(1,1:nm) = YT(1,1:nm)

	call COEFFT(nlayer,rlayer,mu,NN,dept,X)

	C(1:nm,1:ns) = X(1:nm,1:ns)

	d = matmul(Y,C)
	LNT(1:ns) = d(1,1:ns)*Re**2

    end subroutine LOVET

! --------------------------------------------------------

	subroutine propagation(N,nl,nc,y0,lambda,mu,Cg,xi1,xi2,r0,r1,y1)
!  propagate the 6x3 solution from r0 to r1 by using 4th-order Runge-Kutta (RK4) method
	implicit none
	integer :: N
	integer :: nl, nc
	complex*16 :: lambda,mu,Cg,xi1,xi2
	real*8 :: r0,r1
	complex*16 :: y0(nl,nc), y1(nl,nc)
	real*8 :: h    ! step
	complex*16,dimension(nl,nc) :: k1,k2,k3,k4
	complex*16 :: f(nl,nl)
	real*8 :: rt

	h = r1-r0
	
	call fundamatrixA(lambda,mu,Cg,xi1,xi2,N,r0,nl,f)
	k1 = matmul(f,y0)
	k1 = h*k1

	rt = r0+h/2.d0
	y1 = y0+k1/2.d0
	call fundamatrixA(lambda,mu,Cg,xi1,xi2,N,rt,nl,f)
	k2 = matmul(f,y1)
	k2 = h*k2

	y1 = y0+k2/2.d0
	k3 = matmul(f,y1)
	k3 = h*k3

	y1 = y0+k3
	call fundamatrixA(lambda,mu,Cg,xi1,xi2,N,r1,nl,f)
	k4 = matmul(f,y1)
	k4 = h*k4

	y1 = y0 + (k1+2.d0*k2+2.d0*k3+k4)/6.d0

	end subroutine propagation

! --------------------------------------------------------

! --------------------------------------------------------

	subroutine layermatrixT(mu,r1,r2,NN,nm,matr)
! 	layer matrix used to matrix propagation for spheroid
! 	r1 -- upper bound
!   r2 -- lower bound
		implicit none
		complex*16 :: mu
		real*8 :: r1,r2,rd
		integer :: NN,nm
		complex*16 :: matr(nm,nm),a(nm,nm),b(nm,nm),c(nm,nm)

		call fundamatrixT(mu,NN,r1,r1,nm,a)
		call invfundamatrixT(mu,NN,r2,r2,nm,b)
		rd = r1/r2
		c = dcmplx(0.d0,0.d0)
		c(1,1) = dcmplx(rd**NN)
		c(2,2) = dcmplx(rd**(-NN-1))

		matr = matmul(a,c)
		matr = matmul(matr,b)
	end subroutine layermatrixT

! --------------------------------------------------------

	subroutine COEFFS(nlayer,rlayer0,lambdas,mus,rho,g,NN,dept,d)
!
!   solve the boundary value equations
!
		use para
		implicit none
		integer,parameter :: N=3, nm=6, ns=4
		integer :: NN,i,j,k,nds
		integer :: nlayer
		real*8 :: rlayer0(nlayer,2),rho(nlayer),g(nlayer)
		real*8 :: rlayer(nlayer,2)
		complex*16 :: lambdas(nlayer),mus(nlayer)
		complex*16 :: lambda(nlayer),mu(nlayer)
		complex*16 :: Cg(nlayer),xi1(nlayer),xi2(nlayer)
		real*8 :: dept,rdept
		complex*16 :: A(N,N),B(nm,ns),BB(N,ns),X(N,ns),d(N,ns)
		complex*16 :: Yrs(nm,N),Yrs1(nm,N),Yrs2(nm,ns),Yr1(nm,ns)
		REAL*8 :: r1,r2,h
		integer :: IPIV(N), info
		external zgesv


	lambda = lambdas/lambdas(nlayer)
	mu = mus/lambdas(nlayer)
	rlayer = 1.d0-rlayer0/Re
	rdept = 1.d0-dept/Re
	do k=1,nlayer
		Cg(k) = rho(k)*g(k)*Re/lambdas(nlayer)
		xi1(k) = 4.d0*pi*UGC*rho(k)*rho(k)*Re*Re/lambdas(nlayer)
		xi2(k) = 4.d0*pi*UGC*rho(k)*Re/g0
	end do


	! locate source layer
	do k=1,nlayer
		if(rlayer0(k,2)>dept) then
			nds = k
			exit
		endif
	end do

	! initial solution at the Earth's center
	r1 = Re-rlayer0(nlayer,1)
	call Spheroidalsol(lambdas(nlayer),mus(nlayer),rho(nlayer),NN,r1,Yrs)
	Yrs(1,:) = Yrs(1,:)/Re
	Yrs(2,:) = Yrs(2,:)/Re
	Yrs(3,:) = Yrs(3,:)/lambdas(nlayer)
	Yrs(4,:) = Yrs(4,:)/lambdas(nlayer)
	Yrs(5,:) = Yrs(5,:)/Re/g0
	Yrs(6,:) = Yrs(6,:)/g0

	! from Earth's center to surface
	do k=nlayer-1,1,-1
		r2 = rlayer(k,2)
		h = (rlayer(k,1)-rlayer(k,2))/nstep
		h = min(step/Re, h)
		do while(r2 < rlayer(k,1))
			r1 = min(r2+h,rlayer(k,1))
			call propagation(NN,nm,N,Yrs,lambda(k),mu(k),Cg(k),xi1(k),xi2(k),r2,r1,Yrs1)
			call lufact(Yrs1)
			Yrs = Yrs1
			r2 = r1
		enddo
	end do

	! from source's depth to surface
	call sourceS(lambda(nds),mu(nds),NN,rdept,Yrs2)
! 	Yrs2 = Yrs2/Re**3     ! dimensionless

	k = nds
	r2 = rdept
	h = (rlayer(k,1)-r2)/nstep
	h = min(step/Re, h)
	do while(r2 < rlayer(k,1))
		r1 = min(r2+h,rlayer(k,1))
		call propagation(NN,nm,nm,Yrs2,lambda(k),mu(k),Cg(k),xi1(k),xi2(k),r2,r1,Yr1)
! 		call lufact(Yr1)
! 		Yrs2 = Yr1/maxval(zabs(Yr1))
		Yrs2 = Yr1
		r2 = r1
	enddo

	if(nds/=1) then
	do k=nds-1,1,-1
		r2 = rlayer(k,2)
		h = (rlayer(k,1)-rlayer(k,2))/nstep
		h = min(step/Re, h)
		do while(r2 < rlayer(k,1))
			r1 = min(r2+h,rlayer(k,1))
			call propagation(NN,nm,nm,Yrs2,lambda(k),mu(k),Cg(k),xi1(k),xi2(k),r2,r1,Yr1)
! 			call lufact(Yr1)
! 			Yrs2 = Yr1/maxval(zabs(Yr1))
			Yrs2 = Yr1
			r2 = r1
		enddo
	end do
	endif

	B = Yrs2

	A(1:2,1:N) = Yrs(3:4,1:N)
	A(N,1:N) = Yrs(6,1:N)
	BB(1:2,1:ns) = -B(3:4,1:ns)
	BB(3,1:ns) = -B(6,1:ns)

	call zgesv(N,ns,A,N,IPIV,BB,N,info)
    if (info /= 0) then
    	write(*,*) NN
    	write(*,*) A
    	stop "zegsv: info /= 0"
    endif
	X = BB

	BB(1:2,1:ns) = B(1:2,1:ns)
	BB(3,1:ns) = B(5,1:ns)
	A(1:2,1:N) = Yrs(1:2,1:N)
	A(3,1:N) = Yrs(5,1:N)
	d = matmul(A,X)
	d = d+BB
	end subroutine COEFFS

! --------------------------------------------------------

! --------------------------------------------------------

	subroutine COEFFT(nlayer,rlayer,mus,NN,dept,X)
		use para
		implicit none
		INTEGER i,j,k,NN,nds
		integer,parameter :: nm=2,N=3,ns=2
		integer :: nlayer
		real*8 :: rlayer(nlayer,2)
		complex*16 :: mus(nlayer)
		complex*16 :: mu(nlayer)
		real*8 :: dept
		complex*16 :: YT(nm,nm),YT1(nm,nm),YTs1(nm,nm),YTs2(nm,nm)
		complex*16 :: A(N,N),B(nm,ns),BB(N,ns),X(N,ns)    ! work array for linear system
		real*8 :: r1,r2,r_upper
		integer :: IPIV(N), info
		external zgesv


	mu = mus/Gpa2pa

	!  at the Earth's surface
	r1 = Re-rlayer(1,1)
	call fundamatrixT(mu(1),NN,r1,Re,nm,YT)

	!  depth of source
	!  downward propagate solutions from Earth's surface to the source's depth
	YTs1 = dcmplx(0.d0,0.d0)
	do i=1,nm
		YTs1(i,i) = dcmplx(1.d0,0.d0)
	end do
	k = 1
	do 
		if(rlayer(k,2)>dept) then
			r2 = Re-dept
			nds = k
			call fundamatrixT(mu(k),NN,r2,Re,nm,YT1)
			YTs1 = matmul(YT1,YTs1)
			exit
		else
			r2 = Re-rlayer(k,2)
			call fundamatrixT(mu(k),NN,r2,Re,nm,YT1)
			YTs1 = matmul(YT1,YTs1)
		endif
		k = k+1
		r1 = Re-rlayer(k,1)     ! upper boundary of each layer
		call invfundamatrixT(mu(k),NN,r1,Re,nm,YT1)
		YTs1 = matmul(YT1,YTs1)
	end do

	!  upward propagate solutions from Earth's inner layer to the source's depth
	r1 = Re-rlayer(nlayer,1)
	call fundamatrixT(mu(nlayer),NN,r1,r1,nm,YTs2)
	! locate source layer
	do k=1,nlayer
		if(rlayer(k,2)>dept) then
			nds = k
			exit
		endif
	end do

	do k=nlayer-1,nds,-1
		r2 = Re-rlayer(k,2)
		if(k==nds) then
			r_upper = Re-dept
		else
			r_upper = Re-rlayer(k,1)
		endif
		do while(r2 < r_upper)
			r1 = min(r2+step,r_upper)
			call layermatrixT(mu(k),r1,r2,NN,nm,YT1)
			YT1 = YT1/maxval(dabs(real(YT1)))
			YTs2 = matmul(YT1,YTs2)
 			YTs2 = YTs2/maxval(dabs(real(YTs2)))
			r2 = r1
		enddo
	end do

	!  construct array
	do j=1,nm
		A(1,j) = YT(2,j)
	enddo
		A(1,N) = dcmplx(0.d0,0.d0)
	do i=2,N
		do j=1,nm
			A(i,j) =  YTs1(i-1,j)
		enddo
			A(i,N) = -YTs2(i-1,1)
	enddo

	call sourceT(mu(nds),NN,dept,B)
	BB(1,1:ns) = dcmplx(0.d0,0.d0)
	BB(2:N,1:ns) = B(1:nm,1:ns)

	call zgesv(N,ns,A,N,IPIV,BB,N,info)
	X = BB
	end subroutine COEFFT

! --------------------------------------------------------

	subroutine GREENsingle(nmax,theta,HN,LN,LNT,KN,ux,uy,uz,geoid,dgf,dgd)
!
!   geoid(1),dgf(1),dgd(1) -- u12, strike-slip for a vertical fault
!   geoid(2),dgf(2),dgd(2) -- u32, dip-slip for a vertical fault
!   geoid(3),dgf(3),dgd(3) -- u220, horizontal tensile for m=0
!   geoid(4),dgf(4),dgd(4) -- u33, vertical tensile
!
		use para
		implicit none
		integer :: nmax,N,NI(4)
		integer,parameter :: nm=4,ns=4,nt=2
		real*8 :: theta
		real*8 :: geoid(nm),dgf(nm),dgd(nm),ux(nm),uy(nm),uz(nm)
		real*8 :: HN(0:nmax,ns),LN(0:nmax,ns),LNT(0:nmax,nt),KN(0:nmax,ns)
		real*8 :: Pnm(0:nmax,0:2),pdnm(0:nmax,0:2)
		real*8 :: theta1
	
		THETA1=THETA*PI/180.D0

		CALL legendre(NMAX,2,THETA1,PNM,PDNM)

		geoid=0.D0
		dgf=0.D0
		dgd=0.D0
		UX=0.D0
		UY=0.D0
		uz=0.d0

		! the initial degree for the four typical sources, not yet used
		NI(1) = 2
		NI(2) = 1
		NI(3) = 0
		NI(4) = 0

		do N=2,nmax
			geoid(1)=geoid(1)+(-2.D0)*KN(N,1)*PNM(N,2)
			dgf(1)=dgf(1)+(-2.D0)*KN(N,1)*PNM(N,2)*dble(N+1)/Re
			uz(1)=uz(1)+(-2.D0)*HN(N,1)*PNM(N,2)
			UY(1)=UY(1)+(-2.D0)*(LN(N,1)*PDNM(N,2)+2.D0*LNT(N,1)*PNM(N,2)/DSIN(THETA1))
			UX(1)=UX(1)+(-2.D0)*(2.D0*LN(N,1)*PNM(N,2)/DSIN(THETA1)+LNT(N,1)*PDNM(N,2))

			geoid(2)=geoid(2)+(-2.D0)*KN(N,2)*PNM(N,1)
			dgf(2)=dgf(2)+(-2.D0)*KN(N,2)*PNM(N,1)*dble(N+1)/Re
			uz(2)=uz(2)+(-2.D0)*HN(N,2)*PNM(N,1)
			UY(2)=UY(2)+(-2.D0)*(LN(N,2)*PDNM(N,1)+LNT(N,2)*PNM(N,1)/DSIN(THETA1))   ! different from Sun's book where positive of toroidal term is used
			UX(2)=UX(2)+(-2.D0)*(LN(N,2)*PNM(N,1)/DSIN(THETA1)+LNT(N,2)*PDNM(N,1))

			geoid(3)=geoid(3)+KN(N,3)*PNM(N,0)
			dgf(3)=dgf(3)+KN(N,3)*PNM(N,0)*dble(N+1)/Re
			uz(3)=uz(3)+HN(N,3)*PNM(N,0)
			UY(3)=UY(3)+LN(N,3)*PDNM(N,0)

			geoid(4)=geoid(4)+KN(N,4)*PNM(N,0)
			dgf(4)=dgf(4)+KN(N,4)*PNM(N,0)*dble(N+1)/Re
			uz(4)=uz(4)+HN(N,4)*PNM(N,0)
			UY(4)=UY(4)+LN(N,4)*PDNM(N,0)
		end do

		if(abs(theta) < 1.d-15) then
			geoid = 0.d0
			dgf = 0.d0
			UX = 0.d0
			UY = 0.d0
			uz = 0.d0
			do N=2,nmax
				geoid(3) = geoid(3)+KN(N,3)
				geoid(4) = geoid(4)+KN(N,4)
				dgf(3) = dgf(3)+KN(N,3)*dble(N+1)/Re
				dgf(4) = dgf(4)+KN(N,4)*dble(N+1)/Re
				uz(3) = uz(3)+HN(N,3)
				uz(4) = uz(4)+HN(N,4)
			end do
		elseif(abs(theta-180) < 1.d-15) then
			geoid = 0.d0
			dgf = 0.d0
			UX = 0.d0
			UY = 0.d0
			uz = 0.d0
			do N=2,nmax
				geoid(3) = geoid(3)+(-1)**N*HN(N,3)
				geoid(4) = geoid(4)+(-1)**N*HN(N,4)
				dgf(3) = dgf(3)+(-1)**N*KN(N,3)*dble(N+1)/Re
				dgf(4) = dgf(4)+(-1)**N*KN(N,4)*dble(N+1)/Re
				uz(3) = uz(3)+(-1)**N*HN(N,3)
				uz(4) = uz(4)+(-1)**N*HN(N,4)
			end do
		endif

		geoid = geoid/g0
		dgd(1:nm) = dgf(1:nm) - 2.d0*g0*uz(1:nm)/Re

	end subroutine GREENsingle

! --------------------------------------------------------

	subroutine viscoGreen(nlayer,rlayer,lambda,mu,rho,g,tau,t,dept,NAG,angle,nmax,ux,uy,uz,geoid,dgf,dgd)
!
!   Compute the Green's functions of all the angles
!   and the dislocation LOVE numbers at the given source depth
!   input: 
!          dept -- source depth
!          NAG -- dimension of angle
!          angle -- angles of Green's function
!          nmax -- maximum harmonic degrees
!   output:
!          ux, uy, uz -- Green's functions
!
		use para
		implicit none
		integer :: nlayer,nmax,NAG,N,k
		real*8 :: rlayer(nlayer,2),lambda(nlayer),mu(nlayer,2),rho(nlayer),g(nlayer),tau(nlayer,2),t
		real*8 :: dept,angle(NAG)
		integer,parameter :: nm=4,ns=4,nt=2
		real*8 :: ux(NAG,nm),uy(NAG,nm),uz(NAG,nm)
		real*8 :: geoid(NAG,nm),dgf(NAG,nm),dgd(NAG,nm)
		real*8 :: HN(0:nmax,ns),LN(0:nmax,ns),KN(0:nmax,ns),LNT(0:nmax,nt)
		real*8 :: work1(ns),work2(ns),work4(ns),work3(nt)
		real*8 :: tem1(nm),tem2(nm),tem3(nm),tem4(nm),tem5(nm),tem6(nm)

		HN = 0.d0
		LN = 0.d0
		KN = 0.d0
		LNT = 0.d0
		do N=2,nmax
			call LOVESLAP(nlayer,rlayer,lambda,mu,rho,g,tau,t,N,dept,work1,work2,work4)
			HN(N,1:ns) = work1(1:ns)
			LN(N,1:ns) = work2(1:ns)
			KN(N,1:ns) = work4(1:ns)
			call LOVETLAP(nlayer,rlayer,lambda,mu,tau,t,N,dept,work3)
			LNT(N,1:nt) = work3(1:nt)
		end do

		do k=1,NAG
			call GREENsingle(nmax,angle(k),HN,LN,LNT,KN,tem1,tem2,tem3,tem4,tem5,tem6)
			ux(k,1:nm) = tem1(1:nm)
			uy(k,1:nm) = tem2(1:nm)
			uz(k,1:nm) = tem3(1:nm)
			geoid(k,1:nm) = tem4(1:nm)
			dgf(k,1:nm) = tem5(1:nm)
			dgd(k,1:nm) = tem6(1:nm)
		end do
	end subroutine viscoGreen

! --------------------------------------------------------

	subroutine legendre(n,m,theta,pnm,pnmd)
!
!           This subroutine can be used to compute non-normalized legendre function
!           and its differentials
!     Input:  n      -- degree
!             m      -- order, (must smaller than n)
!             theta  -- co-latitude, unit: arc
!     return: pnm    -- legendre function
!             pnmd   -- differential of legendre function
!     porgram by Zhou Xin, 6.21.2011, GUCAS, Beijing
!
	use para
	implicit none
	integer n,m,i,j,kkk
	real*8 pnm(0:n,0:m),pnmd(0:n,0:m),theta,w1,f1,f2,f3,f4,f5,t1,t2,s1,s2

	t1=dsin(theta)
	t2=dcos(theta)
	do i=0,n
		do j=0,m
			pnm(i,j)=0.d0
			pnmd(i,j)=0.d0
		enddo
	enddo
!   compute legendre function
	pnm(0,0)=1.d0
	pnm(1,0)=t2
	pnm(1,1)=t1
	do i=2,n
		f1=2.d0*i-1.d0
		do j=0,m
			if(j.eq.i) then
				pnm(i,j)=f1*t1*pnm(i-1,i-1)
			elseif(j.eq.(i-1)) then
				pnm(i,j)=f1*t2*pnm(i-1,i-1)
			else
				f2=(i-j)*1.d0
				f3=f1/f2
				f4=((i+j)*1.d0-1.d0)/f2
				s1=f3*t2*pnm(i-1,j)
				s2=f4*pnm(i-2,j)
				pnm(i,j)=s1-s2
			endif
		enddo
	enddo

!   compute legendre function's differential
	pnmd(0,0)=0.d0
	pnmd(1,0)=-t1
	pnmd(1,1)=t2
	do i=2,n
	f1=2.d0*i-1.d0
		do j=0,m
			if(j.eq.i) then
				pnmd(i,j)=f1*(t2*pnm(i-1,i-1)+t1*pnmd(i-1,i-1))
			elseif(j.eq.(i-1)) then
				pnmd(i,j)=f1*(-t1*pnm(i-1,i-1)+t2*pnmd(i-1,i-1))
			else
				f2=(i-j)*1.d0
				f3=f1/f2
				f4=((i+j)*1.d0-1.d0)/f2
				s1=f3*(-t1*pnm(i-1,j)+t2*pnmd(i-1,j))
				s2=f4*pnmd(i-2,j)
				pnmd(i,j)=s1-s2
			endif
		enddo
	enddo

	end subroutine legendre
