module para
	implicit none
!	integer,parameter :: Ndeg=100000
	real*8,parameter :: Re=6371.d0, UGC=6.6732d-11, g0=9.8d0, km2m=1.d3
	real*8 :: pi
	real*8,parameter :: step=5.d0  ! step size of matrix propagation

contains

	subroutine fundamatrixS(lambda,mu,rho,N,rr,rn,nm,Yr)
!
!      Spheroidal fundamental matrix, 6x6 matrix consisting of ODE's solutions is obtained analytically.
!      Matrix contains singularity solutions.
!      In the matrix Yr, solutions in the first two column are sigularity at the center.
!      (y1,y3,y2,y4,y5,y6)^T
!      y1 -- vertical displacement
!      y3 -- horizontal displacement
!      y2 -- vertical component of vertical stress
!      y4 -- horizontal component of vertical stress
!      y5 -- gravity potential
!      y6 -- related to gravity
!
	implicit none
	real*8 :: lambda,mu,rho
	integer :: i,nm,N
	real*8 :: rr,rn,NN
	real*8 :: Yr(nm,nm)
	real*8 :: rd          ! rd=rr/Re guarantees r^n can be calculated
	REAL*8 :: fl1,fl2,temp

	rd = rr/rn
	NN = dble(N)

	fl1 = dble(N-2)*lambda+dble(N-4)*mu
	fl2 = dble(N+3)*lambda+dble(N+5)*mu

	Yr(1,1) = -(NN+1.D0)
	Yr(1,2) = -((NN+1.D0)*lambda+(NN+3.D0)*mu)/fl1*NN
	Yr(1,3) =  0.d0
	Yr(1,4) =  (NN*lambda+(NN-2.D0)*mu)/fl2*(NN+1.D0)
	Yr(1,5) =  NN
	Yr(1,6) =  0.d0

	Yr(2,1) = 1.d0
	Yr(2,2) = 1.d0
	Yr(2,3) = 0.d0
	Yr(2,4) = 1.d0
	Yr(2,5) = 1.d0
	Yr(2,6) = 0.d0

	Yr(3,1) = 2.d0*mu*(NN+2.D0)*(1.D0+NN)
	Yr(3,2) = 2.d0*mu*((NN**2+3.D0*NN-1.D0)*lambda+(NN*(NN+3.D0))*mu)/fl1*NN
	Yr(3,3) = 0.d0
	Yr(3,4) = 2.d0*mu*((NN**2-NN-3.D0)*lambda+(NN**2-NN-2.D0)*mu)/fl2*(NN+1.D0)
	Yr(3,5) = 2.d0*mu*NN*(NN-1.D0)
	Yr(3,6) = 0.d0

	Yr(4,1) = -2.d0*mu*(NN+2.D0)
	Yr(4,2) = -2.d0*mu*((NN**2-1.D0)*lambda+(NN**2-2.D0)*mu)/fl1
	Yr(4,3) = 0.d0
	Yr(4,4) = 2.d0*mu*((NN**2+2.D0*NN)*lambda+(NN**2+2.D0*NN-1.D0)*mu)/fl2
	Yr(4,5) = 2.d0*mu*(NN-1.D0)
	Yr(4,6) = 0.d0

	temp = 4.d0*pi*UGC*rho
	Yr(5,1) = 0.d0
	Yr(5,2) = -temp*mu*NN/fl1
	Yr(5,3) = 1.d0
	Yr(5,4) = -temp*mu*(NN+1.d0)/fl2
	Yr(5,5) = 0.d0
	Yr(5,6) = 1.d0

	Yr(6,1) = temp*(NN+1.d0)
	Yr(6,2) = temp*(NN+1.d0)*NN*(lambda+mu)/fl1
	Yr(6,3) = 0.d0
	Yr(6,4) = -temp*dble(NN+1)*(NN*lambda+(3.d0*NN+1.d0)*mu)/fl2
	Yr(6,5) = -temp*NN
	Yr(6,6) = 2.d0*NN+1.d0

	Yr(:,1) = Yr(:,1)*rd**(-N-2)
	Yr(:,2) = Yr(:,2)*rd**(-N)
	Yr(:,3) = Yr(:,3)*rd**(-N-2)
	Yr(:,4) = Yr(:,4)*rd**(N+1)
	Yr(:,5) = Yr(:,5)*rd**(N-1)
	Yr(:,6) = Yr(:,6)*rd**(N-1)

	Yr(3,:) = Yr(3,:)/rr
	Yr(4,:) = Yr(4,:)/rr
	Yr(5,:) = Yr(5,:)*rr

	end subroutine fundamatrixS

! --------------------------------------------------------

	subroutine invfundamatrixS(lambda,mu,rho,N,rr,rn,nm,Yr)
!
!      Invserse of Spheroidal fundamental matrix, 6x6 matrix consisting of ODE's solutions is obtained analytically.
!      Matrix contains singularity solutions.
!
	implicit none
	real*8 :: lambda,mu,rho
	integer :: j,nm,N
	real*8 :: rr,rn
	real*8 :: Yr(nm,nm)
	real*8 :: rd          ! rd=rr/Re guarantees r^n can be calculated
	REAL*8 :: fl1,fl2,temp

	rd = rr/rn

	fl1 = dble(2*N+1)*dble(2*N+3)*(lambda+2*mu)
	fl2 = dble(2*N-1)*dble(2*N+1)*(lambda+2*mu)

	Yr = 0.d0

	Yr(1,1) = (dble(N*N-N-3)*lambda+dble(N*N-N-2)*mu)/fl1
	Yr(1,2) = (-mu+(lambda+mu)*dble(N*N+2*N))/fl1*N
	Yr(1,3) = -(N*(lambda+mu)-2*mu)/fl1/2.d0/mu
	Yr(1,4) = -(dble(N+3)*lambda+dble(N+5)*mu)/fl1*N/2.d0/mu

	temp = dble(N-2)*lambda+dble(N-4)*mu
	Yr(2,1) = -temp*dble(N-1)/fl2
	Yr(2,2) = -temp*(N*N-1.d0)/fl2
	Yr(2,3) = temp/fl2/2.d0/mu
	Yr(2,4) = temp*dble(N+1)/fl2/2.d0/mu

	temp = 4.d0*pi*UGC*rho
	Yr(3,1) = -temp*(dble(2*N+3)*lambda+dble(N*N+3*N+2)*mu)/fl1
	Yr(3,2) = -temp*N*(N+1.d0)**2*mu/fl1
	Yr(3,3) = temp*(N+2.d0)/fl1*0.5d0
	Yr(3,4) = temp*N*(N+1.d0)/fl1*0.5d0
	Yr(3,5) = 1.d0
	Yr(3,6) = -1.d0/(2.d0*N+1.d0)

	temp = dble(N+3)*lambda+dble(N+5)*mu
	Yr(4,1) = -temp*dble(N+2)/fl1
	Yr(4,2) = temp*dble(N+2)/fl1*N
	Yr(4,3) = -temp/fl1/2/mu
	Yr(4,4) = temp*N/fl1/2/mu

	Yr(5,1) = (dble(N*N+3*N-1)*lambda+dble(N+3)*N*mu)/fl2
	Yr(5,2) = -(dble(N*N-1)*lambda+dble(N*N-2)*mu)/fl2*(N+1)
	Yr(5,3) = (dble(N+1)*lambda+dble(N+3)*mu)/fl2/2.d0/mu
	Yr(5,4) = -(dble(N-2)*lambda+dble(N-4)*mu)*(N+1)/fl2/2.d0/mu

	temp = 4.d0*pi*UGC*rho
	Yr(6,1) = -temp*(dble(1-2*N)*lambda+N*dble(N-1)*mu)/fl2
	Yr(6,2) = temp*N*N*dble(N+1)*mu/fl2
	Yr(6,3) = -temp*dble(N-1)/fl2*0.5d0
	Yr(6,4) = temp*dble(N+1)*N/fl2*0.5d0
	Yr(6,6) = 1.d0/(2.d0*N+1.d0)


	Yr(1,:) = Yr(1,:)*rd**(N+2)
	Yr(2,:) = Yr(2,:)*rd**N
	Yr(3,:) = Yr(3,:)*rd**(N+2)
	Yr(4,:) = Yr(4,:)*rd**(-N-1)
	Yr(5,:) = Yr(5,:)*rd**(-N+1)
	Yr(6,:) = Yr(6,:)*rd**(-N+1)

	Yr(:,3) = Yr(:,3)*rr
	Yr(:,4) = Yr(:,4)*rr
	Yr(:,5) = Yr(:,5)/rr

	end subroutine invfundamatrixS

! --------------------------------------------------------

	subroutine fundamatrixT(mu,NN,rr,rn,nm,Yr)
!
!      Toroidal fundamental matrix, 2x2 matrix consisting of ODE's solutions is obtained analytically.
!      Matrix contains singularity solutions.
!      In the matrix Yr, solutions in the last column are sigularity at the center.
!
	implicit none
	real*8 :: mu
	integer :: nm,NN
	real*8 :: rr,rn
	real*8 :: Yr(nm,nm)
	real*8 :: rd          ! rd=rr/Re guarantees r^n can be calculated

	rd = rr/rn

	Yr(1,1) = rd**NN
	Yr(1,2) = rd**(-NN-1)

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
	real*8 :: mu
	integer :: i,nm,N
	real*8 :: rr,rn
	real*8 :: Yr(nm,nm)
	real*8 :: rd          ! rd=rr/Re guarantees r^n can be calculated
	real*8 :: fl1

	rd = rr/rn

	fl1 = dble(2*N+1)
	Yr(1,1) = rd**(-N)*dble(N+2)/fl1
	Yr(1,2) = rd**(-N)/fl1/mu*rr
	Yr(2,1) = rd**(N+1)*dble(N-1)/fl1
	Yr(2,2) = -rd**(N+1)/fl1/mu*rr

	end subroutine invfundamatrixT
! --------------------------------------------------------

	SUBROUTINE sourceS(lambda,mu,N,dept,S)
!
!   spheroidal source function
!   S(1)=y1(rs+)-y1(rs-)      ! vertical displacement
!   S(2)=y3(rs+)-y3(rs-)      ! horizontal displacement
!   S(3)=y2(rs+)-y2(rs-)      ! vertical stress
!   S(4)=y4(rs+)-y4(rs-)      ! horizontal stress
!
	implicit none
	real*8 :: lambda,mu
	integer :: N
	real*8 :: dept,RAD,S(6,4)

	RAD = Re-dept
	
	S = 0.d0

! ************************** Strike-slip **********************************

	S(4,1) = -(2.D0*N+1.D0)/8.D0/PI/N/(N+1.D0)/RAD**3*mu

! ***************************** Dip-slip **********************************

	S(2,2) = (2.D0*N+1.D0)/8.D0/PI/N/(N+1.D0)/RAD**2

! **************************** Horizontal tensile**************************

    S(1,3) = (2.D0*N+1.D0)/4.D0/PI/RAD/RAD*lambda/(lambda+2.D0*mu)
	S(3,3) = -(2.D0*N+1.D0)/2.D0/PI/RAD**3 &
     	*mu*(3.D0*lambda+2.D0*mu)/(lambda+2.D0*mu)
    S(4,3) = (2.D0*N+1.D0)/4.D0/PI/RAD**3 &
     	*mu*(3.D0*lambda+2.D0*mu)/(lambda+2.D0*mu)

! **************************** Vertical tensile***************************

    S(1,4) = (2.D0*N+1.D0)/4.D0/PI/RAD/RAD

! *************************************************************************
	END SUBROUTINE sourceS

! --------------------------------------------------------

	SUBROUTINE sourceT(mu,N,dept,ST)
!
!   toroidal source function
!
		implicit none
		real*8 :: mu
		integer :: N
		real*8 :: dept,RAD,ST(2,2)

	RAD = Re-dept
    ST = 0.d0
! ************************** Strike-slip **********************************

	ST(2,1) = -(2.D0*N+1.D0)*mu/8.D0/PI/N/(N+1.D0)/RAD**3

! ***************************** Dip-slip **********************************

	ST(1,2) = -(2.D0*N+1.D0)/8.D0/PI/N/(N+1.D0)/RAD**2

! *************************************************************************
    end subroutine sourceT

! --------------------------------------------------------


    SUBROUTINE LOVES(nlayer,rlayer,lambda,mu,rho,NN,dept,HN,LN,KN)
! 
!   THIS SUBROUTINE TO CALCULATE spheroidal Love number -- HN,LN,KN
! 
	IMPLICIT NONE
	INTEGER :: i,j,NN
	integer, PARAMETER :: N=9,nm=6,nc=3,ns=4
	integer :: nlayer
	real*8 :: rlayer(nlayer,2),lambda(nlayer),mu(nlayer),rho(nlayer)
	REAL*8 :: dept,HN(ns),LN(ns),KN(ns),r1
	REAL*8 :: X(N,ns),C(nm,ns),d(nc,ns)
	real*8 :: Y(nc,nm),Yr(nm,nm) 

	r1 = Re

	call fundamatrixS(lambda(1),mu(1),rho(1),NN,r1,Re,nm,Yr)

	Y(1:2,1:nm) = Yr(1:2,1:nm)
	Y(3,1:nm) = Yr(5,1:nm)

	if(dept < 2.0) then
		call COEFFSHOMO(nlayer,rlayer,lambda,mu,rho,NN,dept,X)
	else 
		CALL COEFFS(nlayer,rlayer,lambda,mu,rho,NN,dept,X)
	endif
	
	C(1:nm,1:ns) = X(1:nm,1:ns)

	d = matmul(Y,C)

	HN(1:ns) = d(1,1:ns)*Re**2
	LN(1:ns) = d(2,1:ns)*Re**2
	KN(1:ns) = d(3,1:ns)*Re**2

	END SUBROUTINE LOVES

! --------------------------------------------------------

    SUBROUTINE LOVET(nlayer,rlayer,mu,NN,dept,LNT)
! 
!   THIS SUBROUTINE TO CALCULATE Toroidal Love number -- HN,LN
! 
    implicit none
	INTEGER :: i,j,NN
	integer,parameter :: nm=2,N=3,nc=1,ns=2
	integer :: nlayer
	real*8 :: rlayer(nlayer,2),mu(nlayer)
	real*8 :: dept
	REAL*8 :: LNT(ns),d(nc,ns),C(nm,ns)
	real*8 :: X(N,ns)
	real*8 :: YT(nm,nm),Y(nc,nm)
	real*8 :: r1

	!  at the Earth's surface
	r1 = Re
	call fundamatrixT(mu(1),NN,r1,Re,nm,YT)
	Y(1,1:nm) = YT(1,1:nm)
	
	if(dept < 2.0) then
		call COEFFTHOMO(nlayer,rlayer,mu,NN,dept,X)
	else 
		call COEFFT(nlayer,rlayer,mu,NN,dept,X)
	endif

	C(1:nm,1:ns) = X(1:nm,1:ns)

	d = matmul(Y,C)
	LNT(1:ns) = d(1,1:ns)*Re**2

    end subroutine LOVET

! --------------------------------------------------------

	subroutine layermatrixS(lambda,mu,rho,r1,r2,NN,nm,matr)
! 	layer matrix used to matrix propagation for spheroid
! 	r1 -- upper bound
!   r2 -- lower bound
		implicit none
		real*8 :: lambda,mu,rho,r1,r2,rd
		integer :: NN,nm
		real*8 :: matr(nm,nm),a(nm,nm),b(nm,nm),c(nm,nm)

		call fundamatrixS(lambda,mu,rho,NN,r1,r1,nm,a)
		call invfundamatrixS(lambda,mu,rho,NN,r2,r2,nm,b)
		rd = r1/r2
		c = 0.d0
		c(1,1) = rd**(-NN-2)
		c(2,2) = rd**(-NN)
		c(3,3) = rd**(-NN-2)
		c(4,4) = rd**(NN+1)
		c(5,5) = rd**(NN-1)
		c(6,6) = rd**(NN-1)

		matr = matmul(a,c)
		matr = matmul(matr,b)
	end subroutine layermatrixS

! --------------------------------------------------------

	subroutine layermatrixT(mu,r1,r2,NN,nm,matr)
! 	layer matrix used to matrix propagation for spheroid
! 	r1 -- upper bound
!   r2 -- lower bound
		implicit none
		real*8 :: mu,r1,r2,rd
		integer :: NN,nm
		real*8 :: matr(nm,nm),a(nm,nm),b(nm,nm),c(nm,nm)

		call fundamatrixT(mu,NN,r1,r1,nm,a)
		call invfundamatrixT(mu,NN,r2,r2,nm,b)
		rd = r1/r2
		c = 0.d0
		c(1,1) = rd**NN
		c(2,2) = rd**(-NN-1)

		matr = matmul(a,c)
		matr = matmul(matr,b)
	end subroutine layermatrixT

! --------------------------------------------------------

	subroutine COEFFS(nlayer,rlayer,lambda,mu,rho,NN,dept,X)
!
!   solve the boundary value equations
!
		implicit none
		integer,parameter :: N=9, nm=6, ns=4
		integer :: NN,i,j,k,nds
		integer :: nlayer
		real*8 :: rlayer(nlayer,2),lambda(nlayer),mu(nlayer),rho(nlayer)
		real*8 :: dept
		REAL*8 :: A(N,N),B(nm,ns),BB(N,ns),X(N,ns)
		real*8 :: Yr1(nm,nm),Yr(nm,nm),Yrs1(nm,nm),Yrs2(nm,nm)
		REAL*8 :: r1,r2,r_upper
		integer :: IPIV(N), info
		external dgesv

	! Earth's surface
	r1 = Re-rlayer(1,1)
	call fundamatrixS(lambda(1),mu(1),rho(1),NN,r1,Re,nm,Yr)

	!  depth of source
	!  downward propagate solutions from Earth's surface to the source's depth
	Yrs1 = 0.d0
	do i=1,nm
		Yrs1(i,i) = 1.d0
	end do
	k = 1
	do 
		if(rlayer(k,2)>dept) then
			r2 = Re-dept
			nds = k
			call fundamatrixS(lambda(k),mu(k),rho(k),NN,r2,Re,nm,Yr1)
			Yrs1 = matmul(Yr1,Yrs1)
			exit
		else
			r2 = Re-rlayer(k,2)
			call fundamatrixS(lambda(k),mu(k),rho(k),NN,r2,Re,nm,Yr1)
			Yrs1 = matmul(Yr1,Yrs1)
		endif
		k = k+1
		r1 = Re-rlayer(k,1)     ! upper boundary of each layer
		call invfundamatrixS(lambda(k),mu(k),rho(k),NN,r1,Re,nm,Yr1)
		Yrs1 = matmul(Yr1,Yrs1)
	end do

	!  upward propagate solutions from Earth's inner layer to the source's depth
	r1 = Re-rlayer(nlayer,1)
	call fundamatrixS(lambda(nlayer),mu(nlayer),rho(nlayer),NN,r1,r1,nm,Yrs2)
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
			call layermatrixS(lambda(k),mu(k),rho(k),r1,r2,NN,nm,Yr1)
			Yr1 = Yr1/maxval(dabs(Yr1))
			Yrs2 = matmul(Yr1,Yrs2)
 			Yrs2 = Yrs2/maxval(dabs(Yrs2))
			r2 = r1
		enddo
	end do

	!  construct linear system to solve the coefficients
	A(1:2,1:nm) = Yr(3:4,1:nm)
	A(3,1:nm) = Yr(6,1:nm)
	A(1:3,nm+1:N) = 0.d0
	do i=4,N
		do j=1,nm
			A(i,j) =  Yrs1(i-3,j)
		enddo
		do j=nm+1,N
			A(i,j) = -Yrs2(i-3,j-3)
		enddo
	enddo

	call sourceS(lambda(nds),mu(nds),NN,dept,B)
	BB(1:3,1:ns) = 0.d0
	BB(4:N,1:ns) = B(1:nm,1:ns)

	call dgesv(N,ns,A,N,IPIV,BB,N,info)
	X = BB
	end subroutine COEFFS

! --------------------------------------------------------

	subroutine COEFFT(nlayer,rlayer,mu,NN,dept,X)
		implicit none
		INTEGER i,j,k,NN,nds
		integer,parameter :: nm=2,N=3,ns=2
		integer :: nlayer
		real*8 :: rlayer(nlayer,2),mu(nlayer)
		real*8 :: dept
		real*8 :: YT(nm,nm),YT1(nm,nm),YTs1(nm,nm),YTs2(nm,nm)
		real*8 :: A(N,N),B(nm,ns),BB(N,ns),X(N,ns)    ! work array for linear system
		real*8 :: r1,r2,r_upper
		integer :: IPIV(N), info
		external dgesv

	!  at the Earth's surface
	r1 = Re-rlayer(1,1)
	call fundamatrixT(mu(1),NN,r1,Re,nm,YT)

	!  depth of source
	!  downward propagate solutions from Earth's surface to the source's depth
	YTs1 = 0.d0
	do i=1,nm
		YTs1(i,i) = 1.d0
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
			YT1 = YT1/maxval(dabs(YT1))
			YTs2 = matmul(YT1,YTs2)
 			YTs2 = YTs2/maxval(dabs(YTs2))
			r2 = r1
		enddo
	end do

	!  construct array
	do j=1,nm
		A(1,j) = YT(2,j)
	enddo
		A(1,N) = 0.d0
	do i=2,N
		do j=1,nm
			A(i,j) =  YTs1(i-1,j)
		enddo
			A(i,N) = -YTs2(i-1,1)
	enddo

	call sourceT(mu(nds),NN,dept,B)
	BB(1,1:ns) = 0.d0
	BB(2:N,1:ns) = B(1:nm,1:ns)

	call dgesv(N,ns,A,N,IPIV,BB,N,info)
	X = BB
	end subroutine COEFFT

! --------------------------------------------------------

	subroutine COEFFSHOMO(nlayer,rlayer,lambda,mu,rho,NN,dept,X)
!
!   solve the boundary value equations
!
		implicit none
		integer,parameter :: N=9, nm=6, ns=4
		integer :: NN,i,j
		integer :: nlayer
		real*8 :: rlayer(nlayer,2),lambda(nlayer),mu(nlayer),rho(nlayer)
		real*8 :: dept
		REAL*8 :: A(N,N),B(nm,ns),BB(N,ns),X(N,ns)
		real*8 :: Yr(nm,nm),Yrs(nm,nm)
		REAL*8 :: r1,r2,r_upper
		integer :: IPIV(N), info
		external dgesv

	! Earth's surface
	r1 = Re-rlayer(1,1)
	call fundamatrixS(lambda(1),mu(1),rho(1),NN,r1,Re,nm,Yr)

	!  depth of source
	r2 = Re-dept
	call fundamatrixS(lambda(1),mu(1),rho(1),NN,r2,Re,nm,Yrs)

	!  construct linear system to solve the coefficients
	A(1:2,1:nm) = Yr(3:4,1:nm)
	A(3,1:nm) = Yr(6,1:nm)
	A(1:3,nm+1:N) = 0.d0
	do i=4,N
		do j=1,nm
			A(i,j) =  Yrs(i-3,j)
		enddo
		do j=nm+1,N
			A(i,j) = -Yrs(i-3,j-3)
		enddo
	enddo

	call sourceS(lambda(1),mu(1),NN,dept,B)
	BB(1:3,1:ns) = 0.d0
	BB(4:N,1:ns) = B(1:nm,1:ns)

	call dgesv(N,ns,A,N,IPIV,BB,N,info)
	X = BB
	end subroutine COEFFSHOMO

! --------------------------------------------------------

	subroutine COEFFTHOMO(nlayer,rlayer,mu,NN,dept,X)
		implicit none
		INTEGER i,j,k,NN
		integer,parameter :: nm=2,N=3,ns=2
		integer :: nlayer
		real*8 :: rlayer(nlayer,2),mu(nlayer)
		real*8 :: dept
		real*8 :: YT(nm,nm),YTs(nm,nm)
		real*8 :: A(N,N),B(nm,ns),BB(N,ns),X(N,ns)    ! work array for linear system
		real*8 :: r1,r2,r_upper
		integer :: IPIV(N), info
		external dgesv

	!  at the Earth's surface
	r1 = Re-rlayer(1,1)
	call fundamatrixT(mu(1),NN,r1,Re,nm,YT)

	!  depth of source
	r2 = Re-dept
	call fundamatrixT(mu(1),NN,r2,Re,nm,YTs)

	!  construct array
	do j=1,nm
		A(1,j) = YT(2,j)
	enddo
		A(1,N) = 0.d0
	do i=2,N
		do j=1,nm
			A(i,j) =  YTs(i-1,j)
		enddo
			A(i,N) = -YTs(i-1,1)
	enddo

	call sourceT(mu(1),NN,dept,B)
	BB(1,1:ns) = 0.d0
	BB(2:N,1:ns) = B(1:nm,1:ns)

	call dgesv(N,ns,A,N,IPIV,BB,N,info)
	X = BB
	end subroutine COEFFTHOMO

! --------------------------------------------------------

	subroutine GREENsingle(nmax,theta,HN,LN,LNT,KN,ux,uy,uz,geoid,dgf,dgd)
!
!   geoid(1),dgf(1),dgd(1) -- u12, strike-slip for a vertical fault
!   geoid(2),dgf(2),dgd(2) -- u32, dip-slip for a vertical fault
!   geoid(3),dgf(3),dgd(3) -- u220, horizontal tensile for m=0
!   geoid(4),dgf(4),dgd(4) -- u33, vertical tensile
!
		implicit none
		integer :: nmax,N,NI(4)
		integer,parameter :: nm=4,ns=4,nt=2
		real*8 :: theta
		real*8 :: geoid(nm),dgf(nm),dgd(nm),uz(nm),ux(nm),uy(nm)
		real*8 :: HN(0:nmax,ns),LN(0:nmax,ns),LNT(0:nmax,nt),KN(0:nmax,ns)
		real*8 :: Pnm(0:nmax,0:2),pdnm(0:nmax,0:2)
		real*8 :: theta1
	
		THETA1=THETA*PI/180.D0

		CALL legendre(NMAX,2,THETA1,PNM,PDNM)

		ux=0.d0
		uy=0.d0
		uz=0.d0
		geoid=0.D0
		dgf=0.D0
		dgd=0.D0

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
			UY(2)=UY(2)+(-2.D0)*(LN(N,2)*PDNM(N,1)-LNT(N,2)*PNM(N,1)/DSIN(THETA1))   ! different from Sun's book where positive of toroidal term is used
			UX(2)=UX(2)+(-2.D0)*(LN(N,2)*PNM(N,1)/DSIN(THETA1)-LNT(N,2)*PDNM(N,1))

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
			uz = 0.d0
			UX = 0.d0
			UY = 0.d0
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
			uz = 0.d0
			UX = 0.d0
			UY = 0.d0
			do N=2,nmax
				geoid(3) = geoid(3)+(-1)**N*HN(N,3)
				geoid(4) = geoid(4)+(-1)**N*HN(N,4)
				dgf(3) = dgf(3)+(-1)**N*KN(N,3)*dble(N+1)/Re
				dgf(4) = dgf(4)+(-1)**N*KN(N,4)*dble(N+1)/Re
				uz(3) = uz(3)+(-1)**N*HN(N,3)
				uz(4) = uz(4)+(-1)**N*HN(N,4)
			end do
		endif

		geoid = geoid/g0*km2m
		dgd(1:nm) = dgf(1:nm) - 2.d0*g0*uz(1:nm)/Re/km2m

	end subroutine GREENsingle

! --------------------------------------------------------

	subroutine Green(nlayer,rlayer,lambda,mu,rho,dept,NAG,angle,nmax,ux,uy,uz,geoid,dgf,dgd)
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
		implicit none
		integer :: nlayer,nmax,NAG,N,k
		real*8 :: rlayer(nlayer,2),lambda(nlayer),mu(nlayer),rho(nlayer)
		real*8 :: dept,angle(NAG)
		integer,parameter :: nm=4,ns=4,nt=2
		real*8 :: ux(NAG,nm),uy(NAG,nm),uz(NAG,nm)
		real*8 :: geoid(NAG,nm),dgf(NAG,nm),dgd(NAG,nm)
!		real*8 :: HN(0:Ndeg,ns),LN(0:Ndeg,ns),LNT(0:Ndeg,nt)
		real*8 :: HN(0:nmax,ns),LN(0:nmax,ns),KN(0:nmax,ns),LNT(0:nmax,nt)
		real*8 :: work1(ns),work2(ns),work4(ns),work3(nt)
		real*8 :: tem1(nm),tem2(nm),tem3(nm),tem4(nm),tem5(nm),tem6(nm)

! 		nmax = int(2*pi*Re/dept)
! 		nmax = min(nmax,Ndeg)

		HN = 0.d0
		LN = 0.d0
		LNT = 0.d0
		do N=2,nmax
			call LOVES(nlayer,rlayer,lambda,mu,rho,N,dept,work1,work2,work4)
			HN(N,1:ns) = work1(1:ns)
			LN(N,1:ns) = work2(1:ns)
			KN(N,1:ns) = work4(1:ns)
			call LOVET(nlayer,rlayer,mu,N,dept,work3)
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
	end subroutine Green

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

end module
