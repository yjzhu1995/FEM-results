program main
	use para
    use omp_lib
	implicit none
	integer :: i,j,k
	integer :: nlay0, nlay   ! number of Earth's layers
	integer :: NAG    ! THE NUMBER OF ANGLE OF GREEN'S FUNCTION
	integer :: ndep    ! number of source depths
	integer,parameter :: nm=4,ns=4,nt=2
	integer,allocatable :: degmax(:)  ! max degree
	integer :: Ndeg
	real*8,allocatable :: rlayer(:,:),lambda(:),mu(:,:),rho(:),g(:),eta(:,:),tau(:,:)
	real*8,allocatable :: rlayer0(:,:),lambda0(:),mu0(:,:),rho0(:),gg0(:),eta0(:,:)
	real*8,allocatable :: angle(:), sdepth(:), tim(:)
	character*80,allocatable :: outfl(:)
! 	real*8,allocatable,dimension(:,:,:) :: HN,LN,LNT,ux,uy,uz
	real*8,allocatable,dimension(:,:,:,:) :: geoid,dgf,dgd,ux,uy,uz
	integer :: ios,ix
	real*8,allocatable,dimension(:,:) :: work1,work2,work3,work4,work5,work6
! 	real*8,allocatable,dimension(:,:) :: work4,work5,work6
    real :: start,finish

	start = OMP_get_wtime()
! 	call CPU_TIME(start)

	read(*,*) depth_ob
	write(*,*) "The observation depth =",depth_ob

	open(unit=1, file="input.dat", iostat=ios, status="old", action="read")
	if ( ios /= 0 ) stop "Error opening file parameter.txt"
	call skip(1,ix)
	read(1,*) nlay
	call skip(1,ix)
	! read earth model
	allocate(rlayer0(nlay,2))
	allocate(lambda0(nlay))
	allocate(mu0(nlay,2))
	allocate(rho0(nlay))
	allocate(gg0(nlay))
	allocate(eta0(nlay,2))
	do i=1,nlay
		read(1,*) rlayer0(i,1),rlayer0(i,2),lambda0(i),mu0(i,1),mu0(i,2),eta0(i,1),eta0(i,2),rho0(i),gg0(i)
	end do
	! rebuild earth model
	if(depth_ob>1.d-3) then
		nlay0 = nlay+1
		do i=1,nlay0
			if(depth_ob>=rlayer0(i,1) .and. depth_ob<rlayer0(i,2)) then
				n_obdepth = i
				exit				
			end if
		end do
	else
		nlay0 = nlay
		n_obdepth = 0
	end if
	nlay = nlay0
	write(*,*) "the depth_ob is in layer of", n_obdepth
	write(*,*)
	write(*,*) "layers of the Earth model are", nlay
	allocate(rlayer(nlay,2))
	allocate(lambda(nlay))
	allocate(mu(nlay,2))
	allocate(rho(nlay))
	allocate(g(nlay))
	allocate(eta(nlay,2))
	if(depth_ob>1.d-3) then
		rlayer(1:n_obdepth,1) = rlayer0(1:n_obdepth,1)
		rlayer(n_obdepth+1,1) = depth_ob
		rlayer(n_obdepth+2:nlay,1) = rlayer0(n_obdepth+1:nlay-1,1)
		if(n_obdepth/=1) rlayer(1:n_obdepth-1,2) = rlayer0(1:n_obdepth-1,2)
		rlayer(n_obdepth,2) = depth_ob
		rlayer(n_obdepth+1:nlay,2) = rlayer0(n_obdepth:nlay-1,2)
		lambda(1:n_obdepth) = lambda0(1:n_obdepth)
		lambda(n_obdepth+1) = lambda0(n_obdepth)
		lambda(n_obdepth+2:nlay) = lambda0(n_obdepth+1:nlay-1)
		mu(1:n_obdepth,1:2) = mu0(1:n_obdepth,1:2)
		mu(n_obdepth+1,1:2) = mu0(n_obdepth,1:2)
		mu(n_obdepth+2:nlay,1:2) = mu0(n_obdepth+1:nlay-1,1:2)
		eta(1:n_obdepth,1:2) = eta0(1:n_obdepth,1:2)
		eta(n_obdepth+1,1:2) = eta0(n_obdepth,1:2)
		eta(n_obdepth+2:nlay,1:2) = eta0(n_obdepth+1:nlay-1,1:2)
		rho(1:n_obdepth) = rho0(1:n_obdepth)
		rho(n_obdepth+1) = rho0(n_obdepth)
		rho(n_obdepth+2:nlay) = rho0(n_obdepth+1:nlay-1)
		g(1:n_obdepth) = gg0(1:n_obdepth)
		g(n_obdepth+1) = gg0(n_obdepth)
		g(n_obdepth+2:nlay) = gg0(n_obdepth+1:nlay-1)
	else
		rlayer = rlayer0
		lambda = lambda0
		mu = mu0
		rho = rho0
		g = gg0
		eta = eta0
	endif
	deallocate(rlayer0)
	deallocate(lambda0)
	deallocate(mu0)
	deallocate(rho0)
	deallocate(gg0)
	deallocate(eta0)
	call skip(1,ix)
	read(1,*) Ndeg
	call skip(1,ix)
	read(1,*) ndep
	allocate(sdepth(ndep))
	call skip(1,ix)
	do i=1,ndep
		read(1,*) sdepth(i)
	end do
	call skip(1,ix)
	read(1,*) np
	allocate(tim(np))
	do i=1,np
		read(1,*) tim(i)
	end do
	call skip(1,ix)
	read(1,*) NAG
	allocate(angle(NAG))
	call skip(1,ix)
	do i=1,NAG
		read(1,*) angle(i)
	end do
	close(1)

	open(22,file='earthmodel.re')
	do i=1,nlay
		write(22,77) rlayer(i,1),rlayer(i,2),lambda(i),mu(i,1),mu(i,2),eta(i,1),eta(i,2),rho(i),g(i)
	end do
	close(22)
77	format(5f9.3,2E13.5,2f9.3)

	allocate(outfl(ndep))
	open(unit=1, file="depoutfile.txt", iostat=ios, status="old", action="read")
	if ( ios /= 0 ) stop "Error opening file depoutfile.txt"
	call skip(1,ix)
	do i=1,ndep
		read(1,*) outfl(i)
	end do
	close(1)


	pi = 4.d0*datan(1.d0)

	g0 = g(1)

	s2y=1.d0/24.d0/3600.d0/365.25d0      ! year
	allocate(tau(nlay,2))
	do i=1,nlay
		do j=1,2
		tau(i,j)=eta(i,j)/mu(i,2)*s2y                      ! relaxation time, in year
		enddo
	end do

	lambda = lambda*Gpa2pa
	mu = mu*Gpa2pa
	rlayer = rlayer*km2m
	sdepth = sdepth*km2m

	allocate(geoid(np,ndep,NAG,nm))
	allocate(dgf(np,ndep,NAG,nm))
	allocate(dgd(np,ndep,NAG,nm))
	allocate(ux(np,ndep,NAG,nm))
	allocate(uy(np,ndep,NAG,nm))
	allocate(uz(np,ndep,NAG,nm))
	allocate(degmax(ndep))

	allocate(work1(NAG,nm))
	allocate(work2(NAG,nm))
	allocate(work3(NAG,nm))
	allocate(work4(NAG,nm))
	allocate(work5(NAG,nm))
	allocate(work6(NAG,nm))


	do k=1,np
		write(*,*) '*****************************************'
		write(*,*) 'time =',tim(k)
		write(*,*) '-----------------------------------------'
		!$omp parallel do private(i,degmax,work1,work2,work3,work4,work5,work6) &
		!$omp firstprivate(Ndeg,nlay,rlayer,lambda,mu,rho,g,tau,tim,sdepth,NAG,angle)
		do i=1,ndep
			write(*,*) 'source depth =',sdepth(i)
			degmax(i) = int(10*Re/sdepth(i))
			if (Ndeg /= 0) degmax(i) = min(Ndeg,degmax(i))
			call viscoGreen(nlay,rlayer,lambda,mu,rho,g,tau,tim(k),sdepth(i),NAG,angle,degmax(i), &
				work1,work2,work3,work4,work5,work6)
			ux(k,i,1:NAG,1:nm) = work1(1:NAG,1:nm)
			uy(k,i,1:NAG,1:nm) = work2(1:NAG,1:nm)
			uz(k,i,1:NAG,1:nm) = work3(1:NAG,1:nm)
			geoid(k,i,1:NAG,1:nm) = work4(1:NAG,1:nm)
			dgf(k,i,1:NAG,1:nm) = work5(1:NAG,1:nm)
			dgd(k,i,1:NAG,1:nm) = work6(1:NAG,1:nm)
		end do
		!$OMP end parallel do
		write(*,*) '*****************************************'
	end do

	deallocate(work1)
	deallocate(work2)
	deallocate(work3)
	deallocate(work4)
	deallocate(work5)
	deallocate(work6)

	do i=1,ndep
		open(22,file=outfl(i))
		do j=1,NAG
			write(22,66) angle(j),(ux(k,i,j,1),uy(k,i,j,1),uz(k,i,j,1),geoid(k,i,j,1),dgf(k,i,j,1),dgd(k,i,j,1), &
				                   ux(k,i,j,2),uy(k,i,j,2),uz(k,i,j,2),geoid(k,i,j,2),dgf(k,i,j,2),dgd(k,i,j,2), &
				                   ux(k,i,j,3),uy(k,i,j,3),uz(k,i,j,3),geoid(k,i,j,3),dgf(k,i,j,3),dgd(k,i,j,3), &
				                   ux(k,i,j,4),uy(k,i,j,4),uz(k,i,j,4),geoid(k,i,j,4),dgf(k,i,j,4),dgd(k,i,j,4),k=1,np)
		end do
		close(22)
	end do
66	format(f8.4,600E18.8)

	deallocate(outfl)
	deallocate(angle)
	deallocate(sdepth)
	deallocate(ux)
	deallocate(uy)
	deallocate(uz)
	deallocate(geoid)
	deallocate(dgf)
	deallocate(dgd)
	deallocate(degmax)
	deallocate(tim)
	deallocate(rlayer)
	deallocate(lambda)
	deallocate(mu)
	deallocate(eta)
	deallocate(tau)
	deallocate(g)
	deallocate(rho)

	write(*,*) '-------------------'
	write(*,*) 'END PROGRAM'
	write(*,*) '-------------------'
 
    finish = OMP_get_wtime()
!     call CPU_TIME(finish)
    write(*,*) 'it takes',finish-start,'s'

end program	



subroutine skip(nf,ix)
	IMPLICIT REAL*8  (A-H,O-Z)                                        
!*****************************************************************
!  read a line from file nf. if the first character is #, read the
!  next line; if not, backspace
!*****************************************************************
    character*80 :: line

    do 
    	read(nf,'(a80)',end=20) line
    	if(line(1:1).ne.'#') exit
    end do
      backspace(nf)
      ix=0
      return
20    ix=1
      return
end subroutine
