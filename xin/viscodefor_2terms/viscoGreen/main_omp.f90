program main
	use para
    use omp_lib
	implicit none
	integer :: i,j,k,l
	integer :: nlay   ! number of Earth's layers
	integer :: NAG    ! THE NUMBER OF ANGLE OF GREEN'S FUNCTION
	integer :: ndep    ! number of source depths
	integer,parameter :: nm=4,ns=4,nt=2
	integer :: degmax  ! max degree
! 	integer,allocatable :: degmax(:)  ! max degree
	integer :: Ndeg
	real*8,allocatable :: rlayer(:,:),lambda(:),mu(:,:),rho(:),g(:),eta(:,:),tau(:,:)
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
	allocate(rlayer(nlay,2))
	allocate(lambda(nlay))
	allocate(mu(nlay,2))
	allocate(rho(nlay))
	allocate(g(nlay))
	allocate(eta(nlay,2))
	do i=1,nlay
		read(1,*) rlayer(i,1),rlayer(i,2),lambda(i),mu(i,1),mu(i,2),eta(i,1),eta(i,2),rho(i),g(i)
	end do
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

	allocate(work1(NAG,nm))
	allocate(work2(NAG,nm))
	allocate(work3(NAG,nm))
	allocate(work4(NAG,nm))
	allocate(work5(NAG,nm))
	allocate(work6(NAG,nm))


	!!$omp parallel do private(i,k,degmax,work1,work2,work3,work4,work5,work6) &
	!!$omp firstprivate(Ndeg,nlay,rlayer,lambda,mu,rho,g,tau,tim,sdepth,NAG,angle)
	do k=1,np
		write(*,*) '*****************************************'
		write(*,*) 'time =',tim(k)
		write(*,*) '-----------------------------------------'
		!$omp parallel do private(i,degmax,work1,work2,work3,work4,work5,work6) &
		!$omp firstprivate(Ndeg,nlay,rlayer,lambda,mu,rho,g,tau,tim,sdepth,NAG,angle)
		do i=1,ndep
		write(*,*) 'source depth =',sdepth(i)
		degmax = int(10*Re/sdepth(i))
		if (Ndeg /= 0) degmax = min(Ndeg,degmax)
		call viscoGreen(nlay,rlayer,lambda,mu,rho,g,tau,tim(k),sdepth(i),NAG,angle,degmax, &
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
	!!$OMP end parallel do

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
