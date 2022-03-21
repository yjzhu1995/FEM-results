program main
	use para
	implicit none
	integer :: i,j,k
	integer :: nlay   ! number of Earth's layers
	integer :: NAG    ! THE NUMBER OF ANGLE OF GREEN'S FUNCTION
	integer :: ndep    ! number of source depths
	integer,parameter :: nm=4,ns=4,nt=2
	integer,allocatable :: degmax(:)  ! max degree
	real*8,allocatable :: rlayer(:,:),lambda(:),mu(:),rho(:),angle(:), sdepth(:)
	character*80,allocatable :: outfl(:)
! 	real*8,allocatable,dimension(:,:,:) :: HN,LN,LNT,ux,uy,uz
	real*8,allocatable,dimension(:,:,:) :: ux,uy,uz
	real*8,allocatable,dimension(:,:,:) :: geoid,dgf,dgd
	integer :: ios,ix
	real*8,allocatable,dimension(:,:) :: work1,work2,work3,work4,work5,work6
! 	real*8,allocatable,dimension(:,:) :: work4,work5,work6
    real :: start,finish

	call CPU_TIME(start)

	open(unit=1, file="input.dat", iostat=ios, status="old", action="read")
	if ( ios /= 0 ) stop "Error opening file parameter.txt"
	call skip(1,ix)
	read(1,*) nlay
	call skip(1,ix)
	allocate(rlayer(nlay,2))
	allocate(lambda(nlay))
	allocate(mu(nlay))
	allocate(rho(nlay))
	do i=1,nlay
		read(1,*) rlayer(i,1),rlayer(i,2),lambda(i),mu(i),rho(i)
	end do
	call skip(1,ix)
	read(1,*) ndep
	allocate(sdepth(ndep))
	call skip(1,ix)
	do i=1,ndep
		read(1,*) sdepth(i)
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

	allocate(ux(ndep,NAG,nm))
	allocate(uy(ndep,NAG,nm))
	allocate(uz(ndep,NAG,nm))
	allocate(geoid(ndep,NAG,nm))
	allocate(dgf(ndep,NAG,nm))
	allocate(dgd(ndep,NAG,nm))
	allocate(degmax(ndep))

	allocate(work1(NAG,nm))
	allocate(work2(NAG,nm))
	allocate(work3(NAG,nm))
	allocate(work4(NAG,nm))
	allocate(work5(NAG,nm))
	allocate(work6(NAG,nm))


	do i=1,ndep
		write(*,*) "depth =", sdepth(i)
		degmax(i) = int(10*Re/sdepth(i))
		call Green(nlay,rlayer,lambda,mu,rho,sdepth(i),NAG,angle,degmax(i), &
			work1,work2,work3,work4,work5,work6)
		ux(i,1:NAG,1:nm) = work1(1:NAG,1:nm)
		uy(i,1:NAG,1:nm) = work2(1:NAG,1:nm)
		uz(i,1:NAG,1:nm) = work3(1:NAG,1:nm)
		geoid(i,1:NAG,1:nm) = work4(1:NAG,1:nm)
		dgf(i,1:NAG,1:nm) = work5(1:NAG,1:nm)
		dgd(i,1:NAG,1:nm) = work6(1:NAG,1:nm)
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
			write(22,66) angle(j), ux(i,j,1),uy(i,j,1),uz(i,j,1),geoid(i,j,1),dgf(i,j,1),dgd(i,j,1), &
				                   ux(i,j,2),uy(i,j,2),uz(i,j,2),geoid(i,j,2),dgf(i,j,2),dgd(i,j,2), &
				                   ux(i,j,3),uy(i,j,3),uz(i,j,3),geoid(i,j,3),dgf(i,j,3),dgd(i,j,3), &
				                   ux(i,j,4),uy(i,j,4),uz(i,j,4),geoid(i,j,4),dgf(i,j,4),dgd(i,j,4)
		end do
		close(22)
	end do
66	format(f8.4,24E18.8)

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
	deallocate(lambda)
	deallocate(mu)
	deallocate(rho)
	deallocate(rlayer)

	write(*,*) '-------------------'
	write(*,*) 'END PROGRAM'
	write(*,*) '-------------------'
 
    call CPU_TIME(finish)
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
