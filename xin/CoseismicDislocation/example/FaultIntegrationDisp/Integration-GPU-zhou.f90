module class
	real*8, PARAMETER :: PI=3.141592653589793238462643383279502884197
	real*8 :: deg2arc
	real*8, parameter :: Re=6371.d0
!$acc declare create(deg2arc,PI,Re)
	type source
        real*8 :: lat, lon, depth, slip, rake, strike, dip, length, width
	end type source
	type station
		real*8 :: lon, lat, ux, uy, uz
	end type station
end module class



program Integ
	use class
	use openacc
	implicit none
	integer :: i,j,k,m,n,l,ios,stat,ii,jj
	integer :: ndepth, nepdis, ngreenls
	character*128 :: depthfl, epdisfl, greenls
	character*128, allocatable :: greenfl(:)
	real*8,allocatable,dimension(:) :: depth, epdis
	real*8,allocatable :: greenval(:,:,:,:)
	real*8 :: u12x,u12y,u12z,u32x,u32y,u32z,u33x,u33y,u33z,u220x,u220y,u220z
	character*128 :: faultfl, obsfl
	integer :: nfsm, nobs
	type(source),allocatable :: fsm(:),subfsm(:,:,:)
	type(station),allocatable :: obs(:)
	integer :: nsublength, nsubwidth
	real*8 :: dist,u(3),temp,threshold
	real*8 :: sum1,sum2,sum3
    real :: start,finish
    type(source),allocatable :: tt(:,:)
	real*8,allocatable,dimension(:,:) :: green12x,green12y,green12z,green32x, &
	                                     green32y,green32z,green33x,green33y, &
	                                     green33z,green220x,green220y,green220z
	real*8,allocatable,dimension(:,:,:) :: coeffx12,coeffy12,coeffz12,coeffx32, &
	                                     coeffy32,coeffz32,coeffx33,coeffy33, &
	                                     coeffz33,coeffx220,coeffy220,coeffz220
!$acc routine(subfaults) worker
!$acc routine(displ,subfaults,greatcircle) seq
!$acc routine(spline,splint) vector

	call CPU_TIME(start)

	deg2arc = PI/180.d0

!  read parameter file including depth, epicenter distance and green routes
	write(*,*) '-------------------'
	write(*,*) 'reading data files'

	open(unit=1,file='parameter.txt',iostat=ios,status='old')
	if ( ios /= 0 ) stop "Error opening file parameter"
	read(1,*)
	read(1,*) depthfl
	read(1,*)
	read(1,*) epdisfl
	read(1,*)
	read(1,*) greenls
	read(1,*)
	read(1,*) faultfl
	read(1,*)
	read(1,*) obsfl
	read(1,*)
	read(1,*) nsublength, nsubwidth
	close(1)


!  read depth file
	open(unit=2, file=depthfl, iostat=ios, status="old", action="read")
	if ( ios /= 0 ) stop "Error opening file depth"
	ndepth = 0
	do 
		read(2,*,iostat=stat) 
		if (stat<0) exit
		ndepth = ndepth + 1
	end do

	allocate(depth(ndepth))

	rewind(2)
	do k=1,ndepth
		read(2,*) depth(k)
	end do
	close(2)


!   read epicenter distance file
	open(unit=3, file=epdisfl, iostat=ios, status="old", action="read")
	if ( ios /= 0 ) stop "Error opening epicenter distance file"
	nepdis = 0
	do 
		read(3,*,iostat=stat) 
		if (stat<0) exit
		nepdis = nepdis + 1
	end do

	allocate(epdis(nepdis))

	rewind(3)
	do k=1,nepdis
		read(3,*) epdis(k)
	end do
	close(3)

	threshold = epdis(2)*deg2arc

!  read green's function route from list file
	open(unit=4, file=greenls, iostat=ios, status="old", action="read")
	if ( ios /= 0 ) stop "Error opening green list file"
	ngreenls = 0
	do 
		read(4,*,iostat=stat) 
		if (stat<0) exit
		ngreenls = ngreenls + 1
	end do
	close(4)

	! judge ngreenls == ndepth
	if(ndepth /= ngreenls) stop "Error number of green function or number of depth"

	allocate(greenfl(ngreenls))

	open(unit=4, file=greenls, iostat=ios, status="old", action="read")
	if ( ios /= 0 ) stop "Error opening file green funtion list file"
	do k=1,ngreenls
		read(4,*) greenfl(k)
	end do
	close(4)


!  *****

	allocate(greenval(4,ndepth,nepdis,6))

!  read green's function
	do k=1,ndepth
		open(unit=20, file=greenfl(k), iostat=ios, status="old", action="read")
		if ( ios /= 0 ) stop "Error opening file green function file"
		do i=1,nepdis
			read(20,*) temp, greenval(1,k,i,1),greenval(1,k,i,2),greenval(1,k,i,3), &
			                 greenval(1,k,i,4),greenval(1,k,i,5),greenval(1,k,i,6), &
			                 greenval(2,k,i,1),greenval(2,k,i,2),greenval(2,k,i,3), &
			                 greenval(2,k,i,4),greenval(2,k,i,5),greenval(2,k,i,6), &
			                 greenval(4,k,i,1),greenval(4,k,i,2),greenval(4,k,i,3), &
			                 greenval(4,k,i,4),greenval(4,k,i,5),greenval(4,k,i,6), &
			                 greenval(3,k,i,1),greenval(3,k,i,2),greenval(3,k,i,3), &
			                 greenval(3,k,i,4),greenval(3,k,i,5),greenval(3,k,i,6)
		end do
		close(20)
	end do

	deallocate(greenfl)


!  read fault slip model file
	open(unit=5, file=faultfl, iostat=ios, status="old", action="read")
	if ( ios /= 0 ) stop "Error opening fault file"

	read(5,*)
	read(5,*) nfsm

	allocate(fsm(nfsm))

	do k=1,nfsm
		read(5,*) fsm(k)
	end do
	close(5)


!  read observation station file
	open(unit=16, file=obsfl, iostat=ios, status="old", action="read")
	if ( ios /= 0 ) stop "Error opening observation file"
	read(16,*)
	nobs = 0
	do 
		read(16,*,iostat=stat) 
		if (stat<0) exit
		nobs = nobs + 1
	end do

	allocate(obs(nobs))

	rewind(16)
	read(16,*)
	do k=1,nobs
		read(16,*) obs(k)%lat, obs(k)%lon
	end do
	close(16)

	write(*,*) '-------------------'
	write(*,*) 'finish data reading'


!****************************************************
!       inperpolation and convolution
!****************************************************
	allocate(subfsm(nfsm,nsubwidth,nsublength))
	allocate(green12x(ndepth,nepdis))
	allocate(green12y(ndepth,nepdis))
	allocate(green12z(ndepth,nepdis))
	allocate(green32x(ndepth,nepdis))
	allocate(green32y(ndepth,nepdis))
	allocate(green32z(ndepth,nepdis))
	allocate(green33x(ndepth,nepdis))
	allocate(green33y(ndepth,nepdis))
	allocate(green33z(ndepth,nepdis))
	allocate(green220x(ndepth,nepdis))
	allocate(green220y(ndepth,nepdis))
	allocate(green220z(ndepth,nepdis))
	allocate(tt(nsubwidth,nsublength))

	allocate(coeffx12(16,ndepth,nepdis))
	allocate(coeffy12(16,ndepth,nepdis))
	allocate(coeffz12(16,ndepth,nepdis))
	allocate(coeffx32(16,ndepth,nepdis))
	allocate(coeffy32(16,ndepth,nepdis))
	allocate(coeffz32(16,ndepth,nepdis))
	allocate(coeffx33(16,ndepth,nepdis))
	allocate(coeffy33(16,ndepth,nepdis))
	allocate(coeffz33(16,ndepth,nepdis))
	allocate(coeffx220(16,ndepth,nepdis))
	allocate(coeffy220(16,ndepth,nepdis))
	allocate(coeffz220(16,ndepth,nepdis))


!  subfault gridding
	write(*,*) '-----------------------------'
	write(*,*) 'fault disretization'

	!$acc data copy(fsm) copyout(subfsm)
	do k=1,nfsm
		call subfaults(fsm(k), nsubwidth, nsublength, tt)
		!$acc loop collapse(2) independent
		do i=1,nsubwidth
			do j=1,nsublength
				subfsm(k,i,j) = tt(i,j)
			end do
		end do
	end do
	!$acc end data
	

!  interpolation
	write(*,*) '-----------------------------'
	write(*,*) 'interpolation and integration'

	green12x(:,:) = greenval(1,:,:,1)
	green12y(:,:) = greenval(1,:,:,2)
	green12z(:,:) = greenval(1,:,:,3)
	green32x(:,:) = greenval(2,:,:,1)
	green32y(:,:) = greenval(2,:,:,2)
	green32z(:,:) = greenval(2,:,:,3)
	green33x(:,:) = greenval(3,:,:,1)
	green33y(:,:) = greenval(3,:,:,2)
	green33z(:,:) = greenval(3,:,:,3)
	green220x(:,:) = greenval(4,:,:,1)
	green220y(:,:) = greenval(4,:,:,2)
	green220z(:,:) = greenval(4,:,:,3)

	!$acc data copy(ndepth,nepdis,depth,epdis) &
	!$acc& copyin(green12x,green12y,green12z,green32x,green32y,green32z,green33x,green33y,green33z,green220x,green220y,green220z) &
	!$acc& copyout(coeffx12,coeffy12,coeffz12,coeffx32,coeffy32,coeffz32,coeffx33,coeffy33,coeffz33,coeffx220,coeffy220,coeffz220)
	call spline(ndepth,nepdis,depth,epdis,green12x,coeffx12)
	call spline(ndepth,nepdis,depth,epdis,green12y,coeffy12)
	call spline(ndepth,nepdis,depth,epdis,green12z,coeffz12)
	call spline(ndepth,nepdis,depth,epdis,green32x,coeffx32)
	call spline(ndepth,nepdis,depth,epdis,green32y,coeffy32)
	call spline(ndepth,nepdis,depth,epdis,green32z,coeffz32)
	call spline(ndepth,nepdis,depth,epdis,green33x,coeffx33)
	call spline(ndepth,nepdis,depth,epdis,green33y,coeffy33)
	call spline(ndepth,nepdis,depth,epdis,green33z,coeffz33)
	call spline(ndepth,nepdis,depth,epdis,green220x,coeffx220)
	call spline(ndepth,nepdis,depth,epdis,green220y,coeffy220)
	call spline(ndepth,nepdis,depth,epdis,green220z,coeffz220)
	!$acc end data

	!$acc data copy(depth,epdis,subfsm,obs) &
	!$acc& copyin(coeffx12,coeffy12,coeffz12,coeffx32,coeffy32,coeffz32,coeffx33,coeffy33,coeffz33,coeffx220,coeffy220,coeffz220)
	do m=1,nobs
		sum1 = 0.d0
		sum2 = 0.d0
		sum3 = 0.d0
		!$acc loop gang worker vector private(dist,u12x,u12y,u12z,u32x,u32y,u32z,u33x,u33y,u33z,u220x,u220y,u220z) reduction(+:sum1,sum2,sum3)
		do k=1,nfsm
			do i=1,nsubwidth
				do j=1,nsublength
					call greatcircle(subfsm(k,i,j)%lat, subfsm(k,i,j)%lon, obs(m)%lat, obs(m)%lon, dist)
					dist = dacos(dist)/deg2arc
					if(dist < epdis(2)) dist = epdis(2)
					call splint(ndepth,nepdis,depth,epdis,coeffx12,subfsm(k,i,j)%depth,dist,u12x)
					call splint(ndepth,nepdis,depth,epdis,coeffy12,subfsm(k,i,j)%depth,dist,u12y)
					call splint(ndepth,nepdis,depth,epdis,coeffz12,subfsm(k,i,j)%depth,dist,u12z)
					call splint(ndepth,nepdis,depth,epdis,coeffx32,subfsm(k,i,j)%depth,dist,u32x)
					call splint(ndepth,nepdis,depth,epdis,coeffy32,subfsm(k,i,j)%depth,dist,u32y)
					call splint(ndepth,nepdis,depth,epdis,coeffz32,subfsm(k,i,j)%depth,dist,u32z)
					call splint(ndepth,nepdis,depth,epdis,coeffx33,subfsm(k,i,j)%depth,dist,u33x)
					call splint(ndepth,nepdis,depth,epdis,coeffy33,subfsm(k,i,j)%depth,dist,u33y)
					call splint(ndepth,nepdis,depth,epdis,coeffz33,subfsm(k,i,j)%depth,dist,u33z)
					call splint(ndepth,nepdis,depth,epdis,coeffx220,subfsm(k,i,j)%depth,dist,u220x)
					call splint(ndepth,nepdis,depth,epdis,coeffy220,subfsm(k,i,j)%depth,dist,u220y)
					call splint(ndepth,nepdis,depth,epdis,coeffz220,subfsm(k,i,j)%depth,dist,u220z)
					call displ(subfsm(k,i,j),obs(m),u12x,u12y,u12z,u32x,u32y,u32z,u33x,u33y,u33z,u220x,u220y,u220z,u)
					sum1 = sum1+u(1)
					sum2 = sum2+u(2)
					sum3 = sum3+u(3)
				end do
			end do
		end do
		obs(m)%ux = sum1
		obs(m)%uy = sum2
		obs(m)%uz = sum3
	end do
	!$acc end data


!  output subfaults
	open(unit=17, file='subfaults.dat', iostat=ios, action="write")
	if ( ios /= 0 ) stop "Error opening file name"
	do k=1,nfsm
		do i=1,nsubwidth
			do j=1,nsublength
				write(17,*) subfsm(k,i,j)%lon, subfsm(k,i,j)%lat
			end do
		end do
	end do
	close(17)

!  write output file
	open(unit=21, file="displacements.dat", iostat=ios, action="write")
	if ( ios /= 0 ) stop "Error opening file displacements"
	do i=1,nobs
		write(21,33) obs(i)
	end do
	close(21)
33  format(2f10.5,3E20.10)


	deallocate(coeffx12)
	deallocate(coeffy12)
	deallocate(coeffz12)
	deallocate(coeffx32)
	deallocate(coeffy32)
	deallocate(coeffz32)
	deallocate(coeffx33)
	deallocate(coeffy33)
	deallocate(coeffz33)
	deallocate(coeffx220)
	deallocate(coeffy220)
	deallocate(coeffz220)

	! deallocate
	deallocate(tt)
	deallocate(green12x)
	deallocate(green12y)
	deallocate(green12z)
	deallocate(green32x)
	deallocate(green32y)
	deallocate(green32z)
	deallocate(green33x)
	deallocate(green33y)
	deallocate(green33z)
	deallocate(green220x)
	deallocate(green220y)
	deallocate(green220z)
	deallocate(subfsm)
	deallocate(depth)
	deallocate(epdis)
	deallocate(greenval)
	deallocate(fsm)
	deallocate(obs)

	write(*,*) '-------------------'
	write(*,*) 'END PROGRAM'
	write(*,*) '-------------------'
 
    call CPU_TIME(finish)
    write(*,*) 'it takes',finish-start,'s'

end program Integ



subroutine subfaults(ft, ny, nx, subft)
!
!   discretize subfaults into nx by ny superpositions
!
	use class
!$acc routine worker
	implicit none
	integer :: nf, ny, nx
	type(source) :: ft, subft(ny, nx)
	integer :: i,j
	real*8 :: latNW, lonNW, depthNW, lenpatch, widpatch, deppatch
	real*8 :: temp, temp1, temp2, x, y, x1, y1

    temp1 = dsin(ft%strike*deg2arc)
    temp2 = dcos(ft%strike*deg2arc)

	! coordinat of EW corner
    x1 = ft%length*0.5d0
    y1 = -ft%width*dcos(ft%dip*deg2arc)*0.5d0
    x = temp1*x1 - temp2*y1             ! fault cooridnate to E-N coordinate
    y = temp2*x1 + temp1*y1
    latNW = ft%lat - 2.d0*dasin(y/2.d0/Re)/deg2arc
    lonNW = ft%lon - 2.d0*dasin(x/2.d0/Re/dcos(latNW*deg2arc))/deg2arc  ! set NW as reference point
    depthNW = ft%depth - ft%width*dsin(ft%dip*deg2arc)*0.5d0


	! coordinates of subfaults
	lenpatch = ft%length/dble(nx)
	widpatch = ft%width*dcos(ft%dip*deg2arc)/dble(ny)
	deppatch = ft%width*dsin(ft%dip*deg2arc)/dble(ny)

	!$acc loop
    do i=1,nx
    	do j=1,ny
    	    x = temp1*(dble(i)-0.5d0)*lenpatch - temp2*(-(dble(j)-0.5d0)*widpatch)
            y = temp2*(dble(i)-0.5d0)*lenpatch + temp1*(-(dble(j)-0.5d0)*widpatch)
            subft(j,i)%lat = latNW + 2.d0*dasin(y/2.d0/Re)/deg2arc
            subft(j,i)%lon = lonNW + 2.d0*dasin(x/2.d0/Re/dcos(subft(j,i)%lat*deg2arc))/deg2arc
            subft(j,i)%depth = depthNW + (dble(j)-0.5d0)*deppatch
            subft(j,i)%slip = ft%slip
            subft(j,i)%rake = ft%rake
            subft(j,i)%strike = ft%strike
            subft(j,i)%dip = ft%dip
            subft(j,i)%length = lenpatch
            subft(j,i)%width = ft%width/dble(ny)
	    end do
    end do

end subroutine subfaults



subroutine displ(fault,obs,u12x,u12y,u12z,u32x,u32y,u32z,u33x,u33y,u33z,u220x,u220y,u220z,u)
	use class
!$acc routine seq
	implicit none
	type(source) fault
	type(station) obs
	real*8 :: u12x,u12y,u12z,u32x,u32y,u32z,u33x,u33y,u33z,u220x,u220y,u220z
	real*8 :: u(3)       ! output: ux, uy, uz
	real*8 :: v(3)
	integer :: k
	real*8 :: dist, zfault, zobs, cosz, z0
	real*8,dimension(3) :: u12, u32, u13, u33, u22
	real*8 :: rake, dip
!$acc routine(greatcircle,sineangle) seq

	! spherical angles
	call greatcircle(fault%lat, fault%lon, obs%lat, obs%lon, dist)
	dist = dacos(dist)

	call sineangle(fault%lat, fault%lon, obs%lat, obs%lon, dist, zfault, zobs)

	z0 = fault%strike*deg2arc-zfault   ! in arc

	! compute deformation for a single subfault patch
	u12(2) = dsin(2.d0*z0)*u12y
	u32(2) = dsin(z0)*u32y
	u13(2) = dsin(z0+pi/2.d0)*u32y
	u33(2) = u33y
	u12(3) = dsin(2.d0*z0)*u12z
	u32(3) = dsin(z0)*u32z
	u13(3) = dsin(z0+pi/2.d0)*u32z
	u33(3) = u33z
	u22(2) = u220y+dcos(2.d0*z0)*u12y
	u22(3) = u220z+dcos(2.d0*z0)*u12z
	u12(1) = dcos(2.d0*z0)*u12x
	u32(1) = dcos(z0)*u32x
	u13(1) = 0.d0
	u33(1) = 0.d0
	u22(1) = -dsin(2.d0*z0)*u12x


	rake = fault%rake*deg2arc
	dip = fault%dip*deg2arc
	!$acc loop independent
	do k=1,3
		v(k) = dcos(rake)*(u12(k)*dsin(dip)-u13(k)*dcos(dip)) &
		     + dsin(rake)*(0.5d0*(u33(k)-u22(k))*dsin(2*dip)-u32(k)*dcos(2*dip))
	end do

	! transform to the local N-E coordinate
	! rotate a countclockwise angle of pi/2-zobs
	v(2) = -v(2)                ! v(2) is opposite to North 
	u(1) =  v(1)*dcos(zobs) - v(2)*dsin(zobs)
	u(2) =  v(1)*dsin(zobs) + v(2)*dcos(zobs)
	u(3) =  v(3)

	!$acc loop independent
	do k=1,3
		u(k) = u(k) * fault%slip * fault%length * fault%width /Re/Re
	end do

end subroutine displ


subroutine greatcircle(lat1, lon1, lat2, lon2, dist)
! input:
! lat1, lon1 -- latitude and longitude of point 1
! lat2, lon2 -- latitude and longitude of point 1
! output:
! dist -- the great circle length between point 1 and 2, in arc
	use class
!$acc routine seq
	implicit none
	real*8 :: lat1, lon1, lat2, lon2, dist

	dist = dsin(lat1*deg2arc)*dsin(lat2*deg2arc) &
	      +dcos(lat1*deg2arc)*dcos(lat2*deg2arc) &
		  *dcos((lon2-lon1)*deg2arc)
	if(dist > 1.0) dist = 1.d0
	if(dist < -1.0) dist = -1.d0

end subroutine


subroutine sineangle(lat1, lon1, lat2, lon2, dist, z1, z2)
! input:
! lat1, lon1 -- latitude and longitude of point 1, in degree
! lat2, lon2 -- latitude and longitude of point 1, in degree
! dist -- the great circle length between point 1 and 2, in arc
!       dist can be calculated using subroutine greetcircle
! output:
! z1 -- sine of arzimuth angle of point 1, in arc
! z2 -- sine of arzimuth angle of point 2, in arc
! arzimuth angle is the angle between latitude and the direction
!   of dist
	use class
!$acc routine seq
	implicit none
	real*8 :: lat1, lon1, lat2, lon2, dist, z1, z2
	real*8 :: cosz

	z1 = dcos(lat2*deg2arc)*dsin((lon2-lon1)*deg2arc)/dsin(dist)
	if(z1 > 1.0) z1 = 1.d0
	if(z1 < -1.0) z1 = -1.d0
	z1 = dasin(z1)
	cosz = (dsin(lat2*deg2arc)-dsin(lat1*deg2arc)*dcos(dist)) &
	      /dcos(lat1*deg2arc)/dsin(dist)
	if(cosz < 0.0) z1 = PI-z1

	z2 = dcos(lat1*deg2arc)*dsin((lon2-lon1)*deg2arc)/dsin(dist)
	if(z2 > 1.0) z2 = 1.d0
	if(z2 < -1.0) z2 = -1.d0
	z2 = dasin(z2)
	cosz = (dsin(lat1*deg2arc)-dsin(lat2*deg2arc)*dcos(dist)) &
	      /dcos(lat2*deg2arc)/dsin(dist)
	if(cosz < 0.0) z2 = PI-z2

end subroutine

