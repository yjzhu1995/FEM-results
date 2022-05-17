	program cooling
c*********************************************************************
c  program to calculate temperature and heat flow using the plate
c  model with parameters reported by Stein and Stein (1992)
c  or Parson and Sclater parameters (1977)
c
c  Input arguments at commend line
c    age   -- age of plate in years
c    rate  -- half spreading rate in m/yr 
c             (results are very insensitive to this parameter
c              as explained by P&S (1977) in page 806)
c    depth -- maximum depth of t-z profile in kilometers
c    nz    -- number of depth intervals (resulting in nz+1 t-z pairs)
c    model -- model flag. 
c             abs(model) is for parameter set: GDH1 = 1; P&S = 2
c             sign is or model type: plate > 0; halfspace < 0 
c
c  Output
c    t-z pairs in file temp.dat
c    heat flux value printed on screen
c*********************************************************************
	implicit real*8(a-h,o-z)
	character argname*10
	dimension para(2,2)
	data rho,c,cond/3330.,1171.,3.318/
c  the first two in para are GDH1 thinkness and basal temperature
c  the last two in para are P&S thinkness and basal temperature
	data para/95.,1450.,
     .           125.,1350./
	open(2,file='temp.dat')
	year=365.25*24*3600
	call getarg(1,argname)
	read(argname,*) age
	age=age*year
	call getarg(2,argname)
	read(argname,*) rate
	rate=rate/year
	call getarg(3,argname)
	read(argname,*) depth
	call getarg(4,argname)
	read(argname,*) nz
	call getarg(5,argname)
	read(argname,*) mod
	model=abs(mod)
	h=para(1,model)
	tm=para(2,model)
	if(mod.gt.0) call gdh1(age,rate,depth,nz,h,tm)
	if(mod.lt.0) call half(age,depth,nz,tm)
	stop
	end
c
c
	subroutine gdh1(age,rate,depth,nz,h,tm)
	implicit real*8(a-h,o-z)
	dimension t(0:10000)
	data rho,c,cond/3330.,1171.,3.318/
	pi=3.14159265359
	diff=cond/(rho*c)
	dz=depth/nz
	x=rate*age/1000.
	r=rate*h*1000./(2.*diff)
	q=1.
	  do i=0,nz
	  t(i)=i*dz/h
	  end do
	do 200 n=1,50
	pin=n*pi
	bn=sqrt(r*r+pin*pin)-r
	expo=exp(-bn*x/h)
	q=q+2*expo
	cn=2./(pin)
	print*, 'n=',n,' expo=',expo
	  do i=0,nz
	  t(i)=t(i)+cn*expo*sin(pin*i*dz/h)
	  if(expo.lt.1.e-7) then
	    write(2,*) i*dz,tm*t(i)
	    if(i.eq.nz) goto 300
	  end if
	  end do
200	continue
300	q=q*cond*tm/h
	print*, 'heat flux =',q,' mW/m2'
	return
	end
c
c
	subroutine half(age,depth,nz,tm)
	implicit real*8(a-h,o-z)
	data rho,c,cond/3330.,1171.,3.318/
	pi=3.14159265359
	diff=cond/(rho*c)
	dz=depth/nz
	da=2.*sqrt(diff*age)
	do i=0,nz
	  zz=i*dz
          temp=tm*erf(1000.*zz/da)
	  write(2,*) zz,temp
        end do
	q=2.*cond*tm/da/sqrt(pi)*1000.
	print*, 'heat flux =',q,' mW/m2'
	return
	end
	
