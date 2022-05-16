	program brittle
c****************************************************************************
c to calculate "yield" stress at depth points given in the t-z file
c using the Bylerlee's law
c
c Arguments
c   t-z file, pore fluid-pressure ratio, depth of slabface, slabdip, max.thickness
c
c output files
c   b.1 for normal faulting regime, b.2 for thrust faulting
c****************************************************************************
	character fname*15
c fc1 and fc2 are friction coefficients of the Bylerlee's law
	data rho,g,fc1,fc2/3.33,9.8,0.85,0.6/
	call getarg(1,fname)
	open(1,file=fname)
	rewind 1
	call getarg(2,fname)
	read(fname,*) pore
	call getarg(3,fname)
	read(fname,*) deps
	call getarg(4,fname)
	read(fname,*) dip
          dip=dip/180*3.1415926
	call getarg(5,fname)
	read(fname,*) h
	open(2,file='b.1')
	rewind 2
	open(3,file='b.2')
	rewind 3
	fai1=atan(fc1)
	fai2=atan(fc2)
c sigma1 and sigma2 are values of sigma at two breaking pionts in normal
c faulting and thrust faulting, respectively.   
	sigma1=200./(1-sin(fai1))
	sigma2=200./(1+sin(fai2))
c dep_bp1 and dep_bp2 are values of depths of two breaking pionts.
	dep_bp1=(sigma1/(rho*g*(1-pore))-deps)/cos(dip)
	dep_bp2=(sigma2/(rho*g*(1-pore))-deps)/cos(dip)
c dif_str1 and dif_str2 are invariable component of differential stresses
c at dep_bp1 and dep_bp2, respectively.
	dif_str1=(coef1(fai1)-coef1(fai2))*sigma1
	dif_str2=(coef2(fai1)-coef2(fai2))*sigma2
	k=0
10	k=k+1
	read(1,*,end=100) dep,tem
	sigma=rho*g*(deps+dep*cos(dip))*(1-pore)
	if(sigma.lt.sigma2) then
	  write(2,*) dep,coef1(fai1)*sigma
	  write(3,*) dep,coef2(fai1)*sigma
c calculate break point for thrust faulting
	else if(sigma.gt.sigma2.and.sigma0.lt.sigma2) then
	  if(dep_bp2.le.h) then
	    write(2,*) dep_bp2,coef1(fai1)*rho*g*(dep_bp2*cos(dip)+deps)*
     .                 (1-pore)
	    write(3,*) dep_bp2,coef2(fai1)*sigma2
	    k=k+1
	  end if
	else if(sigma.gt.sigma2.and.sigma.lt.sigma1) then
	  write(2,*) dep,coef1(fai1)*sigma
	  write(3,*) dep,dif_str2+coef2(fai2)*sigma
c calculate break point for normal faulting
	else if(sigma.gt.sigma1.and.sigma0.lt.sigma1) then
	  if(dep_bp1.le.h) then
	    write(2,*) dep_bp1,coef1(fai1)*sigma1
	    write(3,*) dep_bp1,dif_str2+coef2(fai2)*rho*g*
     .                           (dep_bp1*cos(dip)+deps)*(1-pore)
	    k=k+1
	  end if
	else if(sigma.gt.sigma1) then
	  write(2,*) dep,dif_str1+coef1(fai2)*sigma
	  write(3,*) dep,dif_str2+coef2(fai2)*sigma
	end if
	sigma0=sigma
	goto 10
100	print*, 'brittle',k-1,' points written into b.1 and b.2'
	stop
	end

c coef1: times of differential stress divided by maximum principal stress
c coef2: times of differential stress divided by minimum principal stress

	function coef1(x)
	coef1=2*sin(x)/(1+sin(x))
	return
	end
	function coef2(x)
	coef2=2*sin(x)/(1-sin(x))
	return
	end
