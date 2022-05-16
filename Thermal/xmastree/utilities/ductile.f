	program ductile
c********************************************************************
c  to calculate "yield" stress for given temperature and strain rate
c  using the power-law (dislocation creep) and dorn law (dislocation 
c  glide)
c
c  Arguments at the commend line
c    fname -- name of file containing depth-temperature pairs
c    rate  -- strain rate
c 
c  Output file
c    d.1   -- depth-stress pairs
c********************************************************************
	character fname*15
	data r,ep,anp,qp,ed,sigd,qd/
     .    8.314,7.e-14,3.,523.e3,5.7e11,8.5e9,549.e3/
	call getarg(1,fname)
	open(1,file=fname)
	rewind 1
	call getarg(2,fname)
	read(fname,*) rate
	open(2,file='d.1')
	rewind 2
	k=0
10	k=k+1
	read(1,*,end=100) dep,tem
	tem=tem+273.16
	if(tem.lt.773.) then
	  power=201.
	else
	  power=(rate/ep*exp(qp/(r*tem)))**(1/anp)
	end if
	if(tem.gt.1173.) then
	  dorn=199
	else
	  dorn=sigd*(1.-sqrt(r*tem/qd*alog(ed/rate)))
	end if
	if(power.lt.200.e6) stress=power
	if(dorn.gt.200.e6) stress=dorn
	write(2,*) dep,1.e-6*stress
	goto 10
100	print*, 'ductile:',k-1,' points written into d.1'
	stop
	end
