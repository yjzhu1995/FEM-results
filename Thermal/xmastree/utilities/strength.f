	program strength
c********************************************************************
c  This program compares two stress profiles over a depth range and 
c  for each depth rejects the high stress value. The remaining low
c  stress values form the new, combined stresss profile 
c
c  Arguments at the commend line
c    fname1  -- name of the first stress profile file
c    fname2  -- name of the second stress profile file
c    fname3  -- name of output file
c
c  Note: The two input files must have identical depth values. But
c    the output file has additional data points at (brittle-ductile)
c    stress transition depths.
c
c********************************************************************
	character fname1*15,fname2*15,fname3*15
	call getarg(1,fname1)
	open(1,file=fname1)
	rewind 1
	call getarg(2,fname2)
	open(2,file=fname2)
	rewind 2
	call getarg(3,fname3)
	open(3,file=fname3)
	rewind 3
	dstr1=0
	bstr1=0
	k=0
10	k=k+1
	read(2,*,end=100) ddep2,dstr2
	read(1,*,end=100) bdep2,bstr2
        if(ddep2.ne.bdep2) then
	  print*, 'stop combd: different ddepth and bdepth at k=',k
	  stop
	end if
	if((dstr2-bstr2)*(dstr1-bstr1).lt.0.and.ddep2.ne.ddep1) then
	  dx=ddep2-ddep1
	  dd=dstr2-dstr1
	  db=bstr2-bstr1
	  dep=ddep1+dx*(bstr1-dstr1)/(dd-db)
	  str=bstr1+db/dx*(dep-ddep1)
	  write(3,*) dep,str
	  print*, 'stress transition at depth ',dep,' stress ',str
	end if
	write(3,*) ddep2,amin1(dstr2,bstr2)
	ddep1=ddep2
	dstr1=dstr2
	bstr1=bstr2
	goto 10
100	print*, fname1,' + ',fname2,'= ',fname3,'(',k-1,' points)'
	stop
	end
