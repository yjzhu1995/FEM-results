	program tstrength
c********************************************************************
c  This program combines two stress profiles over a depth range and 
c  for each depth generates a new stress value with an equation 
c  developed by Shimamoto and Noda (2014,JGR). The new generated
c  stress values form the total stresss profile with a smooth brittle-
c  ductile transition.
c
c  Arguments at the commend line
c    fname1  -- name of the first stress profile file
c    fname2  -- name of the second stress profile file
c    fname3  -- name of output file
c
c  Note: The two input files must have identical depth values. 
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
	k=0
10	k=k+1
	read(1,*,end=100) bdep,bstr
	read(2,*,end=100) ddep,dstr
        if(ddep.ne.bdep) then
	  print*, 'stop combd: different ddepth and bdepth at k=',k
	  stop
	end if
        if(dstr.gt.1.e-10) then
	  tbd_tmp=bstr/dstr
        else
          tbd_tmp=31.
        end if
	if(tbd_tmp.le.1.e-10) then
	  tbd_h=tbd_tmp
	else if (tbd_tmp.ge.30.) then
	  tbd_h=1.0
	else
	  tbd_h=(exp(tbd_tmp)-exp(tbd_tmp*(-1.d0)))/(exp(tbd_tmp)
     .         +exp(tbd_tmp*(-1.d0)))	
	end if
	tbdstr=dstr*tbd_h
	write(3,*) ddep,tbdstr
	goto 10
100	print*, fname1,' + ',fname2,'= ',fname3,'(',k-1,' points)'
	stop
	end
