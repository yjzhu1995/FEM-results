!*******************************************************
!
!                 CUBIC SPLINES
!
!  This code modified from Fu & Sun's FaultIntegration,
!  contains two parts:
!  1. spline, return the coefficients of bicubic spline
!  2. splint, interpolate a value z for a new point
!     (xnew,ynew).
!
!  One can call the subroutine spline first, then call
!  subroutine splint to obtain the interpolation value.
!
!  The code can be running on both CPU and GPU machines.
!
!  Programmed by Xin Zhou.
!
!*******************************************************

subroutine splint(n,m,x,y,a,xnew,ynew,z)
!$acc routine vector
  implicit none
  integer,intent(in) :: n,m
  real*8,intent(in) :: x(n),y(m),a(16,n,m),xnew,ynew
  real*8,intent(out) :: z
  integer :: i,j,k,xloc,yloc

  do i=1,n-1
    if (xnew<=x(i+1)) then
      xloc = i
      exit
    endif
  end do
      
  do j=1,m-1
    if (ynew<=y(j+1)) then
      yloc = j
      exit
    endif
  end do 

  z=0.d0
  k=0
  !$acc loop reduction(+:z)
  do i=1,4
    do j=1,4
      k=k+1
      z=z+a(k,xloc,yloc)*(xnew-x(xloc))**(i-1)*(ynew-y(yloc))**(j-1)
    end do
  end do

end subroutine splint

 
subroutine spline(n,m,x,y,u,a)
!$acc routine vector
  implicit none
  integer,intent(in) :: n,m
  real*8,intent(in) :: x(n),y(m),u(n,m)
  real*8,intent(out) :: a(16,n,m)
  integer :: n1,m1,i,j,k,kk,k1,k2,l1,l2
  real*8 :: p(n,m),q(n,m),r(n,m),b(4,4),c(4,4),d(4,4),e(4,4)
  real*8 :: hx,hy,sum
!$acc routine(mat) seq
!$acc routine(sp1x, sp1y) vector

  n1 = n-1
  m1 = m-1

!  call bounday(n,m,x,y,u,p,q)

  !$acc loop collapse(2) independent
  do i=1,2
    do j=1,4
      b(i,j) = 0.d0
      e(i,j) = 0.d0
    end do
    b(i,i) = 1.d0
    e(i,i) = 1.d0
  end do

  call sp1x(n,m,1,x,u,p)
  call sp1y(n,m,y,u,q)
!  call sp1x(n,m,m1,x,u,p)
  call sp1y(n,m,y,p,r)

  do i=1,n1
    hx=1.d0/(x(i+1)-x(i))
    call mat(hx,b)

    do j=1,m1
      !$acc loop
      do k1=1,3,2
        l1=(k1-1)/2+i
        do k2=1,3,2
          l2=(k2-1)/2+j
          c(k1,k2)=u(l1,l2)
          c(k1,k2+1)=q(l1,l2)
          c(k1+1,k2)=p(l1,l2)
          c(k1+1,k2+1)=r(l1,l2)
        end do
      end do

      do k1=1,4
        do k2=1,4
          sum=0.d0
          !$acc loop reduction(+:sum)
          do k=1,4
            sum=sum+b(k1,k)*c(k,k2)
          end do
          d(k1,k2)=sum
        end do
      end do

      hy=1.d0/(y(j+1)-y(j))  
      call mat(hy,e)

      kk=0
      do k1=1,4
        do k2=1,4
          sum=0.d0
          !$acc loop reduction(+:sum)
          do k=1,4
            sum=sum+d(k1,k)*e(k2,k)
          end do
!          kk=4*(k1-1)+k2
          kk=kk+1
          a(kk,i,j)=sum
        end do
      end do

    end do
  end do

end subroutine spline



subroutine sp1x(n,m,m1,x,u,p)
!$acc routine vector
  implicit none
  integer,intent(in) :: n,m,m1
  real*8,intent(in) :: x(n),u(n,m)
  real*8 :: p(n,m),a(n-1),b(n-1)
  integer :: i,j,n1
  real*8 :: h,hi,af,bt

  !$acc loop independent
  do i=1,m
    p(1,i) = (u(2,i)-u(1,i))/(x(2)-x(1))
    p(n,i) = (u(n,i)-u(n-1,i))/(x(n)-x(n-1))
  end do

  n1=n-1

  !$acc loop
  do j=1,m,m1
    a(1)=0.d0
    b(1)=p(1,j)
    do i=2,n1
      hi=x(i)-x(i-1)
      h=x(i+1)-x(i)
      af=hi/(hi+h)
      bt=3.d0*((1.d0-af)*(u(i,j)-u(i-1,j))/hi+af*(u(i+1,j)-u(i,j))/h)
      a(i)=-af/(2.d0+(1.d0-af)*a(i-1))
      b(i)=(bt-(1.d0-af)*b(i-1))/(2.d0+(1.d0-af)*a(i-1))
    end do

    !$acc loop 
    do i=n1,2,-1
      p(i,j)=a(i)*p(i+1,j)+b(i)
    end do
  end do

end subroutine sp1x


subroutine sp1y(n,m,y,u,q) 
!$acc routine vector
  implicit none
  integer,intent(in) :: m,n
  real*8,intent(in) :: y(m),u(n,m)
  real*8 :: q(n,m),a(m-1),b(m-1)
  integer :: m1,i,j
  real*8 :: h,hj,af,bt

  !$acc loop independent
  do i=1,n
    q(i,1) = (u(i,2)-u(i,1))/(y(2)-y(1))
    q(i,m) = (u(i,m)-u(i,m-1))/(y(m)-y(m-1))
  end do

  m1=m-1

  !$acc loop
  do i=1,n
    a(1)=0.d0
    b(1)=q(i,1)
    do j=2,m1
      hj=y(j)-y(j-1)
      h=y(j+1)-y(j)
      af=hj/(hj+h)
      bt=3.d0*((1.d0-af)*(u(i,j)-u(i,j-1))/hj+af*(u(i,j+1)-u(i,j))/h)
      a(j)=-af/(2.d0+(1.d0-af)*a(j-1))
      b(j)=(bt-(1.d0-af)*b(j-1))/(2.d0+(1.d0-af)*a(j-1))
    end do

    !$acc loop 
    do j=m1,2,-1
      q(i,j) = a(j)*q(i,j+1)+b(j)
    end do
  end do

end subroutine sp1y




subroutine mat(h,b)
!$acc routine seq
  implicit none  
  real*8 :: h,b(4,4)
  real*8 :: h1
  
  h1=h*h
  b(3,1)=-3.d0*h1
  b(3,2)=-2.d0*h
  b(3,3)=-b(3,1)
  b(3,4)=-h
  b(4,1)=2.d0*h*h1
  b(4,2)=h1
  b(4,3)=-b(4,1)
  b(4,4)=h1

end subroutine mat


!subroutine bounday(n,m,x,y,u,p,q)
!!$acc routine vector
!  implicit none
!  integer,intent(in) :: n,m
!  real*8,intent(in) :: x(n),y(m),u(n,m)
!  real*8,intent(out) :: p(n,m),q(n,m)
!  integer :: i
!
!  !$acc parallel present(x,u)
!  !$acc loop independet
!  do i=1,m
!    p(1,i) = (u(2,i)-u(1,i))/(x(2)-x(1))
!    p(n,i) = (u(n,i)-u(n-1,i))/(x(n)-x(n-1))
!  end do
!
!  !$acc loop independet
!  do i=1,n
!    q(i,1) = (u(i,2)-u(i,1))/(y(2)-y(1))
!    q(i,m) = (u(i,m)-u(i,m-1))/(y(m)-y(m-1))
!  end do
!  !$acc end parallel
!
!end subroutine bounday
