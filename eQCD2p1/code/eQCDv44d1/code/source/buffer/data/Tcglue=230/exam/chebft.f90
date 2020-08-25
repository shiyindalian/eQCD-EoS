subroutine chebft(c, n, m)
!Chebyshev fit
  implicit none

  integer n, m, NMAX
!n is the number of fitted points of function£¬m is the order ofchebyshev polynomial
  real(8) c(m), PI
  parameter (NMAX=200, PI=3.141592653589793d0)
  integer j, k
  real(8) fac, f(NMAX)
  DOUBLE PRECISION sum
  common /prefit/ f

  fac=2./n
  do j=1, m
    sum=0.d0
    do k=1, n
      sum=sum+f(k)*cos((PI*(j-1))*((k-0.5d0)/n))
    end do
    c(j)=fac*sum
  end do
  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE chder(a,b,c,cder,n)

  implicit none

  INTEGER n
  real(8) a,b,c(n),cder(n)
  INTEGER j
  real(8) con

  cder(n)=0.
  cder(n-1)=2*(n-1)*c(n)
  if(n.ge.3)then
    do j=n-2,1,-1
      cder(j)=cder(j+2)+2*j*c(j+1)
    end do
  endif
  con=2./(b-a)
  do j=1,n
    cder(j)=cder(j)*con
  end do
  return
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION chebev(a,b,c,m,x)

  implicit none

  INTEGER m
  real(8) chebev,a,b,x,c(m)
  INTEGER j
  real(8) d,dd,sv,y,y2
  if ((x-a)*(x-b).gt.0.) pause 'x not in range in chebev'
  d=0.
  dd=0.
  y=(2.*x-a-b)/(b-a)
  y2=2.*y
  do j=m,2,-1
    sv=d
    d=y2*d-dd+c(j)
    dd=sv
  end do
  chebev=y*d-dd+0.5*c(1)
  return
END
