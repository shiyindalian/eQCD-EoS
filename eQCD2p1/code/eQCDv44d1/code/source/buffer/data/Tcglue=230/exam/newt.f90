subroutine newt(x, n, check, functionn)

  implicit none

  integer n, nn, NP, MAXITS
  logical check
  real(8) x(n), fvec, TOLF, TOLMIN, TOLX, STPMX
  parameter (NP=40, MAXITS=200, TOLF=1.e-4, TOLMIN=1.e-6, TOLX=1.e-7, STPMX=100.)

  common /newtv/ fvec(NP), nn       !Communicates with fmin
  save /newtv/
!Uses fdjac, fmin, lnsrch, lubksb, ludcmp
  integer i, its, j, indx(NP)
  real(8) d, den, f, fold, stpmax, sum, temp, test, fjac(NP,NP), g(NP), p(NP), xold(NP) !, fmin
  external ffmin, functionn



  nn=n
!  f=fmin(x)                         !The vector fvec is also computed by this call
  call ffmin(x, f, functionn)
  test=0.                           !Test for initial guess being a root. Use more stringent
  do i=1, n                         !test than simply TOLF
    if(abs(fvec(i)).gt.test)test=abs(fvec(i))
  end do

  if(test.lt..01*TOLF) then
    check=.false.
	return
  end if

  sum=0.                            !Calculate stpmax for line searches
  do i=1, n
    sum=sum+x(i)**2
  end do
  stpmax=STPMX*max(sqrt(sum),float(n))

  do its=1, MAXITS                  !Start of iteration loop
    call fdjac(n, x, fvec, NP, fjac, functionn)
!If analytic Jacobian is available, you can replace the routine fdjac with your own routine.
    do i=1, n                       !compute Grad f for the line search
	  sum=0.
	  do j=1, n
	    sum=sum+fjac(j,i)*fvec(j)
	  end do
	  g(i)=sum
	end do

	do i=1, n                        !Store x,
	  xold(i)=x(i)
	end do

	fold=f                           !and f
	do i=1, n                        !Right hand side for linear equations
	  p(i)=-fvec(i)
	end do

	call ludcmp(fjac, n, NP, indx, d)  !Solve linear equations by LU decomposition
    call lubksb(fjac, n, NP, indx, p)
    call lnsrch(n, xold, fold, g, p, x, f, stpmax, check, ffmin, functionn)
!lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin

    test=0.                            !Test for convergence on function values
	do i=1, n
	  if(abs(fvec(i)).gt.test)test=abs(fvec(i))
	end do

	if(test.lt.TOLF) then
	  check=.false.
	  return
	end if

	if(check) then                     !Check for gradient of f zero, i.e., spurious convergence
	  test=0.
	  den=max(f,.5*n)
	  do i=1, n
	    temp=abs(g(i))*max(abs(x(i)),1.)/den
		if(temp.gt.test)test=temp
	  end do

	  if(test.lt.TOLMIN)then
	    check=.true.
	  else
	    check=.false.
	  end if
	  return
	endif

	test=0.                            !Test for convergence on delt-x
	do i=1, n
	  temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
	  if(temp.gt.test)test=temp
	end do
	if(test.lt.TOLX)return

  end do

  pause 'MAXITS exceeded in newt'

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ffmin(x, f, functionn)

  implicit none

  integer n, NP
  real(8) f, x(*), fvec
  parameter (NP=40)
  common /newtv/ fvec(NP), n
  save /newtv/
!Uses funcv
  integer i
  real(8) sum
  external functionn

!  call funcv(n, x, fvec)
  call functionn(n, x, fvec)
  sum=0.
  do i=1, n
    sum=sum+fvec(i)**2
  end do

  f=0.5*sum

  return

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fdjac(n, x, fvec, np, df, functionn)

  implicit none

  integer n, np, NMAX
  real(8) df(np,np), fvec(n), x(n), EPS
  parameter (NMAX=40,EPS=1.e-4)
!Uses funcv
  integer i, j
  real(8) h, temp, f(NMAX)
  external functionn

  do j=1, n
    temp=x(j)
	h=EPS*abs(temp)
	if(h.eq.0.)h=EPS
	x(j)=temp+h                       !Trick to reduce finite precision error
	h=x(j)-temp
!	call funcv(n, x, f)
    call functionn(n, x, f)
	x(j)=temp
	do i=1, n                         !Forward difference formula
	  df(i,j)=(f(i)-fvec(i))/h
	end do
  end do


  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine lnsrch(n, xold, fold, g, p, x, f, stpmax, check, func, functionn)

  implicit none

  integer n
  logical check
  real(8) f, fold, stpmax, g(n), p(n), x(n), xold(n), ALF, TOLX
  parameter (ALF=1.e-4,TOLX=1.e-7)
  external func, functionn
!Uses func
  integer i
  real(8) a, alam, alam2, alamin, b, disc, f2, fold2, rhs1, rhs2, slope, sum, temp, &
                 test, tmplam


  check=.false.

  sum=0.
  do i=1, n
    sum=sum+p(i)*p(i)
  end do
  sum=sqrt(sum)

  if(sum.gt.stpmax) then
    do i=1, n
	  p(i)=p(i)*stpmax/sum
	end do
  endif

  slope=0.
  do i=1, n
    slope=slope+g(i)*p(i)
  end do

!  if(slope >= 0.0) pause 'roundoff problem in lnsrch'

  test=0.
  do i=1, n
    temp=abs(p(i))/max(abs(xold(i)),1.)
	if(temp.gt.test)test=temp
  end do

  alamin=TOLX/test
  alam=1.

1 continue
    do i=1, n
	  x(i)=xold(i)+alam*p(i)
	end do

!	f=func(x)
    call func(x, f, functionn)
	if(alam.lt.alamin) then
	  do i=1, n
	    x(i)=xold(i)
	  end do
	  check=.true.
	  return
	else if(f.le.fold+ALF*alam*slope) then
	  return
	else
	  if(alam.eq.1.) then
	    tmplam=-slope/(2.*(f-fold-slope))
	  else
	    rhs1=f-fold-alam*slope
		rhs2=f2-fold2-alam2*slope
		a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
		b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
		if(a.eq.0.) then
		  tmplam=-slope/(2.*b)
		else
		  disc=b*b-3.*a*slope
		  if(disc.lt.0.) pause 'roundoff problem in lnsrch'
		  tmplam=(-b+sqrt(disc))/(3.*a)
		end if

		if(tmplam.gt..5*alam) tmplam=.5*alam

	  endif
	endif
	alam2=alam
	f2=f
	fold2=fold
	alam=max(tmplam,.1*alam)

  goto 1

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lubksb(a, n, np, indx, b)  !forward substitution and backsubstitution

  implicit none

  integer n, np, indx(n)
  real(8) a(np,np), b(n)
  integer i, ii, j, ll
  real(8) sum



  ii=0
  do i=1, n

    ll=indx(i)
	sum=b(ll)
	b(ll)=b(i)

	if (ii.ne.0) then

	  do j=ii, i-1
	    sum=sum-a(i,j)*b(j)
	  end do

	else if (sum.ne.0.) then

	  ii=i

	end if

	b(i)=sum

  end do

  do i=n, 1, -1

    sum=b(i)
	do j=i+1, n
	  sum=sum-a(i,j)*b(j)
	end do

	b(i)=sum/a(i,i)

  end do

  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ludcmp(a, n, np, indx, d)  !¾ØÕóµÄLU·Ö½â

  implicit none

  integer n, np, indx(n), NMAX
  real(8) d, a(np,np), TINY
  parameter (NMAX=10,TINY=1.0e-20)
  integer i, imax, j, k
  real(8) aamax, dum, sum, vv(NMAX)



  d=1.0
  do i=1, n

    aamax=0.0
	do j=1, n
	  if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
    end do
    if (aamax.eq.0.) pause 'singular matrix in ludcmp'
	vv(i)=1.0/aamax

  end do

  do j=1, n

    do i=1, j-1

	  sum=a(i,j)
	  do k=1, i-1

	    sum=sum-a(i,k)*a(k,j)

	  end do
	  a(i,j)=sum

	end do

	aamax=0.0
	do i=j, n

	  sum=a(i,j)
	  do k=1, j-1

	    sum=sum-a(i,k)*a(k,j)

	  end do
	  a(i,j)=sum

	  dum=vv(i)*abs(sum)
	  if (dum.ge.aamax) then
	    imax=i
		aamax=dum
	  end if

	end do

	if (j.ne.imax) then

	  do k=1, n
	    dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
	  end do

	  d=-d
	  vv(imax)=vv(j)

	endif
	indx(j)=imax

	if(a(j,j).eq.0.) a(j,j)=TINY

	if(j.ne.n)then

	  dum=1./a(j,j)
	  do i=j+1,n
	    a(i,j)=a(i,j)*dum
	  end do

	endif

  end do

  return
end


