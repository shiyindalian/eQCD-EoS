SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
!Runge-Kutta driver with adaptive stepsize control. Integrate the starting values ystart(1:nvar)
!from x1 to x2 with accuracy eps, storing intermediate results in the common block /path/.
!h1 should be set as a guessed first stepsize, hmin as the minimum allowed stepsize (can be zero).
!On output nok and nbad are number of good and bad (but retried and fixed) steps taken, and ystart is
!replaced by values at the end of the integration interval. derivs is the user-supplied subroutine for
!calculating the right-hand side derivatives, while rkqs is the name of the stepper routine to be used.
!/path/ contains its own information about how often an intermediate value is to be stored.
  implicit none

  INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
  real(16) eps,h1,hmin,x1,x2,ystart(nvar),TINY
  EXTERNAL derivs,rkqs
  PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=2000,TINY=1.Q-30)
  INTEGER i,kmax,kount,nstp
  real(16) dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),yp(NMAX,KMAXX),yscal(NMAX)
  COMMON /path/ kmax,kount,dxsav,xp,yp
!User storage for intermediate results. Preset dxsav and kmax
  x=x1
  h=sign(h1,x2-x1)
  nok=0
  nbad=0
  kount=0
  do i=1,nvar
    y(i)=ystart(i)
  end do
  if (kmax.gt.0) xsav=x-2.Q+0*dxsav
  do nstp=1,MAXSTP
    call derivs(x,y,dydx)
    do i=1,nvar
      yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
    end do
    if(kmax.gt.0)then
      if(abs(x-xsav).gt.abs(dxsav)) then
        if(kount.lt.kmax-1)then
          kount=kount+1
          xp(kount)=x
          do i=1,nvar
            yp(i,kount)=y(i)
          end do
          xsav=x
        endif
      endif
    endif
    if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
    call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
    if(hdid.eq.h)then
      nok=nok+1
    else
      nbad=nbad+1
    endif
    if((x-x2)*(x2-x1).ge.0.)then
      do i=1,nvar
        ystart(i)=y(i)
      end do
      if(kmax.ne.0)then
        kount=kount+1
        xp(kount)=x
        do i=1,nvar
          yp(i,kount)=y(i)
        end do
      endif
      return
    endif
    if(abs(hnext).lt.hmin) pause 'stepsize smaller than minimum in odeint'
    h=hnext
  end do
  pause 'too many steps in odeint'
  return
END

! (C) Copr. 1986-92 Numerical Recipes Software ,4-#.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
!Given values for n variables y and their derivatives dydx known at x, use the
!fifth-order Cash-Karp Runge-Kutta method to advance the solution over an interval h
!and return the incremented variable as yout. Also return an estimate of the local truncations
!error in yerr using the embeded fourth-order method. The user supplied the subroutine
!derivs(x,y,dydx), which returns derivatives dydx at x.

  implicit none

  INTEGER n,NMAX
  real(16) h,x,dydx(n),y(n),yerr(n),yout(n)
  EXTERNAL derivs
  PARAMETER (NMAX=50)
!U    USES derivs
  INTEGER i
  real(16) ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX), &
       ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53, &
       B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
  PARAMETER (A2=.2Q+0,A3=.3Q+0,A4=.6Q+0,A5=1.Q+0,A6=.875Q+0,B21=.2Q+0,B31=3.Q+0/40.Q+0, &
             B32=9.Q+0/40.Q+0,B41=.3Q+0,B42=-.9Q+0,B43=1.2Q+0,B51=-11.Q+0/54.Q+0,B52=2.5Q+0, &
             B53=-70.Q+0/27.Q+0,B54=35.Q+0/27.Q+0,B61=1631.Q+0/55296.Q+0,B62=175.Q+0/512.Q+0, &
             B63=575.Q+0/13824.Q+0,B64=44275.Q+0/110592.Q+0,B65=253.Q+0/4096.Q+0,C1=37.Q+0/378.Q+0, &
             C3=250.Q+0/621.Q+0,C4=125.Q+0/594.Q+0,C6=512.Q+0/1771.Q+0,DC1=C1-2825.Q+0/27648.Q+0, &
             DC3=C3-18575.Q+0/48384.Q+0,DC4=C4-13525.Q+0/55296.Q+0,DC5=-277.Q+0/14336.Q+0, &
             DC6=C6-.25Q+0)
  do i=1,n
    ytemp(i)=y(i)+B21*h*dydx(i)
  end do
  call derivs(x+A2*h,ytemp,ak2)
  do i=1,n
    ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
  end do
  call derivs(x+A3*h,ytemp,ak3)
  do i=1,n
    ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
  end do
  call derivs(x+A4*h,ytemp,ak4)
  do i=1,n
    ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
  end do
  call derivs(x+A5*h,ytemp,ak5)
  do i=1,n
    ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
  end do
  call derivs(x+A6*h,ytemp,ak6)
  do i=1,n
    yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
  end do
  do i=1,n
    yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
  end do
  return
END
!  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
!Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy
!and adjust stepsize. Input are the dependent variable vector y(1:n) and its derivative dydx(1:n)
!at the starting value of the independent variable x, Also input are the stepsize to be attempted
!htry, the required accuracy eps, and the vector yscal(1:n) against which the error is scaled. On
!output, y and x are replaced by their new values, hdid is the stepsize that was actually accomplished,
!and hnext is the estimated next stepsize. derivs is the user-supplied subroutine that computes the
!right-hand side derivatives.
  implicit none

  INTEGER n,NMAX
  real(16) eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
  EXTERNAL derivs
  PARAMETER (NMAX=50)
!    USES derivs,rkck
  INTEGER i
  real(16) errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,ERRCON
  PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
  h=htry
1 call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
  errmax=0.
  do i=1,n
    errmax=max(errmax,abs(yerr(i)/yscal(i)))
  end do
  errmax=errmax/eps
  if(errmax.gt.1.)then
    htemp=SAFETY*h*(errmax**PSHRNK)
    h=sign(max(abs(htemp),0.1*abs(h)),h)
    xnew=x+h
    if(xnew.eq.x)pause 'stepsize underflow in rkqs'
    goto 1
  else
    if(errmax.gt.ERRCON)then
      hnext=SAFETY*h*(errmax**PGROW)
    else
      hnext=5.*h
    endif
    hdid=h
    x=x+h
    do i=1,n
      y(i)=ytemp(i)
    end do
    return
  endif
END
!  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.

