subroutine intLin(KMAXX,kount,npoint,x_bef,y_bef,x_aft,y_aft)
!This subroutine calculate the kurtosis analytically

  implicit none

  real(16) x_bef(KMAXX),y_bef(KMAXX) !used for interpolation
  real(16) x_aft(npoint),y_aft(npoint) !used for interpolation
  real(16) x
  integer KMAXX,kount,npoint
  integer i,j,j_ini
  integer j1,j2

  j_ini=1

  do i=1, npoint
    x=x_aft(i)
    do j=j_ini, kount
      j1=j
      if(j1<kount)then
        j2=j1+1
        if(x_bef(j2)>x)then
          continue
        else
          y_aft(i)=(x-x_bef(j1))/(x_bef(j2)-x_bef(j1))*(y_bef(j2)-y_bef(j1))+y_bef(j1)
          exit
        end if
      else
        y_aft(i)=y_bef(j1)
        exit
      end if
    end do
    j_ini=j1
  end do

end





