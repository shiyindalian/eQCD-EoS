program main
!This small program is used to fix two bugs in the Fortran codes generated from
!mathematica: replacing "\" with " ", adding "&" at the end of every line.
!change "," to next line

  implicit none
  character*80 buffer
  character compa
  integer sta,i

  open(unit=101,file='dhamilton.f90')
  open(unit=102,file='mass.f90')

  do while(.true.)

    read(101,"(A80)",iostat=sta)buffer
    if(sta/=0)exit

    do i=1, 80
      compa=buffer(i:i)
      if(compa=="\")then
        buffer(i:i)=" "
      end if
      if(compa==",")then
        buffer(i:i)=char(10)
      end if
    end do

    write(102,"(A80, A1)")buffer,"&"

  end do

end
