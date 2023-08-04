program picard
implicit none
integer, parameter :: ndim=2
double precision :: y(ndim), g(ndim), dt, t, tmax, y1(ndim), t1, g1(ndim)
integer :: it, itmax
!
tmax=20.d0; itmax=2000.d0; dt=tmax/itmax
y(1)=0.0001d0; y(2)=1.d0 ; t=0.0001d0 ! Initial Conditions

open(1,file='bessep=1.txt')

do it=1,itmax
call calcG(g,y,t,ndim)
y1(:)=y(:)+g(:)*dt
t1=t+dt
call calcG(g1,y1,t1,ndim)
y(:)=y(:)+0.5d0*(g1(:)+g(:))*dt
t=t+dt
write(1,*)t,y(1)
enddo
pause
end program
!-----------------------------------------------
subroutine calcG(g,y,t,ndim)
implicit none
integer:: ndim
double precision :: g(ndim), y(ndim), t,p
p=1.d0
!
g(1)=y(2)
g(2)=1/(t**2)*(-t*y(2)-(t**2 + p**2)*y(1))
!
return
end subroutine calcG
