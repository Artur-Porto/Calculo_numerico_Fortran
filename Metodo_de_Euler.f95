program euler
implicit none
integer, parameter :: ndim=2
double precision :: y(ndim), g(ndim), dt, t, tmax
integer :: it, itmax
!
tmax=20.; itmax=2000; dt=tmax/itmax
y(1)=1.0d0; y(2)=0.0d0 ; t=0.0d0 ! Initial Conditions
!
open(1,file='pontos.txt')

!
do it=1,itmax
call calcG(g,y,t,ndim)
y(:)=y(:)+g(:)*dt
t=t+dt
write(1,*)t,y(1)
enddo
pause
end program
!------------------------------------------------------------------
subroutine calcG(g,y,t,ndim)
implicit none
integer:: ndim
double precision :: g(ndim), y(ndim), t,k,m
k=1.d0
m=1.d0
!
g(1)=y(2)
g(2)=-y(1)
!
return
end subroutine calcG
