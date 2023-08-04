program runge_kutta
implicit none
integer, parameter :: ndim = 2
double precision :: y(ndim), g(ndim), dt, t, tmax
integer :: it, itmax

tmax=20.; itmax=200; dt=tmax/dfloat(itmax)
!
! condicoes iniciais
y(1)=1.0d0; y(2)=0.0d0 ; t=0.0d0 ! Initial Conditions
!
open(1,file='osc.har.kutta.txt')
do it=2,itmax
call RK2(y,t,dt,ndim)
!print*, t, y(1), y(2)
write(1,*)t,y(1)
enddo
pause
end program
!##################################################################################################################
subroutine RK2(y,t,dt,ndim)
implicit none
integer:: ndim
double precision :: g1(ndim),g2(ndim), y(ndim), y1(ndim), y2(ndim), c1(ndim), c2(ndim)
double precision :: t1,t,t2,dt

real*8,parameter:: v21 = 3.d0/4.d0
real*8,parameter:: a1 = 1.d0/3.d0
real*8,parameter:: a2 = 2.d0/3.d0

t1 = t
y1(:) = y(:)
call calcG(g1,y1,t1,ndim)
c1(:)=dt*g1(:)
!
t2=t+(v21)*dt; y2(:)=y(:)+(v21)*c1(:)
call calcG(g2,y2,t2,ndim)
c2(:)=dt*g2(:)
!
y(:)=y(:)+(a1)*c1(:)+(a2)*c2(:)
t=t+dt
return
end subroutine RK2
!#############################################################################################################
subroutine calcG(g,y,t,ndim)
implicit none
integer:: ndim
double precision :: g(ndim), y(ndim), t
!
g(1)=y(2)
g(2)=-y(1)
!
return
end subroutine calcG
