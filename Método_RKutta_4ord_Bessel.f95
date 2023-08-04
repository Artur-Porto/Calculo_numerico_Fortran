program runge_kutta
implicit none
integer, parameter :: ndim = 2
double precision :: y(ndim), g(ndim), dt, t, tmax
integer :: it, itmax

tmax=20.; itmax=200; dt=tmax/dfloat(itmax)
!
! condicoes iniciais
y(1)=1.0d0; y(2)=0.0001d0 ; t=0.0001d0 ! Initial Conditions
!
open(1,file='bessel.kutta 4.txt')
do it=2,itmax
call RK4(y,t,dt,ndim)
!print*, t, y(1), y(2)
write(1,*)t,y(1)
enddo
pause
end program
!##################################################################################################################
subroutine RK4(y,t,dt,ndim)
implicit none
integer:: ndim
double precision :: g(ndim),g1(ndim),g2(ndim),g3(ndim),g4(ndim)
double precision :: y(ndim),y1(ndim),y2(ndim),y3(ndim),y4(ndim)
double precision :: c1(ndim),c2(ndim),c3(ndim),c4(ndim)
double precision :: t1,t,t2,dt,t3,t4

!c1
t1 = t ; y1(:) = y(:)
call calcG(g1,y1,t1,ndim)
c1(:)=dt*g1(:)
!c2
t2=t+(0.5d0)*dt; y2(:)=y(:)+(0.5d0)*c1(:)
call calcG(g2,y2,t2,ndim)
c2(:)=dt*g2(:)
!c3
t3 = t+dt/2.0d0 ; y3(:) = y(:) + c2(:)/2.0d0
call calcG(g3,y3,t3,ndim)
c3(:) = dt*g3(:)
!c4
t4 = t+dt ; y4(:) = y(:) + c3(:)
call calcG(g4,y4,t4,ndim)
c4(:)= dt*g4(:)
!Y
y(:)=y(:)+(c1(:)+2.0d0*c2(:)+2.0d0*c3(:)+c4(:))/6.0d0
t=t+dt
return
end subroutine RK4
!#############################################################################################################
subroutine calcG(g,y,t,ndim)
implicit none
integer:: ndim,p
double precision :: g(ndim), y(ndim), t
!
p=0
!
g(1)=y(2)
g(2)=(-t*y(2) - (t**2 + p**2)*y(1))/t**2
!
return
end subroutine calcG
