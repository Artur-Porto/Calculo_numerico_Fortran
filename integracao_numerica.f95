program int_num
implicit none

!call trapezio_composto
!call simpson_composto
call tres_oitavos_simpson

pause
end program
!##############################################################################
double precision function f(x)
implicit none
real*8::x,L
real*8,parameter:: pi=acos(-1.)
integer :: n

if (x.lt.0.d0) f=0.d0
if (x.eq.0.d0) f=0.5d0
if (x.gt.0.d0) f=1.d0

n=10
L=1.d0

f=f*dsin(n*pi*x/L)
end function
!##################################################################################

subroutine trapezio_composto
implicit none
real*8,external :: f
real*8::a,b,h,integral
integer:: i,N
a=0.
b=1.
N=304

h=(b-a)/N

integral=(f(a)+f(b))*h/2

do i=1,N-1
    integral=integral+f(a+i*h)*h
enddo

print*,'por trapezio_composto'
print*,integral


end subroutine
!############################################################################################
subroutine simpson_composto
implicit none
real*8,external :: f
real*8::a,b,h,integral
integer:: i,N
a=-1.d0
b=1.d0
N=2000

h=(b-a)/N

integral=(f(a)+f(b))*h/3

do i=1,N-1,2
    integral=integral+f(a+i*h)*h*4/3
enddo
do i=2,N-1,2
    integral=integral+f(a+i*h)*h*2/3
enddo

print*,'por simpson_composto'
print*,integral

end subroutine
!################################################################################################
subroutine tres_oitavos_simpson
implicit none
real*8,external :: f
real*8::a,b,h,integral
integer:: i,N
a=-1.d0
b=1.d0
N=2000

h=(b-a)/(3*N)

do i=1,N
    integral=integral+f(a+(3*i-2)*h)+3*f(a+(3*i-1)*h)+3*f(a+(3*i)*h)+ f(a+(3*i+1)*h)
enddo
integral=integral*3*h/8
print*,'por tres_oitavos_simpson'
print*,integral


end subroutine
