! $Id$
! [compile with: ifort -xHOST -O3 -ip -r8 -fpp billiard.f90 polint.f]
! ODE integration methods tested on a simple anharmonic oscillator

program billiard; implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preheating model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! fields are referred to by their symbolic aliases
integer, parameter :: fields = 2                ! total number of scalar fields being evolved
integer, parameter :: phi = 1, chi = 2          ! symbolic aliases for scalar field components
integer, parameter :: n = 2*(fields+1)          ! total dynamical system dimension

! parameters of the scalar field potential are defined here
real, parameter :: lambda = 1.0, g2 = 1.875     ! scalar field potential
real, parameter :: mu = 1.76274717403907543734  ! Floquet growth per period

! potential and its derivatives are inlined in ...step() routines
#define Vx4(PHI,CHI) (lambda * (PHI)**2 + (2.0*g2) * (CHI)**2) * (PHI)**2
#define M2I(PHI,CHI) (/ lambda*(PHI)**2 + g2*(CHI)**2, g2*(PHI)**2 /)

! state vector packing (implemented as preprocessor macros)
#define $fi$ 1:2
#define $pi$ 3:4
#define $a$  5
#define $p$  6

! initial conditions for homogeneous field components
real, parameter ::  phi0 =  2.339383796213256
real, parameter :: dphi0 = -2.736358272992573


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter :: s = 1

integer i

write (*,'(g)') "# OUTPUT: alpha, t, p, dp;"

do i = 1,1000
        call evolve(-11.0+i/250.0, 1e-3)
end do

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! width of the phase space bundle
function width(bundle)
        real width, bundle(n,-s:s), gii(n), dy(n), a2
        
        ! reference trajectory and displacement
        a2 = bundle($a$,0)**2
        dy = bundle(:,s) - bundle(:,-s)
        
        ! phase space metric
        gii = (/ a2, a2, 1.0/a2, 1.0/a2, 1.0, 1.0 /)
        
        ! return width
        width = sqrt(sum(gii*dy*dy))
end function width

! shrink phase space bundle
subroutine shrink(bundle, width, w)
        real width, w, e; integer i, k
        real bundle(n,-s:s), shrunk(n,-s:s)
        real, parameter :: x(-s:s) = (/ -s:s /)
        
        ! resample the bundle
        do i = 1,n; do k = -s,s
                call polint(x, bundle(i,:), 2*s+1, w*x(k), shrunk(i,k), e)
        end do; end do
        
        ! update the bundle
        bundle = shrunk; width = w*width
end subroutine shrink

! evolve phase space bundle around alpha trajectory
subroutine evolve(alpha, da)
        real, value :: alpha, da
        real w0, bundle(n,-s:s)
        integer l, k
        
        ! machine precision accuracy
        real, parameter :: dt = 0.01
        
        ! initialize phase space bundle
        do k = -s,s; call inity(bundle(:,k), alpha+k*da); end do
        
        ! initial bundle width
        w0 = width(bundle)
        
        ! evolve the bundle, shrinking it if it gets too wide
        do l = 0,5000
                if (mod(l,10) == 0) write (*,'(8g)') alpha, l*dt, &
                        bundle($p$,0), (bundle($p$,s)-bundle($p$,-s))/(2.0*s*da)
                do k = -s,s; call gl10(bundle(:,k), dt); end do
                if (width(bundle) > 10.0*w0) call shrink(bundle, da, 0.05)
        end do
        
        ! flush the frame
        write (*,'(g)') "", ""
end subroutine evolve


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! homogeneous scalar fields in expanding FRW spacetime (EoM in Hamiltonian form)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize dynamical system
subroutine inity(y, alpha)
        real y(n), alpha; real chi0, dchi0, m2(fields), H2
        
        ! initial isocurvature field value
        chi0 = phi0 * exp(mu*alpha)
        
        ! set chi' to terminal velocity, curvature to flat
        H2 = dphi0**2/6.0 + Vx4(phi0,chi0)/12.0
        m2 = M2I(phi0,chi0); dchi0 = -m2(chi)*chi0/sqrt(H2)/2.0
        H2 = H2 + dchi0**2/6.0
        
        ! initial state vector
        y = (/ phi0, chi0, dphi0, dchi0, 1.0, -6.0*sqrt(H2) /)
end subroutine inity

! equations of motion dy/dt = f(y)
subroutine evalf(y, f)
        real, target :: y(n), f(n); real, pointer :: fi(:), pi(:), a, p
        
        ! unpack dynamical system vector
        fi => y($fi$); pi => y($pi$); a => y($a$); p => y($p$)
        
        ! Hamiltonian equations of motion
        f($fi$) = pi/a**2; f($pi$) = -a**4 * M2I(fi(phi),fi(chi)) * fi
        f($a$) = -p/6.0; f($p$) = sum(pi**2)/a**3 - Vx4(fi(phi),fi(chi))*a**3
end subroutine evalf

! Hamiltonian constraint violation
function omegak(y)
        real, target :: y(n); real, pointer :: fi(:), pi(:), a, p
        real omegak, P2, KE, PE
        
        ! unpack dynamical system vector
        fi => y($fi$); pi => y($pi$); a => y($a$); p => y($p$)
        
        ! Hamiltonian pieces
        P2 = p**2/12.0
        KE = sum(pi**2)/a**2/2.0
        PE = Vx4(fi(phi),fi(chi))*a**4/4.0
        
        ! residual curvature
        omegak = (KE+PE)/P2 - 1.0
end function omegak


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! implicit Gauss-Legendre method; symplectic with arbitrary Hamiltonian, A-stable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5
        real y(n), g(n,s), dt; integer i, k
        
        ! Butcher tableau for 8th order Gauss-Legendre method
        real, parameter :: a(s,s) = (/ &
                  0.5923172126404727187856601017997934066Q-1, -1.9570364359076037492643214050884060018Q-2, &
                  1.1254400818642955552716244215090748773Q-2, -0.5593793660812184876817721964475928216Q-2, &
                  1.5881129678659985393652424705934162371Q-3,  1.2815100567004528349616684832951382219Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -2.4592114619642200389318251686004016630Q-2, &
                  1.0318280670683357408953945056355839486Q-2, -2.7689943987696030442826307588795957613Q-3, &
                  1.1377628800422460252874127381536557686Q-1,  2.6000465168064151859240589518757397939Q-1, &
                  1.4222222222222222222222222222222222222Q-1, -2.0690316430958284571760137769754882933Q-2, &
                  4.6871545238699412283907465445931044619Q-3,  1.2123243692686414680141465111883827708Q-1, &
                  2.2899605457899987661169181236146325697Q-1,  3.0903655906408664483376269613044846112Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -0.9687563141950739739034827969555140871Q-2, &
                  1.1687532956022854521776677788936526508Q-1,  2.4490812891049541889746347938229502468Q-1, &
                  2.7319004362580148889172820022935369566Q-1,  2.5888469960875927151328897146870315648Q-1, &
                  0.5923172126404727187856601017997934066Q-1 /)
        real, parameter ::   b(s) = (/ &
                  1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  1.1846344252809454375713202035995868132Q-1 /)
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

end