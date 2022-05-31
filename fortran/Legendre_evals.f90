module legendre
    ! Spectra CALCULATION OF LEGENDRE POLYNOMIALS
    implicit none
    doubleprecision, parameter :: pi = 3.1415926535897932384d0
    doubleprecision, parameter :: sqrt2 = dsqrt(2.d0)
contains

    subroutine gauleg(lower_bound, upper_bound, gauss_pts, gauss_wgths, N)
        ! Generates N gauss points and weights between a lower and upper bound.
        integer, intent(in) :: N
        real*8, intent(in) ::  lower_bound, upper_bound
        real*8, intent(out) :: gauss_pts(N), gauss_wgths(N)

        real*8, parameter :: EPS=1.D-14
        integer :: I, J, M
        real*8 :: PP, PPP, P1, P2, P3, XM, XL, Z, Z1
        logical :: converged

        M=(N+1)/2
        XM=0.5D0*(upper_bound+lower_bound)
        XL=0.5D0*(upper_bound-lower_bound)

        do I=1,M
            converged = .false.
            Z = COS(pi*(I-.25D0)/(N+.5D0))
            do while(.not. converged)
                P1=1.D0
                P2=0.D0

                do J=1,N
                    P3=P2
                    P2=P1
                    P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
                end do

                PP=N*(Z*P1-P2)/( (Z-1.D0)*(Z+1.D0) )
                PPP=N*(Z*P1-P2)
                Z1=Z
                Z=Z1-P1/PP

                if (ABS(Z-Z1) < EPS) converged = .true.
            end do

            gauss_pts(I)=XM-XL*Z
            gauss_pts(N+1-I)=XM+XL*Z
            gauss_wgths(I)=2.D0*XL*(1.D0-Z)*(1.D0+Z)/(PPP*PPP)
            gauss_wgths(N+1-I)=gauss_wgths(I)
        end do

        return
    end subroutine


    subroutine plmbar2(p, dp, dp2, z, lmax, inorm)
        !c
        !c  evaluates normalized associated legendre function p(l,m) as function of
        !c   z=cos(colatitude) using recurrence relation starting with p(l,l)
        !c   and then increasing l keeping m fixed.  normalization is:
        !c   integral(y(l,m)*y(l,m))=4.*pi, where y(l,m) = p(l,m)*exp(i*m*longitude),
        !c   which is incorporated into the recurrence relation. p(k) contains p(l,m)
        !c   with k=(l+1)*l/2+m+1; i.e. m increments through range 0 to l before
        !c   incrementing l. routine is stable in single and double precision to
        !c   l,m = 511 at least; timing proportional to lmax**2
        !c   r.j.o connell 7 sept. 1989

        !c   a.jackson 19 october 1989  code added at end:
        !c!!   (1) choice of normalisation added ---
        !c       if inorm = 1 schmidt normalisation is chosen, where
        !c       p[schmidt](l,m) = 1/sqrt(2*l+1)*p[normalised](l,m)
        !c       inorm = 2 for fully normalised harmonics, otherwise error
        !c   (2) derivatives added and stored in dp(k)
        !c       using same arrangement as for p(k)
        !c

        real*8, intent(in) :: z
        integer, intent(in) :: lmax, inorm

        ! Adding intent(out) makes f2py remove them from the arg and does not work with assumed shapes
        real*8 :: p(:), dp(:), dp2(:)

        real*8 :: pm1, pm2, plm, pmm
        real *8 :: coth, fac, fac1, fac2, f1, f2, fnum, fden, sintsq

        integer :: l, m, k, kstart

        if (lmax < 0) then
            write(6, *) 'Incorrect lmax (should be positive integer, got ', lmax, 'instead'
            stop
        end if
        if (abs(z) > 1.d0) then
            write(6, *) 'Incorrect gauss point (abs value should be lower than 1., got ', abs(z), 'instead'
            stop
        end if
        if(inorm /= 1 .and. inorm /= 2) then
            write(6,*) 'inorm incorrect: use inorm=1 for schmidt normalisation>,/,inorm=2 for full normalisation'
            stop
        end if

        ! --case for p(l,0)
        pm2=1.d0
        p(1)=1.d0
        dp(1)=0.d0
        if (lmax == 0) return

        pm1=z
        p(2)=dsqrt(3.d0)*pm1
        k=2
        do l=2, lmax
            k=k+l
            plm=(dfloat(2*l-1)*z*pm1-dfloat(l-1)*pm2)/dfloat(l)
            p(k)=dsqrt(dfloat(2*l+1))*plm
            pm2=pm1
            pm1=plm
        end do

        ! --case for m > 0
        pmm = 1.d0
        sintsq = (1.d0-z)*(1.d0+z)
        fnum = -1.0d0
        fden = 0.0d0
        kstart = 1
        do m=1 ,lmax
            ! --case for p(m,m)
            kstart = kstart+m+1
            fnum = fnum+2.0d0
            fden = fden+2.0d0
            pmm = pmm*sintsq*fnum/fden
            pm2 = dsqrt(dfloat(4*m+2)*pmm)
            p(kstart) = pm2
            if (m == lmax) exit

            !--case for p(m+1,m)
            pm1=z*dsqrt(dfloat(2*m+3))*pm2
            k = kstart+m+1
            p(k) = pm1

            ! --case for p(l,m) with l > m+1
            if (m < (lmax-1)) then
                do l = m+2, lmax
                    k = k+l
                    f1=dsqrt(dfloat((2*l+1)*(2*l-1))/dfloat((l+m)*(l-m)))
                    f2=dsqrt(dfloat((2*l+1)*(l-m-1)*(l+m-1))/dfloat((2*l-3)*(l+m)*(l-m)))
                    plm=z*f1*pm1-f2*pm2
                    p(k) = plm
                    pm2 = pm1
                    pm1 = plm
                end do
            end if
        end do

        ! choice of normalisation:
        if(inorm == 1) then
            k=1
            do l=1,lmax
                fac=1.d0/dsqrt(dfloat(2*l+1))
                do m=0, l
                    k=k+1
                    p(k)=p(k)*fac
                end do
            end do
        endif

        ! now find derivatives of p(z) wrt theta, where z=cos(theta)
        ! recurrence is same regardless of normalisation, since l=constant
        dp(2)=-p(3)
        dp(3)=p(2)
        k=3
        do l=2, lmax
            k=k+1
            !c!         treat m=0 and m=l separately
            dp(k)=-dsqrt(dfloat(l*(l+1))/2.d0)*p(k+1)
            dp(k+l)=dsqrt(dfloat(l)/2.d0)*p(k+l-1)
            do m=1, l-1
                k=k+1
                fac1=dsqrt( dfloat( (l-m)*(l+m+1) ) )
                fac2=dsqrt( dfloat( (l+m)*(l-m+1) ) )
                if(m == 1) fac2 = fac2*sqrt2
                dp(k)=0.5d0*( fac2*p(k-1) - fac1*p(k+1) )
            end do
            k=k+1
        end do

        ! and now find second deriavtives
        coth = z/dsqrt(sintsq)

        dp2(1)=0.d0
        k=1
        do l=1, lmax
            do m=0, l
                k=k+1
                dp2(k) = ((m*m)/sintsq - l*(l+1))*p(k) - coth*dp(k)
            end do
        end do

        return
    end subroutine

end module
