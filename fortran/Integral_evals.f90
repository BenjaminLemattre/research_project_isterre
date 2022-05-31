module integral
    implicit none
!======================== DEFINITION OF THE LOCAL PARAMETERS OF THE MODULE ==========================!
    doubleprecision, parameter :: R_C = 3485.0d0 ! km
    doubleprecision, parameter :: R_A = 6371.2d0 ! km
    doubleprecision, parameter :: R_SAT = 6671.2d0 ! km
    doubleprecision, parameter :: PI = 3.1415926535897932384d0
    doubleprecision, parameter :: TWOPI = 2.0d0*PI
    complex(kind=8), parameter :: IM = (0.0d0, 1.0d0)

    doubleprecision, parameter :: R_B = R_C+2000.0d0 ! km
    doubleprecision, parameter :: Sigma_c = 3000.d0 ! S/m
    doubleprecision, parameter :: alpha_c = 30.d0
    doubleprecision, parameter :: RCB = R_C/R_B
    doubleprecision, parameter :: RAC = R_A/R_C !6371.2d0/3485.0d0 !== ratio of Earth radius and Earth core R_A/R_C in km

contains
    !============ RETURNS A(b), THE GAUNT-ELSSASSERS INTEGRAL COEFFICIENTS AT CORE SURFACE FROM BR ============!
    subroutine radmats(gauth, gauwt, p1, dp1, d2p1, p2, dp2, d2p2, p3, dp3, d2p3, &
            Ab, g, l1, l2, l3, n1, n2, n3, nn1, nn2, nn3, tmax)
    !=========================================================================
    !     Modified by Loïc Huder from "regular" radmats subroutine
    !=========================================================================
    !
    !     l3max  is the max degree of our magnetic field
    !     l1max  "  "   "     "    "  "   secular variation
    !     l2max   "  "   "     "    "  "  resulting core motion
    !
    !     triangle rule => svinduced = l2 + l3
    !     to be exact 2*tmax-1 > (l1+svinduced)
    !                     pmax > (l1+svinduced)+1 and 2**n
    !
    !     for l3 = l2 = l1 the minimum values for pmax are therefore
    !                      l3  1 ->  2 pmax =  8
    !                      l3  3 ->  4 pmax = 16
    !                      l3  5 -> 10 pmax = 32
    !                      l3 11 -> 20 pmax = 64
    !                      l3 21 -> 42 pmax =128
    !
    !=========================================================================

        integer, intent(in) :: l1,l2,l3,n1,n2,n3,nn1,nn2,nn3,tmax

        ! F2PY: Adding intent(out) to 'Ab' removes it from the list of arg and makes it returned instead !
        real*8, dimension(2*n2, n1) :: Ab
        real*8, dimension(n3), intent(in) :: g
        real*8, dimension(n1) :: r, eqwts
        complex(kind=8), dimension(2*tmax,0:l2) :: eisph

        real*8, dimension(nn1,tmax), intent(in) ::  p1, dp1, d2p1
        real*8, dimension(nn2,tmax), intent(in) ::  p2, dp2, d2p2
        real*8, dimension(nn3,tmax), intent(in) ::  p3, dp3, d2p3

        real*8, dimension(tmax), intent(in) :: gauth, gauwt

        complex(kind=8), dimension(3,3,tmax,2*tmax) :: b

        real*8, dimension(tmax,2*tmax,3,2) :: at, as
        real*8, dimension(2,3) :: ct, cs

        real(8) :: bb, brcmb, phi, theta
        integer :: i, l, m, cmpt, kcnt, l2i, lm, lp1, ltmp, m2i, ph, rsc, rss, th
        integer :: pmax
        pmax = 2*tmax

        ltmp=max(l1,l2)
        ltmp=max(ltmp,l3)
        ltmp=3*ltmp+1

        if(ltmp >= pmax)   stop 'increase pmax in radmats'
        if(ltmp/2 >= tmax) stop 'increase tmax in radmats'

        Ab(:,:)=0.0

        do ph = 1, pmax
             phi = twopi*(ph-1)/pmax
             do m = 0, l2
                eisph(ph,m) = cos(m*phi) + im*sin(m*phi)
             enddo
        enddo

        !   calculate the weights to relate cmb fields to coefficients
        kcnt=0
        bb=RAC*RAC
        do l=1, l1
             bb=bb*RAC
             lp1=l+1
             do m=1, 2*l+1
                kcnt = kcnt + 1
                brcmb = lp1*bb
                eqwts(kcnt) = 1.d0/brcmb
             enddo
        enddo

        !---------

        call calcb(b, g, l3, tmax, pmax, gauth, p3, dp3, d2p3)

        ! loop over each velocity mode

        lm = 1
        do l2i = 1, l2
             rsc = l2i*l2i

             do m2i = 0, l2i
                lm = lm +1

                do th = 1, tmax
                   theta = gauth (th)
                   do ph = 1, pmax
                      call indeqnr(cs, ct, b(1,1,th,ph), p2(lm,th), &
                              dp2(lm,th), d2p2(lm,th), eisph(ph,m2i), l2i, m2i, theta)
                      do i = 1,2
                         do cmpt = 1,3
                            as(th,ph,cmpt,i) = cs(i,cmpt)
                            at(th,ph,cmpt,i) = ct(i,cmpt)
                         enddo
                      enddo
                   enddo
                enddo

                call ast2h(at(1,1,1,1), tmax, pmax, r, l1, gauth, gauwt, p1)
                do i=1, n1
                   Ab(rsc,i)=Ab(rsc,i)+r(i)*eqwts(i)
                enddo
                call ast2h(as(1,1,1,1), tmax, pmax, r, l1, gauth, gauwt, p1)
                do i=1, n1
                   Ab(n2+rsc,i)=Ab(n2+rsc,i)+r(i)*eqwts(i)
                enddo

                if (m2i > 0) then
                   call ast2h(at(1,1,1,2), tmax, pmax, r, l1, gauth, gauwt, p1)
                   do i=1, n1
                      Ab(rss,i)=Ab(rss,i)+r(i)*eqwts(i)
                   enddo
                   call ast2h(as(1,1,1,2), tmax, pmax, r, l1, gauth, gauwt, p1)
                   do i=1, n1
                      Ab(n2+rss,i)=Ab(n2+rss,i)+r(i)*eqwts(i)
                   enddo
                   rsc = rsc + 2
                   rss = rsc + 1
                else
                   rsc = rsc + 1
                   rss = rsc + 1
                endif

             enddo
        enddo

        return
    end subroutine
!====================================================================================================!
!============= DIFFERENTS SUBROUTINES TO CALCULATE GAUNT-ELSSASSER INTEGRAL COEFFICIENTS ============!
!====================================================================================================!
    subroutine indeqnr(cffsr, cfftr, b, p, dp, d2p, eim, l, m, theta)

        complex*16, intent(in) :: b(3,3), eim
        real*8, intent(in) :: p, dp, d2p, theta
        integer, intent(in) :: l, m

        real*8, intent(out) ::  cffsr(2,3), cfftr(2,3)

        complex*16 cffs(3), cfft(3)
        complex*16 y, dydt, dydp, d2ydtt, d2ydtp
        complex*16 t(3,3), s(3,3), dgh1(3), dgh2(3)
        complex*16 b11
        real*8 rr

! ********* warning radius unity used here **************

        y = p*eim
        dydt = dp*eim
        dydp = p*m*im*eim
        d2ydtt = d2p*eim
        d2ydtp = dp*m*im*eim

        rr=1.d0/R_C

        call calct(t, y, dydt, dydp, d2ydtt, d2ydtp, rr, theta, m)
        call calcs(s, y, dydt, dydp, d2ydtt, d2ydtp, rr, theta, m)

        call dgradh(dgh1, b, s, rr, theta)
        call dgradh(dgh2, s, b, rr, theta)

        b11 = b(1,1)

!       ds/dr = -s/r
!       and add the component from the shear equation to 1

        cffs(1) = -b11*rr*s(1,1) + dgh1(1) - dgh2(1) + l*(l+1)*rr*rr*y*b11
        cffs(2) = -b11*rr*s(2,1) + dgh1(2) - dgh2(2)
        cffs(3) = -b11*rr*s(3,1) + dgh1(3) - dgh2(3)

        call dgradh(dgh1, b, t, rr, theta)
        call dgradh(dgh2, t, b, rr, theta)

!       dt/dr = -t/r

        cfft(1) = -b11*rr*t(1,1) + dgh1(1) - dgh2(1)
        cfft(2) = -b11*rr*t(2,1) + dgh1(2) - dgh2(2)
        cfft(3) = -b11*rr*t(3,1) + dgh1(3) - dgh2(3)

!       post-normalise velocities by c, leave shear because we want coefficients
!       to be of form c*v^t (to be dimensionally similar)

        cffs(:)=cffs(:)*R_C
        cfft(:)=cfft(:)*R_C

        cffsr(1,:) = DBLE(cffs(:))
        cffsr(2,:) = aIMAG(cffs(:))
        cfftr(1,:) = DBLE(cfft(:))
        cfftr(2,:) = aIMAG(cfft(:))

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calcb(bv, g, lmax, tmax, pmax, gauth, p, dp, d2p)
!c     JB version of AJ version of dave lloyd calcb
!c     arguments changed to include p dp

        integer, intent(in) :: lmax, tmax, pmax
        real*8, intent(in) ::  g(lmax*(lmax+2)), gauth(tmax)
        real*8, intent(in) :: p((lmax+1)*(lmax+2)/2,tmax)
        real*8, intent(in) :: dp((lmax+1)*(lmax+2)/2,tmax)
        real*8, intent(in) :: d2p((lmax+1)*(lmax+2)/2,tmax)

        complex*16, intent(out) :: bv(3, 3, tmax, pmax)

        real*8 plm, dplm, d2plm
        integer th, ph, lm, lmc, lms, l, m, l12
        real*8 sinth, cs, bl2, phi, cosmp, sinmp, sum1, sum2

        do th = 1, tmax
            sinth = sin (gauth(th))
            cs = cos (gauth(th))/(sinth*sinth)

            do ph = 1, pmax
                bl2 = RAC*RAC

                bv(:, :, th, ph) = (0.,0.)

                 lm = 1
                 do l = 1, lmax
                     bl2 = bl2*RAC
                     lmc = l*l
                     lm  = lm + 1
                     l12 = (l+1)*(l+1)

                     m = 0
                     plm = p(lm,th)
                     dplm = dp(lm,th)
                     d2plm = d2p(lm,th)

                     bv(1,1,th,ph) = bv(1,1,th,ph) + (l+1)*plm*bl2*g(lmc)
                     bv(1,2,th,ph) = bv(1,2,th,ph) + (l+1)*dplm*bl2*g(lmc)

                     bv(2,1,th,ph) = bv(2,1,th,ph) - dplm*bl2*g(lmc)
                     bv(2,2,th,ph) = bv(2,2,th,ph) - d2plm*bl2*g(lmc)

                     do m = 1, l
                         lmc = lmc + 1
                         lms = lmc + 1
                         lm  = lm  + 1
                         plm = p(lm,th)
                         dplm = dp(lm,th)
                         d2plm = d2p(lm,th)

                         phi = m*(ph-1)*twopi/pmax
                         cosmp = cos(phi)
                         sinmp = sin(phi)

                         sum1 = bl2*  ( g(lmc)*cosmp + g(lms)*sinmp)
                         sum2 = bl2*m*(-g(lmc)*sinmp + g(lms)*cosmp)

                         bv(1,1,th,ph) = bv(1,1,th,ph) + (l+1)* plm*sum1
                         bv(1,2,th,ph) = bv(1,2,th,ph) + (l+1)*dplm*sum1
                         bv(1,3,th,ph) = bv(1,3,th,ph) + (l+1)* plm*sum2


                         bv(2,1,th,ph) = bv(2,1,th,ph) - dplm*sum1
                         bv(2,2,th,ph) = bv(2,2,th,ph) - d2plm*sum1
                         bv(2,3,th,ph) = bv(2,3,th,ph) - dplm*sum2

                         bv(3,1,th,ph) = bv(3,1,th,ph) - plm*sum2/sinth
                         bv(3,2,th,ph) = bv(3,2,th,ph)&
                                 + (cs*plm - dplm/sinth)*sum2
                         bv(3,3,th,ph) = bv(3,3,th,ph)&
                                 + plm*m*m*sum1/sinth

                         lmc = lmc + 1

                     end do
                end do
            end do
        end do

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calcs(s, y, dydt, dydp, d2ydtt, d2ydtp, rr, theta, m)

        complex*16, intent(in) :: y, dydt, dydp, d2ydtt, d2ydtp
        real*8, intent(in) :: rr, theta
        integer, intent(in) :: m
        complex*16, intent(out) :: s(3,3)

        real*8 sinth, coth, rsinth

        sinth = sin(theta)
        coth = cos(theta)/sinth
        rsinth = rr/sinth

        s (1,:) = (0.,0.)

        s (2,1) = dydt*rr
        s (3,1) = dydp*rsinth

        s (2,2) = d2ydtt*rr
        s (3,2) = (d2ydtp/sinth - coth/sinth*dydp)*rr

        s (2,3) = d2ydtp*rr
        s (3,3) = -m*m*y*rsinth

    end subroutine
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calct(t, y, dydt, dydp, d2ydtt, d2ydtp, rr, theta, m)
        !c     25 october sign convention changed so that v = curl(rt)

        complex*16, intent(in) :: y, dydt, dydp, d2ydtt, d2ydtp
        real*8, intent(in) :: rr, theta
        integer, intent(in) :: m
        complex*16, intent(out) :: t (3,3)

        real*8 sinth, coth, rsinth

        sinth = sin(theta)
        coth = cos (theta)/sinth
        rsinth = rr/sinth

        t (1, :) = (0.,0.)

        t (2,1) = dydp*rsinth
        t (3,1) = -dydt*rr

        t (2,2) = -(coth/sinth*dydp - d2ydtp/sinth)*rr
        t (3,2) = -d2ydtt*rr

        t (2,3) = -m*m*y*rsinth
        t (3,3) = -d2ydtp*rr
    end subroutine

!c    JB version of routine by dave lloyd: 24 march 88
!c
!c    utilities for manipulating basic vectors. also defined is an extended
!c    vector which includes surface derivatives.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine dgradh(c, a, b, rr, theta)

        complex*16, intent(in) :: a(3), b(3,3)
        real*8, intent(in) :: rr, theta
        complex*16, intent(out) :: c(3)

        real*8 sinth, coth

        sinth = 1.0/sin(theta)
        coth  = cos(theta)*sinth

        c(1) = (a(2)*b(1,2) + a (3)*b(1,3)*sinth - a(2)*b(2,1) - a(3)*b(3,1))*rr
        c(2) = (a(2)*b(2,2) + a (3)*b(2,3)*sinth - a(3)*b(3,1)*coth + a(2)*b(1,1))*rr
        c(3) = (a(2)*b(3,2) + a (3)*b(3,3)*sinth + a(3)*b(1,1) + a(3)*b(2,1)*coth)*rr

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!C   DAVID LLOYD, 13-7-87
!C                18-3-88 VECTOR TRANSFORM ADDED
!C
!C       S P H E R I C A L     T R A N S F O R M S
!C       =========================================


    subroutine AST2H(F, FL, FM, R, RL, GAUTH, GAUWT, P)

        !C     Andy's changed version of Dave Lloyd's ST2H (24 October 1989)
        !C     P,DP Added as argument to save call to legs
        !C     8 November dimension of FTF changed to 128 (see below)
        !C     9 november R dimensioned as RL*(RL+2)
        !C     27 January 1990 dimension of FTF changed to 256 !
        !
        !C!     Calculate (2*l+1)/4PI INTEGRAL { F Y_l^m } dOmega forall l,m

        !C   TAKES A GRID OF POINTS ON A SPHERICAL SURFACE AND TRANSFORMS THEM INTO THEIR
        !C   SPHERICAL HARMONIC SET (UP TO DEGREE AND ORDER RL).
        !
        !C   ASSERT: 2*RL < FL AND RL < FM/2 ==> INTEGRATION EXACT!
        !
        !C   FL, FM MUST BE POWERS OF 2.
        !C   THE THETA DIRECTION GRID POINTS ARE TAKEN TO BE AT THE APPROPRIATE GAUSS-
        !C   QUADRATURE POINTS, WHERE THE THETA VALUES ARE GIVEN IN GAUTH. IF GX ARE
        !C   THE QUADRATURE POINTS FOR (-1, 1) THEN:
        !C             GAUTH (I) = PI/2*(GX (I) +/- 1).
        !C   GAUWT HOLDS THE WEIGHT VALUES FOR (-1, 1) - THIS RANGE BECAUSE IF THE
        !C   SUBST COS (THETA) --> X HOLDS, THEN (0, 2*PI)- --> (1, -1).
        !C   THE PHI DIRECTION GRID POINTS ARE EVENLY SPACED.
        !C   MAXIMUM GRID PHI POINTS IS 64.
        !C   MAXIMUM RESULT DEGREE AND ORDER IS 30.

        implicit none
        integer, intent(in) :: FL, FM, RL
        real*8, intent(in) :: GAUTH(0:FL-1), GAUWT(0:FL-1)
        real*8, intent(in) :: P((RL+1)*(RL+2)/2,0:FL-1)
        real*8, intent(in) :: F (0:FL-1, 0:FM-1)

        real*8, intent(out) :: R(RL*(RL+2))


        real*8 THETA, WEIGHT

        complex*16 FTF(0:255)
        integer TH, L, M, LM, LMC, LMS

        ! INITIALISE R TO BE ZERO.
        R(:) = 0.

        do TH = 0, FL-1
            THETA = GAUTH(TH)
            WEIGHT  = GAUWT(TH)/FM

            ! SLICE F (TH, 0:FM-1) INTO FTF (0:FM-1) FOR FFT
            ! LH: Possible bug here (fails with flag -fcheck=bounds):
            ! We try to make an array of size FM fit into a hard coded size of 256...
            FTF(:) = F(TH, :)

            call FOUR1(FTF, FM, 1)

            LM = 1
            do L = 1, RL
                LMC = L*L
                LM  = LM + 1

                M = 0

                R(LMC) = R(LMC) + P(LM,TH)*DBLE(FTF(M))*WEIGHT

                do M = 1, L
                    LM  = LM + 1
                    LMC = LMC + 1
                    LMS = LMC + 1

                    R(LMC) = R(LMC) + DBLE(FTF(M))*P(LM,TH)*WEIGHT
                    R(LMS) = R(LMS) + aIMAG(FTF(M))*P(LM,TH)*WEIGHT

                    LMC = LMC + 1
                end do

            end do

        end do

        !C     POST-SCALE THE COEFFICIENTS ALLOWING FOR SCHMIDT NORMALISATION
        !C     AND SCALE FACTOR OF FFT.

        do L = 1, RL
            WEIGHT = DBLE(2*L+1)/2
            LMC = L*L

            R(LMC) = R(LMC) * WEIGHT

            do M = 1, L
                LMC = LMC + 1
                LMS = LMC + 1

                R(LMC) = R(LMC) * WEIGHT
                R(LMS) = R(LMS) * WEIGHT

                LMC = LMC + 1
            end do

        end do

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine FOUR1(DATA2, NN, ISIGN)
        integer, intent(in) :: NN, ISIGN
        complex*16, intent(out) :: DATA2(0:255)

        integer :: I, ISTEP, J, K, M, MMAX, N
        doubleprecision :: TEMPI, TEMPR
        real*8 WR,WI,WPR,WPI,WTEMP,THETA
        real*8 DATA(2*NN)

        do I=1, 2*NN, 2
            K=INT(I/2)
            DATA(I)=DBLE(DATA2(K))
            DATA(I+1)=dIMAG(DATA2(K))
        end do

        N=2*NN
        J=1

        do I=1, N, 2
            if(J > I) then
                TEMPR=DATA(J)
                TEMPI=DATA(J+1)
                DATA(J)=DATA(I)
                DATA(J+1)=DATA(I+1)
                DATA(I)=TEMPR
                DATA(I+1)=TEMPI
            end if
            M=N/2
            do while((M >= 2).AND.(J > M))
                J=J-M
                M=M/2
            end do
            J = J+M
        end do
        MMAX=2
        do while(N > MMAX)
            ISTEP=2*MMAX
            THETA=twopi/(ISIGN*MMAX)
            WPR=-2.D0*DSIN(0.5D0*THETA)**2
            WPI=DSIN(THETA)
            WR=1.D0
            WI=0.D0

            do M=1, MMAX, 2
                do I=M, N, ISTEP
                    J=I+MMAX
                    TEMPR=WR*DATA(J)-WI*DATA(J+1)
                    TEMPI=WR*DATA(J+1)+WI*DATA(J)
                    DATA(J)=DATA(I)-TEMPR
                    DATA(J+1)=DATA(I+1)-TEMPI
                    DATA(I)=DATA(I)+TEMPR
                    DATA(I+1)=DATA(I+1)+TEMPI
                end do
                WTEMP=WR
                WR=WR*WPR-WI*WPI+WR
                WI=WI*WPR+WTEMP*WPI+WI
            end do

            MMAX=ISTEP
        end do
        do I=1, 2*NN, 2
            K=INT(I/2)
            DATA2(K) = dCMPLX(DATA(I),DATA(I+1))
        end do

        return
    end subroutine

end module
