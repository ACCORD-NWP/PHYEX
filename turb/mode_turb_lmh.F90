MODULE MODE_TURB_LMH
!     ###############
!
!!****  *MODE_TURB_LMH * Explicit computation of a parameterized
!!                            Horizontal Mixing Length scale (LMH)
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS * None
!!    ------------------
!!
!!    REFERENCE Rogel et al. (2024) 6th ACCORD Newsletter - 3D turbulence
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      L. Rogel       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        27/01/2025
!----------------------------------------------------------------------------
!
  USE MODD_TURB_n,     ONLY: TURB_t
  USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
  USE MODD_CTURB,      ONLY: CSTURB_t
  !USE PARKIND1 , ONLY : JPRB
  USE YOMHOOK ,  ONLY : LHOOK, DR_HOOK, JPHOOK
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: TURB_LMH

  CONTAINS
!------------------------------------------------------------------------------------------------
    SUBROUTINE TURB_LMH(D,TURBN,KGRADIENTSHSP,PDXX,PDYY,PDZX,PDZY,PTHL,PTKEM,PUT,PVT,PHGRAD,PLMV,PLMH)
      TYPE(DIMPHYEX_t),       INTENT(IN)   :: D             ! PHYEX variables dimensions structure
      TYPE(TURB_t),           INTENT(IN)   :: TURBN         ! modn_turbn (turb namelist) structure
      INTEGER,                INTENT(IN)   :: KGRADIENTSHSP
      REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PDXX, PDYY   ! Mesh size (at mass point)
      REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PTHL, PTKEM, PUT, PVT ! thetal, TKE, u-wind, v-wind  on mass point
      REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PDZX, PDZY       ! metric term dzdx, dzdy (flux point)
      REAL, DIMENSION(D%NIJT,D%NKT,KGRADIENTSHSP), INTENT(IN)  :: PHGRAD  ! Horizontal gradients structure
      REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PLMV             ! Vertical mixing length scale
      REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PLMH             ! Horizontal length scale output
      ! Local
      REAL, DIMENSION(D%NIJT,D%NKT) :: ZLWANG
      REAL, DIMENSION(D%NIJT,D%NKT) :: ZLSMAG
      REAL    :: ZDUDX2, ZDUDY2, ZDVDX2, ZDVDY2 ! Local values to avoid 0 division
      REAL    :: ZALPHAX, ZALPHAY ! Slope correction for Smagorinsky (1963) length scale
      ! loops
      INTEGER        :: JIJ,JK
      INTEGER        :: IKTB,IKTE, &
                      & IKT,IKB,IKE,&
                      & IKA,IKU,IKL,&
                      & IIJB,IIJE
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('TURB_LMH',0,ZHOOK_HANDLE)

      IKT=D%NKT
      IKTB=D%NKTB
      IKTE=D%NKTE
      IKB=D%NKB
      IKE=D%NKE
      IKA=D%NKA
      IKU=D%NKU
      IKL=D%NKL
      IIJE=D%NIJE
      IIJB=D%NIJB

      ASSOCIATE(&
      & PDUDX => PHGRAD(:,:,1),&
      & PDUDY => PHGRAD(:,:,2),&
      & PDVDX => PHGRAD(:,:,3),&
      & PDVDY => PHGRAD(:,:,4) &
      & )


      IF ((TURBN%CLMHTURB=='SMAG') .OR. (TURBN%CLMHTURB=='WANG')) THEN
        IF (TURBN%NSMAG == 0) THEN
          DO JK=IKTB,IKTE
            DO JIJ=IIJB,IIJE
              ZLSMAG(JIJ,JK) = TURBN%XCSMAG * SQRT( PDXX(JIJ,JK) * PDYY(JIJ,JK) )
            ENDDO
          ENDDO
        ELSE IF (TURBN%NSMAG == 1) THEN
          DO JK=IKTB,IKTE
            DO JIJ=IIJB,IIJE
              ZALPHAX = COS( ATAN( PDZX(JIJ,JK) / PDXX(JIJ,JK) ) )
              ZALPHAY = COS( ATAN( PDZY(JIJ,JK) / PDYY(JIJ,JK) ) )
              ZLSMAG(JIJ,JK) = TURBN%XCSMAG & 
                    & * SQRT( ZALPHAY *PDXX(JIJ,JK) * ZALPHAX * PDYY(JIJ,JK) )
            ENDDO
          ENDDO
        END IF
      END IF

      SELECT CASE (TURBN%CLMHTURB)

        CASE ('SMAG')

          PLMH(:,:) =  ZLSMAG(:,:)

        CASE ('WANG')

          DO JK=IKTB,IKTE
            DO JIJ=IIJB,IIJE
              ZDUDX2 = MAX(1E-14,PDUDX(JIJ,JK)**2)
              ZDUDY2 = MAX(1E-14,PDUDY(JIJ,JK)**2)
              ZDVDX2 = MAX(1E-14,PDVDX(JIJ,JK)**2)
              ZDVDY2 = MAX(1E-14,PDVDY(JIJ,JK)**2)
              ZLWANG(JIJ,JK) = ( TURBN%XWANG_DELTA & 
                / SQRT(PDXX(JIJ,JK)*PDYY(JIJ,JK)) ) ** TURBN%XWANG_ALPHA &
                & *   SQRT( PUT(JIJ,JK)**2 + PVT(JIJ,JK)**2 ) &
                & / ( (ZDVDX2 + ZDUDY2) * (ZDUDX2 + ZDVDY2) )**0.25
            ENDDO
          ENDDO

          DO JK=IKTB,IKTE
            DO JIJ=IIJB,IIJE
              PLMH(JIJ,JK) = MIN(ZLWANG(JIJ,JK), ZLSMAG(JIJ,JK))
            ENDDO
          ENDDO

        ! CASE ('SHLS') 
          ! A smarter shear-based length scale (to be defined)
          ! NOT IMPLEMENTED
          ! PLMH(:,:) = -999.

        CASE DEFAULT
          ! Use the same mixing length than the 1D vertical
          DO JK=IKTB,IKTE
            DO JIJ=IIJB,IIJE
              PLMH(JIJ,JK) = PLMV(JIJ,JK)
            ENDDO
          ENDDO

      END SELECT
      END ASSOCIATE

      IF (LHOOK) CALL DR_HOOK('TURB_LMH',1,ZHOOK_HANDLE)
    END SUBROUTINE TURB_LMH
END MODULE MODE_TURB_LMH
