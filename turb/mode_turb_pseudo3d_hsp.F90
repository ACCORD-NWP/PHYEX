MODULE MODE_TURB_PSEUDO3D_HSP
!     ###############
!
!!****  *MODE_TURB_PSEUDO3D_HSP * Explicit computation of the parameterized
!!                                Horizontal Shear Production (HSP)
!!
!!    PURPOSE
!!    -------
!!
!!**  INTERFACE *CALL* TURB_PSEUDO3D_HSP
!!    --------
!!
!!**  METHOD See documentation
!!    ------
!!
!!    EXTERNAL * None
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
  USE MODD_CST,        ONLY : CST_t
  USE MODD_CTURB,      ONLY : CSTURB_t
  USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
  USE MODD_CTURB,      ONLY: CSTURB_t
  USE MODD_TURB_n,     ONLY: TURB_t
  !USE PARKIND1 ,       ONLY : JPRB
  USE YOMHOOK ,        ONLY : LHOOK, DR_HOOK, JPHOOK

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: TURB_PSEUDO3D_HSP     ! Explicit computation of the parameterized HSP

  CONTAINS
  !------------------------------------------------------------------------------------------------
    SUBROUTINE TURB_PSEUDO3D_HSP(D,CST,CSTURB,TURBN,KGRADIENTSHSP,PHGRAD,PTKET,PLMH,PHSP,PLMH_OUT)

      TYPE(DIMPHYEX_t), INTENT(IN)    :: D       ! PHYEX variables dimensions structure
      TYPE(CST_t),      INTENT(IN)    :: CST  ! modd_cst phyex constants structure
      TYPE(CSTURB_t),   INTENT(IN)    :: CSTURB  ! modd_csturb turb constant structure
      TYPE(TURB_t),     INTENT(IN)    :: TURBN   ! modn_turbn (turb namelist) structure
      INTEGER,          INTENT(IN)    :: KGRADIENTSHSP ! Horizontal gradients number
      REAL, DIMENSION(D%NIJT,D%NKT,KGRADIENTSHSP), INTENT(IN)  :: PHGRAD ! horizontal gradients
      REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PTKET ! TKE at t on mass level
      REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PLMH  ! horizontal length scale at mass point at t
      REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHSP  ! Horizontal Shear Production
      ! Optional output horizontal scale ! (useful for diags if modified)
      REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT),OPTIONAL   :: PLMH_OUT
      ! Local scalar
      INTEGER        :: JIJ,JK
      INTEGER        :: IKTB,IKTE, &
                      & IKT,IKB,IKE,&
                      & IKA,IKU,IKL,&
                      & IIJB,IIJE
      REAL           :: ZDEF, ZDIV
      !! local arrays
      REAL, DIMENSION(D%NIJT,D%NKT) :: ZCHI      ! horizontal gradients array
      REAL, DIMENSION(D%NIJT,D%NKT) :: ZLMH      ! horizontal length scale  at t
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      !------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('TURB_PSEUDO3D_HSP',0,ZHOOK_HANDLE)

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

      !------------------------------------------------
      ! -- compute pseudo HSP dependence on DEF and DIV 

      ASSOCIATE(&
      & PDUDX => PHGRAD(:,:,1),&
      & PDUDY => PHGRAD(:,:,2),&
      & PDVDX => PHGRAD(:,:,3),&
      & PDVDY => PHGRAD(:,:,4) &
      & )

      DO JK=IKTB,IKTE
        DO JIJ=IIJB,IIJE
           ZDIV=PDUDX(JIJ,JK)+PDVDY(JIJ,JK)                ! div = (du/dx) + (dv/dy)
           ZDEF=SQRT(((PDUDX(JIJ,JK)-PDVDY(JIJ,JK))**2)  & ! Def = sqrt( ((du/dx) - (dv/dy))^2 
               & +   ((PDUDY(JIJ,JK)+PDVDX(JIJ,JK))**2))   !            +((du/dy) + (dv/dx))^2 )
           ! ZCHI(JIJ,JK) = CSTURB%XC0*( (SQRT((ZDEF**2)+(CSTURB%XC1**2)*(ZDIV**2))-CSTURB%XC2*ZDIV)**3 )
           ZCHI(JIJ,JK) = (CSTURB%XC0)**(1./3.) * &
               & (SQRT((ZDEF**2)+(CSTURB%XC1**2)*(ZDIV**2))-CSTURB%XC2*ZDIV)
        ENDDO
      ENDDO

      !--------------------------------------------------------
      ! -- Limiter on horizontal mixing length based on dissipation
      ! -- (from the equilibrium hyopthesis)

      DO JK=IKTB,IKTE
        DO JIJ=IIJB,IIJE
          ZLMH(JIJ,JK) = MIN( PLMH(JIJ,JK), &
            & (TURBN%XCED)**(1./3.)*SQRT(PTKET(JIJ,JK)) / &
            & MAX(ZCHI(JIJ,JK),CST%XMNH_EPSILON ))
        ENDDO
      ENDDO

      !-------------------------------------------
      ! -- Final generic computation of Pseudo-3D HSP
        DO JK=IKTB,IKTE
          DO JIJ=IIJB,IIJE
            PHSP(JIJ,JK) = ZLMH(JIJ,JK)**2 * ZCHI(JIJ,JK)**3
          ENDDO
        ENDDO

      ! Output of horizontal length scale if required
      IF (PRESENT(PLMH_OUT)) THEN
        DO JK=IKTB,IKTE
          DO JIJ=IIJB,IIJE
            PLMH_OUT(JIJ,JK) = ZLMH(JIJ,JK)
          ENDDO
        ENDDO
      ENDIF

      END ASSOCIATE
      IF (LHOOK) CALL DR_HOOK('TURB_PSEUDO3D_HSP',1,ZHOOK_HANDLE)
    END SUBROUTINE TURB_PSEUDO3D_HSP
END MODULE  MODE_TURB_PSEUDO3D_HSP
