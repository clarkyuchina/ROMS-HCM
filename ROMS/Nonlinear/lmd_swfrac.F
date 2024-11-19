#include "cppdefs.h"

#if defined NONLINEAR && (defined LMD_SKPP || defined SOLAR_SOURCE) && \
    defined SOLVE3D


      SUBROUTINE bio_swfrac_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Zscale, Z, Chla, swdk)
!
!svn $Id: lmd_swfrac.F 889 2018-02-10 03:32:52Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine computes the  fraction  of  solar shortwave flux    !
!  penetrating to specified depth (times Zscale) due to exponential    !
!  decay in water chlorophyll-dependent parameters.                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Zscale   scale factor to apply to depth array.                   !
!     Z        vertical height (meters, negative) for                  !
!              desired solar short-wave fraction.                      !
!     Chla     chlorophyll-a concentration .                           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     swdk     shortwave (radiation) fractional decay.                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Ohlmann, J. C., 2003: Ocean Radiant Heating in Climate Models.      !
!     J. Climate,16, 1337–1351,                                        !
!     https://doi.org/10.1175/1520-0442-16.9.1337.                     !
!                                                                      !
!  This routine was adapted from Bill Large 1995 code.                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_mixing
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

      real(r8), intent(in) :: Zscale

      real(r8), intent(in) :: Z(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in) :: Chla(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(out) :: swdk(IminS:ImaxS,JminS:JmaxS)
!
!  Local variable declarations.
!
      integer :: Jindex, i, j

      real(r8) :: A1, A2, B1, B2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Use Ohlmann (2003) two wavelength bands solar absorption model.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          A1 = 0.571_r8 + 0.025_r8* LOG(0.149_r8*Chla(i,j))
          A2 = 0.223_r8 + 0.010_r8* LOG(2.329_r8*Chla(i,j))
          B1 = 0.015_r8 + 0.176_r8*SQRT(0.462_r8*Chla(i,j))
          B2 = 0.688_r8 + 0.060_r8* LOG(0.125_r8*Chla(i,j))
 
          swdk(i,j)=A1*EXP(Zscale*B1*Z(i,j)) + A2*EXP(Zscale*B2*Z(i,j))
        END DO
      END DO
      RETURN
      END SUBROUTINE bio_swfrac_tile


      SUBROUTINE lmd_swfrac_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Zscale, Z, swdk)
!
!svn $Id: lmd_swfrac.F 889 2018-02-10 03:32:52Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine computes the  fraction  of  solar shortwave flux    !
!  penetrating to specified depth (times Zscale) due to exponential    !
!  decay in Jerlov water type.                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Zscale   scale factor to apply to depth array.                   !
!     Z        vertical height (meters, negative) for                  !
!              desired solar short-wave fraction.                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     swdk     shortwave (radiation) fractional decay.                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Paulson, C.A., and J.J. Simpson, 1977: Irradiance meassurements     !
!     in the upper ocean, J. Phys. Oceanogr., 7, 952-956.              !
!                                                                      !
!  This routine was adapted from Bill Large 1995 code.                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_mixing
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

      real(r8), intent(in) :: Zscale

      real(r8), intent(in) :: Z(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(out) :: swdk(IminS:ImaxS,JminS:JmaxS)
!
!  Local variable declarations.
!
      integer :: Jindex, i, j

      real(r8), dimension(IminS:ImaxS) :: fac1, fac2, fac3

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Use Paulson and Simpson (1977) two wavelength bands solar
!  absorption model.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Jindex=INT(MIXING(ng)%Jwtype(i,j))
          fac1(i)=Zscale/lmd_mu1(Jindex)
          fac2(i)=Zscale/lmd_mu2(Jindex)
          fac3(i)=lmd_r1(Jindex)
        END DO
!!DIR$ VECTOR ALWAYS
        DO i=Istr,Iend
          swdk(i,j)=EXP(Z(i,j)*fac1(i))*fac3(i)+                        &
     &              EXP(Z(i,j)*fac2(i))*(1.0_r8-fac3(i))
        END DO
      END DO
      RETURN
      END SUBROUTINE lmd_swfrac_tile
#else
      SUBROUTINE lmd_swfrac
      RETURN
      END SUBROUTINE lmd_swfrac
#endif

