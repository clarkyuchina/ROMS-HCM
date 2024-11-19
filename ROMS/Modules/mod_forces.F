#include "cppdefs.h"
      MODULE mod_forces
!
!svn $Id: mod_forces.F 921 2018-09-06 18:27:34Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Surface momentum stresses.                                          !
!                                                                      !
!  sustr        Kinematic surface momentum flux (wind stress) in       !
!                 the XI-direction (m2/s2) at horizontal U-points.     !
!  sustrG       Latest two-time snapshots of input "sustr" grided      !
!                 data used for interpolation.                         !
!  svstr        Kinematic surface momentum flux (wind stress) in       !
!                 the ETA-direction (m2/s2) at horizontal V-points.    !
!  svstrG       Latest two-time snapshots of input "svstr" grided      !
!                 data used for interpolation.                         !
!  Taux         Surface stress in the XI-direction at rho points       !
!                 from atm model.                                      !
!  Tauy         Surface stress in the ETA-direction at rho points      !
!                 from atm model.                                      !
!                                                                      !
!  Bottom momentum stresses.                                           !
!                                                                      !
!  bustr        Kinematic bottom momentum flux (bottom stress) in      !
!                 the XI-direction (m2/s2) at horizontal U-points.     !
!  bvstr        Kinematic bottom momentum flux (bottom stress) in      !
!                 ETA-direction (m2/s2) at horizontal V-points.        !
!                                                                      !
!  Surface wind induced waves.                                         !
!                                                                      !
!  Hwave        Surface wind induced wave height (m).                  !
!  HwaveG       Latest two-time snapshots of input "Hwave" grided      !
!                 data used for interpolation.                         !
!  Dwave        Surface wind induced mean wave direction (radians).    !
!  DwaveG       Latest two-time snapshots of input "Dwave" grided      !
!                 data used for interpolation.                         !
!  Dwavep       Surface wind induced peak wave direction (radians).    !
!  DwavepG      Latest two-time snapshots of input "Dwavep" grided     !
!                 data used for interpolation.                         !
!  Lwave        Mean surface wavelength read in from swan output       !
!  LwaveG       Latest two-time snapshots of input "Lwave" grided      !
!                 data used for interpolation.                         !
!  Lwavep       Peak surface wavelength read in from swan output       !
!  LwavepG      Latest two-time snapshots of input "Lwavep" grided     !
!                 data used for interpolation.                         !
!  Pwave_top    Wind induced surface wave period (s).                  !
!  Pwave_topG   Latest two-time snapshots of input "Pwave_top" grided  !
!                 data used for interpolation.                         !
!  Pwave_bot    Wind induced bottom wave period (s).                   !
!  Pwave_botG   Latest two-time snapshots of input "Pwave_bot" grided  !
!                 data used for interpolation.                         !
!  Uwave_rms    Bottom orbital velocity read in from swan output       !
!  Uwave_rmsG   Latest two-time snapshots of input "Uwave_rms" grided  !
!                 data used for interpolation.                         !
!  wave_dissip  Wave dissipation                                       !
!  wave_dissipG Latest two-time snapshots of input "wave_dissip"       !
!                 gridded data used for interpolation.                 !
!  Wave_break   Percent of wave breaking for use with roller model.    !
!  Wave_breakG  Latest two-time snapshots of input "wave_break"        !
!                 gridded data used for interpolation.                 !
!  Wave_ds      Wave directional spreading.                            !
!  Wave_qp      Wave spectrum peakedness.                              !
!                                                                      !
!  Solar shortwave radiation flux.                                     !
!                                                                      !
!  srflx        Kinematic surface shortwave solar radiation flux       !
!                 (Celsius m/s) at horizontal RHO-points               !
!  srflxG       Latest two-time snapshots of input "srflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Cloud fraction.                                                     !
!                                                                      !
!  cloud        Cloud fraction (percentage/100).                       !
!  cloudG       Latest two-time snapshots of input "cloud" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface heat fluxes, Atmosphere-Ocean bulk parameterization.        !
!                                                                      !
!  lhflx        Kinematic latent heat flux (degC m/s).                 !
!  lrflx        Kinematic longwave radiation (degC m/s).               !
!  snhflx       Kinematic surface net heat flux (degC m/s).            !
!  shflx        Kinematic sensible heat flux (degC m/s).               !
!                                                                      !
!  Surface air humidity.                                               !
!                                                                      !
!  Hair         Surface air specific (g/kg) or relative humidity       !
!                 (percentage).                                        !
!  HairG        Latest two-time snapshots of input "Hair" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface air pressure.                                               !
!                                                                      !
!  Pair         Surface air pressure (mb).                             !
!  PairG        Latest two-time snapshots of input "Pair" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface air temperature.                                            !
!                                                                      !
!  Tair         Surface air temperature (Celsius)                      !
!  TairG        Latest two-time snapshots of input "Tair" grided       !
!                 data used for interpolation.                         !
!  PotT         Surface air potential temperature (Kelvin)             !
!  Surface Winds.                                                      !
!                                                                      !
!  Uwind        Surface wind in the XI-direction (m/s) at              !
!                 horizontal RHO-points.                               !
!  UwindG       Latest two-time snapshots of input "Uwind" grided      !
!                 data used for interpolation.                         !
!  Vwind        Surface wind in the ETA-direction (m/s) at             !
!                 horizontal RHO-points.                               !
!  VwindG       Latest two-time snapshots of input "Vwind" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Rain fall rate.                                                     !
!                                                                      !
!  evap         Evaporation rate (kg/m2/s).                            !
!  rain         Rain fall rate (kg/m2/s).                              !
!  rainG        Latest two-time snapshots of input "rain" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Snow fall rate.                                                     !
!                                                                      !
!  snow         Snow fall rate (kg/m2/s).                              !
!  snowG        Latest two-time snapshots of input "snow" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface tracer fluxes.                                              !
!                                                                      !
!  stflx        Kinematic surface flux of tracer type variables        !
!                 (temperature: degC m/s; salinity: PSU m/s) at        !
!                 horizontal RHO-points.                               !
!  stflxG       Latest two-time snapshots of input "stflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Bottom tracer fluxes.                                               !
!                                                                      !
!  btflx        Kinematic bottom flux of tracer type variables         !
!                 (temperature: degC m/s; salinity: PSU m/s) at        !
!                horizontal RHO-points.                                !
!  btflxG       Latest two-time snapshots of input "btflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface heat flux correction.                                       !
!                                                                      !
!  dqdt         Kinematic surface net heat flux sensitivity to SST,    !
!                 d(Q)/d(SST), (m/s).                                  !
!  dqdtG        Latest two-time snapshots of input "dqdt" grided       !
!                 data used for interpolation.                         !
!  sst          Sea surface temperature (Celsius).                     !
!  sstG         Latest two-time snapshots of input "sst" grided        !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface freshwater flux correction.                                 !
!                                                                      !
!  sss          Sea surface salinity (PSU).                            !
!  sssG         Latest two-time snapshots of input "sss" grided        !
!                 data used for interpolation.                         !
!  sssflx       Sea surface salinity flux correction.                  !
!  sssflxG      Latest two-time snapshots of input "sssflx" grided     !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface spectral downwelling irradiance.                            !
!                                                                      !
!  SpecIr       Spectral irradiance (NBands) from 400-700 nm at        !
!                 5 nm bandwidth.                                      !
!  avcos        Cosine of average zenith angle of downwelling          !
!                 spectral photons.                                    !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_FORCES
!
!  Nonlinear model state.
!
#ifdef NCEP_FLUXES
          real(r8), pointer :: sustr(:,:)
          real(r8), pointer :: svstr(:,:)
          real(r8), pointer :: nustr(:,:)
          real(r8), pointer :: nvstr(:,:)
          real(r8), pointer :: nustrG(:,:,:)
          real(r8), pointer :: nvstrG(:,:,:)
          real(r8), pointer :: bustr(:,:)
          real(r8), pointer :: bvstr(:,:)
          real(r8), pointer :: srflx(:,:)
          real(r8), pointer :: srflxG(:,:,:)
          real(r8), pointer :: cloud(:,:)
          real(r8), pointer :: cloudG(:,:,:)
          real(r8), pointer :: lhflx(:,:)
          real(r8), pointer :: lhflxG(:,:,:)
          real(r8), pointer :: lrflx(:,:)
          real(r8), pointer :: lrflxG(:,:,:)
          real(r8), pointer :: shflx(:,:)
          real(r8), pointer :: shflxG(:,:,:)
          real(r8), pointer :: Pair(:,:)
          real(r8), pointer :: PairG(:,:,:)
          real(r8), pointer :: rain(:,:)
          real(r8), pointer :: rainG(:,:,:)
# ifdef RUNOFF
          real(r8), pointer :: runoff(:,:)
          real(r8), pointer :: runoffG(:,:,:)
# endif
          real(r8), pointer :: skt(:,:)
          real(r8), pointer :: sktG(:,:,:)
          real(r8), pointer :: icec(:,:)
          real(r8), pointer :: icecG(:,:,:)
          real(r8), pointer :: snow_n(:,:)
          real(r8), pointer :: p_e_n(:,:)
          real(r8), pointer :: wg2_d(:,:)
          real(r8), pointer :: cd_d(:,:)
          real(r8), pointer :: ch_d(:,:)
          real(r8), pointer :: ce_d(:,:)
          real(r8), pointer :: wg2_m(:,:)
          real(r8), pointer :: cd_m(:,:)
          real(r8), pointer :: ch_m(:,:)
          real(r8), pointer :: ce_m(:,:)
          real(r8), pointer :: rhoa_n(:,:)
          real(r8), pointer :: cawdir(:,:)
          real(r8), pointer :: qao_n(:,:)
          real(r8), pointer :: qai_n(:,:)
          real(r8), pointer :: qswi_n(:,:)
          real(r8), pointer :: tau_awx_n(:,:)
          real(r8), pointer :: tau_awy_n(:,:)
          real(r8), pointer :: tau_aix_n(:,:)
          real(r8), pointer :: tau_aiy_n(:,:)
          real(r8), pointer :: sustr_aw(:,:)
          real(r8), pointer :: svstr_aw(:,:)
#else

          real(r8), pointer :: sustr(:,:)
          real(r8), pointer :: svstr(:,:)
# if !defined ANA_SMFLUX && !defined BULK_FLUXES && !defined \
      BULK_FLUXES2D || defined NL_BULK_FLUXES
          real(r8), pointer :: sustrG(:,:,:)
          real(r8), pointer :: svstrG(:,:,:)
# endif
# ifdef ATM2OCN_FLUXES
          real(r8), pointer :: Taux(:,:)
          real(r8), pointer :: Tauy(:,:)
# endif
# ifdef ADJUST_WSTRESS
          real(r8), pointer :: ustr(:,:,:,:)
          real(r8), pointer :: vstr(:,:,:,:)
# endif
          real(r8), pointer :: bustr(:,:)
          real(r8), pointer :: bvstr(:,:)
#endif

#ifdef WAVES_DIR
          real(r8), pointer :: Dwave(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: DwaveG(:,:,:)
# endif
#endif

#ifdef WAVES_DIRP
          real(r8), pointer :: Dwavep(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: DwavepG(:,:,:)
# endif
#endif

#ifdef WAVES_HEIGHT
          real(r8), pointer :: Hwave(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: HwaveG(:,:,:)
# endif
#endif

#ifdef WAVES_LENGTH
          real(r8), pointer :: Lwave(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: LwaveG(:,:,:)
# endif
#endif

#ifdef WAVES_LENGTHP
          real(r8), pointer :: Lwavep(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: LwavepG(:,:,:)
# endif
#endif

#ifdef WAVES_TOP_PERIOD
          real(r8), pointer :: Pwave_top(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Pwave_topG(:,:,:)
# endif
#endif

#ifdef WAVES_BOT_PERIOD
          real(r8), pointer :: Pwave_bot(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Pwave_botG(:,:,:)
# endif
#endif

#if defined BBL_MODEL || defined WAVES_OCEAN || \
    defined SED_BEDLOAD_VANDERA
          real(r8), pointer :: Uwave_rms(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Uwave_rmsG(:,:,:)
# endif
#endif

#if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
    defined WAVES_DISS || defined WDISS_INWAVE
          real(r8), pointer :: Dissip_break(:,:)
          real(r8), pointer :: Dissip_wcap(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Dissip_breakG(:,:,:)
          real(r8), pointer :: Dissip_wcapG(:,:,:)
# endif
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
          real(r8), pointer :: Dissip_fric(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Dissip_fricG(:,:,:)
# endif
#endif

#if defined WEC_ROLLER
          real(r8), pointer :: Dissip_roller(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Dissip_rollerG(:,:,:)
# endif
#endif

#if defined ROLLER_SVENDSEN
          real(r8), pointer :: wave_break(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: wave_breakG(:,:,:)
# endif
#endif

#if defined WAVES_DSPR
          real(r8), pointer :: Wave_ds(:,:)
          real(r8), pointer :: Wave_qp(:,:)
# ifndef ANA_WWAVE
          real(r8), pointer :: Wave_dsG(:,:,:)
          real(r8), pointer :: Wave_qpG(:,:,:)
# endif
#endif

#if defined BULK_FLUXES || defined BULK_FLUXES2D || defined ECOSIM || \
    defined CCSM_FLUXES2D || defined ATM2OCN_FLUXES || \
    defined SPECTRAL_LIGHT
          real(r8), pointer :: Uwind(:,:)
          real(r8), pointer :: Vwind(:,:)
# ifndef ANA_WINDS
          real(r8), pointer :: UwindG(:,:,:)
          real(r8), pointer :: VwindG(:,:,:)
# endif
#endif

#if defined BULK_FLUXES || defined ECOSIM || defined CCSM_FLUXES2D || \
    defined SPECTRAL_LIGHT || (defined SHORTWAVE && defined ANA_SRFLUX)
          real(r8), pointer :: Hair(:,:)
# ifndef ANA_HUMIDITY
          real(r8), pointer :: HairG(:,:,:)
# endif
          real(r8), pointer :: Tair(:,:)
# ifndef ANA_TAIR
          real(r8), pointer :: TairG(:,:,:)
# endif
#endif

#if defined WEC_ROLLER
          real(r8), pointer :: rollA(:,:)
#endif
#ifdef CICE_MODEL
          real(r8), pointer :: PotT(:,:)
          real(r8), pointer :: LW_down(:,:)
          real(r8), pointer :: SW_down(:,:)
#endif
#ifdef ICE_DIAGS
          real(r8), pointer :: LW_down(:,:)
          real(r8), pointer :: SW_down(:,:)
          real(r8), pointer :: lat_ice(:,:)
          real(r8), pointer :: sens_ice(:,:)
          real(r8), pointer :: LW_up_ice(:,:)
          real(r8), pointer :: SW_up_ice(:,:)
          real(r8), pointer :: saltflux_ice(:,:)
          real(r8), pointer :: saltflux_ocean(:,:)
          real(r8), pointer :: snoice(:,:)
#endif
#if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS || \
    defined SPECTRAL_LIGHT  || defined ATM2OCN_FLUXES ||          \
    defined CCSM_FLUXES2D
          real(r8), pointer :: Pair(:,:)
# ifndef ANA_PAIR
          real(r8), pointer :: PairG(:,:,:)
# endif
#endif

#ifdef SOLVE3D

# ifdef SHORTWAVE
          real(r8), pointer :: srflx(:,:)
#  ifndef ANA_SRFLUX
          real(r8), pointer :: srflxG(:,:,:)
#  endif
#  ifdef ALBEDO
          real(r8), pointer :: albedo(:,:)
#   ifdef ICE_MODEL
          real(r8), pointer :: albedo_ice(:,:)
#   endif
#  endif
#  ifdef ALBEDO_FILE
          real(r8), pointer :: albedoG(:,:,:)
#  endif
#  if defined ALBEDO_CLOUD
          real(r8), pointer :: cawdir(:,:)
#  endif
# endif

# if defined RED_TIDE && defined DAILY_SHORTWAVE
          real(r8), pointer :: srflx_avg(:,:)
          real(r8), pointer :: srflxG_avg(:,:,:)
# endif

# if defined ICE_BULK_FLUXES && !defined NCEP_FLUXES
          real(r8), pointer :: sustr_aw(:,:)
          real(r8), pointer :: svstr_aw(:,:)
          real(r8), pointer :: tau_aix_n(:,:)
          real(r8), pointer :: tau_aiy_n(:,:)
          real(r8), pointer :: qai_n(:,:)
          real(r8), pointer :: qi_o_n(:,:)
          real(r8), pointer :: SW_thru_ice(:,:)
          real(r8), pointer :: qswi_n(:,:)
          real(r8), pointer :: qao_n(:,:)
          real(r8), pointer :: qio_n(:,:)
          real(r8), pointer :: qi2_n(:,:)
          real(r8), pointer :: snow_n(:,:)
          real(r8), pointer :: p_e_n(:,:)
# endif 

# ifdef CLOUDS
          real(r8), pointer :: cloud(:,:)
#  ifndef ANA_CLOUD
          real(r8), pointer :: cloudG(:,:,:)
#  endif
# endif
# ifdef BULK_FLUXES
          real(r8), pointer :: lhflx(:,:)
          real(r8), pointer :: lrflx(:,:)
#  ifndef LONGWAVE
          real(r8), pointer :: lrflxG(:,:,:)
#  endif
          real(r8), pointer :: shflx(:,:)
          real(r8), pointer :: snhflx(:,:)
# endif
# if defined ATM2OCN_FLUXES
          real(r8), pointer :: lhflx(:,:)
          real(r8), pointer :: lrflx(:,:)
          real(r8), pointer :: shflx(:,:)
          real(r8), pointer :: snhflx(:,:)
# endif

# if defined BULK_FLUXES || defined ATM2OCN_FLUXES
          real(r8), pointer :: rain(:,:)
#  ifdef SNOWFALL
          real(r8), pointer :: snow(:,:)
#  endif
#  ifndef ANA_RAIN
          real(r8), pointer :: rainG(:,:,:)
#  ifdef SNOWFALL
          real(r8), pointer :: snowG(:,:,:)
#  endif
#  endif
#  ifdef EMINUSP
          real(r8), pointer :: EminusP(:,:)
          real(r8), pointer :: evap(:,:)
#  endif
#  ifdef RUNOFF
          real(r8), pointer :: runoff(:,:)
          real(r8), pointer :: runoffG(:,:,:)
#  endif
# endif

          real(r8), pointer :: stflx(:,:,:)
# if defined PERFECT_RESTART && defined ICE_MODEL
          real(r8), pointer :: stflx_save(:,:,:)
          real(r8), pointer :: sustr_save(:,:)
          real(r8), pointer :: svstr_save(:,:)
# endif
# if !defined ANA_STFLUX || !defined ANA_SSFLUX || \
     !defined ANA_SPFLUX
          real(r8), pointer :: stflxG(:,:,:,:)
# endif
# ifdef ADJUST_STFLUX
          real(r8), pointer :: tflux(:,:,:,:,:)
# endif
          real(r8), pointer :: btflx(:,:,:)
# if !defined ANA_BTFLUX || !defined ANA_BSFLUX || \
     !defined ANA_BPFLUX
          real(r8), pointer :: btflxG(:,:,:,:)
# endif

# if defined QCORRECTION
          real(r8), pointer :: sst(:,:)
          real(r8), pointer :: dqdt(:,:)
#  ifndef ANA_SST
          real(r8), pointer :: sstG(:,:,:)
          real(r8), pointer :: dqdtG(:,:,:)
#  endif
# endif

# ifdef SSH_NUDGING
          real(r8), pointer :: sshC(:,:)
          real(r8), pointer :: sshCG(:,:,:)
# endif

# ifdef HCM_COUPLING
#  ifdef HCM_SST_FILTER
#    ifdef HCM_ROSSBY_FILTER
          real(r8), pointer :: Rrob(:,:)
#    endif
          real(r8), pointer :: sstF(:,:)
          real(r8), pointer :: sstH(:,:)
          real(r8), pointer :: sstL(:,:)
#  endif
#  ifndef HCM_CLIM_SPIN
          real(r8), pointer :: sstC(:,:)
          real(r8), pointer :: sstA(:,:)
#  endif
#  ifdef HCM_WIND_STRESS
          real(r8), pointer :: sustrC(:,:)
          real(r8), pointer :: svstrC(:,:)
          real(r8), pointer :: sustrA(:,:)
          real(r8), pointer :: svstrA(:,:)
#    ifndef HCM_AI_WIND
          real(r8), pointer :: dudtM1(:,:)
          real(r8), pointer :: dudtM2(:,:)
          real(r8), pointer :: dudtM3(:,:)
          real(r8), pointer :: dudtM4(:,:)
          real(r8), pointer :: dudtM5(:,:)
          real(r8), pointer :: dvdtM1(:,:)
          real(r8), pointer :: dvdtM2(:,:)
          real(r8), pointer :: dvdtM3(:,:)
          real(r8), pointer :: dvdtM4(:,:)
          real(r8), pointer :: dvdtM5(:,:)
          real(r8), pointer :: dudtW1(:,:)
          real(r8), pointer :: dudtW2(:,:)
          real(r8), pointer :: dudtW3(:,:)
          real(r8), pointer :: dudtW4(:,:)
          real(r8), pointer :: dudtW5(:,:)
          real(r8), pointer :: dvdtW1(:,:)
          real(r8), pointer :: dvdtW2(:,:)
          real(r8), pointer :: dvdtW3(:,:)
          real(r8), pointer :: dvdtW4(:,:)
          real(r8), pointer :: dvdtW5(:,:)
#    endif
#  endif
#  ifdef HCM_EMINP_FORCE
          real(r8), pointer :: fwfC(:,:)
          real(r8), pointer :: fwfA(:,:)
#    ifndef HCM_AI_FWF
          real(r8), pointer :: dwdtM1(:,:)
          real(r8), pointer :: dwdtM2(:,:)
          real(r8), pointer :: dwdtM3(:,:)
          real(r8), pointer :: dwdtM4(:,:)
          real(r8), pointer :: dwdtM5(:,:)
          real(r8), pointer :: dwdtW1(:,:)
          real(r8), pointer :: dwdtW2(:,:)
          real(r8), pointer :: dwdtW3(:,:)
          real(r8), pointer :: dwdtW4(:,:)
          real(r8), pointer :: dwdtW5(:,:)
#    endif
#  endif
#  ifdef TIW_COUPLING
          real(r8), pointer :: sustrTIW(:,:)
          real(r8), pointer :: svstrTIW(:,:)
#    ifndef HCM_AI_TIW
          real(r8), pointer :: tiwU01(:,:)
          real(r8), pointer :: tiwV01(:,:)
          real(r8), pointer :: tiwT01(:,:)
          real(r8), pointer :: tiwU02(:,:)
          real(r8), pointer :: tiwV02(:,:)
          real(r8), pointer :: tiwT02(:,:)
          real(r8), pointer :: tiwU03(:,:)
          real(r8), pointer :: tiwV03(:,:)
          real(r8), pointer :: tiwT03(:,:)
          real(r8), pointer :: tiwU04(:,:)
          real(r8), pointer :: tiwV04(:,:)
          real(r8), pointer :: tiwT04(:,:)
          real(r8), pointer :: tiwU05(:,:)
          real(r8), pointer :: tiwV05(:,:)
          real(r8), pointer :: tiwT05(:,:)
          real(r8), pointer :: tiwU06(:,:)
          real(r8), pointer :: tiwV06(:,:)
          real(r8), pointer :: tiwT06(:,:)
          real(r8), pointer :: tiwU07(:,:)
          real(r8), pointer :: tiwV07(:,:)
          real(r8), pointer :: tiwT07(:,:)
          real(r8), pointer :: tiwU08(:,:)
          real(r8), pointer :: tiwV08(:,:)
          real(r8), pointer :: tiwT08(:,:)
          real(r8), pointer :: tiwU09(:,:)
          real(r8), pointer :: tiwV09(:,:)
          real(r8), pointer :: tiwT09(:,:)
          real(r8), pointer :: tiwU10(:,:)
          real(r8), pointer :: tiwV10(:,:)
          real(r8), pointer :: tiwT10(:,:)
#    endif
#  endif
#  ifdef HCM_BULK_FLUX
          real(r8), pointer :: TairC(:,:)
          real(r8), pointer :: PairC(:,:)
          real(r8), pointer :: QairC(:,:)
          real(r8), pointer :: SWR(:,:)
          real(r8), pointer :: LWR(:,:)
          real(r8), pointer :: LHF(:,:)
          real(r8), pointer :: SHF(:,:)
#  endif
#  ifdef HCM_SST_FILTER
#    ifdef HCM_ROSSBY_FILTER
          real(r8), pointer :: RrobG(:,:,:)
#    endif
#  endif
#  ifndef HCM_CLIM_SPIN
          real(r8), pointer :: sstCG(:,:,:)
          real(r8), pointer :: sstAG(:,:,:)
#  endif
#  ifdef HCM_WIND_STRESS
          real(r8), pointer :: sustrCG(:,:,:)
          real(r8), pointer :: svstrCG(:,:,:)
          real(r8), pointer :: sustrAG(:,:,:)
          real(r8), pointer :: svstrAG(:,:,:)
#    ifndef HCM_AI_WIND
          real(r8), pointer :: dudtM1G(:,:,:)
          real(r8), pointer :: dudtM2G(:,:,:)
          real(r8), pointer :: dudtM3G(:,:,:)
          real(r8), pointer :: dudtM4G(:,:,:)
          real(r8), pointer :: dudtM5G(:,:,:)
          real(r8), pointer :: dvdtM1G(:,:,:)
          real(r8), pointer :: dvdtM2G(:,:,:)
          real(r8), pointer :: dvdtM3G(:,:,:)
          real(r8), pointer :: dvdtM4G(:,:,:)
          real(r8), pointer :: dvdtM5G(:,:,:)
          real(r8), pointer :: dudtW1G(:,:,:)
          real(r8), pointer :: dudtW2G(:,:,:)
          real(r8), pointer :: dudtW3G(:,:,:)
          real(r8), pointer :: dudtW4G(:,:,:)
          real(r8), pointer :: dudtW5G(:,:,:)
          real(r8), pointer :: dvdtW1G(:,:,:)
          real(r8), pointer :: dvdtW2G(:,:,:)
          real(r8), pointer :: dvdtW3G(:,:,:)
          real(r8), pointer :: dvdtW4G(:,:,:)
          real(r8), pointer :: dvdtW5G(:,:,:)
#    endif
#  endif
#  ifdef HCM_EMINP_FORCE
          real(r8), pointer :: fwfCG(:,:,:)
          real(r8), pointer :: fwfAG(:,:,:)
#    ifndef HCM_AI_FWF
          real(r8), pointer :: dwdtM1G(:,:,:)
          real(r8), pointer :: dwdtM2G(:,:,:)
          real(r8), pointer :: dwdtM3G(:,:,:)
          real(r8), pointer :: dwdtM4G(:,:,:)
          real(r8), pointer :: dwdtM5G(:,:,:)
          real(r8), pointer :: dwdtW1G(:,:,:)
          real(r8), pointer :: dwdtW2G(:,:,:)
          real(r8), pointer :: dwdtW3G(:,:,:)
          real(r8), pointer :: dwdtW4G(:,:,:)
          real(r8), pointer :: dwdtW5G(:,:,:)
#    endif
#  endif
#  ifdef TIW_COUPLING
          real(r8), pointer :: sustrTIWG(:,:,:)
          real(r8), pointer :: svstrTIWG(:,:,:)
#    ifndef HCM_AI_TIW
          real(r8), pointer :: tiwU01G(:,:,:)
          real(r8), pointer :: tiwV01G(:,:,:)
          real(r8), pointer :: tiwT01G(:,:,:)
          real(r8), pointer :: tiwU02G(:,:,:)
          real(r8), pointer :: tiwV02G(:,:,:)
          real(r8), pointer :: tiwT02G(:,:,:)
          real(r8), pointer :: tiwU03G(:,:,:)
          real(r8), pointer :: tiwV03G(:,:,:)
          real(r8), pointer :: tiwT03G(:,:,:)
          real(r8), pointer :: tiwU04G(:,:,:)
          real(r8), pointer :: tiwV04G(:,:,:)
          real(r8), pointer :: tiwT04G(:,:,:)
          real(r8), pointer :: tiwU05G(:,:,:)
          real(r8), pointer :: tiwV05G(:,:,:)
          real(r8), pointer :: tiwT05G(:,:,:)
          real(r8), pointer :: tiwU06G(:,:,:)
          real(r8), pointer :: tiwV06G(:,:,:)
          real(r8), pointer :: tiwT06G(:,:,:)
          real(r8), pointer :: tiwU07G(:,:,:)
          real(r8), pointer :: tiwV07G(:,:,:)
          real(r8), pointer :: tiwT07G(:,:,:)
          real(r8), pointer :: tiwU08G(:,:,:)
          real(r8), pointer :: tiwV08G(:,:,:)
          real(r8), pointer :: tiwT08G(:,:,:)
          real(r8), pointer :: tiwU09G(:,:,:)
          real(r8), pointer :: tiwV09G(:,:,:)
          real(r8), pointer :: tiwT09G(:,:,:)
          real(r8), pointer :: tiwU10G(:,:,:)
          real(r8), pointer :: tiwV10G(:,:,:)
          real(r8), pointer :: tiwT10G(:,:,:)
#    endif
#  endif
#  ifdef HCM_BULK_FLUX
          real(r8), pointer :: TairCG(:,:,:)
          real(r8), pointer :: PairCG(:,:,:)
          real(r8), pointer :: QairCG(:,:,:)
          real(r8), pointer :: SWRG(:,:,:)
          real(r8), pointer :: LWRG(:,:,:)
          real(r8), pointer :: LHFG(:,:,:)
          real(r8), pointer :: SHFG(:,:,:)
#  endif
# endif

# if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
          real(r8), pointer :: sss(:,:)
#  ifndef ANA_SSS
          real(r8), pointer :: sssG(:,:,:)
#  endif
# endif

# if defined SSSFLX
          real(r8), pointer :: sssflx(:,:)
#  ifndef ANA_SSSFLX
          real(r8), pointer :: sssflxG(:,:,:)
#  endif
# endif

# ifdef FASTICE_CLIMATOLOGY
          real(r8), pointer :: fastice_clm(:,:)
#  ifndef ANA_FASTICE
          real(r8), pointer :: fastice_clmG(:,:,:)
#  endif
# endif

# if defined ECOSIM  || defined SPECTRAL_LIGHT
          real(r8), pointer :: SpecIr(:,:,:)
          real(r8), pointer :: avcos(:,:,:)
# endif
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
          real(r8), pointer :: tl_sustr(:,:)
          real(r8), pointer :: tl_svstr(:,:)
# ifdef ADJUST_WSTRESS
          real(r8), pointer :: tl_ustr(:,:,:,:)
          real(r8), pointer :: tl_vstr(:,:,:,:)
# endif
          real(r8), pointer :: tl_bustr(:,:)
          real(r8), pointer :: tl_bvstr(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: tl_stflx(:,:,:)
          real(r8), pointer :: tl_btflx(:,:,:)
#  ifdef ADJUST_STFLUX
          real(r8), pointer :: tl_tflux(:,:,:,:,:)
#  endif
#  ifdef SHORTWAVE
          real(r8), pointer :: tl_srflx(:,:)
#  endif
#  ifdef BULK_FLUXES
          real(r8), pointer :: tl_lhflx(:,:)
          real(r8), pointer :: tl_lrflx(:,:)
          real(r8), pointer :: tl_shflx(:,:)
#   ifdef EMINUSP
          real(r8), pointer :: tl_evap(:,:)
#   endif
#  endif
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
          real(r8), pointer :: ad_sustr(:,:)
          real(r8), pointer :: ad_svstr(:,:)
# ifdef ADJUST_WSTRESS
          real(r8), pointer :: ad_ustr(:,:,:,:)
          real(r8), pointer :: ad_vstr(:,:,:,:)
# endif
          real(r8), pointer :: ad_bustr(:,:)
          real(r8), pointer :: ad_bvstr(:,:)
          real(r8), pointer :: ad_bustr_sol(:,:)
          real(r8), pointer :: ad_bvstr_sol(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: ad_stflx(:,:,:)
          real(r8), pointer :: ad_btflx(:,:,:)
#  ifdef ADJUST_STFLUX
          real(r8), pointer :: ad_tflux(:,:,:,:,:)
#  endif
#  ifdef SHORTWAVE
          real(r8), pointer :: ad_srflx(:,:)
#  endif
#  ifdef BULK_FLUXES
          real(r8), pointer :: ad_lhflx(:,:)
          real(r8), pointer :: ad_lrflx(:,:)
          real(r8), pointer :: ad_shflx(:,:)
#   ifdef EMINUSP
          real(r8), pointer :: ad_evap(:,:)
#   endif
#  endif
# endif
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Working arrays to store adjoint impulse forcing, background error
!  covariance, background-error standard deviations, or descent
!  conjugate vectors (directions).
!
# if defined FOUR_DVAR || defined IMPULSE
#  ifdef ADJUST_WSTRESS
          real(r8), pointer :: b_sustr(:,:)
          real(r8), pointer :: b_svstr(:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
          real(r8), pointer :: b_stflx(:,:,:)
#  endif
#  ifdef FOUR_DVAR
#   ifdef ADJUST_WSTRESS
          real(r8), pointer :: d_sustr(:,:,:)
          real(r8), pointer :: d_svstr(:,:,:)
          real(r8), pointer :: e_sustr(:,:)
          real(r8), pointer :: e_svstr(:,:)
#   endif
#   if defined ADJUST_STFLUX && defined SOLVE3D
          real(r8), pointer :: d_stflx(:,:,:,:)
          real(r8), pointer :: e_stflx(:,:,:)
#   endif
#  endif
# endif
#endif

        END TYPE T_FORCES

        TYPE (T_FORCES), allocatable :: FORCES(:)

      CONTAINS

      SUBROUTINE allocate_forces (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
#if defined BIOLOGY
      USE mod_biology
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      real(r8) :: size2d
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( FORCES(Ngrids) )
!
!  Set horizontal array size.
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Nonlinear model state
!
#ifdef NCEP_FLUXES
      allocate ( FORCES(ng) % sustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % svstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % nustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % nvstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % nustrG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % nvstrG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % bustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % bvstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d
      allocate ( FORCES(ng) % btflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d
      allocate ( FORCES(ng) % srflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % srflxG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % cloud(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % cloudG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % lhflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % lhflxG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % lrflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % lrflxG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % shflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % shflxG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % Pair(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % PairG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % rain(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % rainG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# ifdef RUNOFF
      allocate ( FORCES(ng) % runoff(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % runoffG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
      allocate ( FORCES(ng) % skt(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % sktG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % icec(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % icecG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % wg2_d(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % cd_d(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % ch_d(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % ce_d(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % wg2_m(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % cd_m(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % ch_m(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % ce_m(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % rhoa_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % cawdir(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % snow_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % p_e_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % qao_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % qio_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % qi2_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % qai_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % qswi_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tau_awx_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tau_awy_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tau_aix_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tau_aiy_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % sustr_aw(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % svstr_aw(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#else
      allocate ( FORCES(ng) % sustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % svstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifdef ATM2OCN_FLUXES
      allocate ( FORCES(ng) % Taux(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % Tauy(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif
# if !defined ANA_SMFLUX && !defined BULK_FLUXES && !defined \
     BULK_FLUXES2D || defined NL_BULK_FLUXES
      allocate ( FORCES(ng) % sustrG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

      allocate ( FORCES(ng) % svstrG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
# ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d

      allocate ( FORCES(ng) % vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d
# endif
      allocate ( FORCES(ng) % bustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % bvstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#endif

#ifdef WAVES_DIR
      allocate ( FORCES(ng) % Dwave(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % DwaveG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_DIRP
      allocate ( FORCES(ng) % Dwavep(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % DwavepG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_HEIGHT
      allocate ( FORCES(ng) % Hwave(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % HwaveG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_LENGTH
      allocate ( FORCES(ng) % Lwave(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % LwaveG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_LENGTHP
      allocate ( FORCES(ng) % Lwavep(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % LwavepG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_TOP_PERIOD
      allocate ( FORCES(ng) % Pwave_top(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Pwave_topG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#ifdef WAVES_BOT_PERIOD
      allocate ( FORCES(ng) % Pwave_bot(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Pwave_botG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined BBL_MODEL || defined WAVES_OCEAN
      allocate ( FORCES(ng) % Uwave_rms(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Uwave_rmsG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
    defined WAVES_DISS || defined WDISS_INWAVE
      allocate ( FORCES(ng) % Dissip_break(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % Dissip_wcap(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Dissip_breakG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % Dissip_wcapG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
      allocate ( FORCES(ng) % Dissip_fric(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Dissip_fricG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined WEC_ROLLER
      allocate ( FORCES(ng) % Dissip_roller(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Dissip_rollerG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined ROLLER_SVENDSEN
      allocate ( FORCES(ng) % wave_break(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Wave_breakG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined WAVES_DSPR
      allocate ( FORCES(ng) % Wave_ds(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % Wave_qp(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WWAVE
      allocate ( FORCES(ng) % Wave_dsG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % Wave_qpG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined WEC_ROLLER
      allocate ( FORCES(ng) % rollA(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#endif
#if defined BULK_FLUXES || defined BULK_FLUXES2D || defined ECOSIM || \
    defined SPECTRAL_LIGHT || defined ATM2OCN_FLUXES
      allocate ( FORCES(ng) % Uwind(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % Vwind(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_WINDS
      allocate ( FORCES(ng) % UwindG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % VwindG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif

#if defined BULK_FLUXES || defined ECOSIM || defined CCSM_FLUXES2D ||  \
    defined SPECTRAL_LIGHT || (defined SHORTWAVE && defined ANA_SRFLUX)
      allocate ( FORCES(ng) % Hair(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_HUMIDITY
      allocate ( FORCES(ng) % HairG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif

      allocate ( FORCES(ng) % Tair(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_TAIR
      allocate ( FORCES(ng) % TairG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif
#ifdef CICE_MODEL
      allocate ( FORCES(ng) % PotT(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % LW_down(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % SW_down(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#endif

# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS || \
     defined SPECTRAL_LIGHT || defined ATM2OCN_FLUXES || defined CCSM_FLUXES2D
      allocate ( FORCES(ng) % Pair(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifndef ANA_PAIR
      allocate ( FORCES(ng) % PairG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
#endif
#ifdef CICE_MODEL
      allocate ( FORCES(ng) % PotT(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % LW_down(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % SW_down(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#endif
#ifdef ICE_DIAGS
      allocate ( FORCES(ng) % LW_down(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % SW_down(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % lat_ice(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % sens_ice(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % LW_up_ice(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % SW_up_ice(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % saltflux_ice(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % saltflux_ocean(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % snoice(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#endif

#ifdef SOLVE3D

# ifdef SHORTWAVE
      allocate ( FORCES(ng) % srflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_SRFLUX
      allocate ( FORCES(ng) % srflxG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
#  ifdef ALBEDO
      allocate ( FORCES(ng) % albedo(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#   ifdef ICE_MODEL
      allocate ( FORCES(ng) % albedo_ice(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#   endif
#  endif
#  ifdef ALBEDO_FILE
      allocate ( FORCES(ng) % albedoG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
#  if defined ALBEDO_CLOUD
      allocate ( FORCES(ng) % cawdir(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
# endif
# if defined RED_TIDE && defined DAILY_SHORTWAVE
      allocate ( FORCES(ng) % srflx_avg(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % srflxG_avg(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif
# if defined ICE_BULK_FLUXES && !defined NCEP_FLUXES
      allocate ( FORCES(ng) % sustr_aw(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % svstr_aw(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tau_aix_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tau_aiy_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % qai_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % qi_o_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % SW_thru_ice(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % qswi_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % qao_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % qio_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % qi2_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % snow_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % p_e_n(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif

# ifdef CLOUDS
      allocate ( FORCES(ng) % cloud(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_CLOUD
      allocate ( FORCES(ng) % cloudG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
# endif

# ifdef BULK_FLUXES
      allocate ( FORCES(ng) % lhflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % lrflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef LONGWAVE
      allocate ( FORCES(ng) % lrflxG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif

      allocate ( FORCES(ng) % shflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % snhflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif
# ifdef ATM2OCN_FLUXES
      allocate ( FORCES(ng) % lhflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % lrflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % shflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % snhflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif

# if defined BULK_FLUXES || defined ATM2OCN_FLUXES
      allocate ( FORCES(ng) % rain(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  ifdef SNOWFALL
      allocate ( FORCES(ng) % snow(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  ifndef ANA_RAIN
      allocate ( FORCES(ng) % rainG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  ifdef SNOWFALL
      allocate ( FORCES(ng) % snowG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
#  endif
#  ifdef EMINUSP
      allocate ( FORCES(ng) % EminusP(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % evap(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  ifdef RUNOFF
      allocate ( FORCES(ng) % runoff(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % runoffG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
# endif

      allocate ( FORCES(ng) % stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d
# if defined PERFECT_RESTART && defined ICE_MODEL
      allocate ( FORCES(ng) % stflx_save(LBi:UBi,LBj:UBj,NAT) )
      Dmem(ng)=Dmem(ng)+REAL(NAT,r8)*size2d
      allocate ( FORCES(ng) % sustr_save(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % svstr_save(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif

# if !defined ANA_STFLUX || !defined ANA_SSFLUX || \
     !defined ANA_SPFLUX
      allocate ( FORCES(ng) % stflxG(LBi:UBi,LBj:UBj,2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(NT(ng),r8)*size2d

# endif
# ifdef ADJUST_STFLUX
      allocate ( FORCES(ng) % tflux(LBi:UBi,LBj:UBj,nfrec(ng),          &
     &                              2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(nfrec(ng)*NT(ng),r8)*size2d

# endif
      allocate ( FORCES(ng) % btflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d

# if !defined ANA_BTFLUX || !defined ANA_BSFLUX || \
     !defined ANA_BPFLUX
      allocate ( FORCES(ng) % btflxG(LBi:UBi,LBj:UBj,2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(NT(ng),r8)*size2d
# endif

# ifdef QCORRECTION
      allocate ( FORCES(ng) % sst(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dqdt(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  ifndef ANA_SST
      allocate ( FORCES(ng) % sstG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dqdtG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
# endif

# ifdef SSH_NUDGING
      allocate ( FORCES(ng) % sshC(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % sshCG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
# endif

# ifdef HCM_COUPLING
#  ifdef HCM_SST_FILTER
#   ifdef HCM_ROSSBY_FILTER
      allocate ( FORCES(ng) % Rrob(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#   endif
      allocate ( FORCES(ng) % sstF(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % sstH(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % sstL(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  ifndef HCM_CLIM_SPIN
      allocate ( FORCES(ng) % sstC(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % sstA(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  ifdef HCM_WIND_STRESS
      allocate ( FORCES(ng) % sustrC(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % svstrC(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % sustrA(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % svstrA(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#    ifndef HCM_AI_WIND
      allocate ( FORCES(ng) % dudtM1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dudtM2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dudtM3(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dudtM4(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dudtM5(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dvdtM1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dvdtM2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dvdtM3(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dvdtM4(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dvdtM5(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dudtW1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dudtW2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dudtW3(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dudtW4(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dudtW5(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dvdtW1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dvdtW2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dvdtW3(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dvdtW4(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dvdtW5(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#    endif
#  endif
#  ifdef HCM_EMINP_FORCE  
      allocate ( FORCES(ng) % fwfC(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % fwfA(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#    ifndef HCM_AI_FWF
      allocate ( FORCES(ng) % dwdtM1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dwdtM2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dwdtM3(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dwdtM4(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dwdtM5(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dwdtW1(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dwdtW2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dwdtW3(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dwdtW4(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % dwdtW5(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#    endif
#  endif
#  ifdef TIW_COUPLING
      allocate ( FORCES(ng) % sustrTIW(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % svstrTIW(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#    ifndef HCM_AI_FWF
      allocate ( FORCES(ng) % tiwU01(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwV01(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwT01(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwU02(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwV02(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwT02(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwU03(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwV03(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwT03(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwU04(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwV04(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwT04(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwU05(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwV05(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwT05(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwU06(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwV06(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwT06(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwU07(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwV07(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwT07(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwU08(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwV08(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwT08(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwU09(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwV09(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwT09(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwU10(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwV10(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % tiwT10(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#    endif
#  endif
#  ifdef HCM_BULK_FLUX
      allocate ( FORCES(ng) % TairC(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % PairC(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % QairC(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % SWR(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % LWR(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % LHF(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( FORCES(ng) % SHF(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif

#  ifdef HCM_SST_FILTER
#   ifdef HCM_ROSSBY_FILTER
      allocate ( FORCES(ng) % RrobG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#   endif
#  endif
#  ifndef HCM_CLIM_SPIN
      allocate ( FORCES(ng) % sstCG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % sstAG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
#  ifdef HCM_WIND_STRESS
      allocate ( FORCES(ng) % sustrCG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % svstrCG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % sustrAG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % svstrAG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#    ifndef HCM_AI_WIND
      allocate ( FORCES(ng) % dudtM1G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dudtM2G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dudtM3G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dudtM4G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dudtM5G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dvdtM1G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dvdtM2G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dvdtM3G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dvdtM4G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dvdtM5G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dudtW1G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dudtW2G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dudtW3G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dudtW4G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dudtW5G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dvdtW1G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dvdtW2G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dvdtW3G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dvdtW4G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dvdtW5G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#    endif
#  endif
#  ifdef HCM_EMINP_FORCE
      allocate ( FORCES(ng) % fwfCG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % fwfAG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#    ifndef HCM_AI_FWF
      allocate ( FORCES(ng) % dwdtM1G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dwdtM2G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dwdtM3G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dwdtM4G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dwdtM5G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dwdtW1G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dwdtW2G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dwdtW3G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dwdtW4G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % dwdtW5G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#    endif
#  endif
#  ifdef TIW_COUPLING
      allocate ( FORCES(ng) % sustrTIWG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % svstrTIWG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#    ifndef HCM_AI_TIW
      allocate ( FORCES(ng) % tiwU01G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwV01G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwT01G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwU02G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwV02G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwT02G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwU03G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwV03G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwT03G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwU04G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwV04G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwT04G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwU05G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwV05G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwT05G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwU06G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwV06G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwT06G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwU07G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwV07G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwT07G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwU08G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwV08G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwT08G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwU09G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwV09G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwT09G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwU10G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwV10G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % tiwT10G(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#    endif
#  endif
#  ifdef HCM_BULK_FLUX
      allocate ( FORCES(ng) % TairCG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % PairCG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % QairCG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % SWRG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % LWRG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % LHFG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
      allocate ( FORCES(ng) % SHFG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif

# endif

# if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
      allocate ( FORCES(ng) % sss(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_SSS
      allocate ( FORCES(ng) % sssG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
# endif

# if defined SSSFLX
      allocate ( FORCES(ng) % sssflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_SSSFLX
      allocate ( FORCES(ng) % sssflxG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
# endif

# if defined FASTICE_CLIMATOLOGY
      allocate ( FORCES(ng) % fastice_clm(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  ifndef ANA_FASTICE
      allocate ( FORCES(ng) % fastice_clmG(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d
#  endif
# endif

# if defined ECOSIM || defined SPECTRAL_LIGHT 
      allocate ( FORCES(ng) % SpecIr(LBi:UBi,LBj:UBj,NBands) )
      Dmem(ng)=Dmem(ng)+REAL(NBands,r8)*size2d

      allocate ( FORCES(ng) % avcos(LBi:UBi,LBj:UBj,NBands) )
      Dmem(ng)=Dmem(ng)+REAL(NBands,r8)*size2d
# endif

#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state
!
      allocate ( FORCES(ng) % tl_sustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % tl_svstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d

      allocate ( FORCES(ng) % tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d
# endif
      allocate ( FORCES(ng) % tl_bustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % tl_bvstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifdef SOLVE3D
      allocate ( FORCES(ng) % tl_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d

      allocate ( FORCES(ng) % tl_btflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d
#  ifdef ADJUST_STFLUX
      allocate ( FORCES(ng) % tl_tflux(LBi:UBi,LBj:UBj,Nfrec(ng),       &
     &                                 2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng)*NT(ng),r8)*size2d
#  endif
#  ifdef SHORTWAVE
      allocate ( FORCES(ng) % tl_srflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  ifdef BULK_FLUXES
      allocate ( FORCES(ng) % tl_lhflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % tl_lrflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % tl_shflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#   ifdef EMINUSP
      allocate ( FORCES(ng) % tl_evap(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#   endif
#  endif
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state
!
      allocate ( FORCES(ng) % ad_sustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_svstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

# ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d

      allocate ( FORCES(ng) % ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng),r8)*size2d
# endif
      allocate ( FORCES(ng) % ad_bustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_bvstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_bustr_sol(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_bvstr_sol(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# ifdef SOLVE3D
      allocate ( FORCES(ng) % ad_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d

      allocate ( FORCES(ng) % ad_btflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d
#  ifdef ADJUST_STFLUX
      allocate ( FORCES(ng) % ad_tflux(LBi:UBi,LBj:UBj,Nfrec(ng),       &
     &                                 2,NT(ng)) )
      Dmem(ng)=Dmem(ng)+2.0_r8*REAL(Nfrec(ng)*NT(ng),r8)*size2d
#  endif
#  ifdef SHORTWAVE
      allocate ( FORCES(ng) % ad_srflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  ifdef BULK_FLUXES
      allocate ( FORCES(ng) % ad_lhflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_lrflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % ad_shflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#   ifdef EMINUSP
      allocate ( FORCES(ng) % ad_evap(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#   endif
#  endif
# endif
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Working arrays to store adjoint impulse forcing, background error
!  covariance, background-error standard deviations, or descent
!  conjugate vectors (directions).
!
# if defined FOUR_DVAR || defined IMPULSE
#  ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % b_sustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % b_svstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      allocate ( FORCES(ng) % b_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d
#  endif
# endif
# ifdef FOUR_DVAR
#  ifdef ADJUST_WSTRESS
      allocate ( FORCES(ng) % d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nfrec(ng),r8)*size2d

      allocate ( FORCES(ng) % d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nfrec(ng),r8)*size2d

      allocate ( FORCES(ng) % e_sustr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( FORCES(ng) % e_svstr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      allocate ( FORCES(ng) % d_stflx(LBi:UBi,LBj:UBj,                  &
     &                                Nfrec(ng),NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(Nfrec(ng)*NT(ng),r8)*size2d

      allocate ( FORCES(ng) % e_stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng),r8)*size2d
#  endif
# endif
#endif

      RETURN
      END SUBROUTINE allocate_forces

      SUBROUTINE initialize_forces (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
#ifdef BIOLOGY
      USE mod_biology
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k
#ifdef SOLVE3D
      integer :: itrc
#endif

      real(r8), parameter :: IniVal = 0.0_r8

#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#else
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
#ifdef NCEP_FLUXES
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            FORCES(ng) % sustr(i,j) = IniVal
            FORCES(ng) % svstr(i,j) = IniVal
            FORCES(ng) % nustr(i,j) = IniVal
            FORCES(ng) % nvstr(i,j) = IniVal
            FORCES(ng) % nustrG(i,j,1) = IniVal
            FORCES(ng) % nustrG(i,j,2) = IniVal
            FORCES(ng) % nvstrG(i,j,1) = IniVal
            FORCES(ng) % nvstrG(i,j,2) = IniVal
            FORCES(ng) % bustr(i,j) = IniVal
            FORCES(ng) % bvstr(i,j) = IniVal
            FORCES(ng) % srflx(i,j) = IniVal
            FORCES(ng) % srflxG(i,j,1) = IniVal
            FORCES(ng) % srflxG(i,j,2) = IniVal
            FORCES(ng) % cloud(i,j) = IniVal
            FORCES(ng) % cloudG(i,j,1) = IniVal
            FORCES(ng) % cloudG(i,j,2) = IniVal
            FORCES(ng) % lhflx(i,j) = IniVal
            FORCES(ng) % lhflxG(i,j,1) = IniVal
            FORCES(ng) % lhflxG(i,j,2) = IniVal
            FORCES(ng) % lrflx(i,j) = IniVal
            FORCES(ng) % lrflxG(i,j,1) = IniVal
            FORCES(ng) % lrflxG(i,j,2) = IniVal
            FORCES(ng) % shflx(i,j) = IniVal
            FORCES(ng) % shflxG(i,j,1) = IniVal
            FORCES(ng) % shflxG(i,j,2) = IniVal
            FORCES(ng) % Pair(i,j) = IniVal
            FORCES(ng) % PairG(i,j,1) = IniVal
            FORCES(ng) % PairG(i,j,2) = IniVal
            FORCES(ng) % rain(i,j) = IniVal
            FORCES(ng) % rainG(i,j,1) = IniVal
            FORCES(ng) % rainG(i,j,2) = IniVal
# ifdef RUNOFF
            FORCES(ng) % runoff(i,j) = IniVal
            FORCES(ng) % runoffG(i,j,1) = IniVal
            FORCES(ng) % runoffG(i,j,2) = IniVal
# endif
            FORCES(ng) % skt(i,j) = IniVal
            FORCES(ng) % sktG(i,j,1) = IniVal
            FORCES(ng) % sktG(i,j,2) = IniVal
            FORCES(ng) % icec(i,j) = IniVal
            FORCES(ng) % icecG(i,j,1) = IniVal
            FORCES(ng) % icecG(i,j,2) = IniVal
            DO itrc=1,NT(ng)
              FORCES(ng) % stflx(i,j,itrc) = IniVal
              FORCES(ng) % btflx(i,j,itrc) = IniVal
            END DO
            FORCES(ng) % wg2_d(i,j) = IniVal
            FORCES(ng) % cd_d(i,j) = IniVal
            FORCES(ng) % ch_d(i,j) = IniVal
            FORCES(ng) % ce_d(i,j) = IniVal
            FORCES(ng) % wg2_m(i,j) = IniVal
            FORCES(ng) % cd_m(i,j) = IniVal
            FORCES(ng) % ch_m(i,j) = IniVal
            FORCES(ng) % ce_m(i,j) = IniVal
            FORCES(ng) % rhoa_n(i,j) = IniVal
            FORCES(ng) % cawdir(i,j) = IniVal
            FORCES(ng) % snow_n(i,j) = IniVal
            FORCES(ng) % p_e_n(i,j) = IniVal
            FORCES(ng) % qao_n(i,j) = IniVal
            FORCES(ng) % qio_n(i,j) = IniVal
            FORCES(ng) % qi2_n(i,j) = IniVal
            FORCES(ng) % qai_n(i,j) = IniVal
            FORCES(ng) % qswi_n(i,j) = IniVal
            FORCES(ng) % tau_awx_n(i,j) = IniVal
            FORCES(ng) % tau_awy_n(i,j) = IniVal
            FORCES(ng) % tau_aix_n(i,j) = IniVal
            FORCES(ng) % tau_aiy_n(i,j) = IniVal
            FORCES(ng) % sustr_aw(i,j) = IniVal
            FORCES(ng) % svstr_aw(i,j) = IniVal
          END DO
        END DO
      END IF

       IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
#else
# ifndef NCEP
       IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
# endif
# ifdef ADJUST_WSTRESS
            DO k=1,Nfrec(ng)
              FORCES(ng) % ustr(i,j,k,1) = IniVal
              FORCES(ng) % ustr(i,j,k,2) = IniVal
              FORCES(ng) % vstr(i,j,k,1) = IniVal
              FORCES(ng) % vstr(i,j,k,2) = IniVal
            END DO
# endif
            FORCES(ng) % sustr(i,j) = IniVal
            FORCES(ng) % svstr(i,j) = IniVal
# ifdef ATM2OCN_FLUXES
            FORCES(ng) % Taux(i,j) = IniVal
            FORCES(ng) % Tauy(i,j) = IniVal
# endif
# if !defined ANA_SMFLUX     && !defined BULK_FLUXES || \
     defined NL_BULK_FLUXES
            FORCES(ng) % sustrG(i,j,1) = IniVal
            FORCES(ng) % sustrG(i,j,2) = IniVal
            FORCES(ng) % svstrG(i,j,1) = IniVal
            FORCES(ng) % svstrG(i,j,2) = IniVal
# endif
#endif
            FORCES(ng) % bustr(i,j) = IniVal
            FORCES(ng) % bvstr(i,j) = IniVal
#ifdef WAVES_DIR
            FORCES(ng) % Dwave(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % DwaveG(i,j,1) = IniVal
            FORCES(ng) % DwaveG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_DIRP
            FORCES(ng) % Dwavep(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % DwavepG(i,j,1) = IniVal
            FORCES(ng) % DwavepG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_HEIGHT
            FORCES(ng) % Hwave(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % HwaveG(i,j,1) = IniVal
            FORCES(ng) % HwaveG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_LENGTH
            FORCES(ng) % Lwave(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % LwaveG(i,j,1) = IniVal
            FORCES(ng) % LwaveG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_LENGTHP
            FORCES(ng) % Lwavep(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % LwavepG(i,j,1) = IniVal
            FORCES(ng) % LwavepG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_TOP_PERIOD
            FORCES(ng) % Pwave_top(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Pwave_topG(i,j,1) = IniVal
            FORCES(ng) % Pwave_topG(i,j,2) = IniVal
# endif
#endif
#ifdef WAVES_BOT_PERIOD
            FORCES(ng) % Pwave_bot(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Pwave_botG(i,j,1) = IniVal
            FORCES(ng) % Pwave_botG(i,j,2) = IniVal
# endif
#endif
#if defined BBL_MODEL || defined WAVES_OCEAN
            FORCES(ng) % Uwave_rms(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Uwave_rmsG(i,j,1) = IniVal
            FORCES(ng) % Uwave_rmsG(i,j,2) = IniVal
# endif
#endif

#if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
    defined WAVES_DISS || defined WDISS_INWAVE
            FORCES(ng) % Dissip_break(i,j) = IniVal
            FORCES(ng) % Dissip_wcap(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Dissip_breakG(i,j,1) = IniVal
            FORCES(ng) % Dissip_breakG(i,j,2) = IniVal
            FORCES(ng) % Dissip_wcapG(i,j,1) = IniVal
            FORCES(ng) % Dissip_wcapG(i,j,2) = IniVal
# endif
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
            FORCES(ng) % Dissip_fric(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Dissip_fricG(i,j,1) = IniVal
            FORCES(ng) % Dissip_fricG(i,j,2) = IniVal
# endif
#endif

#if defined WEC_ROLLER
            FORCES(ng) % Dissip_roller(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Dissip_rollerG(i,j,1) = IniVal
            FORCES(ng) % Dissip_rollerG(i,j,2) = IniVal
# endif
#endif

#if defined ROLLER_SVENDSEN
            FORCES(ng) % Wave_break(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Wave_breakG(i,j,1) = IniVal
            FORCES(ng) % Wave_breakG(i,j,2) = IniVal
# endif
#endif

#if defined WAVES_DSPR
            FORCES(ng) % Wave_ds(i,j) = IniVal
            FORCES(ng) % Wave_qp(i,j) = IniVal
# ifndef ANA_WWAVE
            FORCES(ng) % Wave_dsG(i,j,1) = IniVal
            FORCES(ng) % Wave_dsG(i,j,2) = IniVal
            FORCES(ng) % Wave_qpG(i,j,1) = IniVal
            FORCES(ng) % Wave_qpG(i,j,2) = IniVal
# endif
#endif

#if defined WEC_ROLLER
            FORCES(ng) % rollA(i,j) = IniVal
#endif
#if defined BULK_FLUXES || defined BULK_FLUXES2D || defined ECOSIM || \
    defined SPECTRAL_LIGHT || defined ATM2OCN_FLUXES
            FORCES(ng) % Uwind(i,j) = IniVal
            FORCES(ng) % Vwind(i,j) = IniVal
# ifndef ANA_WINDS
            FORCES(ng) % UwindG(i,j,1) = IniVal
            FORCES(ng) % UwindG(i,j,2) = IniVal
            FORCES(ng) % VwindG(i,j,1) = IniVal
            FORCES(ng) % VwindG(i,j,2) = IniVal
# endif
#endif
#if defined BULK_FLUXES || defined ECOSIM || defined CCSM_FLUXES2D || \
    defined SPECTRAL_LIGHT || (defined SHORTWAVE && defined ANA_SRFLUX)
            FORCES(ng) % Hair(i,j) = IniVal
            FORCES(ng) % Tair(i,j) = IniVal
# ifndef ANA_HUMIDITY
            FORCES(ng) % HairG(i,j,1) = IniVal
            FORCES(ng) % HairG(i,j,2) = IniVal
# endif
# ifndef ANA_TAIR
            FORCES(ng) % TairG(i,j,1) = IniVal
            FORCES(ng) % TairG(i,j,2) = IniVal
# endif
#endif
#ifdef CICE_MODEL
            FORCES(ng) % PotT(i,j) = IniVal
            FORCES(ng) % LW_down(i,j) = IniVal
            FORCES(ng) % SW_down(i,j) = IniVal
#endif
#ifdef ICE_DIAGS
            FORCES(ng) % LW_down(i,j) = IniVal
            FORCES(ng) % SW_down(i,j) = IniVal
            FORCES(ng) % lat_ice(i,j) = IniVal
            FORCES(ng) % sens_ice(i,j) = IniVal
            FORCES(ng) % LW_up_ice(i,j) = IniVal
            FORCES(ng) % SW_up_ice(i,j) = IniVal
            FORCES(ng) % saltflux_ice(i,j) = IniVal
            FORCES(ng) % saltflux_ocean(i,j) = IniVal
            FORCES(ng) % snoice(i,j) = IniVal
#endif
#if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS || \
    defined SPECTRAL_LIGHT || defined ATM2OCN_FLUXES || defined CCSM_FLUXES2D
            FORCES(ng) % Pair(i,j) = IniVal
# ifndef ANA_PAIR
            FORCES(ng) % PairG(i,j,1) = IniVal
            FORCES(ng) % PairG(i,j,2) = IniVal
# endif
#endif

#ifdef SOLVE3D
# ifdef SHORTWAVE
            FORCES(ng) % srflx(i,j) = IniVal
#  ifndef ANA_SRFLUX
            FORCES(ng) % srflxG(i,j,1) = IniVal
            FORCES(ng) % srflxG(i,j,2) = IniVal
#  endif
#  ifdef ALBEDO
            FORCES(ng) % albedo(i,j) = IniVal
#   ifdef ICE_MODEL
            FORCES(ng) % albedo_ice(i,j) = IniVal
#   endif
#  endif
#  ifdef ALBEDO_FILE
            FORCES(ng) % albedoG(i,j,1) = IniVal
            FORCES(ng) % albedoG(i,j,2) = IniVal
#  endif
#  if defined ALBEDO_CLOUD
            FORCES(ng) % cawdir(i,j) = IniVal
#  endif
# endif
# if defined RED_TIDE && defined DAILY_SHORTWAVE
            FORCES(ng) % srflx_avg(i,j) = IniVal
            FORCES(ng) % srflxG_avg(i,j,1) = IniVal
            FORCES(ng) % srflxG_avg(i,j,2) = IniVal
# endif
# if defined ICE_BULK_FLUXES && !defined NCEP_FLUXES
            FORCES(ng) % sustr_aw(i,j) = IniVal
            FORCES(ng) % svstr_aw(i,j) = IniVal
            FORCES(ng) % tau_aix_n(i,j) = IniVal
            FORCES(ng) % tau_aiy_n(i,j) = IniVal
            FORCES(ng) % qai_n(i,j) = IniVal
            FORCES(ng) % qi_o_n(i,j) = IniVal
            FORCES(ng) % SW_thru_ice(i,j) = IniVal
            FORCES(ng) % qswi_n(i,j) = IniVal
            FORCES(ng) % qao_n(i,j) = IniVal
            FORCES(ng) % qio_n(i,j) = IniVal
            FORCES(ng) % qi2_n(i,j) = IniVal
            FORCES(ng) % snow_n(i,j) = IniVal
            FORCES(ng) % p_e_n(i,j) = IniVal
# endif
# ifdef CLOUDS
            FORCES(ng) % cloud(i,j) = IniVal
#  ifndef ANA_CLOUD
            FORCES(ng) % cloudG(i,j,1) = IniVal
            FORCES(ng) % cloudG(i,j,2) = IniVal
#  endif
# endif
# if defined BULK_FLUXES || defined ATM2OCN_FLUXES
            FORCES(ng) % lhflx(i,j) = IniVal
            FORCES(ng) % lrflx(i,j) = IniVal
            FORCES(ng) % shflx(i,j) = IniVal
            FORCES(ng) % snhflx(i,j) = IniVal
# endif
# if defined BULK_FLUXES || defined ATM2OCN_FLUXES
            FORCES(ng) % rain(i,j) = IniVal
#  ifdef SNOWFALL
            FORCES(ng) % snow(i,j) = IniVal
#  endif
#  ifndef ANA_RAIN
            FORCES(ng) % rainG(i,j,1) = IniVal
            FORCES(ng) % rainG(i,j,2) = IniVal
#   ifdef SNOWFALL
            FORCES(ng) % snowG(i,j,1) = IniVal
            FORCES(ng) % snowG(i,j,2) = IniVal
#   endif
#  endif
#  ifdef EMINUSP
            FORCES(ng) % EminusP(i,j) = IniVal
            FORCES(ng) % evap(i,j) = IniVal
#  endif
#  ifdef RUNOFF
            FORCES(ng) % runoff(i,j) = IniVal
            FORCES(ng) % runoffG(i,j,1) = IniVal
            FORCES(ng) % runoffG(i,j,2) = IniVal
#  endif
# endif
# if !defined LONGWAVE && defined BULK_FLUXES
            FORCES(ng) % lrflxG(i,j,1) = IniVal
            FORCES(ng) % lrflxG(i,j,2) = IniVal
# endif
# if !defined ANA_RAIN && defined BULK_FLUXES
            FORCES(ng) % rainG(i,j,1) = IniVal
            FORCES(ng) % rainG(i,j,2) = IniVal
#  ifdef SNOWFALL
            FORCES(ng) % snowG(i,j,1) = IniVal
            FORCES(ng) % snowG(i,j,2) = IniVal
#  endif
# endif
# ifdef QCORRECTION
            FORCES(ng) % sst(i,j) = IniVal
            FORCES(ng) % dqdt(i,j) = IniVal
#  ifndef ANA_SST
            FORCES(ng) % sstG(i,j,1) = IniVal
            FORCES(ng) % sstG(i,j,2) = IniVal
            FORCES(ng) % dqdtG(i,j,1) = IniVal
            FORCES(ng) % dqdtG(i,j,2) = IniVal
#  endif
# endif

# ifdef SSH_NUDGING
            FORCES(ng) % sshC(i,j) = IniVal
            FORCES(ng) % sshCG(i,j,1) = IniVal
            FORCES(ng) % sshCG(i,j,2) = IniVal
# endif

# ifdef HCM_COUPLING
#  ifdef HCM_SST_FILTER
#   ifdef HCM_ROSSBY_FILTER
            FORCES(ng) % Rrob(i,j) = IniVal
#   endif
            FORCES(ng) % sstF(i,j) = IniVal
            FORCES(ng) % sstH(i,j) = IniVal
            FORCES(ng) % sstL(i,j) = IniVal
#  endif
#  ifndef HCM_CLIM_SPIN
            FORCES(ng) % sstC(i,j) = IniVal
            FORCES(ng) % sstA(i,j) = IniVal
#  endif
#  ifdef HCM_WIND_STRESS
            FORCES(ng) % sustrC(i,j) = IniVal
            FORCES(ng) % svstrC(i,j) = IniVal
            FORCES(ng) % sustrA(i,j) = IniVal
            FORCES(ng) % svstrA(i,j) = IniVal
#    ifndef HCM_AI_WIND
            FORCES(ng) % dudtM1(i,j) = IniVal
            FORCES(ng) % dudtM2(i,j) = IniVal
            FORCES(ng) % dudtM3(i,j) = IniVal
            FORCES(ng) % dudtM4(i,j) = IniVal
            FORCES(ng) % dudtM5(i,j) = IniVal
            FORCES(ng) % dvdtM1(i,j) = IniVal
            FORCES(ng) % dvdtM2(i,j) = IniVal
            FORCES(ng) % dvdtM3(i,j) = IniVal
            FORCES(ng) % dvdtM4(i,j) = IniVal
            FORCES(ng) % dvdtM5(i,j) = IniVal
            FORCES(ng) % dudtW1(i,j) = IniVal
            FORCES(ng) % dudtW2(i,j) = IniVal
            FORCES(ng) % dudtW3(i,j) = IniVal
            FORCES(ng) % dudtW4(i,j) = IniVal
            FORCES(ng) % dudtW5(i,j) = IniVal
            FORCES(ng) % dvdtW1(i,j) = IniVal
            FORCES(ng) % dvdtW2(i,j) = IniVal
            FORCES(ng) % dvdtW3(i,j) = IniVal
            FORCES(ng) % dvdtW4(i,j) = IniVal
            FORCES(ng) % dvdtW5(i,j) = IniVal
#    endif
#  endif
#  ifdef HCM_EMINP_FORCE
            FORCES(ng) % fwfC(i,j) = IniVal
            FORCES(ng) % fwfA(i,j) = IniVal
#    ifndef HCM_AI_FWF
            FORCES(ng) % dwdtM1(i,j) = IniVal
            FORCES(ng) % dwdtM2(i,j) = IniVal
            FORCES(ng) % dwdtM3(i,j) = IniVal
            FORCES(ng) % dwdtM4(i,j) = IniVal
            FORCES(ng) % dwdtM5(i,j) = IniVal
            FORCES(ng) % dwdtW1(i,j) = IniVal
            FORCES(ng) % dwdtW2(i,j) = IniVal
            FORCES(ng) % dwdtW3(i,j) = IniVal
            FORCES(ng) % dwdtW4(i,j) = IniVal
            FORCES(ng) % dwdtW5(i,j) = IniVal
#    endif
#  endif
#  ifdef TIW_COUPLING
            FORCES(ng) % sustrTIW(i,j) = IniVal
            FORCES(ng) % svstrTIW(i,j) = IniVal
#    ifndef HCM_AI_TIW
            FORCES(ng) % tiwU01(i,j) = IniVal
            FORCES(ng) % tiwV01(i,j) = IniVal
            FORCES(ng) % tiwT01(i,j) = IniVal
            FORCES(ng) % tiwU02(i,j) = IniVal
            FORCES(ng) % tiwV02(i,j) = IniVal
            FORCES(ng) % tiwT02(i,j) = IniVal
            FORCES(ng) % tiwU03(i,j) = IniVal
            FORCES(ng) % tiwV03(i,j) = IniVal
            FORCES(ng) % tiwT03(i,j) = IniVal
            FORCES(ng) % tiwU04(i,j) = IniVal
            FORCES(ng) % tiwV04(i,j) = IniVal
            FORCES(ng) % tiwT04(i,j) = IniVal
            FORCES(ng) % tiwU05(i,j) = IniVal
            FORCES(ng) % tiwV05(i,j) = IniVal
            FORCES(ng) % tiwT05(i,j) = IniVal
            FORCES(ng) % tiwU06(i,j) = IniVal
            FORCES(ng) % tiwV06(i,j) = IniVal
            FORCES(ng) % tiwT06(i,j) = IniVal
            FORCES(ng) % tiwU07(i,j) = IniVal
            FORCES(ng) % tiwV07(i,j) = IniVal
            FORCES(ng) % tiwT07(i,j) = IniVal
            FORCES(ng) % tiwU08(i,j) = IniVal
            FORCES(ng) % tiwV08(i,j) = IniVal
            FORCES(ng) % tiwT08(i,j) = IniVal
            FORCES(ng) % tiwU09(i,j) = IniVal
            FORCES(ng) % tiwV09(i,j) = IniVal
            FORCES(ng) % tiwT09(i,j) = IniVal
            FORCES(ng) % tiwU10(i,j) = IniVal
            FORCES(ng) % tiwV10(i,j) = IniVal
            FORCES(ng) % tiwT10(i,j) = IniVal
#    endif
#  endif
#  ifdef HCM_BULK_FLUX
            FORCES(ng) % TairC(i,j) = IniVal
            FORCES(ng) % PairC(i,j) = IniVal
            FORCES(ng) % QairC(i,j) = IniVal
            FORCES(ng) % SWR(i,j) = IniVal
            FORCES(ng) % LWR(i,j) = IniVal
            FORCES(ng) % LHF(i,j) = IniVal
            FORCES(ng) % SHF(i,j) = IniVal
#  endif

#  ifdef HCM_SST_FILTER
#   ifdef HCM_ROSSBY_FILTER
            FORCES(ng) % RrobG(i,j,1) = IniVal
            FORCES(ng) % RrobG(i,j,2) = IniVal
#   endif
#  endif
#  ifndef HCM_CLIM_SPIN
            FORCES(ng) % sstCG(i,j,1) = IniVal
            FORCES(ng) % sstCG(i,j,2) = IniVal
            FORCES(ng) % sstAG(i,j,1) = IniVal
            FORCES(ng) % sstAG(i,j,2) = IniVal
#  endif
#  ifdef HCM_WIND_STRESS
            FORCES(ng) % sustrCG(i,j,1) = IniVal
            FORCES(ng) % sustrCG(i,j,2) = IniVal
            FORCES(ng) % svstrCG(i,j,1) = IniVal
            FORCES(ng) % svstrCG(i,j,2) = IniVal
            FORCES(ng) % sustrAG(i,j,1) = IniVal
            FORCES(ng) % sustrAG(i,j,2) = IniVal
            FORCES(ng) % svstrAG(i,j,1) = IniVal
            FORCES(ng) % svstrAG(i,j,2) = IniVal
#    ifndef HCM_AI_WIND
            FORCES(ng) % dudtM1G(i,j,1) = IniVal
            FORCES(ng) % dudtM1G(i,j,2) = IniVal
            FORCES(ng) % dudtM2G(i,j,1) = IniVal
            FORCES(ng) % dudtM2G(i,j,2) = IniVal
            FORCES(ng) % dudtM3G(i,j,1) = IniVal
            FORCES(ng) % dudtM3G(i,j,2) = IniVal
            FORCES(ng) % dudtM4G(i,j,1) = IniVal
            FORCES(ng) % dudtM4G(i,j,2) = IniVal
            FORCES(ng) % dudtM5G(i,j,1) = IniVal
            FORCES(ng) % dudtM5G(i,j,2) = IniVal
            FORCES(ng) % dvdtM1G(i,j,1) = IniVal
            FORCES(ng) % dvdtM1G(i,j,2) = IniVal
            FORCES(ng) % dvdtM2G(i,j,1) = IniVal
            FORCES(ng) % dvdtM2G(i,j,2) = IniVal
            FORCES(ng) % dvdtM3G(i,j,1) = IniVal
            FORCES(ng) % dvdtM3G(i,j,2) = IniVal
            FORCES(ng) % dvdtM4G(i,j,1) = IniVal
            FORCES(ng) % dvdtM4G(i,j,2) = IniVal
            FORCES(ng) % dvdtM5G(i,j,1) = IniVal
            FORCES(ng) % dvdtM5G(i,j,2) = IniVal
            FORCES(ng) % dudtW1G(i,j,1) = IniVal
            FORCES(ng) % dudtW1G(i,j,2) = IniVal
            FORCES(ng) % dudtW2G(i,j,1) = IniVal
            FORCES(ng) % dudtW2G(i,j,2) = IniVal
            FORCES(ng) % dudtW3G(i,j,1) = IniVal
            FORCES(ng) % dudtW3G(i,j,2) = IniVal
            FORCES(ng) % dudtW4G(i,j,1) = IniVal
            FORCES(ng) % dudtW4G(i,j,2) = IniVal
            FORCES(ng) % dudtW5G(i,j,1) = IniVal
            FORCES(ng) % dudtW5G(i,j,2) = IniVal
            FORCES(ng) % dvdtW1G(i,j,1) = IniVal
            FORCES(ng) % dvdtW1G(i,j,2) = IniVal
            FORCES(ng) % dvdtW2G(i,j,1) = IniVal
            FORCES(ng) % dvdtW2G(i,j,2) = IniVal
            FORCES(ng) % dvdtW3G(i,j,1) = IniVal
            FORCES(ng) % dvdtW3G(i,j,2) = IniVal
            FORCES(ng) % dvdtW4G(i,j,1) = IniVal
            FORCES(ng) % dvdtW4G(i,j,2) = IniVal
            FORCES(ng) % dvdtW5G(i,j,1) = IniVal
            FORCES(ng) % dvdtW5G(i,j,2) = IniVal
#    endif
#  endif
#  ifdef HCM_EMINP_FORCE
            FORCES(ng) % fwfCG(i,j,1) = IniVal
            FORCES(ng) % fwfCG(i,j,2) = IniVal
            FORCES(ng) % fwfAG(i,j,1) = IniVal
            FORCES(ng) % fwfAG(i,j,2) = IniVal
#    ifndef HCM_AI_FWF
            FORCES(ng) % dwdtM1G(i,j,1) = IniVal
            FORCES(ng) % dwdtM1G(i,j,2) = IniVal
            FORCES(ng) % dwdtM2G(i,j,1) = IniVal
            FORCES(ng) % dwdtM2G(i,j,2) = IniVal
            FORCES(ng) % dwdtM3G(i,j,1) = IniVal
            FORCES(ng) % dwdtM3G(i,j,2) = IniVal
            FORCES(ng) % dwdtM4G(i,j,1) = IniVal
            FORCES(ng) % dwdtM4G(i,j,2) = IniVal
            FORCES(ng) % dwdtM5G(i,j,1) = IniVal
            FORCES(ng) % dwdtM5G(i,j,2) = IniVal
            FORCES(ng) % dwdtW1G(i,j,1) = IniVal
            FORCES(ng) % dwdtW1G(i,j,2) = IniVal
            FORCES(ng) % dwdtW2G(i,j,1) = IniVal
            FORCES(ng) % dwdtW2G(i,j,2) = IniVal
            FORCES(ng) % dwdtW3G(i,j,1) = IniVal
            FORCES(ng) % dwdtW3G(i,j,2) = IniVal
            FORCES(ng) % dwdtW4G(i,j,1) = IniVal
            FORCES(ng) % dwdtW4G(i,j,2) = IniVal
            FORCES(ng) % dwdtW5G(i,j,1) = IniVal
            FORCES(ng) % dwdtW5G(i,j,2) = IniVal
#    endif
#  endif
#  ifdef TIW_COUPLING
            FORCES(ng) % sustrTIWG(i,j,1) = IniVal
            FORCES(ng) % sustrTIWG(i,j,2) = IniVal
            FORCES(ng) % svstrTIWG(i,j,1) = IniVal
            FORCES(ng) % svstrTIWG(i,j,2) = IniVal
#    ifndef HCM_AI_TIW
            FORCES(ng) % tiwU01G(i,j,1) = IniVal
            FORCES(ng) % tiwU01G(i,j,2) = IniVal
            FORCES(ng) % tiwV01G(i,j,1) = IniVal
            FORCES(ng) % tiwV01G(i,j,2) = IniVal
            FORCES(ng) % tiwT01G(i,j,1) = IniVal
            FORCES(ng) % tiwT01G(i,j,2) = IniVal
            FORCES(ng) % tiwU02G(i,j,1) = IniVal
            FORCES(ng) % tiwU02G(i,j,2) = IniVal
            FORCES(ng) % tiwV02G(i,j,1) = IniVal
            FORCES(ng) % tiwV02G(i,j,2) = IniVal
            FORCES(ng) % tiwT02G(i,j,1) = IniVal
            FORCES(ng) % tiwT02G(i,j,2) = IniVal
            FORCES(ng) % tiwU03G(i,j,1) = IniVal
            FORCES(ng) % tiwU03G(i,j,2) = IniVal
            FORCES(ng) % tiwV03G(i,j,1) = IniVal
            FORCES(ng) % tiwV03G(i,j,2) = IniVal
            FORCES(ng) % tiwT03G(i,j,1) = IniVal
            FORCES(ng) % tiwT03G(i,j,2) = IniVal
            FORCES(ng) % tiwU04G(i,j,1) = IniVal
            FORCES(ng) % tiwU04G(i,j,2) = IniVal
            FORCES(ng) % tiwV04G(i,j,1) = IniVal
            FORCES(ng) % tiwV04G(i,j,2) = IniVal
            FORCES(ng) % tiwT04G(i,j,1) = IniVal
            FORCES(ng) % tiwT04G(i,j,2) = IniVal
            FORCES(ng) % tiwU05G(i,j,1) = IniVal
            FORCES(ng) % tiwU05G(i,j,2) = IniVal
            FORCES(ng) % tiwV05G(i,j,1) = IniVal
            FORCES(ng) % tiwV05G(i,j,2) = IniVal
            FORCES(ng) % tiwT05G(i,j,1) = IniVal
            FORCES(ng) % tiwT05G(i,j,2) = IniVal
            FORCES(ng) % tiwU06G(i,j,1) = IniVal
            FORCES(ng) % tiwU06G(i,j,2) = IniVal
            FORCES(ng) % tiwV06G(i,j,1) = IniVal
            FORCES(ng) % tiwV06G(i,j,2) = IniVal
            FORCES(ng) % tiwT06G(i,j,1) = IniVal
            FORCES(ng) % tiwT06G(i,j,2) = IniVal
            FORCES(ng) % tiwU07G(i,j,1) = IniVal
            FORCES(ng) % tiwU07G(i,j,2) = IniVal
            FORCES(ng) % tiwV07G(i,j,1) = IniVal
            FORCES(ng) % tiwV07G(i,j,2) = IniVal
            FORCES(ng) % tiwT07G(i,j,1) = IniVal
            FORCES(ng) % tiwT07G(i,j,2) = IniVal
            FORCES(ng) % tiwU08G(i,j,1) = IniVal
            FORCES(ng) % tiwU08G(i,j,2) = IniVal
            FORCES(ng) % tiwV08G(i,j,1) = IniVal
            FORCES(ng) % tiwV08G(i,j,2) = IniVal
            FORCES(ng) % tiwT08G(i,j,1) = IniVal
            FORCES(ng) % tiwT08G(i,j,2) = IniVal
            FORCES(ng) % tiwU09G(i,j,1) = IniVal
            FORCES(ng) % tiwU09G(i,j,2) = IniVal
            FORCES(ng) % tiwV09G(i,j,1) = IniVal
            FORCES(ng) % tiwV09G(i,j,2) = IniVal
            FORCES(ng) % tiwT09G(i,j,1) = IniVal
            FORCES(ng) % tiwT09G(i,j,2) = IniVal
            FORCES(ng) % tiwU10G(i,j,1) = IniVal
            FORCES(ng) % tiwU10G(i,j,2) = IniVal
            FORCES(ng) % tiwV10G(i,j,1) = IniVal
            FORCES(ng) % tiwV10G(i,j,2) = IniVal
            FORCES(ng) % tiwT10G(i,j,1) = IniVal
            FORCES(ng) % tiwT10G(i,j,2) = IniVal
#    endif
#  endif
#  ifdef HCM_BULK_FLUX
            FORCES(ng) % TairCG(i,j,1) = IniVal
            FORCES(ng) % TairCG(i,j,2) = IniVal
            FORCES(ng) % PairCG(i,j,1) = IniVal
            FORCES(ng) % PairCG(i,j,2) = IniVal
            FORCES(ng) % QairCG(i,j,1) = IniVal
            FORCES(ng) % QairCG(i,j,2) = IniVal
            FORCES(ng) % SWRG(i,j,1) = IniVal
            FORCES(ng) % SWRG(i,j,2) = IniVal
            FORCES(ng) % LWRG(i,j,1) = IniVal
            FORCES(ng) % LWRG(i,j,2) = IniVal
            FORCES(ng) % LHFG(i,j,1) = IniVal
            FORCES(ng) % LHFG(i,j,2) = IniVal
            FORCES(ng) % SHFG(i,j,1) = IniVal
            FORCES(ng) % SHFG(i,j,2) = IniVal
#  endif

# endif


# if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
            FORCES(ng) % sss(i,j) = IniVal
#  ifndef ANA_SSS
            FORCES(ng) % sssG(i,j,1) = IniVal
            FORCES(ng) % sssG(i,j,2) = IniVal
#  endif
# endif
# if defined SSSFLX
            FORCES(ng) % sssflx(i,j) = IniVal
#  ifndef ANA_SSSFLX
            FORCES(ng) % sssflxG(i,j,1) = IniVal
            FORCES(ng) % sssflxG(i,j,2) = IniVal
#  endif
# endif
# if defined FASTICE_CLIMATOLOGY
            FORCES(ng) % fastice_clm(i,j) = IniVal
#  ifndef ANA_FASTICE
            FORCES(ng) % fastice_clmG(i,j,1) = IniVal
            FORCES(ng) % fastice_clmG(i,j,2) = IniVal
#  endif
# endif
# if defined PERFECT_RESTART && defined ICE_MODEL
            FORCES(ng) % stflx_save(i,j,1) = IniVal
            FORCES(ng) % stflx_save(i,j,2) = IniVal
            FORCES(ng) % sustr_save(i,j) = IniVal
            FORCES(ng) % svstr_save(i,j) = IniVal
# endif
            DO itrc=1,NT(ng)
# ifdef ADJUST_STFLUX
              DO k=1,Nfrec(ng)
                FORCES(ng) % tflux(i,j,k,1,itrc) = IniVal
                FORCES(ng) % tflux(i,j,k,2,itrc) = IniVal
              END DO
# endif
              FORCES(ng) % stflx(i,j,itrc) = IniVal
              FORCES(ng) % btflx(i,j,itrc) = IniVal
# if !defined ANA_STFLUX || !defined ANA_SSFLUX || \
     !defined ANA_SPFLUX
              FORCES(ng) % stflxG(i,j,1,itrc) = IniVal
              FORCES(ng) % stflxG(i,j,2,itrc) = IniVal
# endif
# if !defined ANA_BTFLUX || !defined ANA_BSFLUX || \
     !defined ANA_BPFLUX
              FORCES(ng) % btflxG(i,j,1,itrc) = IniVal
              FORCES(ng) % btflxG(i,j,2,itrc) = IniVal
# endif
            END DO
# if defined ECOSIM || defined SPECTRAL_LIGHT 
            DO itrc=1,NBands
              FORCES(ng) % SpecIr(i,j,itrc) = IniVal
              FORCES(ng) % avcos(i,j,itrc) = IniVal
            END DO
# endif
#endif
          END DO
        END DO
      END IF

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      IF ((model.eq.0).or.(model.eq.iTLM).or.(model.eq.iRPM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
# ifdef ADJUST_WSTRESS
            DO k=1,Nfrec(ng)
              FORCES(ng) % tl_ustr(i,j,k,1) = IniVal
              FORCES(ng) % tl_ustr(i,j,k,2) = IniVal
              FORCES(ng) % tl_vstr(i,j,k,1) = IniVal
              FORCES(ng) % tl_vstr(i,j,k,2) = IniVal
            END DO
# endif
            FORCES(ng) % tl_sustr(i,j) = IniVal
            FORCES(ng) % tl_svstr(i,j) = IniVal
            FORCES(ng) % tl_bustr(i,j) = IniVal
            FORCES(ng) % tl_bvstr(i,j) = IniVal
          END DO
# ifdef SOLVE3D
#  ifdef SHORTWAVE
          DO i=Imin,Imax
            FORCES(ng) % tl_srflx(i,j) = IniVal
          END DO
#  endif
#  ifdef BULK_FLUXES
          DO i=Imin,Imax
            FORCES(ng) % tl_lhflx(i,j) = IniVal
            FORCES(ng) % tl_lrflx(i,j) = IniVal
            FORCES(ng) % tl_shflx(i,j) = IniVal
#   ifdef EMINUSP
            FORCES(ng) % tl_evap(i,j) = IniVal
#   endif
          END DO
#  endif
          DO itrc=1,NT(ng)
            DO i=Imin,Imax
#  ifdef ADJUST_STFLUX
              DO k=1,Nfrec(ng)
                FORCES(ng) % tl_tflux(i,j,k,1,itrc) = IniVal
                FORCES(ng) % tl_tflux(i,j,k,2,itrc) = IniVal
              END DO
#  endif
              FORCES(ng) % tl_stflx(i,j,itrc) = IniVal
              FORCES(ng) % tl_btflx(i,j,itrc) = IniVal
            END DO
          END DO
# endif
        END DO
      END IF
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      IF ((model.eq.0).or.(model.eq.iADM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
# ifdef ADJUST_WSTRESS
            DO k=1,Nfrec(ng)
              FORCES(ng) % ad_ustr(i,j,k,1) = IniVal
              FORCES(ng) % ad_ustr(i,j,k,2) = IniVal
              FORCES(ng) % ad_vstr(i,j,k,1) = IniVal
              FORCES(ng) % ad_vstr(i,j,k,2) = IniVal
            END DO
# endif
            FORCES(ng) % ad_sustr(i,j) = IniVal
            FORCES(ng) % ad_svstr(i,j) = IniVal
            FORCES(ng) % ad_bustr(i,j) = IniVal
            FORCES(ng) % ad_bvstr(i,j) = IniVal
            FORCES(ng) % ad_bustr_sol(i,j) = IniVal
            FORCES(ng) % ad_bvstr_sol(i,j) = IniVal
          END DO
# ifdef SOLVE3D
#  ifdef SHORTWAVE
          DO i=Imin,Imax
            FORCES(ng) % ad_srflx(i,j) = IniVal
          END DO
#  endif
#  ifdef BULK_FLUXES
          DO i=Imin,Imax
            FORCES(ng) % ad_lhflx(i,j) = IniVal
            FORCES(ng) % ad_lrflx(i,j) = IniVal
            FORCES(ng) % ad_shflx(i,j) = IniVal
#   ifdef EMINUSP
            FORCES(ng) % ad_evap(i,j) = IniVal
#   endif
          END DO
#  endif
          DO itrc=1,NT(ng)
            DO i=Imin,Imax
#  ifdef ADJUST_STFLUX
              DO k=1,Nfrec(ng)
                FORCES(ng) % ad_tflux(i,j,k,1,itrc) = IniVal
                FORCES(ng) % ad_tflux(i,j,k,2,itrc) = IniVal
              END DO
#  endif
              FORCES(ng) % ad_stflx(i,j,itrc) = IniVal
              FORCES(ng) % ad_btflx(i,j,itrc) = IniVal
            END DO
          END DO
# endif
        END DO
      END IF
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Working arrays to store adjoint impulse forcing, background error
!  covariance, background-error standard deviations, or descent
!  conjugate vectors (directions).
!
# if defined FOUR_DVAR || defined IMPULSE
      IF (model.eq.0) THEN
#  ifdef ADJUST_WSTRESS
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            FORCES(ng) % b_sustr(i,j) = IniVal
            FORCES(ng) % b_svstr(i,j) = IniVal
#   ifdef FOUR_DVAR
            FORCES(ng) % e_sustr(i,j) = IniVal
            FORCES(ng) % e_svstr(i,j) = IniVal
            DO k=1,Nfrec(ng)
              FORCES(ng) % d_sustr(i,j,k) = IniVal
              FORCES(ng) % d_svstr(i,j,k) = IniVal
            END DO
#   endif
          END DO
        END DO
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
        DO itrc=1,NT(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FORCES(ng) % b_stflx(i,j,itrc) = IniVal
#   ifdef FOUR_DVAR
              FORCES(ng) % e_stflx(i,j,itrc) = IniVal
              DO k=1,Nfrec(ng)
                FORCES(ng) % d_stflx(i,j,k,itrc) = IniVal
              END DO
#   endif
            END DO
          END DO
        END DO
#  endif
      END IF
# endif
#endif

      RETURN
      END SUBROUTINE initialize_forces

      END MODULE mod_forces
