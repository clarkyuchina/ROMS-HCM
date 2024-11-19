#include "cppdefs.h"
      MODULE hcm_forc_mod

#if defined HCM_COUPLING
#    ifndef HCM_CLIM_SPIN 
      implicit none
!
      PRIVATE
      PUBLIC  :: hcm_forcing
!
      CONTAINS
!
!***********************************************************************      
      SUBROUTINE hcm_forcing(ng)
!***********************************************************************
!             
      USE mod_param
      USE mod_scalars
      USE mod_grid
      USE mod_parallel
      USE mod_forces
      USE mod_ocean
      USE mod_stepping
      USE dateclock_mod,        ONLY : time_string
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!       
      integer :: tile
      integer :: Ntiles

!--MPI-related parameters--
      integer :: len_buffer,request,ierr
      integer,dimension(mpi_status_size) :: rstatus
      real(r8),allocatable,dimension(:) ::send_buffer
      real(r8),allocatable,dimension(:) ::recv_buffer
 
#     ifdef HCM_INITIAL_KICK
      real(r8):: lonr, latr
      real(r8):: cffX, cffY, cffL
#     endif
      real(r8):: cffT
      real(r8):: sstC, sstM
      real(r8):: cffM, cffR
#     ifdef HCM_SSTA_REINFORCE
      real(r8):: PI_half
#     endif

      integer :: i, j
      integer :: LBi,UBi,LBj,UBj,dim1,dim2
      integer :: ILB,IUB,JLB,JUB
      integer :: IstrR, IendR, JstrR, JendR

      integer :: tagSSTA,tagSSTN
      real(r8),allocatable,dimension(:,:) :: sstA
#ifdef HCM_AI_COUPLING
      real(r8):: COUPnum
      real(r8):: AI_NUM
      logical :: AI_FINISH
      integer :: tagLON,tagLAT
      character (len=22) :: fileT
      character (len=:),allocatable :: fileN
      real(r8),allocatable,dimension(:,:) :: lonM, latM
      real(r8),allocatable,dimension(:,:) :: AIDT
#endif
#ifdef HCM_WIND_STRESS
      integer :: tagU,tagV,tagUc,tagVc,tagUa,tagVa
      real(r8),allocatable,dimension(:,:) :: sustr, svstr
      real(r8),allocatable,dimension(:,:) :: sustrC, svstrC
      real(r8),allocatable,dimension(:,:) :: sustrA, svstrA
#  ifndef HCM_AI_WIND
      integer :: tagUw, tagVw
      real(r8):: cffU1, cffU2, cffU3, cffU4, cffU5
      real(r8):: cffV1, cffV2, cffV3, cffV4, cffV5
      integer :: tagUm1,tagUm2,tagUm3,tagUm4,tagUm5
      integer :: tagVm1,tagVm2,tagVm3,tagVm4,tagVm5
      real(r8),allocatable,dimension(:,:) :: Um1, Um2, Um3, Um4, Um5
      real(r8),allocatable,dimension(:,:) :: Vm1, Vm2, Vm3, Vm4, Vm5
#  endif
#endif
#ifdef HCM_EMINP_FORCE
      integer :: tagW,tagWc,tagWa
      real(r8),allocatable,dimension(:,:) :: fwf, fwfC, fwfA
#  ifndef HCM_AI_FWF
      integer :: tagWw
      real(r8):: cffW1, cffW2, cffW3, cffW4, cffW5
      integer :: tagWm1,tagWm2,tagWm3,tagWm4,tagWm5
      real(r8),allocatable,dimension(:,:) :: Wm1, Wm2, Wm3, Wm4, Wm5
#  endif
#endif
#ifdef TIW_COUPLING
      real(r8):: cffTIW
      integer :: tagTIWu,  tagTIWv
      real(r8),allocatable,dimension(:,:) :: sustrTIW, svstrTIW
#  ifndef HCM_AI_TIW
      integer :: tagTIWw
      integer :: tagTIWu01,tagTIWu02,tagTIWu03,tagTIWu04,tagTIWu05
      integer :: tagTIWu06,tagTIWu07,tagTIWu08,tagTIWu09,tagTIWu10
      integer :: tagTIWv01,tagTIWv02,tagTIWv03,tagTIWv04,tagTIWv05
      integer :: tagTIWv06,tagTIWv07,tagTIWv08,tagTIWv09,tagTIWv10
      real(r8):: cffTIW01, cffTIW02, cffTIW03, cffTIW04, cffTIW05
      real(r8):: cffTIW06, cffTIW07, cffTIW08, cffTIW09, cffTIW10
      real(r8),allocatable,dimension(:,:) :: TIWu01, TIWu02, TIWu03, TIWu04, TIWu05
      real(r8),allocatable,dimension(:,:) :: TIWu06, TIWu07, TIWu08, TIWu09, TIWu10
      real(r8),allocatable,dimension(:,:) :: TIWv01, TIWv02, TIWv03, TIWv04, TIWv05
      real(r8),allocatable,dimension(:,:) :: TIWv06, TIWv07, TIWv08, TIWv09, TIWv10
#  else
      integer :: tagSSTH
      real(r8),allocatable,dimension(:,:) :: sstH
#  endif
#endif

#ifdef HCM_SSTA_REINFORCE
      PI_half=2*ATAN(1.0_r8);
#endif

      Ntiles=NtileI(ng)*NtileJ(ng)-1

      tagSSTA=101
      tagSSTN=102
#ifdef HCM_AI_COUPLING
      tagLON=103
      tagLAT=104
#endif
#ifdef HCM_WIND_STRESS
      tagU=201
      tagV=202
      tagUc=203
      tagVc=204
      tagUa=205
      tagVa=206
#  ifndef HCM_AI_WIND
      tagUw=211
      tagVw=212
      tagUm1=221
      tagUm2=222
      tagUm3=223
      tagUm4=224
      tagUm5=225
      tagVm1=231
      tagVm2=232
      tagVm3=233
      tagVm4=234
      tagVm5=235
#  endif
#endif
#ifdef HCM_EMINP_FORCE
      tagW=301
      tagWc=302
      tagWa=303
#  ifndef HCM_AI_FWF
      tagWw=311
      tagWm1=321
      tagWm2=322
      tagWm3=323
      tagWm4=324
      tagWm5=325
#  endif
#endif
#ifdef TIW_COUPLING
      tagTIWu=401
      tagTIWv=402
#  ifndef HCM_AI_TIW
      tagTIWw=411
      tagTIWu01=421
      tagTIWu02=422
      tagTIWu03=423
      tagTIWu04=424
      tagTIWu05=425
      tagTIWu06=426
      tagTIWu07=427
      tagTIWu08=428
      tagTIWu09=429
      tagTIWu10=430
      tagTIWv01=431
      tagTIWv02=432
      tagTIWv03=433
      tagTIWv04=434
      tagTIWv05=435
      tagTIWv06=436
      tagTIWv07=437
      tagTIWv08=438
      tagTIWv09=439
      tagTIWv10=440
#  else
      tagSSTH=451
#  endif
#endif

!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(tile)
      IUB=BOUNDS(ng)%UBi(tile)
      JLB=BOUNDS(ng)%LBj(tile)
      JUB=BOUNDS(ng)%UBj(tile)

      allocate(sstA(ILB:IUB,JLB:JUB))
#ifdef HCM_AI_COUPLING
      allocate(lonM(ILB:IUB,JLB:JUB))
      allocate(latM(ILB:IUB,JLB:JUB))
      allocate(AIDT(ILB:IUB,JLB:JUB))
#endif
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
      allocate(Um1(ILB:IUB,JLB:JUB))
      allocate(Um2(ILB:IUB,JLB:JUB))
      allocate(Um3(ILB:IUB,JLB:JUB))
      allocate(Um4(ILB:IUB,JLB:JUB))
      allocate(Um5(ILB:IUB,JLB:JUB))
      allocate(Vm1(ILB:IUB,JLB:JUB))
      allocate(Vm2(ILB:IUB,JLB:JUB))
      allocate(Vm3(ILB:IUB,JLB:JUB))
      allocate(Vm4(ILB:IUB,JLB:JUB))
      allocate(Vm5(ILB:IUB,JLB:JUB))
#  endif
      allocate(sustr(ILB:IUB,JLB:JUB))
      allocate(svstr(ILB:IUB,JLB:JUB))
      allocate(sustrC(ILB:IUB,JLB:JUB))
      allocate(svstrC(ILB:IUB,JLB:JUB))
      allocate(sustrA(ILB:IUB,JLB:JUB))
      allocate(svstrA(ILB:IUB,JLB:JUB))
#endif        
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
      allocate(Wm1(ILB:IUB,JLB:JUB))
      allocate(Wm2(ILB:IUB,JLB:JUB))
      allocate(Wm3(ILB:IUB,JLB:JUB))
      allocate(Wm4(ILB:IUB,JLB:JUB))
      allocate(Wm5(ILB:IUB,JLB:JUB))
#  endif
      allocate(fwf(ILB:IUB,JLB:JUB))
      allocate(fwfC(ILB:IUB,JLB:JUB))
      allocate(fwfA(ILB:IUB,JLB:JUB))
#endif  
#ifdef TIW_COUPLING
#  ifndef HCM_AI_TIW
      allocate(TIWu01(ILB:IUB,JLB:JUB))
      allocate(TIWu02(ILB:IUB,JLB:JUB))
      allocate(TIWu03(ILB:IUB,JLB:JUB))
      allocate(TIWu04(ILB:IUB,JLB:JUB))
      allocate(TIWu05(ILB:IUB,JLB:JUB))
      allocate(TIWu06(ILB:IUB,JLB:JUB))
      allocate(TIWu07(ILB:IUB,JLB:JUB))
      allocate(TIWu08(ILB:IUB,JLB:JUB))
      allocate(TIWu09(ILB:IUB,JLB:JUB))
      allocate(TIWu10(ILB:IUB,JLB:JUB))
      allocate(TIWv01(ILB:IUB,JLB:JUB))
      allocate(TIWv02(ILB:IUB,JLB:JUB))
      allocate(TIWv03(ILB:IUB,JLB:JUB))
      allocate(TIWv04(ILB:IUB,JLB:JUB))
      allocate(TIWv05(ILB:IUB,JLB:JUB))
      allocate(TIWv06(ILB:IUB,JLB:JUB))
      allocate(TIWv07(ILB:IUB,JLB:JUB))
      allocate(TIWv08(ILB:IUB,JLB:JUB))
      allocate(TIWv09(ILB:IUB,JLB:JUB))
      allocate(TIWv10(ILB:IUB,JLB:JUB))
#  else
      allocate(sstH(ILB:IUB,JLB:JUB))
#  endif
      allocate(sustrTIW(ILB:IUB,JLB:JUB))
      allocate(svstrTIW(ILB:IUB,JLB:JUB))
#endif  


!-----------------------------------------------------------------------
! SEND To MPI if MyRank/=MyMaster
!-----------------------------------------------------------------------
      if (MyRank/=MyMaster) then
       !WRITE (*,*) ' MyRank/=MyMaster '
        tile=MyRank
        LBi=BOUNDS(ng) % LBi(tile)
        UBi=BOUNDS(ng) % UBi(tile)
        LBj=BOUNDS(ng) % LBj(tile)
        UBj=BOUNDS(ng) % UBj(tile)
        IstrR=LBi+NghostPoints
        IendR=UBi-NghostPoints
        JstrR=LBj+NghostPoints
        JendR=UBj-NghostPoints
        if (LBi<NghostPoints) then
        IstrR=LBi
        endif
        if (UBi>Lm(ng)) then
        IendR=UBi
        endif
        if (LBj<NghostPoints) then
        JstrR=LBj
        endif
        if (UBj>=Mm(ng)) then
        JendR=UBj
        endif
!-----------------------------------------------------------------------
!beg weight
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
        cffU1=0.0_r8
        cffU2=0.0_r8
        cffU3=0.0_r8
        cffU4=0.0_r8
        cffU5=0.0_r8
        cffV1=0.0_r8
        cffV2=0.0_r8
        cffV3=0.0_r8
        cffV4=0.0_r8
        cffV5=0.0_r8
#  endif
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
        cffW1=0.0_r8
        cffW2=0.0_r8
        cffW3=0.0_r8
        cffW4=0.0_r8
        cffW5=0.0_r8
#  endif
#endif
#ifdef TIW_COUPLING
#  ifndef HCM_AI_TIW
        cffTIW01=0.0_r8
        cffTIW02=0.0_r8
        cffTIW03=0.0_r8
        cffTIW04=0.0_r8
        cffTIW05=0.0_r8
        cffTIW06=0.0_r8
        cffTIW07=0.0_r8
        cffTIW08=0.0_r8
        cffTIW09=0.0_r8
        cffTIW10=0.0_r8
#  endif
#endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            sstC=FORCES(ng)%sstC(i,j)
#    ifdef HCM_SST_FILTER
            sstM=FORCES(ng)%sstA(i,j)
#    else
            sstM=OCEAN(ng)%t(i,j,N(ng),nrhs(ng),itemp)-sstC
#    endif
#    ifdef HCM_INITIAL_KICK
            lonr=GRID(ng)%lonr(i,j)
            latr=GRID(ng)%latr(i,j)
            IF (tdays(ng).lt.hcm_kick_time(ng)) THEN
              cffX=lonr-hcm_kick_locx(ng)
              cffY=latr-hcm_kick_locy(ng)
              cffX=cffX/hcm_kick_lenx(ng)
              cffY=cffY/hcm_kick_leny(ng)
              cffL=EXP(-cffX*cffX)*EXP(-cffY*cffY)
              cffT=cffL*hcm_kick_temp(ng)
            ELSE
              cffT=sstM
            ENDIF
#    else
            cffT=sstM
#    endif
#    ifdef WET_DRY
            cffT=cffT*GRID(ng)%rmask_wet(i,j)
#    endif
#    ifdef MASKING
            cffT=cffT*GRID(ng)%rmask(i,j)
#    endif
            sstA(i,j)=cffT

#ifdef HCM_AI_COUPLING
            lonM(i,j)=GRID(ng)%lonr(i,j)
            latM(i,j)=GRID(ng)%latr(i,j)
#endif
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
            cffU1=cffU1+cffT*FORCES(ng)%dudtW1(i,j)
            cffU2=cffU2+cffT*FORCES(ng)%dudtW2(i,j)
            cffU3=cffU3+cffT*FORCES(ng)%dudtW3(i,j)
            cffU4=cffU4+cffT*FORCES(ng)%dudtW4(i,j)
            cffU5=cffU5+cffT*FORCES(ng)%dudtW5(i,j)
            cffV1=cffV1+cffT*FORCES(ng)%dvdtW1(i,j)
            cffV2=cffV2+cffT*FORCES(ng)%dvdtW2(i,j)
            cffV3=cffV3+cffT*FORCES(ng)%dvdtW3(i,j)
            cffV4=cffV4+cffT*FORCES(ng)%dvdtW4(i,j)
            cffV5=cffV5+cffT*FORCES(ng)%dvdtW5(i,j)
#  endif
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
            cffW1=cffW1+cffT*FORCES(ng)%dwdtW1(i,j)
            cffW2=cffW2+cffT*FORCES(ng)%dwdtW2(i,j)
            cffW3=cffW3+cffT*FORCES(ng)%dwdtW3(i,j)
            cffW4=cffW4+cffT*FORCES(ng)%dwdtW4(i,j)
            cffW5=cffW5+cffT*FORCES(ng)%dwdtW5(i,j)
#  endif
#endif
#ifdef TIW_COUPLING
            cffTIW=FORCES(ng)%sstH(i,j)
#    ifdef WET_DRY
            cffTIW=cffTIW*GRID(ng)%rmask_wet(i,j)
#    endif
#    ifdef MASKING
            cffTIW=cffTIW*GRID(ng)%rmask(i,j)
#    endif
#  ifndef HCM_AI_TIW
            cffTIW01=cffTIW01+cffTIW*FORCES(ng)%tiwT01(i,j)
            cffTIW02=cffTIW02+cffTIW*FORCES(ng)%tiwT02(i,j)
            cffTIW03=cffTIW03+cffTIW*FORCES(ng)%tiwT03(i,j)
            cffTIW04=cffTIW04+cffTIW*FORCES(ng)%tiwT04(i,j)
            cffTIW05=cffTIW05+cffTIW*FORCES(ng)%tiwT05(i,j)
            cffTIW06=cffTIW06+cffTIW*FORCES(ng)%tiwT06(i,j)
            cffTIW07=cffTIW07+cffTIW*FORCES(ng)%tiwT07(i,j)
            cffTIW08=cffTIW08+cffTIW*FORCES(ng)%tiwT08(i,j)
            cffTIW09=cffTIW09+cffTIW*FORCES(ng)%tiwT09(i,j)
            cffTIW10=cffTIW10+cffTIW*FORCES(ng)%tiwT10(i,j)
#  else
            sstH(i,j)=cffTIW
#  endif
#endif
          ENDDO
        ENDDO

        dim1=IendR-IstrR+1
        dim2=JendR-JstrR+1
        len_buffer=dim1*dim2
        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(sstA(IstrR:IendR,JstrR:JendR),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagSSTA,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#ifdef HCM_AI_COUPLING
        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(lonM(IstrR:IendR,JstrR:JendR),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagLON,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(latM(IstrR:IendR,JstrR:JendR),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagLAT,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#endif

#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
        len_buffer=5
        allocate(send_buffer(len_buffer))
        send_buffer(1)=cffU1
        send_buffer(2)=cffU2
        send_buffer(3)=cffU3
        send_buffer(4)=cffU4
        send_buffer(5)=cffU5
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagUw,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        len_buffer=5
        allocate(send_buffer(len_buffer))
        send_buffer(1)=cffV1
        send_buffer(2)=cffV2
        send_buffer(3)=cffV3
        send_buffer(4)=cffV4
        send_buffer(5)=cffV5
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagVw,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#  endif
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
        len_buffer=5
        allocate(send_buffer(len_buffer))
        send_buffer(1)=cffW1
        send_buffer(2)=cffW2
        send_buffer(3)=cffW3
        send_buffer(4)=cffW4
        send_buffer(5)=cffW5
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagWw,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#  endif
#endif
#ifdef TIW_COUPLING
#  ifndef HCM_AI_TIW
        len_buffer=10
        allocate(send_buffer(len_buffer))
        send_buffer(1)=cffTIW01
        send_buffer(2)=cffTIW02
        send_buffer(3)=cffTIW03
        send_buffer(4)=cffTIW04
        send_buffer(5)=cffTIW05
        send_buffer(6)=cffTIW06
        send_buffer(7)=cffTIW07
        send_buffer(8)=cffTIW08
        send_buffer(9)=cffTIW09
        send_buffer(10)=cffTIW10
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWw,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#  else
        len_buffer=dim1*dim2
        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(sstH(IstrR:IendR,JstrR:JendR),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagSSTH,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#  endif
#endif
!-----------------------------------------------------------------------
!beg pattern
        do i=LBi,UBi
        do j=LBj,UBj
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
         Um1(i,j)=FORCES(ng)%dudtM1(i,j)
         Um2(i,j)=FORCES(ng)%dudtM2(i,j)
         Um3(i,j)=FORCES(ng)%dudtM3(i,j)
         Um4(i,j)=FORCES(ng)%dudtM4(i,j)
         Um5(i,j)=FORCES(ng)%dudtM5(i,j)
         Vm1(i,j)=FORCES(ng)%dvdtM1(i,j)
         Vm2(i,j)=FORCES(ng)%dvdtM2(i,j)
         Vm3(i,j)=FORCES(ng)%dvdtM3(i,j)
         Vm4(i,j)=FORCES(ng)%dvdtM4(i,j)
         Vm5(i,j)=FORCES(ng)%dvdtM5(i,j)
#  endif
         sustrC(i,j)=FORCES(ng)%sustrC(i,j)
         svstrC(i,j)=FORCES(ng)%svstrC(i,j)
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
         Wm1(i,j)=FORCES(ng)%dwdtM1(i,j)
         Wm2(i,j)=FORCES(ng)%dwdtM2(i,j)
         Wm3(i,j)=FORCES(ng)%dwdtM3(i,j)
         Wm4(i,j)=FORCES(ng)%dwdtM4(i,j)
         Wm5(i,j)=FORCES(ng)%dwdtM5(i,j)
#  endif
         fwfC(i,j)=FORCES(ng)%fwfC(i,j)
#endif
#ifdef TIW_COUPLING
#  ifndef HCM_AI_TIW
         TIWu01(i,j)=FORCES(ng)%tiwU01(i,j)
         TIWu02(i,j)=FORCES(ng)%tiwU02(i,j)
         TIWu03(i,j)=FORCES(ng)%tiwU03(i,j)
         TIWu04(i,j)=FORCES(ng)%tiwU04(i,j)
         TIWu05(i,j)=FORCES(ng)%tiwU05(i,j)
         TIWu06(i,j)=FORCES(ng)%tiwU06(i,j)
         TIWu07(i,j)=FORCES(ng)%tiwU07(i,j)
         TIWu08(i,j)=FORCES(ng)%tiwU08(i,j)
         TIWu09(i,j)=FORCES(ng)%tiwU09(i,j)
         TIWu10(i,j)=FORCES(ng)%tiwU10(i,j)
         TIWv01(i,j)=FORCES(ng)%tiwV01(i,j)
         TIWv02(i,j)=FORCES(ng)%tiwV02(i,j)
         TIWv03(i,j)=FORCES(ng)%tiwV03(i,j)
         TIWv04(i,j)=FORCES(ng)%tiwV04(i,j)
         TIWv05(i,j)=FORCES(ng)%tiwV05(i,j)
         TIWv06(i,j)=FORCES(ng)%tiwV06(i,j)
         TIWv07(i,j)=FORCES(ng)%tiwV07(i,j)
         TIWv08(i,j)=FORCES(ng)%tiwV08(i,j)
         TIWv09(i,j)=FORCES(ng)%tiwV09(i,j)
         TIWv10(i,j)=FORCES(ng)%tiwV10(i,j)
#  endif
#endif
        enddo
        enddo
        dim1=UBi-LBi+1
        dim2=UBj-LBj+1
        len_buffer=dim1*dim2
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Um1(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagUm1,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Um2(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagUm2,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Um3(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagUm3,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Um4(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagUm4,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Um5(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagUm5,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Vm1(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagVm1,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Vm2(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagVm2,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Vm3(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagVm3,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Vm4(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagVm4,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Vm5(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagVm5,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#  endif
        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(sustrC(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagUc,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(svstrC(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagVc,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Wm1(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagWm1,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Wm2(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagWm2,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Wm3(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagWm3,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Wm4(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagWm4,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(Wm5(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagWm5,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#  endif
        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(fwfC(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagWc,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#endif
#ifdef TIW_COUPLING
#  ifndef HCM_AI_TIW
        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWu01(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWu01,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWu02(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWu02,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWu03(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWu03,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWu04(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWu04,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWu05(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWu05,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWu06(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWu06,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWu07(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWu07,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWu08(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWu08,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWu09(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWu09,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWu10(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWu10,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWv01(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWv01,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWv02(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWv02,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWv03(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWv03,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWv04(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWv04,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWv05(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWv05,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWv06(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWv06,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWv07(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWv07,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWv08(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWv08,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWv09(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWv09,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(TIWv10(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagTIWv10,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#  endif
#endif

!-----------------------------------------------------------------------
! END MyRank/=MyMaster
      endif


!-----------------------------------------------------------------------
! RECEIVE From MPI if MyRank==MyMaster
!-----------------------------------------------------------------------
      if (MyRank==MyMaster) then
       !WRITE (*,*) ' MyRank==MyMaster '
        LBi=BOUNDS(ng) % LBi(MyRank)
        UBi=BOUNDS(ng) % UBi(MyRank)
        LBj=BOUNDS(ng) % LBj(MyRank)
        UBj=BOUNDS(ng) % UBj(MyRank)
        IstrR=LBi+NghostPoints
        IendR=UBi-NghostPoints
        JstrR=LBj+NghostPoints
        JendR=UBj-NghostPoints
        if (LBi<NghostPoints) then
        IstrR=LBi
        endif
        if (UBi>Lm(ng)) then
        IendR=UBi
        endif
        if (LBj<NghostPoints) then
        JstrR=LBj
        endif
        if (UBj>=Mm(ng)) then
        JendR=UBj
        endif
!-----------------------------------------------------------------------
!beg weight for dudt 
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
        cffU1=0.0_r8
        cffU2=0.0_r8
        cffU3=0.0_r8
        cffU4=0.0_r8
        cffU5=0.0_r8
        cffV1=0.0_r8
        cffV2=0.0_r8
        cffV3=0.0_r8
        cffV4=0.0_r8
        cffV5=0.0_r8
#  endif
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
        cffW1=0.0_r8
        cffW2=0.0_r8
        cffW3=0.0_r8
        cffW4=0.0_r8
        cffW5=0.0_r8
#  endif
#endif
#ifdef TIW_COUPLING
#  ifndef HCM_AI_TIW
        cffTIW01=0.0_r8
        cffTIW02=0.0_r8
        cffTIW03=0.0_r8
        cffTIW04=0.0_r8
        cffTIW05=0.0_r8
        cffTIW06=0.0_r8
        cffTIW07=0.0_r8
        cffTIW08=0.0_r8
        cffTIW09=0.0_r8
        cffTIW10=0.0_r8
#  endif
#endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            sstC=FORCES(ng)%sstC(i,j)
#    ifdef HCM_SST_FILTER
            sstM=FORCES(ng)%sstA(i,j)
#    else
            sstM=OCEAN(ng)%t(i,j,N(ng),nrhs(ng),itemp)-sstC
#    endif
#    ifdef HCM_INITIAL_KICK
            lonr=GRID(ng)%lonr(i,j)
            latr=GRID(ng)%latr(i,j)
            IF (tdays(ng).lt.hcm_kick_time(ng)) THEN
              cffX=lonr-hcm_kick_locx(ng)
              cffY=latr-hcm_kick_locy(ng)
              cffX=cffX/hcm_kick_lenx(ng)
              cffY=cffY/hcm_kick_leny(ng)
              cffL=EXP(-cffX*cffX)*EXP(-cffY*cffY)
              cffT=cffL*hcm_kick_temp(ng)
            ELSE
              cffT=sstM
            ENDIF
#    else
            cffT=sstM
#    endif
#    ifdef WET_DRY
            cffT=cffT*GRID(ng)%rmask_wet(i,j)
#    endif
#    ifdef MASKING
            cffT=cffT*GRID(ng)%rmask(i,j)
#    endif
            sstA(i,j)=cffT

#ifdef HCM_AI_COUPLING
            lonM(i,j)=GRID(ng)%lonr(i,j)
            latM(i,j)=GRID(ng)%latr(i,j)
#endif
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
            cffU1=cffU1+cffT*FORCES(ng)%dudtW1(i,j)
            cffU2=cffU2+cffT*FORCES(ng)%dudtW2(i,j)
            cffU3=cffU3+cffT*FORCES(ng)%dudtW3(i,j)
            cffU4=cffU4+cffT*FORCES(ng)%dudtW4(i,j)
            cffU5=cffU5+cffT*FORCES(ng)%dudtW5(i,j)
            cffV1=cffV1+cffT*FORCES(ng)%dvdtW1(i,j)
            cffV2=cffV2+cffT*FORCES(ng)%dvdtW2(i,j)
            cffV3=cffV3+cffT*FORCES(ng)%dvdtW3(i,j)
            cffV4=cffV4+cffT*FORCES(ng)%dvdtW4(i,j)
            cffV5=cffV5+cffT*FORCES(ng)%dvdtW5(i,j)
#  endif
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
            cffW1=cffW1+cffT*FORCES(ng)%dwdtW1(i,j)
            cffW2=cffW2+cffT*FORCES(ng)%dwdtW2(i,j)
            cffW3=cffW3+cffT*FORCES(ng)%dwdtW3(i,j)
            cffW4=cffW4+cffT*FORCES(ng)%dwdtW4(i,j)
            cffW5=cffW5+cffT*FORCES(ng)%dwdtW5(i,j)
#  endif
#endif
#ifdef TIW_COUPLING
            cffTIW=FORCES(ng)%sstH(i,j)
#    ifdef WET_DRY
            cffTIW=cffTIW*GRID(ng)%rmask_wet(i,j)
#    endif
#    ifdef MASKING
            cffTIW=cffTIW*GRID(ng)%rmask(i,j)
#    endif
#  ifndef HCM_AI_TIW
            cffTIW01=cffTIW01+cffTIW*FORCES(ng)%tiwT01(i,j)
            cffTIW02=cffTIW02+cffTIW*FORCES(ng)%tiwT02(i,j)
            cffTIW03=cffTIW03+cffTIW*FORCES(ng)%tiwT03(i,j)
            cffTIW04=cffTIW04+cffTIW*FORCES(ng)%tiwT04(i,j)
            cffTIW05=cffTIW05+cffTIW*FORCES(ng)%tiwT05(i,j)
            cffTIW06=cffTIW06+cffTIW*FORCES(ng)%tiwT06(i,j)
            cffTIW07=cffTIW07+cffTIW*FORCES(ng)%tiwT07(i,j)
            cffTIW08=cffTIW08+cffTIW*FORCES(ng)%tiwT08(i,j)
            cffTIW09=cffTIW09+cffTIW*FORCES(ng)%tiwT09(i,j)
            cffTIW10=cffTIW10+cffTIW*FORCES(ng)%tiwT10(i,j)
#  else
            sstH(i,j)=cffTIW
#  endif
#endif
          ENDDO
        ENDDO

        do tile=1,Ntiles
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
          len_buffer=5
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagUw,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          cffU1=cffU1+recv_buffer(1)
          cffU2=cffU2+recv_buffer(2)
          cffU3=cffU3+recv_buffer(3)
          cffU4=cffU4+recv_buffer(4)
          cffU5=cffU5+recv_buffer(5)
          deallocate(recv_buffer)

          len_buffer=5
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagVw,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          cffV1=cffV1+recv_buffer(1)
          cffV2=cffV2+recv_buffer(2)
          cffV3=cffV3+recv_buffer(3)
          cffV4=cffV4+recv_buffer(4)
          cffV5=cffV5+recv_buffer(5)
          deallocate(recv_buffer)
#  endif
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
          len_buffer=5
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagWw,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          cffW1=cffW1+recv_buffer(1)
          cffW2=cffW2+recv_buffer(2)
          cffW3=cffW3+recv_buffer(3)
          cffW4=cffW4+recv_buffer(4)
          cffW5=cffW5+recv_buffer(5)
          deallocate(recv_buffer)
#  endif
#endif
#ifdef TIW_COUPLING
#  ifndef HCM_AI_TIW
          len_buffer=10
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWw,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          cffTIW01=cffTIW01+recv_buffer(1)
          cffTIW02=cffTIW02+recv_buffer(2)
          cffTIW03=cffTIW03+recv_buffer(3)
          cffTIW04=cffTIW04+recv_buffer(4)
          cffTIW05=cffTIW05+recv_buffer(5)
          cffTIW06=cffTIW06+recv_buffer(6)
          cffTIW07=cffTIW07+recv_buffer(7)
          cffTIW08=cffTIW08+recv_buffer(8)
          cffTIW09=cffTIW09+recv_buffer(9)
          cffTIW10=cffTIW10+recv_buffer(10)
          deallocate(recv_buffer)
#  endif
#endif
        enddo
!-----------------------------------------------------------------------
!beg pattern
        do i=LBi,UBi
        do j=LBj,UBj
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
         Um1(i,j)=FORCES(ng)%dudtM1(i,j)
         Um2(i,j)=FORCES(ng)%dudtM2(i,j)
         Um3(i,j)=FORCES(ng)%dudtM3(i,j)
         Um4(i,j)=FORCES(ng)%dudtM4(i,j)
         Um5(i,j)=FORCES(ng)%dudtM5(i,j)
         Vm1(i,j)=FORCES(ng)%dvdtM1(i,j)
         Vm2(i,j)=FORCES(ng)%dvdtM2(i,j)
         Vm3(i,j)=FORCES(ng)%dvdtM3(i,j)
         Vm4(i,j)=FORCES(ng)%dvdtM4(i,j)
         Vm5(i,j)=FORCES(ng)%dvdtM5(i,j)
#  endif
         sustrC(i,j)=FORCES(ng)%sustrC(i,j)
         svstrC(i,j)=FORCES(ng)%svstrC(i,j)
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
         Wm1(i,j)=FORCES(ng)%dwdtM1(i,j)
         Wm2(i,j)=FORCES(ng)%dwdtM2(i,j)
         Wm3(i,j)=FORCES(ng)%dwdtM3(i,j)
         Wm4(i,j)=FORCES(ng)%dwdtM4(i,j)
         Wm5(i,j)=FORCES(ng)%dwdtM5(i,j)
#  endif
         fwfC(i,j)=FORCES(ng)%fwfC(i,j)
#endif
#ifdef TIW_COUPLING
#  ifndef HCM_AI_TIW
         TIWu01(i,j)=FORCES(ng)%tiwU01(i,j)
         TIWu02(i,j)=FORCES(ng)%tiwU02(i,j)
         TIWu03(i,j)=FORCES(ng)%tiwU03(i,j)
         TIWu04(i,j)=FORCES(ng)%tiwU04(i,j)
         TIWu05(i,j)=FORCES(ng)%tiwU05(i,j)
         TIWu06(i,j)=FORCES(ng)%tiwU06(i,j)
         TIWu07(i,j)=FORCES(ng)%tiwU07(i,j)
         TIWu08(i,j)=FORCES(ng)%tiwU08(i,j)
         TIWu09(i,j)=FORCES(ng)%tiwU09(i,j)
         TIWu10(i,j)=FORCES(ng)%tiwU10(i,j)
         TIWv01(i,j)=FORCES(ng)%tiwV01(i,j)
         TIWv02(i,j)=FORCES(ng)%tiwV02(i,j)
         TIWv03(i,j)=FORCES(ng)%tiwV03(i,j)
         TIWv04(i,j)=FORCES(ng)%tiwV04(i,j)
         TIWv05(i,j)=FORCES(ng)%tiwV05(i,j)
         TIWv06(i,j)=FORCES(ng)%tiwV06(i,j)
         TIWv07(i,j)=FORCES(ng)%tiwV07(i,j)
         TIWv08(i,j)=FORCES(ng)%tiwV08(i,j)
         TIWv09(i,j)=FORCES(ng)%tiwV09(i,j)
         TIWv10(i,j)=FORCES(ng)%tiwV10(i,j)
#  endif
#endif
        enddo
        enddo

        do tile=1,Ntiles
          LBi=BOUNDS(ng) % LBi(tile)
          UBi=BOUNDS(ng) % UBi(tile)
          LBj=BOUNDS(ng) % LBj(tile)
          UBj=BOUNDS(ng) % UBj(tile)
          IstrR=LBi+NghostPoints
          IendR=UBi-NghostPoints
          JstrR=LBj+NghostPoints
          JendR=UBj-NghostPoints
          if (LBi<NghostPoints) then
          IstrR=LBi
          endif
          if (UBi>Lm(ng)) then
          IendR=UBi
          endif
          if (LBj<NghostPoints) then
          JstrR=LBj
          endif
          if (UBj>=Mm(ng)) then
          JendR=UBj
          endif
          dim1=IendR-IstrR+1
          dim2=JendR-JstrR+1
          len_buffer=dim1*dim2
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagSSTA,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          sstA(IstrR:IendR,JstrR:JendR)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
#ifdef HCM_AI_COUPLING
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagLON,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          lonM(IstrR:IendR,JstrR:JendR)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagLAT,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          latM(IstrR:IendR,JstrR:JendR)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
#endif

          dim1=UBi-LBi+1
          dim2=UBj-LBj+1
          len_buffer=dim1*dim2
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagUm1,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Um1(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagUm2,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Um2(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagUm3,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Um3(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagUm4,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Um4(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagUm5,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Um5(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagVm1,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Vm1(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagVm2,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Vm2(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagVm3,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Vm3(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagVm4,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Vm4(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagVm5,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Vm5(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
#  endif
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagUc,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          sustrC(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagVc,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          svstrC(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagWm1,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Wm1(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagWm2,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Wm2(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagWm3,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Wm3(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagWm4,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Wm4(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagWm5,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          Wm5(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
#  endif
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagWc,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          fwfC(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
#endif
#ifdef TIW_COUPLING
#  ifndef HCM_AI_TIW
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWu01,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWu01(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWu02,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWu02(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWu03,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWu03(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWu04,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWu04(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWu05,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWu05(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWu06,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWu06(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWu07,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWu07(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWu08,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWu08(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWu09,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWu09(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWu10,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWu10(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWv01,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWv01(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWv02,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWv02(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWv03,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWv03(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWv04,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWv04(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWv05,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWv05(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWv06,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWv06(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWv07,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWv07(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWv08,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWv08(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWv09,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWv09(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagTIWv10,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          TIWv10(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
#  else
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagSSTH,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          sstH(IstrR:IendR,JstrR:JendR)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
#  endif
#endif
        enddo


! END MyRank==MyMaster
      endif
!-----------------------------------------------------------------------
      call  mpi_barrier(MPI_COMM_WORLD, ierr)
!-----------------------------------------------------------------------




! Calculate HCM forcing and scatter to different ranks
      if (MyRank==MyMaster) then
        LBi=BOUNDS(ng) % LBi(MyRank)
        UBi=BOUNDS(ng) % UBi(MyRank)
        LBj=BOUNDS(ng) % LBj(MyRank)
        UBj=BOUNDS(ng) % UBj(MyRank)
!=======================================================================
        IF (MOD(time(ng),hcm_coup_time(ng)).eq.0) THEN
          WRITE (*,*) ' Hybrid Coupled Model (HCM) coupling '
#ifdef HCM_AI_COUPLING
          COUPnum=time(ng)/hcm_coup_time(ng)
          CALL time_string (time(ng),fileT)
          OPEN  (1101,file='./COUP/COUPnum.txt')
          WRITE (1101,*) COUPnum
          CLOSE (1101)
          OPEN  (1102,file='./COUP/LON_'//fileT//'.txt')
          WRITE (1102,"(F15.6)") lonM 
          CLOSE (1102)
          OPEN  (1103,file='./COUP/LAT_'//fileT//'.txt')
          WRITE (1103,"(F15.6)") latM 
          CLOSE (1103)
          OPEN  (1104,file='./COUP/SSTA_'//fileT//'.txt')
          WRITE (1104,"(F15.6)") sstA 
          CLOSE (1104)
#  ifdef HCM_AI_TIW
          OPEN  (1105,file='./COUP/SSTH_'//fileT//'.txt')
          WRITE (1105,"(F15.6)") sstH
          CLOSE (1105)
#  endif
          OPEN  (1199,file='./COUP/TIME_'//fileT//'.txt')
          WRITE (1199,*) time(ng)
          CLOSE (1199)
#endif
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
          OPEN  (1201,file='./Build/hcm_dudt_weight.dat')
          WRITE (1201,*) cffU1, cffU2, cffU3, cffU4, cffU5
          CLOSE (1201)
          WRITE (*,*) ' HCM weights for U wind stress'
          WRITE (*,*) ' Mod 1:', cffU1, ' Mod 2:', cffU2
          WRITE (*,*) ' Mod 3:', cffU3, ' Mod 4:', cffU4
          WRITE (*,*) ' Mod 5:', cffU5

          OPEN  (1202,file='./Build/hcm_dvdt_weight.dat')
          WRITE (1202,*) cffV1, cffV2, cffV3, cffV4, cffV5
          CLOSE (1202)
          WRITE (*,*) ' HCM weights for V wind stress'
          WRITE (*,*) ' Mod 1:', cffV1, ' Mod 2:', cffV2
          WRITE (*,*) ' Mod 3:', cffV3, ' Mod 4:', cffV4
          WRITE (*,*) ' Mod 5:', cffV5
#  endif
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
          OPEN  (1203,file='./Build/hcm_dwdt_weight.dat')
          WRITE (1203,*) cffW1, cffW2, cffW3, cffW4, cffW5
          CLOSE (1203)
          WRITE (*,*) ' HCM weights for fresh water flux (E-P)'
          WRITE (*,*) ' Mod 1:', cffW1, ' Mod 2:', cffW2
          WRITE (*,*) ' Mod 3:', cffW3, ' Mod 4:', cffW4
          WRITE (*,*) ' Mod 5:', cffW5
#  endif
#endif
#ifdef TIW_COUPLING
#  ifndef HCM_AI_TIW
          OPEN  (1204,file='./Build/hcm_tiw_weight.dat')
          WRITE (1204,*) cffTIW01, cffTIW02, cffTIW03, cffTIW04, cffTIW05,  &
     &                   cffTIW06, cffTIW07, cffTIW08, cffTIW09, cffTIW10
          CLOSE (1204)
          WRITE (*,*) ' HCM weights for TIW wind stress'
          WRITE (*,*) ' Mod  1:', cffTIW01, ' Mod  2:', cffTIW02
          WRITE (*,*) ' Mod  3:', cffTIW03, ' Mod  4:', cffTIW04
          WRITE (*,*) ' Mod  5:', cffTIW05, ' Mod  6:', cffTIW06
          WRITE (*,*) ' Mod  7:', cffTIW07, ' Mod  8:', cffTIW08
          WRITE (*,*) ' Mod  9:', cffTIW09, ' Mod 10:', cffTIW10
#  endif
#endif
        ELSE
#ifdef HCM_AI_COUPLING
          OPEN  (1101,file='./COUP/COUPnum.txt')
          READ (1101,*) COUPnum
          CLOSE (1101)
          CALL time_string (hcm_coup_time(ng)*COUPnum,fileT)
#endif
#ifdef HCM_WIND_STRESS
#  ifndef HCM_AI_WIND
          cffU1=0.0_r8
          cffU2=0.0_r8
          cffU3=0.0_r8
          cffU4=0.0_r8
          cffU5=0.0_r8
          OPEN  (2101,file='./Build/hcm_dudt_weight.dat')
          READ  (2101,*) cffU1, cffU2, cffU3, cffU4, cffU5
          CLOSE (2101)

          cffV1=0.0_r8
          cffV2=0.0_r8
          cffV3=0.0_r8
          cffV4=0.0_r8
          cffV5=0.0_r8
          OPEN  (2102,file='./Build/hcm_dvdt_weight.dat')
          READ  (2102,*) cffV1, cffV2, cffV3, cffV4, cffV5
          CLOSE (2102)
#  endif
#endif
#ifdef HCM_EMINP_FORCE
#  ifndef HCM_AI_FWF
          cffW1=0.0_r8
          cffW2=0.0_r8
          cffW3=0.0_r8
          cffW4=0.0_r8
          cffW5=0.0_r8
          OPEN  (2103,file='./Build/hcm_dwdt_weight.dat')
          READ  (2103,*) cffW1, cffW2, cffW3, cffW4, cffW5
          CLOSE (2103)
#  endif
#endif
#ifdef TIW_COUPLING
#  ifndef HCM_AI_TIW
          cffTIW01=0.0_r8
          cffTIW02=0.0_r8
          cffTIW03=0.0_r8
          cffTIW04=0.0_r8
          cffTIW05=0.0_r8
          cffTIW06=0.0_r8
          cffTIW07=0.0_r8
          cffTIW08=0.0_r8
          cffTIW09=0.0_r8
          cffTIW10=0.0_r8
          OPEN  (2104,file='./Build/hcm_tiw_weight.dat')
          READ  (2104,*) cffTIW01, cffTIW02, cffTIW03, cffTIW04, cffTIW05,  &
     &                   cffTIW06, cffTIW07, cffTIW08, cffTIW09, cffTIW10
          CLOSE (2104)
#  endif
#endif
        ENDIF
!-----------------------------------------------------------------------
#ifdef HCM_AI_COUPLING
        AI_NUM=1.0_r8
        inquire(file='./COUP/FINISH_'//fileT//'.txt',exist=AI_FINISH)
        DO WHILE (AI_FINISH == 0)
        inquire(file='./COUP/FINISH_'//fileT//'.txt',exist=AI_FINISH)
        AI_NUM=AI_NUM+1.0_r8
        ENDDO
#endif

#ifdef HCM_WIND_STRESS
#  ifdef HCM_AI_WIND
        fileN='./COUP/sustrA_'//fileT//'.txt'
        WRITE(*,*) ' Load Data From: ',fileN
        OPEN  (3101,file=fileN)
        do j=JLB,JUB
          do i=ILB,IUB
          AIDT(i,j)=0.0_r8
          READ(3101,"(F15.6)") AIDT(i,j)
          enddo
        enddo
        CLOSE (3101)
#  endif
        do i=ILB+1,IUB
          do j=JLB,JUB
#  ifndef HCM_AI_WIND
            cffM=cffU1*Um1(i-1,j)+cffU2*Um2(i-1,j)+cffU3*Um3(i-1,j)     &
     &          +cffU4*Um4(i-1,j)+cffU5*Um5(i-1,j)
            cffR=cffU1*Um1(i,j)+cffU2*Um2(i,j)+cffU3*Um3(i,j)           &
     &          +cffU4*Um4(i,j)+cffU5*Um5(i,j)
#     ifdef HCM_INITIAL_KICK
        IF (tdays(ng).lt.hcm_kick_time(ng)) THEN
          sustrA(i,j)=0.5_r8*hcm_alpha_tau(ng)*(cffM+cffR)
        ELSE
#       ifdef HCM_SST_FILTER
#         ifdef HCM_SSTA_REINFORCE
          sustrA(i,j)=0.5_r8*PI_half*hcm_alpha_tau(ng)*(cffM+cffR)
#         else
          sustrA(i,j)=0.5_r8*hcm_alpha_tau(ng)*(cffM+cffR)
#         endif
#       else
          sustrA(i,j)=0.5_r8*hcm_alpha_tau(ng)*(cffM+cffR)
#       endif
        ENDIF
#     else
#       ifdef HCM_SST_FILTER
#         ifdef HCM_SSTA_REINFORCE
          sustrA(i,j)=0.5_r8*PI_half*hcm_alpha_tau(ng)*(cffM+cffR)
#         else
          sustrA(i,j)=0.5_r8*hcm_alpha_tau(ng)*(cffM+cffR)
#         endif
#       else
          sustrA(i,j)=0.5_r8*hcm_alpha_tau(ng)*(cffM+cffR)
#       endif
#     endif
#  else
          cffM=AIDT(i-1,j)
          cffR=AIDT(i,j)
          sustrA(i,j)=0.5_r8*(cffM+cffR)
#  endif
          sustr(i,j)=sustrC(i,j)+sustrA(i,j)
          sustr(i,j)=sustr(i,j)/rho0
          enddo
        enddo
#  ifdef HCM_AI_WIND
        fileN='./COUP/svstrA_'//fileT//'.txt'
        WRITE(*,*) ' Load Data From: ',fileN
        OPEN  (3102,file=fileN)
        do j=JLB,JUB
          do i=ILB,IUB
          AIDT(i,j)=0.0_r8
          READ(3102,"(F15.6)") AIDT(i,j)
          enddo
        enddo
        CLOSE (3102)
#  endif
        do i=ILB,IUB
          do j=JLB+1,JUB
#  ifndef HCM_AI_WIND
            cffM=cffV1*Vm1(i,j-1)+cffV2*Vm2(i,j-1)+cffV3*Vm3(i,j-1)     &
     &          +cffV4*Vm4(i,j-1)+cffV5*Vm5(i,j-1)
            cffR=cffV1*Vm1(i,j)+cffV2*Vm2(i,j)+cffV3*Vm3(i,j)           &
     &          +cffV4*Vm4(i,j)+cffV5*Vm5(i,j)
#     ifdef HCM_INITIAL_KICK
        IF (tdays(ng).lt.hcm_kick_time(ng)) THEN
          svstrA(i,j)=0.5_r8*hcm_alpha_tau(ng)*(cffM+cffR)
        ELSE
#       ifdef HCM_SST_FILTER
#         ifdef HCM_SSTA_REINFORCE
          svstrA(i,j)=0.5_r8*PI_half*hcm_alpha_tau(ng)*(cffM+cffR)
#         else
          svstrA(i,j)=0.5_r8*hcm_alpha_tau(ng)*(cffM+cffR)
#         endif
#       else
          svstrA(i,j)=0.5_r8*hcm_alpha_tau(ng)*(cffM+cffR)
#       endif
        ENDIF
#     else
#       ifdef HCM_SST_FILTER
#         ifdef HCM_SSTA_REINFORCE
          svstrA(i,j)=0.5_r8*PI_half*hcm_alpha_tau(ng)*(cffM+cffR)
#         else
          svstrA(i,j)=0.5_r8*hcm_alpha_tau(ng)*(cffM+cffR)
#         endif
#       else
          svstrA(i,j)=0.5_r8*hcm_alpha_tau(ng)*(cffM+cffR)
#       endif
#     endif
#  else
          cffM=AIDT(i,j-1)
          cffR=AIDT(i,j)
          svstrA(i,j)=0.5_r8*(cffM+cffR)
#  endif
          svstr(i,j)=svstrC(i,j)+svstrA(i,j)
          svstr(i,j)=svstr(i,j)/rho0
          enddo
        enddo
        FORCES(ng)%sustrA(LBi:UBi,LBj:UBj)=sustrA(LBi:UBi,LBj:UBj)
        FORCES(ng)%svstrA(LBi:UBi,LBj:UBj)=svstrA(LBi:UBi,LBj:UBj)
        FORCES(ng)%sustr(LBi:UBi,LBj:UBj)=sustr(LBi:UBi,LBj:UBj)
        FORCES(ng)%svstr(LBi:UBi,LBj:UBj)=svstr(LBi:UBi,LBj:UBj)
#endif
#ifdef HCM_EMINP_FORCE
#  ifdef HCM_AI_FWF
        fileN='./COUP/FWFA_'//fileT//'.txt'
        WRITE(*,*) ' Load Data From: ',fileN
        OPEN  (3103,file=fileN)
        do j=JLB,JUB
          do i=ILB,IUB
          AIDT(i,j)=0.0_r8
          READ(3103,"(F15.6)") AIDT(i,j)
          enddo
        enddo
        CLOSE (3103)
#  endif
        do i=ILB,IUB
          do j=JLB,JUB
#  ifndef HCM_AI_FWF
            cffR=cffW1*Wm1(i,j)+cffW2*Wm2(i,j)+cffW3*Wm3(i,j)           &
     &          +cffW4*Wm4(i,j)+cffW5*Wm5(i,j)
#     ifdef HCM_INITIAL_KICK
        IF (tdays(ng).lt.hcm_kick_time(ng)) THEN
          fwfA(i,j)=hcm_alpha_fwf(ng)*cffR
        ELSE
#       ifdef HCM_SST_FILTER
#         ifdef HCM_SSTA_REINFORCE
          fwfA(i,j)=PI_half*hcm_alpha_fwf(ng)*cffR
#         else
          fwfA(i,j)=hcm_alpha_fwf(ng)*cffR
#         endif
#       else
          fwfA(i,j)=hcm_alpha_fwf(ng)*cffR
#       endif
        ENDIF
#     else
#       ifdef HCM_SST_FILTER
#         ifdef HCM_SSTA_REINFORCE
          fwfA(i,j)=PI_half*hcm_alpha_fwf(ng)*cffR
#         else
          fwfA(i,j)=hcm_alpha_fwf(ng)*cffR
#         endif
#       else
          fwfA(i,j)=hcm_alpha_fwf(ng)*cffR
#       endif
#     endif
#  else
          fwfA(i,j)=AIDT(i,j)
#  endif
          fwf(i,j)=fwfC(i,j)+fwfA(i,j)
          enddo
        enddo
        FORCES(ng)%fwfA(LBi:UBi,LBj:UBj)=fwfA(LBi:UBi,LBj:UBj)
        FORCES(ng)%stflx(LBi:UBi,LBj:UBj,isalt)=fwf(LBi:UBi,LBj:UBj)
#endif
#ifdef TIW_COUPLING
#  ifdef HCM_AI_TIW
        fileN='./COUP/sustrTIW_'//fileT//'.txt'
        WRITE(*,*) ' Load Data From: ',fileN
        OPEN  (3104,file=fileN)
        do j=JLB,JUB
          do i=ILB,IUB
          AIDT(i,j)=0.0_r8
          READ(3104,"(F15.6)") AIDT(i,j)
          enddo
        enddo
        CLOSE (3104)
#  endif
        do i=ILB+1,IUB
          do j=JLB,JUB
#  ifndef HCM_AI_TIW
            cffM=cffTIW01*TIWu01(i-1,j)+cffTIW02*TIWu02(i-1,j)          &
     &          +cffTIW03*TIWu03(i-1,j)+cffTIW04*TIWu04(i-1,j)          &
     &          +cffTIW05*TIWu05(i-1,j)+cffTIW06*TIWu06(i-1,j)          &
     &          +cffTIW07*TIWu07(i-1,j)+cffTIW08*TIWu08(i-1,j)          &
     &          +cffTIW09*TIWu09(i-1,j)+cffTIW10*TIWu10(i-1,j)
            cffR=cffTIW01*TIWu01(i,j)+cffTIW02*TIWu02(i,j)              &
     &          +cffTIW03*TIWu03(i,j)+cffTIW04*TIWu04(i,j)              &
     &          +cffTIW05*TIWu05(i,j)+cffTIW06*TIWu06(i,j)              &
     &          +cffTIW07*TIWu07(i,j)+cffTIW08*TIWu08(i,j)              &
     &          +cffTIW09*TIWu09(i,j)+cffTIW10*TIWu10(i,j)
          sustrTIW(i,j)=0.5_r8*hcm_alpha_tiw(ng)*(cffM+cffR)
#  else
          cffM=AIDT(i-1,j)
          cffR=AIDT(i,j)
          sustrTIW(i,j)=0.5_r8*(cffM+cffR)
#  endif
          enddo
        enddo
#  ifdef HCM_AI_TIW
        fileN='./COUP/svstrTIW_'//fileT//'.txt'
        WRITE(*,*) ' Load Data From: ',fileN
        OPEN  (3105,file=fileN)
        do j=JLB,JUB
          do i=ILB,IUB
          AIDT(i,j)=0.0_r8
          READ(3105,"(F15.6)") AIDT(i,j)
          enddo
        enddo
        CLOSE (3105)
#  endif
        do i=ILB,IUB
          do j=JLB+1,JUB
#  ifndef HCM_AI_TIW
            cffM=cffTIW01*TIWv01(i,j-1)+cffTIW02*TIWv02(i,j-1)          &
     &          +cffTIW03*TIWv03(i,j-1)+cffTIW04*TIWv04(i,j-1)          &
     &          +cffTIW05*TIWv05(i,j-1)+cffTIW06*TIWv06(i,j-1)          &
     &          +cffTIW07*TIWv07(i,j-1)+cffTIW08*TIWv08(i,j-1)          &
     &          +cffTIW09*TIWv09(i,j-1)+cffTIW10*TIWv10(i,j-1)
            cffR=cffTIW01*TIWv01(i,j)+cffTIW02*TIWv02(i,j)              &
     &          +cffTIW03*TIWv03(i,j)+cffTIW04*TIWv04(i,j)              &
     &          +cffTIW05*TIWv05(i,j)+cffTIW06*TIWv06(i,j)              &
     &          +cffTIW07*TIWv07(i,j)+cffTIW08*TIWv08(i,j)              &
     &          +cffTIW09*TIWv09(i,j)+cffTIW10*TIWv10(i,j)
          svstrTIW(i,j)=0.5_r8*hcm_alpha_tiw(ng)*(cffM+cffR)
#  else
          cffM=AIDT(i,j-1)
          cffR=AIDT(i,j)
          svstrTIW(i,j)=0.5_r8*(cffM+cffR)
#  endif
          enddo
        enddo
        FORCES(ng)%sustrTIW(LBi:UBi,LBj:UBj)=sustrTIW(LBi:UBi,LBj:UBj)
        FORCES(ng)%svstrTIW(LBi:UBi,LBj:UBj)=svstrTIW(LBi:UBi,LBj:UBj)
#endif
!-----------------------------------------------------------------------
        FORCES(ng)%sstA(LBi:UBi,LBj:UBj)=sstA(LBi:UBi,LBj:UBj)
!-----------------------------------------------------------------------
        do tile=1,Ntiles
          LBi=BOUNDS(ng) % LBi(tile)
          UBi=BOUNDS(ng) % UBi(tile)
          LBj=BOUNDS(ng) % LBj(tile)
          UBj=BOUNDS(ng) % UBj(tile)
          dim1=UBi-LBi+1
          dim2=UBj-LBj+1
          len_buffer=dim1*dim2
          allocate(send_buffer(len_buffer))
          send_buffer=reshape(sstA(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagSSTN,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)
#ifdef HCM_WIND_STRESS
          allocate(send_buffer(len_buffer))
          send_buffer=reshape(sustr(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagU,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)

          allocate(send_buffer(len_buffer))
          send_buffer=reshape(svstr(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagV,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)

          allocate(send_buffer(len_buffer))
          send_buffer=reshape(sustrA(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagUa,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)

          allocate(send_buffer(len_buffer))
          send_buffer=reshape(svstrA(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagVa,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)
#endif
#ifdef HCM_EMINP_FORCE
          allocate(send_buffer(len_buffer))
          send_buffer=reshape(fwf(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagW,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)

          allocate(send_buffer(len_buffer))
          send_buffer=reshape(fwfA(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagWa,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)
#endif
#ifdef TIW_COUPLING
          allocate(send_buffer(len_buffer))
          send_buffer=reshape(sustrTIW(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagTIWu,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)

          allocate(send_buffer(len_buffer))
          send_buffer=reshape(svstrTIW(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagTIWv,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)
#endif
        enddo
!=======================================================================
      else
        LBi=BOUNDS(ng) % LBi(MyRank)
        UBi=BOUNDS(ng) % UBi(MyRank)
        LBj=BOUNDS(ng) % LBj(MyRank)
        UBj=BOUNDS(ng) % UBj(MyRank)
        dim1=UBi-LBi+1
        dim2=UBj-LBj+1
        len_buffer=dim1*dim2
        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagSSTN,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%sstA(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)
#ifdef HCM_WIND_STRESS
        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagU,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%sustr(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)

        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagV,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%svstr(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)

        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagUa,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%sustrA(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)

        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagVa,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%svstrA(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)
#endif
#ifdef HCM_EMINP_FORCE
        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagW,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%stflx(LBi:UBi,LBj:UBj,isalt)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)

        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagWa,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%fwfA(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)
#endif
#ifdef TIW_COUPLING
        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagTIWu,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%sustrTIW(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)

        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagTIWv,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%svstrTIW(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)
#endif

      endif
!-----------------------------------------------------------------------
      call  mpi_barrier(MPI_COMM_WORLD, ierr)
!-----------------------------------------------------------------------
!


      END SUBROUTINE hcm_forcing
#    endif
#endif

      END MODULE hcm_forc_mod
