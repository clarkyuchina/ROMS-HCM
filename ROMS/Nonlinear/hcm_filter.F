#include "cppdefs.h"
      MODULE hcm_filt_mod
#ifdef HCM_COUPLING
#  ifdef HCM_SST_FILTER
!
!svn $Id: set_vbc.F 889 2018-02-10 03:32:52Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module sets vertical boundary conditons for momentum and       !
!  tracers.                                                            !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: hcm_filter
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE hcm_filter (ng)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_grid
      USE mod_parallel
      USE mod_forces
      USE mod_ocean
      USE mod_stepping
      USE mod_boundary
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
      integer :: tagT,tagF,tagH,tagL,tagX,tagY,tagM
      integer,dimension(mpi_status_size) :: rstatus
      real(r8),allocatable,dimension(:) ::send_buffer
      real(r8),allocatable,dimension(:) ::recv_buffer
      real(r8),allocatable,dimension(:,:) :: sstF, sstH, sstL, mask
      real(r8),allocatable,dimension(:,:) :: lonM, latM, sstM, M 
      integer :: i, j, loopN, LN, BN
      integer :: LBi,UBi,LBj,UBj,dim1,dim2
      integer :: ILB,IUB,JLB,JUB

#  ifdef HCM_ROSSBY_FILTER
      integer :: tagDX,tagDY,tagRb
      integer :: ir,jr,RX,RY
      real(r8):: RaT,Vdd,Cdd,Ndd
      real(r8),allocatable,dimension(:,:) :: disX, disY, disR 
#  endif
#  ifdef HCM_LOESS_FILTER
      integer ::nx,ny
      real(r8)::span_x,span_y
      real(r8),allocatable,dimension(:) ::xin,yin
      real(r8):: Vdd,Cdd
#  endif
#  ifdef HCM_RUNNING_FILTER
      integer :: ir,jr,RX,RY
      real(r8):: Vdd,Cdd,Ndd
#  endif
#  ifndef HCM_CLIM_SPIN
      integer :: tagC,tagA
      real(r8),allocatable,dimension(:,:) :: sstC, sstA
#    ifdef HCM_INITIAL_KICK
      real(r8):: cffX, cffY, cffL
#    endif
#    ifdef SSTA_FILTER_OPT1
      real(r8):: VddA,CddA
      real(r8),allocatable,dimension(:,:) :: sstP
#    endif
#  endif      


      Ntiles=NtileI(ng)*NtileJ(ng)-1

      tagT=10
      tagF=11
      tagH=12
      tagL=13
      tagX=97
      tagY=98
      tagM=99
#  ifdef HCM_ROSSBY_FILTER
      tagDX=21
      tagDY=22
      tagRb=23
#  endif
#  ifndef HCM_CLIM_SPIN
      tagC=31
      tagA=32
#  endif
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(tile)
      IUB=BOUNDS(ng)%UBi(tile)
      JLB=BOUNDS(ng)%LBj(tile)
      JUB=BOUNDS(ng)%UBj(tile)

      allocate(sstF(ILB:IUB,JLB:JUB))
      allocate(sstH(ILB:IUB,JLB:JUB))
      allocate(sstL(ILB:IUB,JLB:JUB))
      allocate(mask(ILB:IUB,JLB:JUB))

      allocate(lonM(ILB:IUB,JLB:JUB))
      allocate(latM(ILB:IUB,JLB:JUB))
      allocate(sstM(ILB:IUB,JLB:JUB))
      allocate(   M(ILB:IUB,JLB:JUB))
        
#  ifdef HCM_ROSSBY_FILTER
      allocate(disX(ILB:IUB,JLB:JUB))
      allocate(disY(ILB:IUB,JLB:JUB))
      allocate(disR(ILB:IUB,JLB:JUB))
#  endif

#  ifndef HCM_CLIM_SPIN
      allocate(sstC(ILB:IUB,JLB:JUB))
      allocate(sstA(ILB:IUB,JLB:JUB))
#    ifdef SSTA_FILTER_OPT1
      allocate(sstP(ILB:IUB,JLB:JUB))
#    endif
#  endif

!=======================================================================
! SEND SST To MPI if MyRank/=MyMaster
!=======================================================================
      if (MyRank/=MyMaster) then
       !WRITE (*,*) ' MyRank/=MyMaster '
        tile=MyRank
        LBi=BOUNDS(ng) % LBi(tile)
        UBi=BOUNDS(ng) % UBi(tile)
        LBj=BOUNDS(ng) % LBj(tile)
        UBj=BOUNDS(ng) % UBj(tile)
        do i=LBi,UBi
        do j=LBj,UBj
#    ifdef HCM_ROSSBY_FILTER
        disX(i,j)=1.0_r8/GRID(ng)%pm(i,j)
        disY(i,j)=1.0_r8/GRID(ng)%pn(i,j)
        disR(i,j)=FORCES(ng)%Rrob(i,j)
#    endif
        lonM(i,j)=GRID(ng)%lonr(i,j)
        latM(i,j)=GRID(ng)%latr(i,j)
        sstM(i,j)=OCEAN(ng)%t(i,j,N(ng),nrhs(ng),itemp)
#    ifndef HCM_CLIM_SPIN
        sstC(i,j)=FORCES(ng)%sstC(i,j)
#    endif
           M(i,j)=1.0_r8
#    ifdef WET_DRY
           M(i,j)=M(i,j)*GRID(ng)%rmask_wet(i,j)
#    endif
#    ifdef MASKING
           M(i,j)=M(i,j)*GRID(ng)%rmask(i,j)
#    endif
        enddo
        enddo

        dim1=UBi-LBi+1
        dim2=UBj-LBj+1
        len_buffer=dim1*dim2
#    ifdef HCM_ROSSBY_FILTER
        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(disX(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagDX,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(disY(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagDY,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(disR(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagRb,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#    endif

#    ifndef HCM_CLIM_SPIN
        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(sstC(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagC,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
#    endif

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(sstM(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagT,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(lonM(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagX,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(latM(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagY,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)

        allocate(send_buffer(len_buffer))
        send_buffer(len_buffer)=0.0
        send_buffer=reshape(M(LBi:UBi,LBj:UBj),(/dim1*dim2/))
        call mpi_isend(send_buffer,len_buffer,MPI_REAL8,0,tagM,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        deallocate(send_buffer)
       !END MyRank/=MyMaster
      endif
!=======================================================================
! RECEIVE SST From MPI if MyRank==MyMaster
!=======================================================================
      if (MyRank==MyMaster) then
       !WRITE (*,*) ' MyRank==MyMaster '
        LBi=BOUNDS(ng) % LBi(MyRank)
        UBi=BOUNDS(ng) % UBi(MyRank)
        LBj=BOUNDS(ng) % LBj(MyRank)
        UBj=BOUNDS(ng) % UBj(MyRank)
        do i=LBi,UBi
        do j=LBj,UBj
#    ifdef HCM_ROSSBY_FILTER
        disX(i,j)=1.0_r8/GRID(ng)%pm(i,j)
        disY(i,j)=1.0_r8/GRID(ng)%pn(i,j)
        disR(i,j)=FORCES(ng)%Rrob(i,j)
#    endif
        lonM(i,j)=GRID(ng)%lonr(i,j)
        latM(i,j)=GRID(ng)%latr(i,j)
        sstM(i,j)=OCEAN(ng)%t(i,j,N(ng),nrhs(ng),itemp)
#    ifndef HCM_CLIM_SPIN
        sstC(i,j)=FORCES(ng)%sstC(i,j)
#    endif
           M(i,j)=1.0_r8
#    ifdef WET_DRY
           M(i,j)=M(i,j)*GRID(ng)%rmask_wet(i,j)
#    endif
#    ifdef MASKING
           M(i,j)=M(i,j)*GRID(ng)%rmask(i,j)
#    endif
        enddo
        enddo

        do tile=1,Ntiles
          LBi=BOUNDS(ng) % LBi(tile)
          UBi=BOUNDS(ng) % UBi(tile)
          LBj=BOUNDS(ng) % LBj(tile)
          UBj=BOUNDS(ng) % UBj(tile)
          dim1=UBi-LBi+1
          dim2=UBj-LBj+1
          len_buffer=dim1*dim2
#    ifdef HCM_ROSSBY_FILTER
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagDX,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          disX(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagDY,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          disY(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagRb,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          disR(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
#    endif

#    ifndef HCM_CLIM_SPIN
          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagC,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          sstC(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
#    endif

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagT,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          sstM(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagX,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          lonM(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagY,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          latM(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)

          allocate(recv_buffer(len_buffer))
          recv_buffer(len_buffer)=0.0
          call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,tile,tagM,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
             M(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
          deallocate(recv_buffer)
        enddo
       ! END MyRank==MyMaster
      endif
!=======================================================================
      call  mpi_barrier(MPI_COMM_WORLD, ierr)
!=======================================================================
!


! Apply SST filter to SST if MyRank==MyMaster
!=======================================================================
      if (MyRank==MyMaster) then
        LN = floor(hcm_filt_loop(ng))
        BN = floor(hcm_filt_buff(ng))

        do i=ILB,IUB
          do j=JLB,JUB
            sstF(i,j)=sstM(i,j)
            mask(i,j)=   M(i,j)
          enddo
        enddo

        do i=ILB,IUB
          do j=JLB,JUB
            sstL(i,j)=sstM(i,j)
          enddo
        enddo

#ifndef HCM_CLIM_SPIN
        do i=ILB,IUB
          do j=JLB,JUB
            sstA(i,j)=sstM(i,j)-sstC(i,j)
#         ifdef SSTA_FILTER_OPT1
            sstP(i,j)=sstM(i,j)-sstC(i,j)
#         endif
          enddo
        enddo
#endif

#ifdef HCM_ROSSBY_FILTER 
        do loopN=1,LN
          do i=ILB+BN,IUB-BN
          do j=JLB+BN,JUB-BN
            RX = nint(hcm_filt_lenx(ng)*disR(i,j)/disX(i,j))
            RY = nint(hcm_filt_leny(ng)*disR(i,j)/disY(i,j))
            if(RX>0.or.RY>0)then
              sstM(i,j)=0.0_r8
              Vdd=0.0_r8
              Cdd=0.0_r8
              Ndd=0.0_r8
#         ifdef SSTA_FILTER_OPT1
              sstP(i,j)=0.0_r8
              VddA=0.0_r8
              CddA=0.0_r8
#         endif      
              do ir=i-RX,i+RX
              do jr=j-RY,j+RY
                Ndd=Ndd+1.0_r8
                if(RX>0.and.RY>0)then
                  RaT = sqrt(((ir-i)/RX*1.0_r8)**2 + ((jr-j)/RY*1.0_r8)**2)
                else
                  if(RX>0)then
                    RaT = (ir-i)/RX*1.0_r8
                  endif
                  if(RY>0)then
                    RaT = (jr-j)/RY*1.0_r8
                  endif
                endif
                if(ir<ILB.or.ir>IUB.or.jr<JLB.or.jr>JUB.or.RaT>1.0_r8)then
                  Vdd=Vdd+0.0_r8
                  Cdd=Cdd+0.0_r8
#         ifdef SSTA_FILTER_OPT1
                  VddA=VddA+0.0_r8
                  CddA=CddA+0.0_r8
#         endif      
                else
                  if(sstL(ir,jr) .ne. -999)then
                    Vdd=Vdd+sstL(ir,jr)*mask(ir,jr)
                    Cdd=Cdd+mask(ir,jr)
#         ifdef SSTA_FILTER_OPT1
                    VddA=VddA+sstA(ir,jr)*mask(ir,jr)
                    CddA=CddA+mask(ir,jr)
#         endif 
                  else
                    Vdd=Vdd+0.0_r8
                    Cdd=Cdd+0.0_r8
#         ifdef SSTA_FILTER_OPT1
                    VddA=VddA+0.0_r8
                    CddA=CddA+0.0_r8
#         endif      
                  endif     
                endif
              enddo
              enddo
              if(Cdd.gt.0.0_r8) then
                sstM(i,j)=Vdd/Cdd
#         ifdef SSTA_FILTER_OPT1
                sstP(i,j)=VddA/CddA
#         endif
              endif
              if(Cdd*5.0 .lt. Ndd*4.0)then
                sstM(i,j)=-999
              endif
            endif
          enddo
          enddo
          do i=ILB,IUB
            do j=JLB,JUB
              sstL(i,j)=sstM(i,j)
#         ifdef SSTA_FILTER_OPT1
              sstA(i,j)=sstP(i,j)
#         endif      
            enddo
          enddo
        enddo
#elif defined HCM_LOESS_FILTER
        nx=IUB-ILB;
        ny=JUB-JLB;
        span_x=hcm_filt_lenx(ng)
        span_y=hcm_filt_leny(ng)
        allocate(xin(ILB:IUB))
        allocate(yin(JLB:JUB))
        do i=ILB,IUB
          Vdd=0.0_r8
          Cdd=0.0_r8
          do j=JLB,JUB
            Vdd=Vdd+lonM(i,j)
            Cdd=Cdd+1.0_r8
          enddo
          xin(i)=Vdd/Cdd
        enddo
        do j=JLB,JUB
          Vdd=0.0_r8
          Cdd=0.0_r8
          do i=ILB,IUB
            Vdd=Vdd+latM(i,j)
            Cdd=Cdd+1.0_r8
          enddo
          yin(j)=Vdd/Cdd
        enddo
        do loopN=1,LN
          call loess2(sstL,xin,yin,nx,ny,span_x,span_y,sstM)
#       ifdef SSTA_FILTER_OPT1
          call loess2(sstA,xin,yin,nx,ny,span_x,span_y,sstP)
#       endif      
          do i=ILB+BN,IUB-BN
            do j=JLB+BN,JUB-BN
              sstL(i,j)=sstM(i,j)
#       ifdef SSTA_FILTER_OPT1
              sstA(i,j)=sstP(i,j)
#       endif      
            enddo
          enddo
        enddo
#elif defined HCM_RUNNING_FILTER
        RX = floor(hcm_filt_lenx(ng)/2.0_r8)
        RY = floor(hcm_filt_leny(ng)/2.0_r8)
        do loopN=1,LN
          do i=ILB+BN,IUB-BN
          do j=JLB+BN,JUB-BN
            sstM(i,j)=0.0_r8
            Vdd=0.0_r8
            Cdd=0.0_r8
            Ndd=0.0_r8
#         ifdef SSTA_FILTER_OPT1
            sstP(i,j)=0.0_r8
            VddA=0.0_r8
            CddA=0.0_r8
#         endif      
            do ir=i-RX,i+RX
            do jr=j-RY,j+RY
              Ndd=Ndd+1.0_r8
              if(ir<ILB.or.ir>IUB.or.jr<JLB.or.jr>JUB)then
                Vdd=Vdd+0.0_r8
                Cdd=Cdd+0.0_r8
#         ifdef SSTA_FILTER_OPT1
                VddA=VddA+0.0_r8
                CddA=CddA+0.0_r8
#         endif      
              else
                if(sstL(ir,jr) .ne. -999)then
                  Vdd=Vdd+sstL(ir,jr)*mask(ir,jr)
                  Cdd=Cdd+mask(ir,jr)
#         ifdef SSTA_FILTER_OPT1
                  VddA=VddA+sstA(ir,jr)*mask(ir,jr)
                  CddA=CddA+mask(ir,jr)
#         endif      
                else
                  Vdd=Vdd+0.0_r8
                  Cdd=Cdd+0.0_r8
#         ifdef SSTA_FILTER_OPT1
                  VddA=VddA+0.0_r8
                  CddA=CddA+0.0_r8
#         endif      
                endif
              endif
            enddo
            enddo
            if(Cdd.gt.0.0_r8) then
              sstM(i,j)=Vdd/Cdd
#         ifdef SSTA_FILTER_OPT1
              sstP(i,j)=VddA/CddA
#         endif      
            endif
            if(Cdd*5.0 .lt. Ndd*4.0)then
              sstM(i,j)=-999
            endif
          enddo
          enddo
          do i=ILB,IUB
            do j=JLB,JUB
              sstL(i,j)=sstM(i,j)
#         ifdef SSTA_FILTER_OPT1
              sstA(i,j)=sstP(i,j)
#         endif      
            enddo
          enddo
        enddo
#endif


        do i=ILB,IUB
          do j=JLB,JUB
            sstH(i,j)=sstF(i,j)-sstL(i,j)
          enddo
        enddo


        do i=ILB,IUB
          do j=JLB,JUB
            if(sstL(i,j) .eq. -999) then
              sstL(i,j)=sstF(i,j)
              sstH(i,j)=0.0_r8
#ifdef SSTA_FILTER_OPT1
              sstA(i,j)=0.0_r8
#endif
            endif
          enddo
        enddo


#ifdef SSTA_FILTER_OPT2
        do i=ILB,IUB
          do j=JLB,JUB
            sstA(i,j)=sstL(i,j)-sstC(i,j)
          enddo
        enddo
#endif


#ifdef HCM_INITIAL_KICK
        IF (tdays(ng).lt.hcm_kick_time(ng)) THEN
          do i=ILB,IUB
            do j=JLB,JUB
              cffX=lonM(i,j)-hcm_kick_locx(ng)
              cffY=latM(i,j)-hcm_kick_locy(ng)
              cffX=cffX/hcm_kick_lenx(ng)
              cffY=cffY/hcm_kick_leny(ng)
              cffL=EXP(-cffX*cffX)*EXP(-cffY*cffY)
              sstA(i,j)=cffL*hcm_kick_temp(ng)
            enddo
          enddo
        ENDIF
#endif


      endif
!=======================================================================
!


! Scatter SST to different ranks
!=======================================================================
      if (MyRank==MyMaster) then
        LBi=BOUNDS(ng) % LBi(MyRank)
        UBi=BOUNDS(ng) % UBi(MyRank)
        LBj=BOUNDS(ng) % LBj(MyRank)
        UBj=BOUNDS(ng) % UBj(MyRank)
        FORCES(ng)%sstF(LBi:UBi,LBj:UBj)=sstF(LBi:UBi,LBj:UBj)
        FORCES(ng)%sstH(LBi:UBi,LBj:UBj)=sstH(LBi:UBi,LBj:UBj)
        FORCES(ng)%sstL(LBi:UBi,LBj:UBj)=sstL(LBi:UBi,LBj:UBj)
#    ifndef HCM_CLIM_SPIN
        FORCES(ng)%sstA(LBi:UBi,LBj:UBj)=sstA(LBi:UBi,LBj:UBj)
#    endif
        do tile=1,Ntiles
          LBi=BOUNDS(ng) % LBi(tile)
          UBi=BOUNDS(ng) % UBi(tile)
          LBj=BOUNDS(ng) % LBj(tile)
          UBj=BOUNDS(ng) % UBj(tile)
          dim1=UBi-LBi+1
          dim2=UBj-LBj+1
          len_buffer=dim1*dim2
          allocate(send_buffer(len_buffer))
          send_buffer=reshape(sstF(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagF,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)

          allocate(send_buffer(len_buffer))
          send_buffer=reshape(sstH(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagH,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)

          allocate(send_buffer(len_buffer))
          send_buffer=reshape(sstL(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagL,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)

#    ifndef HCM_CLIM_SPIN
          allocate(send_buffer(len_buffer))
          send_buffer=reshape(sstA(LBi:UBi,LBj:UBj),(/len_buffer/))
          call mpi_isend(send_buffer,len_buffer,MPI_REAL8,tile,tagA,MPI_COMM_WORLD,request,ierr)
          call mpi_wait(request,rstatus,ierr)
          deallocate(send_buffer)
#    endif
        enddo
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
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagF,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%sstF(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)

        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagH,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%sstH(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)

        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagL,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%sstL(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)

#    ifndef HCM_CLIM_SPIN
        allocate(recv_buffer(len_buffer))
        recv_buffer(len_buffer)=0.0_r8
        call mpi_irecv(recv_buffer,len_buffer,MPI_REAL8,MyMaster,tagA,MPI_COMM_WORLD,request,ierr)
        call mpi_wait(request,rstatus,ierr)
        FORCES(ng)%sstA(LBi:UBi,LBj:UBj)=reshape(recv_buffer,(/dim1,dim2/))
        deallocate(recv_buffer)
#    endif
      endif

! scatter finished
      call  mpi_barrier(MPI_COMM_WORLD, ierr)
!=======================================================================

      END SUBROUTINE hcm_filter
!***********************************************************************
#  endif
#endif


#ifdef HCM_LOESS_FILTER
!=======================================================================
! LOESS FILTER BLOCK 
!=======================================================================
       subroutine loess2(z,xin,yin,nx,ny,span_x,span_y,zout)
!-----------------------------------------------------------------------
!
! Usage: CALL LOESS2(Z,XIN,YIN,NX,NY,SPAN_X,SPAN_Y,ZOUT)
!
! Z is the input data array to be smoothed with NX rows and NY columns
! XIN is the input x-grid vector with NX elements
! YIN is the input y-grid vector with NY elements
! NX is an input integer with the number of x grid points
! NY is an input integer with the number of y grid points
! SPAN_X is the smoothing half-span in the x-direction
! SPAN_Y is the smoothing half-span in the y-direction
! ZOUT is the output double-precision data array the same size as Z
!
! Fortran 90 version
!
! Written by Larry O'Neill, based on an earlier version by
! Michael Schlax.
!
!-----------------------------------------------------------------------
!     Input/output variables
      integer, intent(in) :: nx, ny
      real*8, intent(in) :: xin(nx), yin(ny)
      real*8, dimension(nx,ny), intent(in) :: z
      real*8, intent(in) :: span_x, span_y
      real*8, dimension(nx,ny), intent(out) :: zout
!     Local variables
      integer i, j, ii, jj, ilo, ihi, jlo, jhi
      integer, parameter :: maxreg=20000
      real*8, dimension(maxreg,2) :: des
      real*8, parameter  :: amiss=-1.e34
      real*8, dimension(maxreg)   :: zsel, w
      real*8, dimension(nx,ny)    :: zhat
      real*8 :: dist, tem
      real*8 :: dx(nx)
      real*8 :: dy(ny)
      integer nreg
!-----------------------------------------------------------------------
      zhat = amiss

        do j=1,ny
         dy = yin - yin(j)
         jlo = max(1,minval(maxloc(dy,MASK=dy<-span_y)))
         jhi =       maxval(minloc(dy,MASK=dy> span_y))
         if(jhi==1) jhi=ny
         do i=1,nx
            if (z(i,j) .ne. amiss) then
            dx = xin - xin(i)
               ilo = max(1,minval(maxloc(dx,MASK=dx<-span_x)))
               ihi =       maxval(minloc(dx,MASK=dx> span_x))
               if(ihi==1) ihi=nx
               nreg = 0
               do 50 jj=jlo,jhi
                  do 40 ii=ilo,ihi
                     if(z(ii,jj) .ne. amiss)then
                        dist = dsqrt((dx(ii)/span_x)**2+                &
     &                               (dy(jj)/span_y)**2)
                        if(dist .le. 1.D0)then
                           nreg=nreg+1
                           if(nreg .gt. maxreg) then
                              write(*,*)'Dim error in loess2',i,j,nreg
                              return
                           endif
                           des(nreg,1)=dx(ii)
                           des(nreg,2)=dy(jj)
                           zsel(nreg) =z(ii,jj)
                           tem        =(1.D0-dist**3)
                           w(nreg)    =tem*tem*tem
                        endif
                     endif
   40             continue
   50          continue

               if(nreg .lt. 3)then
                  zhat(i,j)=amiss
               else
                  call regsm2(nreg,zsel,des,w,zhat(i,j),amiss,maxreg)
               endif
            end if
         end do
      end do
      zout = zhat
      return
      end
!=======================================================================
      subroutine regsm2(n,z,des,w,zhat,amiss,lddes)
!-----------------------------------------------------------------------
      integer, parameter :: np_iocas=6
      integer lddes
      real*8 z(lddes),des(lddes,2),w(lddes)
      real*8 a(np_iocas),xtx(np_iocas,np_iocas),tem(np_iocas)
      real*8 d, zmax, zmin, zi, zhat, amiss
      integer indx(np_iocas)
      integer i, j, k
      integer n, ier
!
! make the x-matrix
!
      do 20 j=1,np_iocas
         a(j)=0.D0
         do 10 i=1,np_iocas
            xtx(i,j)=0.D0
   10    continue
   20 continue
      zmax=-1.D35
      zmin=1.D35
      do 100 i=1,n
         if(z(i) .gt. zmax)zmax=z(i)
         if(z(i) .lt. zmin)zmin=z(i)
         zi=z(i)*w(i)
         tem(1)=w(i)
         tem(2)=des(i,1)*w(i)
         tem(3)=des(i,1)*des(i,1)*w(i)
         tem(4)=des(i,2)*w(i)
         tem(5)=des(i,2)*des(i,2)*w(i)
         tem(6)=des(i,1)*des(i,2)*w(i)
         do 40 k=1,np_iocas
            a(k)=a(k)+tem(k)*zi
            do 30 j=1,np_iocas
               xtx(j,k)=xtx(j,k)+tem(j)*tem(k)
   30       continue
   40    continue
  100 continue
!
! now solve it
!
      call ludcmp(xtx,np_iocas,np_iocas,indx,d,ier)
      if(ier .ne. 0)then
         zhat=amiss
         return
      endif
      call lubksb(xtx,np_iocas,np_iocas,indx,a)
      zhat=a(1)
!      ibad=0
      if(zhat.ge.zmax .or. zhat.le.zmin)then
!         write(*,*)'loess1 bad ',zmin,zhat,zmax
!         zhat=amiss
!         ibad=1
      endif
      return
      end
!=======================================================================
      SUBROUTINE LUDCMP(A,N,NP_IOCAS,INDX,D,ier)
!-----------------------------------------------------------------------
!      IMPLICIT REAL*8(A-H,O-Z)
      integer, PARAMETER :: NMAX=100
      integer N, NP_IOCAS, ier, I, J, K, IMAX
      real*8 AAMAX, D, SUM,DUM
      INTEGER INDX(N)
      REAL*8 A(NP_IOCAS,NP_IOCAS),VV(NMAX)
      REAL*8, PARAMETER :: TINY=1.D-20
      ier=0
      D=1.D0
      DO 12 I=1,N
        AAMAX=0.D0
        DO 11 J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.D0)then
            ier=1
!            write(*,*)'Singular matrix.'
            return
        endif
        VV(I)=1.D0/AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.D0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J) = IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.D0)A(J,J)=TINY
          DUM=1.D0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.D0)A(N,N)=TINY
      RETURN
      END
!=======================================================================
      SUBROUTINE LUBKSB(A,N,NP_IOCAS,INDX,B)
!-----------------------------------------------------------------------
!      IMPLICIT REAL*8(A-H,O-Z)
      integer N, NP_IOCAS, I, II, LL, J
      real*8 SUM
      real*8 A(NP_IOCAS,NP_IOCAS),B(N)
      integer INDX(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.D0) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
!=======================================================================
! END LOESS FILTER BLOCK 
!=======================================================================
#endif


      END MODULE hcm_filt_mod
