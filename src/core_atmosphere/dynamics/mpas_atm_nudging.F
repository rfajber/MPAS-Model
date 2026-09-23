! Copyright (c) 2026,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

#ifdef MPAS_OPENACC
#define MPAS_ACC_TIMER_START(X) call mpas_timer_start(X)
#define MPAS_ACC_TIMER_STOP(X) call mpas_timer_stop(X)
#else
#define MPAS_ACC_TIMER_START(X)
#define MPAS_ACC_TIMER_STOP(X)
#endif

module mpas_atm_nudging

 use mpas_derived_types
 use mpas_pool_routines
 use mpas_kind_types
 use mpas_constants, only : rvord
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop
 use mpas_atm_boundaries, only : mpas_atm_get_bdy_state
 use mpas_log, only : mpas_log_write

 private

 public :: atm_add_tend_nudging, &
           pre_add_tend_nudging, &
           post_add_tend_nudging

 contains


!==================================================================================================
 pure real (kind=RKIND) function nudging_taper(zagl, zbot, ztop) result(taper)
!==================================================================================================
! Fraction of the full nudging strength to apply at height zagl above ground level. Ramps linearly
! from zero at zbot to one at ztop; a non-positive depth (ztop <= zbot) disables the taper.

      implicit none

      real (kind=RKIND), intent(in) :: zagl, zbot, ztop

      if (ztop <= zbot) then
         taper = 1.0_RKIND
      else if (zagl <= zbot) then
         taper = 0.0_RKIND
      else if (zagl >= ztop) then
         taper = 1.0_RKIND
      else
         taper = (zagl - zbot) / (ztop - zbot)
      end if

 end function nudging_taper


!==================================================================================================
 pure subroutine nudging_filter_update(f, d_in, a)
!==================================================================================================
! Advance one column of the low-pass filter cascade by a single timestep and leave the filtered
! value in the last stage. Each stage is a first-order section, f <- f + a*(input - f), with the
! output of one stage feeding the next; cascading sharpens the rolloff and flattens the passband
! relative to a single section.
!
! Note that the cascade is an IIR filter: its entire memory is the stage values themselves. A true
! boxcar running mean would instead require the whole window of past samples to be retained, since
! the sample leaving the window changes every step.

      implicit none

      real (kind=RKIND), dimension(:), intent(inout) :: f
      real (kind=RKIND), intent(in) :: d_in
      real (kind=RKIND), intent(in) :: a

      integer :: s
      real (kind=RKIND) :: x

      x = d_in
      do s = 1, size(f)
         f(s) = f(s) + a * (x - f(s))
         x = f(s)
      end do

 end subroutine nudging_filter_update


!==================================================================================================
 subroutine pre_add_tend_nudging(configs, structs)
!==================================================================================================

      implicit none

      type (mpas_pool_type), intent(in) :: configs
      type (mpas_pool_type), intent(inout) :: structs

#ifdef MPAS_OPENACC
      type (mpas_pool_type), pointer :: tend
      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: diag
      type (mpas_pool_type), pointer :: nudging

      real (kind=RKIND), pointer :: config_nudging_filter_tau
      real (kind=RKIND), pointer :: config_nudging_tau_nwfa
      real (kind=RKIND), pointer :: config_nudging_tau_nifa
      real (kind=RKIND), dimension(:,:), pointer :: rho_edge, rho_zz, theta_m, u
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, tend_scalars
      real (kind=RKIND), dimension(:,:,:), pointer :: filt_u, filt_theta, filt_qv, filt_rho
      real (kind=RKIND), dimension(:,:,:), pointer :: filt_nwfa, filt_nifa

      call mpas_pool_get_subpool(structs, 'tend', tend)
      call mpas_pool_get_subpool(structs, 'state', state)
      call mpas_pool_get_subpool(structs, 'diag', diag)

      MPAS_ACC_TIMER_START('atm_srk3: nudging ACC_data_xfer')
      call mpas_pool_get_array(state, 'u', u, 1)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)
      !$acc update self(u, theta_m, scalars, rho_zz, rho_edge)

      call mpas_pool_get_array(tend, 'scalars_tend', tend_scalars)
      !$acc update self(tend_scalars)

      call mpas_pool_get_config(configs, 'config_nudging_filter_tau', config_nudging_filter_tau)
      if (config_nudging_filter_tau > 0.0_RKIND) then
         call mpas_pool_get_subpool(structs, 'nudging', nudging)
         call mpas_pool_get_array(nudging, 'nudging_filter_u', filt_u)
         call mpas_pool_get_array(nudging, 'nudging_filter_theta', filt_theta)
         call mpas_pool_get_array(nudging, 'nudging_filter_qv', filt_qv)
         call mpas_pool_get_array(nudging, 'nudging_filter_rho_zz', filt_rho)
         !$acc update self(filt_u, filt_theta, filt_qv, filt_rho)

         call mpas_pool_get_config(configs, 'config_nudging_tau_nwfa', config_nudging_tau_nwfa)
         if (config_nudging_tau_nwfa > 0.0_RKIND) then
            call mpas_pool_get_array(nudging, 'nudging_filter_nwfa', filt_nwfa)
            !$acc update self(filt_nwfa)
         end if

         call mpas_pool_get_config(configs, 'config_nudging_tau_nifa', config_nudging_tau_nifa)
         if (config_nudging_tau_nifa > 0.0_RKIND) then
            call mpas_pool_get_array(nudging, 'nudging_filter_nifa', filt_nifa)
            !$acc update self(filt_nifa)
         end if
      end if
      MPAS_ACC_TIMER_STOP('atm_srk3: nudging ACC_data_xfer')
#endif

 end subroutine pre_add_tend_nudging


!==================================================================================================
 subroutine post_add_tend_nudging(configs, structs, tend_ru, tend_rtheta, tend_rho)
!==================================================================================================

      implicit none

      type (mpas_pool_type), intent(in) :: configs
      type (mpas_pool_type), intent(inout) :: structs
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_ru
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_rtheta
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_rho

#ifdef MPAS_OPENACC
      type (mpas_pool_type), pointer :: tend
      type (mpas_pool_type), pointer :: nudging

      real (kind=RKIND), pointer :: config_nudging_filter_tau
      real (kind=RKIND), pointer :: config_nudging_tau_nwfa
      real (kind=RKIND), pointer :: config_nudging_tau_nifa
      real (kind=RKIND), dimension(:,:,:), pointer :: tend_scalars
      real (kind=RKIND), dimension(:,:,:), pointer :: filt_u, filt_theta, filt_qv, filt_rho
      real (kind=RKIND), dimension(:,:,:), pointer :: filt_nwfa, filt_nifa

      call mpas_pool_get_subpool(structs, 'tend', tend)

      MPAS_ACC_TIMER_START('atm_srk3: nudging ACC_data_xfer')
      !$acc update device(tend_ru, tend_rtheta, tend_rho)

      call mpas_pool_get_array(tend, 'scalars_tend', tend_scalars)
      !$acc update device(tend_scalars)

      call mpas_pool_get_config(configs, 'config_nudging_filter_tau', config_nudging_filter_tau)
      if (config_nudging_filter_tau > 0.0_RKIND) then
         call mpas_pool_get_subpool(structs, 'nudging', nudging)
         call mpas_pool_get_array(nudging, 'nudging_filter_u', filt_u)
         call mpas_pool_get_array(nudging, 'nudging_filter_theta', filt_theta)
         call mpas_pool_get_array(nudging, 'nudging_filter_qv', filt_qv)
         call mpas_pool_get_array(nudging, 'nudging_filter_rho_zz', filt_rho)
         !$acc update device(filt_u, filt_theta, filt_qv, filt_rho)

         call mpas_pool_get_config(configs, 'config_nudging_tau_nwfa', config_nudging_tau_nwfa)
         if (config_nudging_tau_nwfa > 0.0_RKIND) then
            call mpas_pool_get_array(nudging, 'nudging_filter_nwfa', filt_nwfa)
            !$acc update device(filt_nwfa)
         end if

         call mpas_pool_get_config(configs, 'config_nudging_tau_nifa', config_nudging_tau_nifa)
         if (config_nudging_tau_nifa > 0.0_RKIND) then
            call mpas_pool_get_array(nudging, 'nudging_filter_nifa', filt_nifa)
            !$acc update device(filt_nifa)
         end if
      end if
      MPAS_ACC_TIMER_STOP('atm_srk3: nudging ACC_data_xfer')
#endif

 end subroutine post_add_tend_nudging


!==================================================================================================
 subroutine atm_add_tend_nudging(clock, block, itimestep, dt, tend_ru, tend_rtheta, tend_rho)
!==================================================================================================
! Relax the model state toward the driving state held in the 'lbc' pool, which is read from the
! 'lbc_in' stream and time-interpolated to the current model time.
!
! The relaxation is formulated on the uncoupled fields so that each variable has an independent,
! physically meaningful timescale:
!
!     d(u)/dt      = (1/tau_u)     * (u_drv - u)
!     d(theta)/dt  = (1/tau_theta) * (theta_drv - theta)
!     d(qv)/dt     = (1/tau_qv)    * (qv_drv - qv)
!     d(rho_zz)/dt = (1/tau_rho)   * (rho_zz_drv - rho_zz)
!     d(nwfa)/dt   = (1/tau_nwfa)  * (nwfa_drv - nwfa)
!     d(nifa)/dt   = (1/tau_nifa)  * (nifa_drv - nifa)
!
! The two aerosol number concentrations are those of the aerosol-aware Thompson microphysics, and
! are nudged only when that scheme is active. They take no part in the moist coupling of theta_m;
! like any other scalar they are coupled to the mass field alone.
!
! and the resulting increments are then mass-coupled for the dynamics, following the same chain
! rule used by the IAU scheme in mpas_atm_iau.F. Because surface pressure is not prognostic in
! MPAS, it is nudged indirectly through the dry density rho_zz.
!
! When config_nudging_filter_tau is positive, the departure from the driving state is passed
! through a temporal low-pass filter before being relaxed, so that only the slowly evolving part
! of the departure is forced. Because atmospheric scales satisfy tau ~ L/U, damping the fast part
! of the departure leaves the model's own small-scale variability largely free while still holding
! the large-scale flow to the driving data.
!
! The filter is applied to the departure rather than to the model and driving states separately.
! For a linear filter LP[x_model] - LP[x_drv] == LP[x_model - x_drv], so the two are identical,
! but filtering the departure needs one set of filter states instead of two. It also means that
! any phase lag the filter introduces is common to both terms and cancels exactly, so the result
! does not depend on the filter's phase response.

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      type (block_type), intent(inout) :: block
      integer, intent(in) :: itimestep
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_ru
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_rtheta
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_rho

      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: diag
      type (mpas_pool_type), pointer :: mesh
      type (mpas_pool_type), pointer :: tend
      type (mpas_pool_type), pointer :: nudging

      real (kind=RKIND), pointer :: config_nudging_tau_u
      real (kind=RKIND), pointer :: config_nudging_tau_theta
      real (kind=RKIND), pointer :: config_nudging_tau_qv
      real (kind=RKIND), pointer :: config_nudging_tau_rho
      real (kind=RKIND), pointer :: config_nudging_tau_nwfa
      real (kind=RKIND), pointer :: config_nudging_tau_nifa
      real (kind=RKIND), pointer :: config_nudging_zbot
      real (kind=RKIND), pointer :: config_nudging_ztop
      real (kind=RKIND), pointer :: config_nudging_filter_tau
      logical, pointer :: config_do_restart

      integer, pointer :: nCells_ptr, nEdges_ptr, nCellsSolve_ptr, nEdgesSolve_ptr
      integer, pointer :: nVertLevels_ptr, index_qv_ptr, index_nwfa_ptr, index_nifa_ptr
      integer :: nCells, nEdges, nCellsSolve, nEdgesSolve, nVertLevels, index_qv
      integer :: index_nwfa, index_nifa

      real (kind=RKIND), dimension(:,:), pointer :: u, theta_m, rho_zz, rho_edge, zgrid
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, tend_scalars
      real (kind=RKIND), dimension(:,:,:), pointer :: filt_u, filt_theta, filt_qv, filt_rho
      real (kind=RKIND), dimension(:,:,:), pointer :: filt_nwfa, filt_nifa
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:,:), allocatable :: u_drv, theta_drv, rho_zz_drv, qv_drv
      real (kind=RKIND), dimension(:,:), allocatable :: nwfa_drv, nifa_drv

      real (kind=RKIND) :: alpha_u, alpha_theta, alpha_qv, alpha_rho
      real (kind=RKIND) :: alpha_nwfa, alpha_nifa
      real (kind=RKIND) :: zbot, ztop, zagl, taper
      real (kind=RKIND) :: qv, theta, d_rho_zz, d_theta, d_qv, d_rqv, d_rtheta, d_u
      real (kind=RKIND) :: aer, d_aer
      real (kind=RKIND) :: filter_coef
      logical :: filter_active, filter_init
      integer :: nStages
      integer :: iCell, iEdge, k, cell1, cell2

      call mpas_pool_get_config(block % configs, 'config_nudging_tau_u', config_nudging_tau_u)
      call mpas_pool_get_config(block % configs, 'config_nudging_tau_theta', config_nudging_tau_theta)
      call mpas_pool_get_config(block % configs, 'config_nudging_tau_qv', config_nudging_tau_qv)
      call mpas_pool_get_config(block % configs, 'config_nudging_tau_rho', config_nudging_tau_rho)
      call mpas_pool_get_config(block % configs, 'config_nudging_tau_nwfa', config_nudging_tau_nwfa)
      call mpas_pool_get_config(block % configs, 'config_nudging_tau_nifa', config_nudging_tau_nifa)
      call mpas_pool_get_config(block % configs, 'config_nudging_zbot', config_nudging_zbot)
      call mpas_pool_get_config(block % configs, 'config_nudging_ztop', config_nudging_ztop)
      call mpas_pool_get_config(block % configs, 'config_nudging_filter_tau', config_nudging_filter_tau)
      call mpas_pool_get_config(block % configs, 'config_do_restart', config_do_restart)

      !
      ! A non-positive relaxation timescale disables nudging for that variable
      !
      alpha_u     = 0.0_RKIND
      alpha_theta = 0.0_RKIND
      alpha_qv    = 0.0_RKIND
      alpha_rho   = 0.0_RKIND
      alpha_nwfa  = 0.0_RKIND
      alpha_nifa  = 0.0_RKIND
      if (config_nudging_tau_u     > 0.0_RKIND) alpha_u     = 1.0_RKIND / config_nudging_tau_u
      if (config_nudging_tau_theta > 0.0_RKIND) alpha_theta = 1.0_RKIND / config_nudging_tau_theta
      if (config_nudging_tau_qv    > 0.0_RKIND) alpha_qv    = 1.0_RKIND / config_nudging_tau_qv
      if (config_nudging_tau_rho   > 0.0_RKIND) alpha_rho   = 1.0_RKIND / config_nudging_tau_rho
      if (config_nudging_tau_nwfa  > 0.0_RKIND) alpha_nwfa  = 1.0_RKIND / config_nudging_tau_nwfa
      if (config_nudging_tau_nifa  > 0.0_RKIND) alpha_nifa  = 1.0_RKIND / config_nudging_tau_nifa

      if (alpha_u <= 0.0_RKIND .and. alpha_theta <= 0.0_RKIND .and. &
          alpha_qv <= 0.0_RKIND .and. alpha_rho <= 0.0_RKIND .and. &
          alpha_nwfa <= 0.0_RKIND .and. alpha_nifa <= 0.0_RKIND) then
         return
      end if

      zbot = config_nudging_zbot
      ztop = config_nudging_ztop

      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'tend', tend)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells_ptr)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges_ptr)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve_ptr)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve_ptr)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels_ptr)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv_ptr)
      call mpas_pool_get_dimension(state, 'index_nwfa', index_nwfa_ptr)
      call mpas_pool_get_dimension(state, 'index_nifa', index_nifa_ptr)

      nCells = nCells_ptr
      nEdges = nEdges_ptr
      nCellsSolve = nCellsSolve_ptr
      nEdgesSolve = nEdgesSolve_ptr
      nVertLevels = nVertLevels_ptr
      index_qv = index_qv_ptr

      !
      ! nwfa and nifa exist only with the aerosol-aware Thompson microphysics; the Registry gives
      ! their constituent index as -1 otherwise. Asking to nudge a scalar the run does not carry
      ! is a configuration error, so say so rather than silently doing nothing.
      !
      index_nwfa = -1
      index_nifa = -1
      if (associated(index_nwfa_ptr)) index_nwfa = index_nwfa_ptr
      if (associated(index_nifa_ptr)) index_nifa = index_nifa_ptr

      if (alpha_nwfa > 0.0_RKIND .and. index_nwfa <= 0) then
         call mpas_log_write('config_nudging_tau_nwfa is positive, but the water-friendly aerosol', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('number nwfa is not a scalar in this run. Aerosol nudging requires the', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('aerosol-aware Thompson microphysics (config_microp_scheme = ''mp_thompson_aerosols'').', &
                             messageType=MPAS_LOG_CRIT)
      end if

      if (alpha_nifa > 0.0_RKIND .and. index_nifa <= 0) then
         call mpas_log_write('config_nudging_tau_nifa is positive, but the ice-friendly aerosol', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('number nifa is not a scalar in this run. Aerosol nudging requires the', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('aerosol-aware Thompson microphysics (config_microp_scheme = ''mp_thompson_aerosols'').', &
                             messageType=MPAS_LOG_CRIT)
      end if

      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(state, 'u', u, 1)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)
      call mpas_pool_get_array(tend, 'scalars_tend', tend_scalars)

      !
      ! Temporal low-pass filter of the departure, if enabled
      !
      filter_active = (config_nudging_filter_tau > 0.0_RKIND)
      filter_coef = 0.0_RKIND
      nStages = 0
      filter_init = .false.

      if (filter_active) then
         call mpas_pool_get_subpool(block % structs, 'nudging', nudging)
         call mpas_pool_get_array(nudging, 'nudging_filter_u', filt_u)
         call mpas_pool_get_array(nudging, 'nudging_filter_theta', filt_theta)
         call mpas_pool_get_array(nudging, 'nudging_filter_qv', filt_qv)
         call mpas_pool_get_array(nudging, 'nudging_filter_rho_zz', filt_rho)
         if (alpha_nwfa > 0.0_RKIND) call mpas_pool_get_array(nudging, 'nudging_filter_nwfa', filt_nwfa)
         if (alpha_nifa > 0.0_RKIND) call mpas_pool_get_array(nudging, 'nudging_filter_nifa', filt_nifa)

         nStages = size(filt_theta, 1)

         ! Splitting the requested timescale across the cascade keeps the centroid of the filter
         ! kernel at config_nudging_filter_tau as stages are added, so changing the number of
         ! stages sharpens the filter without also changing how far back it reaches.
         filter_coef = min(1.0_RKIND, dt * real(nStages, RKIND) / config_nudging_filter_tau)

         ! On a cold start the filter states are meaningless, so seed every stage with the current
         ! departure; the cascade update is then a no-op for this step and the filter starts in
         ! equilibrium instead of ramping up from zero. On a restart they come from the restart file.
         filter_init = (itimestep == 1 .and. .not. config_do_restart)
      end if

      !
      ! Driving state, linearly interpolated in time to the current model time
      !
      allocate(u_drv(nVertLevels,nEdges+1))
      allocate(theta_drv(nVertLevels,nCells+1))
      allocate(rho_zz_drv(nVertLevels,nCells+1))
      allocate(qv_drv(nVertLevels,nCells+1))

      call mpas_atm_get_bdy_state(clock, block, nVertLevels, nEdges, 'u', 0.0_RKIND, u_drv)
      call mpas_atm_get_bdy_state(clock, block, nVertLevels, nCells, 'theta', 0.0_RKIND, theta_drv)
      call mpas_atm_get_bdy_state(clock, block, nVertLevels, nCells, 'rho_zz', 0.0_RKIND, rho_zz_drv)
      call mpas_atm_get_bdy_state(clock, block, nVertLevels, nCells, 'qv', 0.0_RKIND, qv_drv)

      if (alpha_nwfa > 0.0_RKIND) then
         allocate(nwfa_drv(nVertLevels,nCells+1))
         call mpas_atm_get_bdy_state(clock, block, nVertLevels, nCells, 'nwfa', 0.0_RKIND, nwfa_drv)
      end if

      if (alpha_nifa > 0.0_RKIND) then
         allocate(nifa_drv(nVertLevels,nCells+1))
         call mpas_atm_get_bdy_state(clock, block, nVertLevels, nCells, 'nifa', 0.0_RKIND, nifa_drv)
      end if

      !
      ! Horizontal momentum
      !
      if (alpha_u > 0.0_RKIND) then
         do iEdge = 1, nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k = 1, nVertLevels
               zagl = 0.25_RKIND * (zgrid(k,cell1) + zgrid(k+1,cell1) + zgrid(k,cell2) + zgrid(k+1,cell2)) &
                      - 0.5_RKIND * (zgrid(1,cell1) + zgrid(1,cell2))
               taper = nudging_taper(zagl, zbot, ztop)

               d_u = u_drv(k,iEdge) - u(k,iEdge)
               if (filter_active) then
                  if (filter_init) filt_u(:,k,iEdge) = d_u
                  call nudging_filter_update(filt_u(:,k,iEdge), d_u, filter_coef)
                  d_u = filt_u(nStages,k,iEdge)
               end if

               tend_ru(k,iEdge) = tend_ru(k,iEdge) + taper * alpha_u * rho_edge(k,iEdge) * d_u
            end do
         end do
      end if

      !
      ! Mass, potential temperature, and water vapor. The coupled tendencies follow from
      ! differentiating rho_zz*qv and rho_zz*theta*(1 + rvord*qv).
      !
      ! Note that the last term below is rvord*rho_zz*theta*d(qv)/dt and not, as in the
      ! analogous conversion in mpas_atm_iau.F, rvord*theta*d(rho_zz*qv)/dt. The two differ by
      ! rvord*theta*qv*d(rho_zz)/dt, which vanishes only when the density is left alone; since
      ! the mass field is nudged here, the extra term would be a real error.
      !
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            zagl = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell)) - zgrid(1,iCell)
            taper = nudging_taper(zagl, zbot, ztop)

            qv = scalars(index_qv,k,iCell)
            theta = theta_m(k,iCell) / (1.0_RKIND + rvord * qv)

            d_rho_zz = rho_zz_drv(k,iCell) - rho_zz(k,iCell)
            d_theta = theta_drv(k,iCell) - theta
            d_qv = qv_drv(k,iCell) - qv

            if (filter_active) then
               if (filter_init) then
                  filt_rho(:,k,iCell) = d_rho_zz
                  filt_theta(:,k,iCell) = d_theta
                  filt_qv(:,k,iCell) = d_qv
               end if
               call nudging_filter_update(filt_rho(:,k,iCell), d_rho_zz, filter_coef)
               call nudging_filter_update(filt_theta(:,k,iCell), d_theta, filter_coef)
               call nudging_filter_update(filt_qv(:,k,iCell), d_qv, filter_coef)
               d_rho_zz = filt_rho(nStages,k,iCell)
               d_theta = filt_theta(nStages,k,iCell)
               d_qv = filt_qv(nStages,k,iCell)
            end if

            d_rho_zz = taper * alpha_rho * d_rho_zz
            d_theta = taper * alpha_theta * d_theta
            d_qv = taper * alpha_qv * d_qv

            d_rqv = d_qv * rho_zz(k,iCell) + qv * d_rho_zz
            d_rtheta = d_theta * rho_zz(k,iCell) + theta * d_rho_zz

            tend_rho(k,iCell) = tend_rho(k,iCell) + d_rho_zz
            tend_scalars(index_qv,k,iCell) = tend_scalars(index_qv,k,iCell) + d_rqv
            tend_rtheta(k,iCell) = tend_rtheta(k,iCell) &
                                   + (1.0_RKIND + rvord * qv) * d_rtheta &
                                   + rvord * rho_zz(k,iCell) * theta * d_qv

            !
            ! Aerosol number concentrations. These are ordinary scalars: they are coupled to the
            ! mass field exactly as qv is, but contribute nothing to theta_m, so there is no
            ! counterpart to the rtheta correction above. d_rho_zz has already been tapered and
            ! scaled by alpha_rho at this point, so the coupling term matches the one used for qv.
            !
            if (alpha_nwfa > 0.0_RKIND) then
               aer = scalars(index_nwfa,k,iCell)
               d_aer = nwfa_drv(k,iCell) - aer

               if (filter_active) then
                  if (filter_init) filt_nwfa(:,k,iCell) = d_aer
                  call nudging_filter_update(filt_nwfa(:,k,iCell), d_aer, filter_coef)
                  d_aer = filt_nwfa(nStages,k,iCell)
               end if

               d_aer = taper * alpha_nwfa * d_aer
               tend_scalars(index_nwfa,k,iCell) = tend_scalars(index_nwfa,k,iCell) &
                                                  + d_aer * rho_zz(k,iCell) + aer * d_rho_zz
            end if

            if (alpha_nifa > 0.0_RKIND) then
               aer = scalars(index_nifa,k,iCell)
               d_aer = nifa_drv(k,iCell) - aer

               if (filter_active) then
                  if (filter_init) filt_nifa(:,k,iCell) = d_aer
                  call nudging_filter_update(filt_nifa(:,k,iCell), d_aer, filter_coef)
                  d_aer = filt_nifa(nStages,k,iCell)
               end if

               d_aer = taper * alpha_nifa * d_aer
               tend_scalars(index_nifa,k,iCell) = tend_scalars(index_nifa,k,iCell) &
                                                  + d_aer * rho_zz(k,iCell) + aer * d_rho_zz
            end if
         end do
      end do

      deallocate(u_drv)
      deallocate(theta_drv)
      deallocate(rho_zz_drv)
      deallocate(qv_drv)
      if (allocated(nwfa_drv)) deallocate(nwfa_drv)
      if (allocated(nifa_drv)) deallocate(nifa_drv)

 end subroutine atm_add_tend_nudging
!==================================================================================================

end module mpas_atm_nudging
