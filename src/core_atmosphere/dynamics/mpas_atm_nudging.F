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
 subroutine pre_add_tend_nudging(structs)
!==================================================================================================

      implicit none

      type (mpas_pool_type), intent(inout) :: structs

#ifdef MPAS_OPENACC
      type (mpas_pool_type), pointer :: tend
      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: diag

      real (kind=RKIND), dimension(:,:), pointer :: rho_edge, rho_zz, theta_m, u
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, tend_scalars

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
      MPAS_ACC_TIMER_STOP('atm_srk3: nudging ACC_data_xfer')
#endif

 end subroutine pre_add_tend_nudging


!==================================================================================================
 subroutine post_add_tend_nudging(structs, tend_ru, tend_rtheta, tend_rho)
!==================================================================================================

      implicit none

      type (mpas_pool_type), intent(inout) :: structs
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_ru
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_rtheta
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_rho

#ifdef MPAS_OPENACC
      type (mpas_pool_type), pointer :: tend

      real (kind=RKIND), dimension(:,:,:), pointer :: tend_scalars

      call mpas_pool_get_subpool(structs, 'tend', tend)

      MPAS_ACC_TIMER_START('atm_srk3: nudging ACC_data_xfer')
      !$acc update device(tend_ru, tend_rtheta, tend_rho)

      call mpas_pool_get_array(tend, 'scalars_tend', tend_scalars)
      !$acc update device(tend_scalars)
      MPAS_ACC_TIMER_STOP('atm_srk3: nudging ACC_data_xfer')
#endif

 end subroutine post_add_tend_nudging


!==================================================================================================
 subroutine atm_add_tend_nudging(clock, block, tend_ru, tend_rtheta, tend_rho)
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
!
! and the resulting increments are then mass-coupled for the dynamics, following the same chain
! rule used by the IAU scheme in mpas_atm_iau.F. Because surface pressure is not prognostic in
! MPAS, it is nudged indirectly through the dry density rho_zz.

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      type (block_type), intent(inout) :: block
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_ru
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_rtheta
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_rho

      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: diag
      type (mpas_pool_type), pointer :: mesh
      type (mpas_pool_type), pointer :: tend

      real (kind=RKIND), pointer :: config_nudging_tau_u
      real (kind=RKIND), pointer :: config_nudging_tau_theta
      real (kind=RKIND), pointer :: config_nudging_tau_qv
      real (kind=RKIND), pointer :: config_nudging_tau_rho
      real (kind=RKIND), pointer :: config_nudging_zbot
      real (kind=RKIND), pointer :: config_nudging_ztop

      integer, pointer :: nCells_ptr, nEdges_ptr, nCellsSolve_ptr, nEdgesSolve_ptr
      integer, pointer :: nVertLevels_ptr, index_qv_ptr
      integer :: nCells, nEdges, nCellsSolve, nEdgesSolve, nVertLevels, index_qv

      real (kind=RKIND), dimension(:,:), pointer :: u, theta_m, rho_zz, rho_edge, zgrid
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, tend_scalars
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:,:), allocatable :: u_drv, theta_drv, rho_zz_drv, qv_drv

      real (kind=RKIND) :: alpha_u, alpha_theta, alpha_qv, alpha_rho
      real (kind=RKIND) :: zbot, ztop, zagl, taper
      real (kind=RKIND) :: qv, theta, d_rho_zz, d_theta, d_qv, d_rqv, d_rtheta
      integer :: iCell, iEdge, k, cell1, cell2

      call mpas_pool_get_config(block % configs, 'config_nudging_tau_u', config_nudging_tau_u)
      call mpas_pool_get_config(block % configs, 'config_nudging_tau_theta', config_nudging_tau_theta)
      call mpas_pool_get_config(block % configs, 'config_nudging_tau_qv', config_nudging_tau_qv)
      call mpas_pool_get_config(block % configs, 'config_nudging_tau_rho', config_nudging_tau_rho)
      call mpas_pool_get_config(block % configs, 'config_nudging_zbot', config_nudging_zbot)
      call mpas_pool_get_config(block % configs, 'config_nudging_ztop', config_nudging_ztop)

      !
      ! A non-positive relaxation timescale disables nudging for that variable
      !
      alpha_u     = 0.0_RKIND
      alpha_theta = 0.0_RKIND
      alpha_qv    = 0.0_RKIND
      alpha_rho   = 0.0_RKIND
      if (config_nudging_tau_u     > 0.0_RKIND) alpha_u     = 1.0_RKIND / config_nudging_tau_u
      if (config_nudging_tau_theta > 0.0_RKIND) alpha_theta = 1.0_RKIND / config_nudging_tau_theta
      if (config_nudging_tau_qv    > 0.0_RKIND) alpha_qv    = 1.0_RKIND / config_nudging_tau_qv
      if (config_nudging_tau_rho   > 0.0_RKIND) alpha_rho   = 1.0_RKIND / config_nudging_tau_rho

      if (alpha_u <= 0.0_RKIND .and. alpha_theta <= 0.0_RKIND .and. &
          alpha_qv <= 0.0_RKIND .and. alpha_rho <= 0.0_RKIND) then
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

      nCells = nCells_ptr
      nEdges = nEdges_ptr
      nCellsSolve = nCellsSolve_ptr
      nEdgesSolve = nEdgesSolve_ptr
      nVertLevels = nVertLevels_ptr
      index_qv = index_qv_ptr

      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(state, 'u', u, 1)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)
      call mpas_pool_get_array(tend, 'scalars_tend', tend_scalars)

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

               tend_ru(k,iEdge) = tend_ru(k,iEdge) &
                                  + taper * alpha_u * rho_edge(k,iEdge) * (u_drv(k,iEdge) - u(k,iEdge))
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

            d_rho_zz = taper * alpha_rho * (rho_zz_drv(k,iCell) - rho_zz(k,iCell))
            d_theta = taper * alpha_theta * (theta_drv(k,iCell) - theta)
            d_qv = taper * alpha_qv * (qv_drv(k,iCell) - qv)

            d_rqv = d_qv * rho_zz(k,iCell) + qv * d_rho_zz
            d_rtheta = d_theta * rho_zz(k,iCell) + theta * d_rho_zz

            tend_rho(k,iCell) = tend_rho(k,iCell) + d_rho_zz
            tend_scalars(index_qv,k,iCell) = tend_scalars(index_qv,k,iCell) + d_rqv
            tend_rtheta(k,iCell) = tend_rtheta(k,iCell) &
                                   + (1.0_RKIND + rvord * qv) * d_rtheta &
                                   + rvord * rho_zz(k,iCell) * theta * d_qv
         end do
      end do

      deallocate(u_drv)
      deallocate(theta_drv)
      deallocate(rho_zz_drv)
      deallocate(qv_drv)

 end subroutine atm_add_tend_nudging
!==================================================================================================

end module mpas_atm_nudging
