! Copyright (c) 2026,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

module mpas_atm_held_suarez

 use mpas_derived_types
 use mpas_pool_routines
 use mpas_kind_types
 use mpas_constants, only : rvord, rgas, cp, p0
 use mpas_log, only : mpas_log_write

 private

 public :: atm_add_tend_held_suarez

 contains


!==================================================================================================
 pure real (kind=RKIND) function hs_theta_equilibrium(p, sinlat2, coslat2, &
                                                      delta_T_y, delta_theta_z, T_strat) result(theta_eq)
!==================================================================================================
! Equilibrium potential temperature of Held and Suarez (1994), expressed in theta rather than T.
!
! The published form is
!
!     T_eq = max{ T_strat, [315 - dT_y sin^2(lat) - d(theta)_z log(p/p0) cos^2(lat)] (p/p0)^kappa }
!
! Dividing through by the Exner function (p/p0)^kappa turns this into a relation for theta and
! cancels the leading factor, leaving only the stratospheric floor carrying a pressure dependence.
! Nothing is approximated here; this is the same field expressed in the variable MPAS prognoses.

      implicit none

      real (kind=RKIND), intent(in) :: p
      real (kind=RKIND), intent(in) :: sinlat2, coslat2
      real (kind=RKIND), intent(in) :: delta_T_y, delta_theta_z, T_strat

      real (kind=RKIND) :: kappa, exner, theta_troposphere, theta_floor

      kappa = rgas / cp
      exner = (p / p0)**kappa

      theta_troposphere = 315.0_RKIND - delta_T_y * sinlat2 &
                          - delta_theta_z * log(p / p0) * coslat2
      theta_floor = T_strat / exner

      theta_eq = max(theta_floor, theta_troposphere)

 end function hs_theta_equilibrium


!==================================================================================================
 subroutine atm_held_suarez_setup(block)
!==================================================================================================
! Precompute the equilibrium potential temperature and the two relaxation rate profiles.
!
! Held and Suarez define these against sigma = p/p_s, which in a model with a moving pressure
! coordinate would have to be re-evaluated every step. MPAS uses a fixed height coordinate, so
! instead they are evaluated once against the static reference hydrostatic profile held in
! pressure_base. The resulting fields are constant in time, which is exact for the benchmark's
! flat, topography-free sphere and turns the forcing into a stored multiply-add.
!
! Because everything here derives from static data, the fields are simply rebuilt at the start of
! every run rather than carried in the restart file.

      implicit none

      type (block_type), intent(inout) :: block

      type (mpas_pool_type), pointer :: mesh
      type (mpas_pool_type), pointer :: diag
      type (mpas_pool_type), pointer :: held_suarez

      real (kind=RKIND), pointer :: config_hs_delta_T_y
      real (kind=RKIND), pointer :: config_hs_delta_theta_z
      real (kind=RKIND), pointer :: config_hs_T_strat
      real (kind=RKIND), pointer :: config_hs_sigma_b
      real (kind=RKIND), pointer :: config_hs_tau_a
      real (kind=RKIND), pointer :: config_hs_tau_s
      real (kind=RKIND), pointer :: config_hs_tau_f

      integer, pointer :: nCells_ptr, nEdges_ptr, nVertLevels_ptr
      integer :: nCells, nEdges, nVertLevels

      real (kind=RKIND), dimension(:), pointer :: latCell
      real (kind=RKIND), dimension(:,:), pointer :: pressure_base
      real (kind=RKIND), dimension(:,:), pointer :: hs_theta_eq, hs_k_t, hs_k_v
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:,:), allocatable :: k_v_cell

      real (kind=RKIND) :: k_a, k_s, k_f, sigma_b
      real (kind=RKIND) :: sinlat2, coslat2, p, sigma, damp
      integer :: iCell, iEdge, k, cell1, cell2

#ifdef MPAS_OPENACC
      ! The forcing is computed on the host and the dynamics tendency arrays it writes into live
      ! on the device, so without data movement its contribution would be silently discarded.
      ! Fail here rather than produce a plausible but unforced simulation.
      call mpas_log_write('Held-Suarez forcing is not supported in OpenACC builds.', &
                          messageType=MPAS_LOG_CRIT)
#endif

      call mpas_pool_get_config(block % configs, 'config_hs_delta_T_y', config_hs_delta_T_y)
      call mpas_pool_get_config(block % configs, 'config_hs_delta_theta_z', config_hs_delta_theta_z)
      call mpas_pool_get_config(block % configs, 'config_hs_T_strat', config_hs_T_strat)
      call mpas_pool_get_config(block % configs, 'config_hs_sigma_b', config_hs_sigma_b)
      call mpas_pool_get_config(block % configs, 'config_hs_tau_a', config_hs_tau_a)
      call mpas_pool_get_config(block % configs, 'config_hs_tau_s', config_hs_tau_s)
      call mpas_pool_get_config(block % configs, 'config_hs_tau_f', config_hs_tau_f)

      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_subpool(block % structs, 'held_suarez', held_suarez)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells_ptr)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges_ptr)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels_ptr)
      nCells = nCells_ptr
      nEdges = nEdges_ptr
      nVertLevels = nVertLevels_ptr

      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_base)

      call mpas_pool_get_array(held_suarez, 'hs_theta_eq', hs_theta_eq)
      call mpas_pool_get_array(held_suarez, 'hs_k_t', hs_k_t)
      call mpas_pool_get_array(held_suarez, 'hs_k_v', hs_k_v)

      sigma_b = config_hs_sigma_b
      k_a = 1.0_RKIND / config_hs_tau_a
      k_s = 1.0_RKIND / config_hs_tau_s
      k_f = 1.0_RKIND / config_hs_tau_f

      allocate(k_v_cell(nVertLevels,nCells+1))
      k_v_cell(:,:) = 0.0_RKIND

      do iCell = 1, nCells
         sinlat2 = sin(latCell(iCell))**2
         coslat2 = cos(latCell(iCell))**2

         do k = 1, nVertLevels
            p = pressure_base(k,iCell)

            ! The benchmark has no topography, so the reference surface pressure is p0
            sigma = p / p0

            hs_theta_eq(k,iCell) = hs_theta_equilibrium(p, sinlat2, coslat2, &
                                                        config_hs_delta_T_y, &
                                                        config_hs_delta_theta_z, &
                                                        config_hs_T_strat)

            damp = max(0.0_RKIND, (sigma - sigma_b) / (1.0_RKIND - sigma_b))

            hs_k_t(k,iCell) = k_a + (k_s - k_a) * damp * coslat2**2
            k_v_cell(k,iCell) = k_f * damp
         end do
      end do

      ! The Rayleigh drag acts on the edge-normal velocity, so average the cell-based rate to edges
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1, nVertLevels
            hs_k_v(k,iEdge) = 0.5_RKIND * (k_v_cell(k,cell1) + k_v_cell(k,cell2))
         end do
      end do

      deallocate(k_v_cell)

      call mpas_log_write('')
      call mpas_log_write('Held-Suarez forcing is active')
      call mpas_log_write('  min/max equilibrium potential temperature = $r $r', &
                          realArgs=(/minval(hs_theta_eq(:,1:nCells)), maxval(hs_theta_eq(:,1:nCells))/))
      call mpas_log_write('  max thermal relaxation rate = $r s^-1', &
                          realArgs=(/maxval(hs_k_t(:,1:nCells))/))
      call mpas_log_write('  max Rayleigh damping rate   = $r s^-1', &
                          realArgs=(/maxval(hs_k_v(:,1:nEdges))/))
      call mpas_log_write('')

 end subroutine atm_held_suarez_setup


!==================================================================================================
 subroutine atm_add_tend_held_suarez(block, itimestep, tend_ru, tend_rtheta)
!==================================================================================================
! Apply the Held and Suarez (1994) idealized forcing: Newtonian relaxation of the temperature
! field toward a zonally symmetric equilibrium, and Rayleigh damping of the low-level winds.
!
!     d(theta)/dt = -k_t(lat,sigma) * (theta - theta_eq(lat,sigma))
!     d(u)/dt     = -k_v(sigma) * u
!
! This replaces the physics parameterizations rather than supplementing them, so it is intended to
! be run with config_physics_suite = 'none'. The benchmark is dry; no moisture or mass forcing is
! applied, and the coupling below reduces to a single term when the water vapor mixing ratio is
! zero.

      implicit none

      type (block_type), intent(inout) :: block
      integer, intent(in) :: itimestep
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_ru
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_rtheta

      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: diag
      type (mpas_pool_type), pointer :: mesh
      type (mpas_pool_type), pointer :: held_suarez

      integer, pointer :: nCellsSolve_ptr, nEdgesSolve_ptr, nVertLevels_ptr, index_qv_ptr
      integer :: nCellsSolve, nEdgesSolve, nVertLevels, index_qv

      real (kind=RKIND), dimension(:,:), pointer :: u, theta_m, rho_zz, rho_edge
      real (kind=RKIND), dimension(:,:), pointer :: hs_theta_eq, hs_k_t, hs_k_v
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      real (kind=RKIND) :: qv, theta, moist_coef
      integer :: iCell, iEdge, k

      !
      ! The equilibrium state and both rate profiles derive only from static data, so they are
      ! rebuilt once at the start of the run rather than recomputed every step or stored in the
      ! restart file.
      !
      if (itimestep == 1) then
         call atm_held_suarez_setup(block)
      end if

      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'held_suarez', held_suarez)

      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve_ptr)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve_ptr)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels_ptr)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv_ptr)
      nCellsSolve = nCellsSolve_ptr
      nEdgesSolve = nEdgesSolve_ptr
      nVertLevels = nVertLevels_ptr
      index_qv = index_qv_ptr

      call mpas_pool_get_array(state, 'u', u, 1)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)

      call mpas_pool_get_array(held_suarez, 'hs_theta_eq', hs_theta_eq)
      call mpas_pool_get_array(held_suarez, 'hs_k_t', hs_k_t)
      call mpas_pool_get_array(held_suarez, 'hs_k_v', hs_k_v)

      !
      ! Rayleigh damping of the low-level winds
      !
      do iEdge = 1, nEdgesSolve
         do k = 1, nVertLevels
            tend_ru(k,iEdge) = tend_ru(k,iEdge) &
                               - hs_k_v(k,iEdge) * rho_edge(k,iEdge) * u(k,iEdge)
         end do
      end do

      !
      ! Newtonian relaxation toward the equilibrium potential temperature. Only theta is forced,
      ! so differentiating rho_zz*theta*(1 + rvord*qv) leaves a single surviving term.
      !
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            qv = scalars(index_qv,k,iCell)
            moist_coef = 1.0_RKIND + rvord * qv
            theta = theta_m(k,iCell) / moist_coef

            tend_rtheta(k,iCell) = tend_rtheta(k,iCell) &
                                   + moist_coef * rho_zz(k,iCell) * hs_k_t(k,iCell) &
                                     * (hs_theta_eq(k,iCell) - theta)
         end do
      end do

 end subroutine atm_add_tend_held_suarez
!==================================================================================================

end module mpas_atm_held_suarez
