!    Copyright (C) 2020 Imperial College London and others.
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Affero General Public License
!    as published by the Free Software Foundation,
!    version 3.0 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without seven the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module multi_pipes

    use fldebug

    use fields

    use reference_counting
    use memory_diagnostics

    use global_parameters

    use Petsc_Tools
    use Sparse_tools
    use sparse_tools_petsc

    use surface_integrals, only :integrate_over_surface_element
    use shape_functions_Linear_Quadratic
    use shape_functions_NDim
    use shape_functions_prototype
    use matrix_operations
    use Copy_Outof_State
    use boundary_conditions

    use multi_tools
    use multi_data_types

    use write_state_module, only: write_state

    implicit none

    private

    public  :: ASSEMBLE_PIPE_TRANSPORT_AND_CTY, MOD_1D_FORCE_BAL_C, retrieve_pipes_coords, pipe_coords, initialize_pipes_package_and_gamma

    !Parameters for boundary conditions copied from cv_adv-diff.F90
    INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1, WIC_T_BC_ROBIN = 2, &
        WIC_T_BC_DIRI_ADV_AND_ROBIN = 3, WIC_D_BC_DIRICHLET = 1, &
        WIC_U_BC_DIRICHLET = 1, &
        WIC_U_BC_ROBIN = 2, &
        WIC_U_BC_DIRI_ADV_AND_ROBIN = 3, &
        WIC_U_BC_DIRICHLET_INOUT = 2, &
        WIC_P_BC_DIRICHLET = 1, &
        WIC_P_BC_FREE = 2

    real:: tolerancePipe = 1d-2!> tolerancePipe has to be around 1e-2 because that is the precision of the nastran input file

contains

  !>@brief: This sub modifies either Mmat%CT or the Advection-diffusion equation for 1D pipe modelling
  SUBROUTINE MOD_1D_CT_AND_ADV( state, packed_state, final_phase, wells_first_phase, Mdims, ndgln, WIC_T_BC_ALL,WIC_D_BC_ALL, WIC_U_BC_ALL, SUF_T_BC_ALL,SUF_D_BC_ALL,SUF_U_BC_ALL, &
                  getcv_disc, getct, Mmat, Mspars, DT, MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE, mass_pipe, MASS_PIPE_FOR_COUP, &
                  INV_SIGMA, OPT_VEL_UPWIND_COEFS_NEW, eles_with_pipe, thermal, CV_BETA, bcs_outfluxes, outfluxes, assemble_collapsed_to_one_phase )
      type(state_type), intent(inout) :: packed_state
      type(state_type), dimension(:), intent(in) :: state
      type(multi_dimensions), intent(in) :: Mdims
      type(multi_ndgln), intent(in) :: ndgln
      type (multi_matrices), intent(inout) :: Mmat
      type (multi_sparsities), intent(in) :: Mspars
      integer, dimension(:,:,:), intent( in ) :: WIC_T_BC_ALL, WIC_D_BC_ALL, WIC_U_BC_ALL
      real, dimension(:,:,:), intent( in ) :: SUF_T_BC_ALL,SUF_D_BC_ALL,SUF_U_BC_ALL
      real, dimension(:,:,:,:), intent( in ) :: OPT_VEL_UPWIND_COEFS_NEW
      real, dimension(:,:),intent( inout ) :: INV_SIGMA
      real, dimension(:),intent( inout ) :: MASS_CVFEM2PIPE, MASS_PIPE2CVFEM, MASS_CVFEM2PIPE_TRUE ! of length NCMC
      real, dimension(:),intent( inout ) :: mass_pipe, MASS_PIPE_FOR_COUP ! of length Mdims%cv_nonods
      logical, intent( in ) :: getcv_disc, getct, thermal, assemble_collapsed_to_one_phase
      real, intent(in) :: DT, CV_BETA
      !Variables that are used to define the pipe pos.
      type(pipe_coords), dimension(:), intent(in):: eles_with_pipe
      integer, intent(in) :: final_phase, wells_first_phase
      ! Local variables
      real, dimension(:,:, :), allocatable, intent(inout):: bcs_outfluxes!<= if allocated then calculate outfluxes
      type (multi_outfluxes), intent(inout) :: outfluxes
      ! Local variables
      INTEGER :: CV_NODI, CV_NODJ, IPHASE, COUNT, CV_SILOC, SELE, cv_iloc, cv_jloc, jphase, assembly_phase
      INTEGER :: cv_ncorner, cv_lnloc, u_lnloc, i_indx, j_indx, ele, cv_gi, iloop, ICORNER, NPIPES, i
      integer, dimension(:), pointer :: cv_neigh_ptr
      integer, dimension(:), allocatable:: CV_GL_LOC, CV_GL_GL, X_GL_GL, MAT_GL_GL, u_GL_LOC, u_GL_GL
      integer, dimension(Mdims%ndim+1) :: CV_LOC_CORNER, U_LOC_CORNER!Allocate with number of corners; NCORNER
      integer, dimension(Mdims%ndim+1, Mdims%ndim+1) :: CV_MID_SIDE, U_MID_SIDE
      real, dimension(:), allocatable:: cvweigh, cv_nodpos, u_nodpos, lcv_b, &
          PIPE_DIAM_GI_vol, PIPE_DIAM_GI_suf, suf_detwei, vol_detwei, vol_detwei2
              !Variables to speed up memory allocations
      real, dimension( max(Mdims%nphase, 100) *6 ) :: memory_limiters!Get biggest between maximum nunber of fields (100), or phases
      real, dimension(Mdims%ndim) :: direction, direction_NORM
      real, dimension(final_phase*2) :: INV_SIGMA_GI, TUPWIND_OUT,  DUPWIND_OUT, TUPWIND_in,  DUPWIND_in, ndotq, income, income_j,&
          FEMTGI, FEMdGI, T_CV_NODI, T_CV_NODJ, D_CV_NODI, D_CV_NODJ, limt, limd, fvt, limdt, bczero, LOC_CT_RHS_U_ILOC, xi_limit
      real, dimension(Mdims%ndim, final_phase*2) :: ct_con, UGI_ALL
      !Variables for boundary conditions
      INTEGER, dimension(Mdims%cv_nonods):: WIC_B_BC_ALL_NODS
      INTEGER, dimension(final_phase*2,Mdims%cv_nonods):: WIC_T_BC_ALL_NODS, WIC_D_BC_ALL_NODS, WIC_U_BC_ALL_NODS
      real, dimension(final_phase*2,Mdims%cv_nonods):: SUF_T_BC_ALL_NODS, SUF_D_BC_ALL_NODS, RVEC_SUM_T, RVEC_SUM_D, RVEC_SUM_U
      real, dimension(Mdims%ndim, final_phase*2) :: MEAN_U
      real, dimension(Mdims%ndim,final_phase*2,Mdims%cv_nonods):: SUF_U_BC_ALL_NODS

      !Variables for shape functions
      real, dimension(:), allocatable :: CVN_VOL_ADJ
      real, dimension(:,:), allocatable:: cvn, n, nlx, un, unlx, sbcvfen, sbcvfenslx, sbufen, L_CVFENX_ALL, L_UFENX_ALL,L_UFEN_REVERSED,L_UFEN,&
          x_all_corn
      real, dimension(:,:), allocatable:: cvn_fem, cvn_femlx, cvnn_fem, cvnn_femlx
      real, dimension(:,:,:), allocatable:: L_CVFENX_ALL_REVERSED
      logical :: CV_QUADRATIC, U_QUADRATIC, ndiff, diff, ELE_HAS_PIPE, U_P0DG
      logical :: CALC_SIGMA_PIPE

      real, dimension(Mdims%ndim) :: T1, T2, TT1, TT2, NN1
      real, dimension(final_phase*2) :: LOC_CV_RHS_I, LOC_MAT_II, LOC_MAT_JJ, LOC_MAT_IJ
      real :: T1TT1, TT1TT1, T1TT2, T2TT1, T2TT2, DET_SQRT, INV_SIGMA_ND, N1NN1, INV_SIGMA_NANO_ND, MIN_INV_SIG
      real :: cv_ldx, u_ldx, dx, ele_angle, cv_m, sigma_gi, M_CVFEM2PIPE, M_PIPE2CVFEM, rnorm_sign, suf_area, PIPE_DIAM_END, MIN_DIAM
      real, dimension(1) :: R1, R2, RZ
      real, dimension(final_phase*2) :: TMAX, TMIN, DENMAX, DENMIN
      integer :: ierr, PIPE_NOD_COUNT, NPIPES_IN_ELE, ipipe, CV_LILOC, CV_LJLOC, U_LILOC, &
          u_iloc, x_iloc, cv_knod, idim, cv_lkloc, u_lkloc, u_knod, gi, ncorner, cv_lngi, u_lngi, cv_bngi, bgi, &
          icorner1, icorner2, icorner3, icorner4, JCV_NOD1, JCV_NOD2, CV_NOD, JCV_NOD, JU_NOD, &
          U_NOD, U_SILOC, COUNT2, MAT_KNOD, MAT_NODI, COUNT3, IPRES, k, iofluxes, compact_phase, global_phase
      real, dimension(final_phase*2, Mdims%cv_nonods):: tmax_all, tmin_all, denmax_all, denmin_all
      type(tensor_field), pointer :: t_all, den_all, u_all, aux_tensor_pointer, tfield, tfield2, t2_all, only_den_all
      type(scalar_field), pointer :: pipe_diameter, sigma1_pipes, sfield
      type(vector_field), pointer :: X
      !Logical to check if we using a conservative method or not, to save cpu time
      logical :: conservative_advection
      !Parameters of the simulation
      logical, parameter :: UPWIND_PIPES = .false.! Used for testing...
      logical, parameter :: PIPE_MIN_DIAM=.true. ! Use the pipe min diamter along a pipe element edge and min inv_sigma (max. drag reflcting min pipe diameter)
      logical, parameter :: LUMP_COUPLING_RES_PIPES = .true. ! Lump the coupling term which couples the pressure between the pipe and reservior.
      real, parameter :: INFINY=1.0E+20
      integer, parameter :: WIC_B_BC_DIRICHLET = 1
      !Variables extra for outfluxes
      type(tensor_field), pointer ::temp_field, salt_field
      logical :: compute_outfluxes


      conservative_advection = abs(cv_beta) > 0.99


      !If we are going to calculate the outfluxes (this is done when GETCT=.true.)
      compute_outfluxes = GETCT
      IF ( compute_outfluxes ) THEN
        !We may need to retrieve some extra fields
        if (has_temperature) then
            temp_field => extract_tensor_field( packed_state, "PackedTemperature" )
            if (outfluxes%calculate_flux)outfluxes%totout(2, :,:) = 0
        end if
        !Arash
        if (has_salt) then
            salt_field => extract_tensor_field( packed_state, "PackedSoluteMassFraction" )
            if (outfluxes%calculate_flux)outfluxes%totout(3, :,:) = 0
        end if


      end if


      CALC_SIGMA_PIPE = have_option("/porous_media/wells_and_pipes/well_options/calculate_sigma_pipe") ! Calculate sigma based on friction factors...
      NCORNER = Mdims%ndim + 1
      ! default limiting NVD diagram...
      XI_LIMIT = 2.0 ; ndiff=.false. ; diff=.true.
      CV_QUADRATIC = (Mdims%cv_nloc==6.AND.Mdims%ndim==2) .OR. (Mdims%cv_nloc==10.AND.Mdims%ndim==3)
      U_QUADRATIC = (Mdims%u_nloc==6.AND.Mdims%ndim==2) .OR. (Mdims%u_nloc==10.AND.Mdims%ndim==3)
      U_P0DG = Mdims%u_nloc == 1
      ! Calculate the local corner nodes...
      CALL CALC_CORNER_NODS( CV_LOC_CORNER, Mdims%ndim, Mdims%cv_nloc, CV_QUADRATIC, CV_MID_SIDE )
      CALL CALC_CORNER_NODS( U_LOC_CORNER, Mdims%ndim, Mdims%u_nloc, U_QUADRATIC, U_MID_SIDE )
      if ( CV_QUADRATIC ) then
          cv_lngi = 9; cv_lnloc = 3; cv_bngi = 2
      else
          cv_lngi = 6; cv_lnloc = 2; cv_bngi = 1
      end if
      if ( U_QUADRATIC ) then
          u_lnloc = 3
      else
          u_lnloc = 2
          !For P0DG only one node exists
          if (U_P0DG) u_lnloc = 1
      end if

      u_lngi = cv_lngi
      allocate( cvweigh(cv_lngi), cvn(cv_lnloc, cv_lngi), n(cv_lnloc, cv_lngi), &
          nlx(cv_lnloc, cv_lngi), un(u_lnloc, u_lngi), unlx(u_lnloc, u_lngi), &
          cv_nodpos(cv_lnloc), u_nodpos(u_lnloc), lcv_b(cv_lnloc), sbcvfen(cv_lnloc, cv_bngi), sbufen(u_lnloc, cv_bngi) )
      allocate(cvn_fem(cv_lnloc, cv_lngi), cvn_femlx(cv_lnloc, cv_lngi), cvnn_fem(cv_lnloc, cv_lngi), cvnn_femlx(cv_lnloc, cv_lngi))
      ! calculate shape functions...
      call quad_1d_shape( cv_lngi, cv_lnloc, u_lnloc, cvn, cvweigh, n, nlx, un, unlx )
      call quad_1d_shape( cv_lngi, cv_lnloc, cv_lnloc, cvn, cvweigh, cvn_fem, cvn_femlx, cvnn_fem, cvnn_femlx )
      ! calculate sbcvfen.
      ! determine the CV node positions...
      cv_nodpos(1) = -1.0
      cv_ldx = 2.0 / real(cv_lnloc-1)
      do i = 2, cv_lnloc
          cv_nodpos(i) = cv_nodpos(i-1) + cv_ldx
      end do
      ! for calculating sbufen.
      ! determine the U node positions...
      u_nodpos(1) = -1.0
      u_ldx = 2.0 / real(u_lnloc-1)
      do i = 2, u_lnloc
          u_nodpos(i) = u_nodpos(i-1) + u_ldx
      end do
      do i = 2, cv_lnloc
          lcv_b(i-1) = 0.5 * ( cv_nodpos(i-1) + cv_nodpos(i) )
      end do
      do CV_Liloc = 1, cv_lnloc
          do bgi = 1, cv_bngi
              SBCVFEN( cv_liloc, bgi ) = lagran( ndiff, lcv_b(bgi), cv_liloc, cv_lnloc, cv_nodpos )
          end do
      end do
      do U_Liloc = 1, U_lnloc
          do bgi = 1, cv_bngi
              SBUFEN( u_liloc, bgi ) = lagran( ndiff, lcv_b(bgi), U_liloc, u_lnloc, u_nodpos )
          end do
      end do

      ! Adjust DX for volume integrations...
      allocate(CVN_VOL_ADJ(CV_LNLOC))
      CVN_VOL_ADJ=0.5
      IF(CV_LNLOC==3) THEN
          CVN_VOL_ADJ(1)=0.25
          CVN_VOL_ADJ(2)=0.5
          CVN_VOL_ADJ(3)=0.25
      ENDIF

      ! SET UP THE SURFACE BC'S !SPRINT_TO_DO I DON'T THINK WE ACTUALLY NEED TO DO THIS
      WIC_B_BC_ALL_NODS=0; WIC_T_BC_ALL_NODS=0; WIC_D_BC_ALL_NODS=0
      WIC_U_BC_ALL_NODS=0; RVEC_SUM_T=0.0; RVEC_SUM_D=0.0
      RVEC_SUM_U=0.0; SUF_T_BC_ALL_NODS=0.0
      SUF_D_BC_ALL_NODS=0.0; SUF_U_BC_ALL_NODS=0.0
      DO SELE = 1, Mdims%stotel
          DO CV_SILOC = 1, Mdims%cv_snloc
              CV_NOD = ndgln%suf_cv( (SELE-1)*Mdims%cv_snloc + CV_SILOC )
              WIC_B_BC_ALL_NODS( CV_NOD ) = WIC_B_BC_DIRICHLET
              DO IPHASE = 1, final_phase
                global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
                compact_phase = iphase + (Mdims%npres - 1)*final_phase
                IF ( WIC_T_BC_ALL( 1,global_phase, SELE ) == WIC_T_BC_DIRICHLET ) THEN
                    WIC_T_BC_ALL_NODS( compact_phase, CV_NOD ) = WIC_T_BC_DIRICHLET
                    SUF_T_BC_ALL_NODS( compact_phase, CV_NOD ) = SUF_T_BC_ALL_NODS( compact_phase, CV_NOD ) +&
                            SUF_T_BC_ALL( 1,global_phase, (SELE-1)*Mdims%cv_snloc + CV_SILOC )
                    RVEC_SUM_T( compact_phase, CV_NOD ) = RVEC_SUM_T( compact_phase, CV_NOD ) + 1.0
                END IF
                IF ( WIC_D_BC_ALL( 1,global_phase, SELE ) == WIC_D_BC_DIRICHLET ) THEN
                    WIC_D_BC_ALL_NODS( compact_phase, CV_NOD ) = WIC_D_BC_DIRICHLET
                    SUF_D_BC_ALL_NODS( compact_phase, CV_NOD ) = SUF_D_BC_ALL_NODS( compact_phase, CV_NOD ) + &
                          SUF_D_BC_ALL( 1,global_phase, (SELE-1)*Mdims%cv_snloc + CV_SILOC )
                    RVEC_SUM_D( compact_phase, CV_NOD ) = RVEC_SUM_D( compact_phase, CV_NOD ) + 1.0
                END IF
              end do
          END DO ! DO CV_SILOC=1,Mdims%cv_snloc
          MEAN_U = 0.0!sprint_to_do why don't we do directly SUF_U_BC_ALL_NODS without MEAN_U???
          DO U_SILOC = 1, Mdims%u_snloc
            U_NOD = ndgln%suf_u( (SELE-1)*Mdims%u_snloc + U_SILOC )
            DO IPHASE=1, final_phase
              global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
              compact_phase = iphase + (Mdims%npres - 1)*final_phase
              IF ( WIC_U_BC_ALL( 1, global_phase, SELE ) == WIC_U_BC_DIRICHLET ) THEN
                  MEAN_U(:,compact_phase) = MEAN_U(:,compact_phase) + &
                        SUF_U_BC_ALL( :,global_phase, (SELE-1)*Mdims%u_snloc + U_SILOC ) / REAL(Mdims%u_snloc)
              END IF
            end do
          END DO ! DO U_SILOC=1,Mdims%u_snloc
          DO CV_SILOC = 1, Mdims%cv_snloc
            CV_NOD = ndgln%suf_cv( (SELE-1)*Mdims%cv_snloc + CV_SILOC )

            DO IPHASE=1, final_phase
              global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
              compact_phase = iphase + (Mdims%npres - 1)*final_phase
              IF ( WIC_U_BC_ALL( 1, global_phase, SELE ) == WIC_U_BC_DIRICHLET ) THEN
                  WIC_U_BC_ALL_NODS( compact_phase, CV_NOD ) = WIC_U_BC_DIRICHLET
                  SUF_U_BC_ALL_NODS( :,compact_phase, CV_NOD ) = SUF_U_BC_ALL_NODS( :,compact_phase, CV_NOD ) +&
                            MEAN_U(:,compact_phase)
                  RVEC_SUM_U( compact_phase, CV_NOD ) = RVEC_SUM_U( compact_phase, CV_NOD ) + 1.0
              END IF
            end do

          END DO ! DO CV_SILOC=1,Mdims%cv_snloc
      END DO
      DO CV_NOD = 1, Mdims%cv_nonods
        DO IPHASE = 1, final_phase
          compact_phase = iphase + (Mdims%npres - 1)*final_phase
          IF ( WIC_T_BC_ALL_NODS( compact_phase, CV_NOD ) /= 0 ) SUF_T_BC_ALL_NODS( compact_phase,CV_NOD ) = &
                                      SUF_T_BC_ALL_NODS( compact_phase,CV_NOD ) / RVEC_SUM_T( compact_phase,CV_NOD )
          IF ( WIC_D_BC_ALL_NODS( compact_phase, CV_NOD ) /= 0 ) SUF_D_BC_ALL_NODS( compact_phase,CV_NOD ) = &
                                      SUF_D_BC_ALL_NODS( compact_phase,CV_NOD ) / RVEC_SUM_D( compact_phase,CV_NOD )
          IF ( WIC_U_BC_ALL_NODS( compact_phase, CV_NOD ) /= 0 ) SUF_U_BC_ALL_NODS( :,compact_phase,CV_NOD ) =&
                                      SUF_U_BC_ALL_NODS( :,compact_phase,CV_NOD ) / RVEC_SUM_U( compact_phase,CV_NOD )
        END DO
      END DO
      ! END OF SET UP THE SURFACE B.C'S



      T_ALL => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
      DEN_ALL => extract_tensor_field( packed_state, "PackedDensity" )
      U_ALL => extract_tensor_field( packed_state, "PackedVelocity" )
      PIPE_DIAMETER => extract_scalar_field( state(1), "DiameterPipe" )
      X => EXTRACT_VECTOR_FIELD( PACKED_STATE, "PressureCoordinate" )
      sigma1_pipes => extract_scalar_field( state(1), "Sigma" )
      TMAX_ALL=-INFINY; TMIN_ALL=+INFINY; DENMAX_ALL=-INFINY; DENMIN_ALL = +INFINY

      if (thermal) then
          !Change pointers
          nullify(T_ALL); nullify(DEN_ALL); nullify(U_ALL);
          T_ALL => extract_tensor_field( packed_state, "PackedTemperature" )
          !this is rho * Cp. This is to make it consistent with the advection term in cv-adv-diff
          DEN_ALL => extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )!this is Cp * Rho.
          U_ALL => extract_tensor_field( packed_state, "PackedVelocity" )!for consistency with cv_assemb!Used to be =>PackedNonlinearVelocity like in cv_assemb
          T2_ALL => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )!for multiphase
          only_den_all => extract_tensor_field( packed_state, "PackedDensity" )
      end if

      DO CV_NODI = 1, Mdims%cv_nonods
          IF ( PIPE_DIAMETER%VAL(CV_NODI) > 1e-8 ) THEN
              do count = Mspars%small_acv%fin(cv_nodi), Mspars%small_acv%fin(cv_nodi+1)-1
                  cv_nodj = Mspars%small_acv%col(count)
                  IF ( PIPE_DIAMETER%VAL(CV_NODJ) /= 0.0 ) THEN
                    DO IPHASE=1, final_phase
                      global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
                      compact_phase = iphase + (Mdims%npres - 1)*final_phase
                        TMAX_ALL( compact_phase, CV_NODI ) = max( TMAX_ALL( compact_phase, CV_NODI ), &
                                                            T_ALL%val( 1, global_phase, cv_nodj ) )
                        TMIN_ALL( compact_phase, CV_NODI ) = min( TMIN_ALL( compact_phase, CV_NODI ),&
                                                            T_ALL%val( 1, global_phase, cv_nodj ) )
                        DENMAX_ALL( compact_phase, CV_NODI ) = max( DENMAX_ALL( compact_phase, CV_NODI ),&
                                                            DEN_ALL%val( 1, global_phase, cv_nodj ) )
                        DENMIN_ALL( compact_phase, CV_NODI ) = min( DENMIN_ALL( compact_phase, CV_NODI ),&
                                                            DEN_ALL%val( 1, global_phase, cv_nodj ) )
                    END DO
                  END IF
              END DO
          END IF
      END DO

      allocate( CV_GL_LOC( cv_lnloc ), CV_GL_GL( cv_lnloc ), X_GL_GL( cv_lnloc ), MAT_GL_GL( cv_lnloc ), &
          U_GL_LOC( u_lnloc ), U_GL_GL( u_lnloc ) )
      allocate( L_CVFENX_ALL(cv_lnloc, cv_lngi), L_UFENX_ALL(u_lnloc, cv_lngi) , PIPE_DIAM_GI_vol(cv_lngi) )
      allocate( PIPE_DIAM_GI_suf(cv_bngi) )
      allocate( L_CVFENX_ALL_REVERSED(Mdims%ndim, cv_lnloc, cv_lngi), L_UFEN_REVERSED(cv_lngi, u_lnloc), L_UFEN(u_lnloc, cv_lngi) )
      allocate( suf_detwei(cv_bngi), vol_detwei(cv_lngi), vol_detwei2(cv_lngi) )
      allocate( x_all_corn(Mdims%ndim, NCORNER) )
      mass_pipe = 0.0; MASS_PIPE_FOR_COUP = 0.0
      MASS_CVFEM2PIPE = 0.0; MASS_PIPE2CVFEM = 0.0; MASS_CVFEM2PIPE_TRUE = 0.0
      INV_SIGMA = 0.0
      !Populate INV_SIGMA inside the pipes to ensure that flux from pipes to the domain can happen
      DO IPHASE = wells_first_phase, final_phase*2
          call assign_val(INV_SIGMA(iphase, : ),1./sigma1_pipes%val)
      end do


      DO k = 1, size(eles_with_pipe)
          ELE = eles_with_pipe(k)%ele!Element with pipe
          X_ALL_CORN(:,1:NCORNER) = x%val(:, ndgln%x( ( ELE - 1 ) * Mdims%cv_nloc + CV_LOC_CORNER(1:NCORNER)) )

          DO IPIPE = 1, eles_with_pipe(k)%npipes
              ! DEFINE CV_LILOC:
              CV_LILOC = 1
              ICORNER = eles_with_pipe(k)%pipe_corner_nds1( IPIPE )
              ICORNER1 = ICORNER
              CV_ILOC = CV_LOC_CORNER( ICORNER )
              CV_GL_LOC( CV_LILOC ) = CV_ILOC
              CV_LILOC = CV_LNLOC
              ICORNER = eles_with_pipe(k)%pipe_corner_nds2( IPIPE )
              ICORNER2 = ICORNER
              CV_ILOC = CV_LOC_CORNER( ICORNER )
              CV_GL_LOC( CV_LILOC ) = CV_ILOC
              IF ( CV_QUADRATIC ) CV_GL_LOC( 2 ) = CV_MID_SIDE( ICORNER1, ICORNER2 )
              DO CV_LILOC = 1, CV_LNLOC
                  CV_ILOC = CV_GL_LOC( CV_LILOC )
                  X_ILOC = CV_ILOC
                  CV_GL_GL( CV_LILOC ) = ndgln%cv( (ELE-1)*Mdims%cv_nloc + CV_ILOC )
                  X_GL_GL( CV_LILOC ) = ndgln%x( (ELE-1)*Mdims%cv_nloc + X_ILOC )
                  MAT_GL_GL( CV_LILOC ) = ndgln%mat( (ELE-1)*Mdims%cv_nloc + CV_ILOC )
              END DO
              ! DEFINE U_LILOC:
              if (U_P0DG) then
                  U_GL_LOC = 1
              else
                  U_LILOC = 1
                  ICORNER = eles_with_pipe(k)%pipe_corner_nds1( IPIPE )
                  U_ILOC = U_LOC_CORNER( ICORNER )
                  U_GL_LOC( U_LILOC ) = U_ILOC
                  U_LILOC = U_LNLOC
                  ICORNER = eles_with_pipe(k)%pipe_corner_nds2( IPIPE )
                  U_ILOC = U_LOC_CORNER( ICORNER )
                  U_GL_LOC( U_LILOC ) = U_ILOC
                  IF ( U_QUADRATIC ) U_GL_LOC( 2 ) = U_MID_SIDE( ICORNER1, ICORNER2 )
              end if
              DO U_LILOC = 1, U_LNLOC
                  U_ILOC = U_GL_LOC( U_LILOC )
                  U_GL_GL( U_LILOC ) = ndgln%u( (ELE-1)*Mdims%u_nloc + U_ILOC )
              END DO
              DIRECTION = X_ALL_CORN( :, ICORNER2 ) - X_ALL_CORN( :, ICORNER1 )
              DX = SQRT( SUM( DIRECTION**2 ) )
              DIRECTION = DIRECTION / DX
              ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
              L_CVFENX_ALL(:,:) = 2.0 * cvn_femlx(:,:) / DX
              L_UFENX_ALL(:,:) = 2.0 * UNLX(:,:) / DX

              IF(PIPE_MIN_DIAM) THEN
                  MIN_DIAM = MINVAL( PIPE_diameter%val( CV_GL_GL( : ) ) )
                  PIPE_DIAM_GI_VOL = MIN_DIAM
                  PIPE_DIAM_GI_SUF = MIN_DIAM
              ELSE
                  PIPE_DIAM_GI_VOL = 0.0
                  PIPE_DIAM_GI_SUF = 0.0
                  DO CV_LILOC = 1, CV_LNLOC
                      CV_KNOD = CV_GL_GL( CV_LILOC )
                      PIPE_DIAM_GI_VOL = PIPE_DIAM_GI_VOL + PIPE_diameter%val( CV_KNOD ) * cvn_fem( CV_LILOC, : )
                      PIPE_DIAM_GI_SUF = PIPE_DIAM_GI_SUF + PIPE_diameter%val( CV_KNOD ) * SBCVFEN( CV_LILOC, : )
                  END DO
                  PIPE_DIAM_GI_VOL = max( 0.0, PIPE_DIAM_GI_VOL )
                  PIPE_DIAM_GI_SUF = max( 0.0, PIPE_DIAM_GI_SUF )
              ENDIF
              DO IDIM = 1, Mdims%ndim
                  L_CVFENX_ALL_REVERSED(IDIM,:,:) = cvn_femlx * DIRECTION(IDIM)
              END DO
              DO U_LILOC = 1, U_LNLOC
                  L_UFEN_REVERSED( :, U_LILOC ) = un( U_LILOC, : )
              END DO
              ! Calculate element angle sweeped out by element and pipe
              IF ( Mdims%ndim == 2 ) THEN
                  ELE_ANGLE = PI
              ELSE
                  ! find the nodes other than the pipe end corner nodes...
                  ICORNER3 = 0
                  DO ICORNER = 1, NCORNER
                      IF ( ICORNER /= ICORNER1 ) THEN
                          IF ( ICORNER /= ICORNER2 ) THEN
                              IF ( ICORNER3 == 0 ) THEN
                                  ICORNER3 = ICORNER
                              ELSE
                                  ICORNER4 = ICORNER
                              END IF
                          END IF
                      END IF
                  END DO
                  ELE_ANGLE = CALC_ELE_ANGLE_3D( X_ALL_CORN(:, ICORNER1 ), X_ALL_CORN(:, ICORNER2) , &
                      &                         X_ALL_CORN(:, ICORNER3 ), X_ALL_CORN(:, ICORNER4) )
              END IF
              ! Adjust according to the volume of the pipe...
              SUF_DETWEI( : )  = 1.0               * PI * ( ( 0.5*PIPE_DIAM_GI_suf )**2 ) * ELE_ANGLE / ( 2.0 * PI )
              VOL_DETWEI( : )  = cvweigh( : ) * DX * PI * ( ( 0.5*PIPE_DIAM_GI_vol )**2 ) * ELE_ANGLE / ( 2.0 * PI )
              ! We do NOT multiply by r**2 here because we have not divided by r**2 in the mass transfer correlation (in cv_assemb)
              VOL_DETWEI2( : ) = cvweigh( : ) * DX * PI *                                      ELE_ANGLE / ( 2.0 * PI )
              DO IDIM = 1, Mdims%ndim
                  L_CVFENX_ALL_REVERSED( IDIM, :, : ) = L_CVFENX_ALL( :, : ) * DIRECTION( IDIM )
              END DO
              DO U_LILOC = 1, U_LNLOC
                  L_UFEN_REVERSED( :, U_LILOC ) =  UN( U_LILOC, : )
              END DO
              ! Add contributions from the volume...
              IF( LUMP_COUPLING_RES_PIPES ) THEN
                  DO CV_LILOC = 1, CV_LNLOC
                      CV_NODI = CV_GL_GL(CV_LILOC)
                      MASS_PIPE(CV_NODI) = MASS_PIPE(CV_NODI) + SUM( CVN(CV_LILOC,:) *  CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI( : ) )
                      MASS_PIPE_FOR_COUP(CV_NODI) = MASS_PIPE_FOR_COUP(CV_NODI) + SUM( CVN(CV_LILOC,:) *  CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI2( : ) )
                      DO COUNT2 = Mspars%CMC%fin(CV_NODI), Mspars%CMC%fin(CV_NODI+1)-1
                          IF ( CV_NODI==Mspars%CMC%col(count2)) count=count2
                      end do
                      MASS_CVFEM2PIPE(COUNT) = MASS_CVFEM2PIPE(COUNT) + SUM( CVN(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI2( : ) )
                      MASS_PIPE2CVFEM(COUNT) = MASS_PIPE2CVFEM(COUNT) + SUM( CVN(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI2( : ) ) ! this is possibly bugged...
                      DO CV_LJLOC = 1, CV_LNLOC
                          CV_NODJ = CV_GL_GL(CV_LJLOC)
                          DO COUNT = Mspars%CMC%fin(CV_NODI), Mspars%CMC%fin(CV_NODI+1)-1
                              IF ( CV_NODI==CV_NODJ ) THEN
                                  MASS_CVFEM2PIPE_TRUE(COUNT) = MASS_CVFEM2PIPE_TRUE(COUNT) + SUM( CVN(CV_LILOC,:) * CVN_FEM(CV_LJLOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI( : ) )
                              END IF
                          END DO
                      END DO
                  END DO
              ELSE
                  DO CV_LILOC = 1, CV_LNLOC
                      CV_NODI = CV_GL_GL(CV_LILOC)
                      MASS_PIPE(CV_NODI) = MASS_PIPE(CV_NODI) + SUM( CVN(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI( : ) )
                      DO CV_LJLOC = 1, CV_LNLOC
                          CV_NODJ = CV_GL_GL(CV_LJLOC)
                          DO COUNT = Mspars%CMC%fin(CV_NODI), Mspars%CMC%fin(CV_NODI+1)-1
                              IF ( CV_NODI==CV_NODJ ) THEN
                                  MASS_CVFEM2PIPE(COUNT) = MASS_CVFEM2PIPE(COUNT) + SUM( CVN(CV_LILOC,:) * CVN_FEM(CV_LJLOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI2( : ) )
                                  MASS_PIPE2CVFEM(COUNT) = MASS_PIPE2CVFEM(COUNT) + SUM( CVN(CV_LJLOC,:) * CVN_FEM(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI2( : ) )
                                  MASS_CVFEM2PIPE_TRUE(COUNT) = MASS_CVFEM2PIPE_TRUE(COUNT) + SUM( CVN(CV_LJLOC,:) * CVN_FEM(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI( : ) )
                              END IF
                          END DO
                      END DO
                  END DO
              END IF
              ! Calculate INV_SIGMA_APPROX for the 1st n_in_pres phases...
              ! Determine the tangent and bi-normal vectors from the normal NormX, NormY, NormZ:
              IF(Mdims%ndim==2) THEN
                  RZ=0.0
                  CALL GET_TANG_BINORM( DIRECTION(1), DIRECTION(2), RZ,           T1(1), T1(2), R1,    T2(1), T2(2), R2, 1 )
              ELSE
                  CALL GET_TANG_BINORM( DIRECTION(1), DIRECTION(2), DIRECTION(3), T1(1), T1(2), T1(3), T2(1), T2(2), T2(3), 1 )
              END IF
              DO CV_LILOC = 1, CV_LNLOC
                  MAT_NODI = MAT_GL_GL(CV_LILOC)
                  CV_NODI = CV_GL_GL(CV_LILOC)
                  DO IPHASE = 1, final_phase
                      if (is_porous_media) then
                              TT1 = MATMUL( OPT_VEL_UPWIND_COEFS_NEW(:,:,IPHASE, MAT_NODI), T1 )
                              T1TT1 = SUM(T1*TT1)
                          IF ( Mdims%ndim==3 ) THEN
                              TT2 = MATMUL( OPT_VEL_UPWIND_COEFS_NEW(:,:,IPHASE, MAT_NODI), T2 )
                              T1TT2 = SUM( T1*TT2 )
                              T2TT1 = SUM( T2*TT1 )
                              T2TT2 = SUM( T2*TT2 )
                          ELSE
                              T1TT2 = 0.0
                              T2TT1 = 0.0
                              T2TT2 = T1TT1
                          END IF
                      else
                          T1TT1 =SUM(T1*T1)
                          IF ( Mdims%ndim==3 ) THEN
                              T1TT2 = SUM( T1*T2 )
                              T2TT1 = SUM( T2*T1 )
                              T2TT2 = SUM( T2*T2 )
                          else
                              T1TT2 = 0.0
                              T2TT1 = 0.0
                              T2TT2 = T1TT1
                          end if
                      end if
                      DET_SQRT = SQRT( ABS( T1TT1*T2TT2 - T1TT2*T2TT1 ) )
                      INV_SIGMA_ND = 1.0 / MAX( 1.E-25, DET_SQRT)
                      INV_SIGMA(IPHASE,CV_NODI) = INV_SIGMA(IPHASE,CV_NODI) + INV_SIGMA_ND * SUM( CVN(CV_LILOC,:) * CVN_VOL_ADJ(CV_LILOC) * VOL_DETWEI( : ) )
                      ! For the nano laterals...
                      if (is_porous_media) then
                          NN1 = MATMUL( OPT_VEL_UPWIND_COEFS_NEW(:,:,IPHASE, MAT_NODI), DIRECTION )
                          N1NN1 = SUM( DIRECTION*NN1 )
                      else
                          N1NN1 = SUM( DIRECTION )
                      end if
                      INV_SIGMA_NANO_ND = 1.0/MAX(1.E-25, N1NN1)
                  END DO
              END DO
              ! The number of CV basis functions...
              DO BGI = 1, cv_bngi
                  DO ILOOP = 1, 2
                      IF ( ILOOP==1 ) THEN
                          CV_LILOC = BGI
                          CV_LJLOC = BGI + 1
                          CV_NODI = CV_GL_GL(CV_LILOC)
                          CV_NODJ = CV_GL_GL(CV_LJLOC)
                          direction_norm = + direction
                      ELSE
                          CV_LILOC = BGI + 1
                          CV_LJLOC = BGI
                          CV_NODI = CV_GL_GL(CV_LILOC)
                          CV_NODJ = CV_GL_GL(CV_LJLOC)
                          direction_norm = - direction
                      END IF
                      TUPWIND_OUT=0.0; DUPWIND_OUT=0.0
                      TUPWIND_IN=0.0; DUPWIND_IN=0.0

                      !Only for the wells domain
                      DO IPHASE=1, final_phase
                        global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
                        compact_phase = iphase + (Mdims%npres - 1)*final_phase
                      ! DO IPHASE = n_in_pres+1, final_phase
                          ! CV incomming T:
                          IF ( T_ALL%val( 1, global_phase, CV_NODI ) > T_ALL%val( 1, global_phase, CV_NODJ ) ) THEN
                              TUPWIND_OUT( compact_phase ) = TMAX_ALL( compact_phase, CV_NODI )
                          ELSE
                              TUPWIND_OUT( compact_phase ) = TMIN_ALL( compact_phase, CV_NODI )
                          END IF
                          IF ( DEN_ALL%val( 1, global_phase, CV_NODI ) > DEN_ALL%val( 1, global_phase, CV_NODJ ) ) THEN
                              DUPWIND_OUT( compact_phase ) = DENMAX_ALL( compact_phase, CV_NODI )
                          ELSE
                              DUPWIND_OUT( compact_phase ) = DENMIN_ALL( compact_phase, CV_NODI )
                          END IF
                          ! CV outgoing T:
                          IF ( T_ALL%val( 1, global_phase, CV_NODI ) < T_ALL%val( 1, global_phase, CV_NODJ ) ) THEN
                              TUPWIND_IN( compact_phase ) = TMAX_ALL( compact_phase, CV_NODJ )
                          ELSE
                              TUPWIND_IN( compact_phase ) = TMIN_ALL( compact_phase, CV_NODJ )
                          END IF
                          IF ( DEN_ALL%val( 1, global_phase, CV_NODI ) < DEN_ALL%val( 1, global_phase, CV_NODJ ) ) THEN
                              DUPWIND_IN( compact_phase ) = DENMAX_ALL( compact_phase, CV_NODJ )
                          ELSE
                              DUPWIND_IN( compact_phase ) = DENMIN_ALL( compact_phase, CV_NODJ )
                          END IF
                      END DO
                      ! Value of sigma in the force balance eqn...
                      INV_SIGMA_GI = 1.0!sprint_to_do remove INV_SIGMA_GI??
                      ! Velocity in the pipe
                      UGI_ALL(:,:) = 0.0
                      DO U_LKLOC = 1, U_LNLOC
                          U_KNOD = U_GL_GL(U_LKLOC)
                          DO IPHASE=1, final_phase
                            global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
                            compact_phase = iphase + (Mdims%npres - 1)*final_phase
                            UGI_ALL(:,compact_phase) = UGI_ALL(:,compact_phase) + SBUFEN( u_lkloc, bgi ) * U_ALL%val(:,global_phase,u_kNOD)
                          end do
                      END DO
                      NDOTQ = matmul(DIRECTION_norm, UGI_all)

                      ! When NDOTQ == 0, INCOME_J has to be 1 as well, not 0
                      WHERE ( NDOTQ <= 0.0 )
                          INCOME_J = 0.0
                      ELSE WHERE
                          INCOME_J = 1.0
                      END WHERE
                      INCOME = 1.0 - INCOME_J
                      ! high order values...
                      FEMTGI=0.0 ; FEMDGI=0.0
                      DO CV_LKLOC = 1, CV_LNLOC
                          CV_KNOD = CV_GL_GL(CV_LKLOC)
                          DO IPHASE=1, final_phase
                            global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
                            compact_phase = iphase + (Mdims%npres - 1)*final_phase
                            FEMTGI(compact_phase) = FEMTGI(compact_phase) + &
                                    CVN_FEM(CV_LJLOC,BGI) * T_ALL%val( 1, global_phase, CV_KNOD)
                            FEMDGI(compact_phase) = FEMDGI(compact_phase) + &
                                    CVN_FEM(CV_LJLOC,BGI) * DEN_ALL%val( 1, global_phase, CV_KNOD)
                          end do
                      END DO
                      FEMDGI = max( 0.0,FEMDGI )
                      DO IPHASE=1, final_phase
                        global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
                        compact_phase = iphase + (Mdims%npres - 1)*final_phase
                        T_CV_NODI(compact_phase) = T_ALL%val( 1, global_phase, CV_NODI)
                        T_CV_NODJ(compact_phase) = T_ALL%val( 1, global_phase, CV_NODJ)
                        D_CV_NODI(compact_phase) = DEN_ALL%val( 1, global_phase, CV_NODI)
                        D_CV_NODJ(compact_phase) = DEN_ALL%val( 1, global_phase, CV_NODJ)
                      end do
                      IF ( UPWIND_PIPES ) THEN ! Used for testing...
                          LIMT = T_CV_NODI*(1.0-INCOME) + T_CV_NODJ*INCOME
                          LIMD = D_CV_NODI*(1.0-INCOME) + D_CV_NODJ*INCOME
                      ELSE
                          ! Call the limiter for T...
                          CALL ONVDLIM_ANO_MANY( final_phase*2, &!Only do it for the phases within the pipe
                              LIMT, FEMTGI, INCOME, &
                              T_CV_NODI, T_CV_NODJ, XI_LIMIT, &
                              TUPWIND_IN, TUPWIND_OUT , &
                              memory_limiters(1:Mdims%nphase), memory_limiters(Mdims%nphase + 1:Mdims%nphase*2),&
                              memory_limiters(2*Mdims%nphase + 1:Mdims%nphase*3), memory_limiters(3*Mdims%nphase + 1:Mdims%nphase*4),&
                              memory_limiters(4*Mdims%nphase + 1:Mdims%nphase*5), memory_limiters(5*Mdims%nphase + 1:Mdims%nphase*6) )
                          ! Call the limiter for D...
                          CALL ONVDLIM_ANO_MANY( final_phase*2, &
                              LIMD, FEMDGI, INCOME, &
                              D_CV_NODI, D_CV_NODJ, XI_LIMIT, &
                              DUPWIND_IN, DUPWIND_OUT, &
                              memory_limiters(1:Mdims%nphase), memory_limiters(Mdims%nphase + 1:Mdims%nphase*2),&
                              memory_limiters(2*Mdims%nphase + 1:Mdims%nphase*3), memory_limiters(3*Mdims%nphase + 1:Mdims%nphase*4),&
                              memory_limiters(4*Mdims%nphase + 1:Mdims%nphase*5), memory_limiters(5*Mdims%nphase + 1:Mdims%nphase*6) )
                      END IF
                      LIMDT = LIMD * LIMT
                      IF ( GETCT ) THEN ! Obtain the CV discretised Mmat%CT eqations plus RHS
                          DO U_LKLOC = 1, U_LNLOC
                              U_KNOD = U_GL_GL(U_LKLOC)
                              DO IDIM = 1, Mdims%ndim
                                  CT_CON(IDIM,:) = SBUFEN( U_LKLOC, BGI ) * LIMDT * suf_DETWEI( BGI ) * DIRECTION_norm(IDIM) * INV_SIGMA_GI / D_CV_NODI
                              END DO
                              ! Put into Mmat%CT matrix...
                              DO COUNT = Mspars%CT%fin(CV_NODI), Mspars%CT%fin(CV_NODI+1)-1
                                  IF ( Mspars%CT%col(COUNT)==U_KNOD ) then
                                      Mmat%CT( :, wells_first_phase:final_phase*2, COUNT ) = &
                                          Mmat%CT( :, wells_first_phase:final_phase*2, COUNT ) + CT_CON( :, wells_first_phase:final_phase*2 )
                                      exit
                                  end if
                              END DO
                          END DO
                      END IF
                      IF ( GETCV_DISC ) THEN

                          FVT = T_CV_NODI*(1.0-INCOME) + T_CV_NODJ*INCOME
                          ! Put results into the RHS vector
                          LOC_CV_RHS_I = 0.0; LOC_MAT_II = 0.; LOC_MAT_IJ = 0.
                          do iphase = wells_first_phase, final_phase*2
                              LOC_CV_RHS_I( IPHASE ) =  LOC_CV_RHS_I( IPHASE ) &
                                    ! subtract 1st order adv. soln.
                                  + suf_DETWEI( bGI ) * NDOTQ(IPHASE) * LIMD(IPHASE) * FVT(IPHASE) &
                                  - suf_DETWEI( bGI ) * NDOTQ(IPHASE) * LIMDT(IPHASE) ! hi order adv
                              if (.not.conservative_advection)  LOC_CV_RHS_I( IPHASE ) = LOC_CV_RHS_I( IPHASE ) &
                                  - suf_DETWEI( bGI ) * NDOTQ(IPHASE) * LIMD(IPHASE) * T_CV_NODI(IPHASE) &
                                  + suf_DETWEI( bGI ) * NDOTQ(IPHASE) * LIMD(IPHASE) * T_CV_NODI(IPHASE)
                          end do
                          ! Put into matrix...
                          do iphase = wells_first_phase, final_phase*2
                              LOC_MAT_II(iphase) = LOC_MAT_II(iphase) +suf_DETWEI( bGI ) * NDOTQ(iphase) * ( 1. - INCOME(iphase) ) * LIMD(iphase)
                              LOC_MAT_IJ(iphase) = LOC_MAT_IJ(iphase) +suf_DETWEI( bGI ) * NDOTQ(iphase) * INCOME(iphase) * LIMD(iphase)
                              if (.not.conservative_advection) LOC_MAT_II(iphase) = LOC_MAT_II(iphase)  - suf_DETWEI( bGI ) * NDOTQ(iphase) * LIMD(iphase)

                          end do
                          do iphase = wells_first_phase, final_phase*2
                            assembly_phase = iphase
                            !For the RHS collapsing to assemble into phase 2 can be done just here
                            if (assemble_collapsed_to_one_phase) assembly_phase = 2
                            call addto(Mmat%CV_RHS,assembly_phase, CV_NODI,LOC_CV_RHS_I(IPHASE))
                            !Introduce the information into the petsc_ACV matrix
                            call addto(Mmat%petsc_ACV,assembly_phase,assembly_phase,cv_nodi,cv_nodi, LOC_MAT_II(iphase) )
                            call addto(Mmat%petsc_ACV,assembly_phase,assembly_phase,cv_nodi,cv_nodj, LOC_MAT_IJ(iphase) )
                          end do
                      END IF
                  END DO ! DO ILOOP = 1, 2
              END DO ! DO BGI = 1, CV_BNGI
              ! Add velocity and saturation b.c's to matrix and rhs...
              JCV_NOD1 = ndgln%cv( (ELE-1)*Mdims%cv_nloc + CV_LOC_CORNER( ICORNER1 ) )
              JCV_NOD2 = ndgln%cv( (ELE-1)*Mdims%cv_nloc + CV_LOC_CORNER( ICORNER2 ) )
              JCV_NOD=0
              IF ( WIC_B_BC_ALL_NODS( JCV_NOD1 ) == WIC_B_BC_DIRICHLET ) THEN
                  CV_LILOC = 1
                  JCV_NOD = JCV_NOD1
                  U_LILOC = 1
                  JU_NOD = U_GL_GL( U_LILOC )
                  direction_norm = - direction ! for the b.c it must be -ve at the bottom of element
              end if
              IF ( WIC_B_BC_ALL_NODS( JCV_NOD2 ) == WIC_B_BC_DIRICHLET ) THEN
                  CV_LILOC = CV_LNLOC
                  JCV_NOD = JCV_NOD2
                  U_LILOC = U_LNLOC
                  JU_NOD = U_GL_GL( U_LILOC )
                  direction_norm = + direction ! for the b.c it must be +ve at the top of element
              END IF

              IF ( JCV_NOD /= 0 ) THEN
                  !We are in the boundary of the domain
                  PIPE_DIAM_END = PIPE_diameter%val( JCV_NOD )
                  NDOTQ = 0.0
                  DO iphase= 1,final_phase
                    compact_phase = iphase + (Mdims%npres - 1)*final_phase
                    if (WIC_U_BC_ALL_NODS( compact_phase, JCV_NOD ) == WIC_U_BC_DIRICHLET ) then
                        NDOTQ(compact_phase) = dot_product( direction_norm, SUF_U_BC_ALL_NODS(:,compact_phase,JCV_NOD) )
                    else
                        global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
                        NDOTQ(compact_phase) = dot_product( direction_norm, U_ALL%val(:,global_phase,JU_NOD) )
                    end if
                  END DO
                  INCOME = 0.5 * ( 1. + SIGN( 1.0, -NDOTQ ) )
                  LIMT=0.0
                  LIMD=0.0
                  FVT =0.0
                  DO IPHASE=1,final_phase
                    compact_phase = iphase + (Mdims%npres - 1)*final_phase
                    global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
                    IF ( WIC_T_BC_ALL_NODS( compact_phase, JCV_NOD ) == WIC_T_BC_DIRICHLET ) THEN
                        LIMT(compact_phase)=T_ALL%val(1,global_phase,JCV_NOD)*(1.0-INCOME(compact_phase)) + SUF_T_BC_ALL_NODS(compact_phase,JCV_NOD)*INCOME(compact_phase)
                    ELSE
                        LIMT(compact_phase)=T_ALL%val(1,global_phase,JCV_NOD)
                    END IF
                    FVT(compact_phase) = LIMT(compact_phase)
                  END DO
                  DO IPHASE=1,final_phase
                    compact_phase = iphase + (Mdims%npres - 1)*final_phase
                    global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
                    IF ( WIC_D_BC_ALL_NODS( compact_phase, JCV_NOD ) == WIC_D_BC_DIRICHLET ) THEN
                        LIMD(compact_phase)=DEN_ALL%val(1,global_phase,JCV_NOD)*(1.0-INCOME(compact_phase)) + SUF_D_BC_ALL_NODS(compact_phase,JCV_NOD)*INCOME(compact_phase)
                    ELSE
                        LIMD(compact_phase)=DEN_ALL%val(1,global_phase,JCV_NOD)
                    END IF
                  END DO
                  LIMDT = LIMD * LIMT


                  ! Add in Mmat%C matrix contribution: (DG velocities)
                  ! In this section we multiply the shape functions over the GI points. i.e: we perform the integration
                  ! over the element of the pressure like source term.
                  ! Put into matrix
                  ! Prepare aid variable NMX_ALL to improve the speed of the calculations
                  suf_area = PI * ( (0.5*PIPE_DIAM_END)**2 ) * ELE_ANGLE / ( 2.0 * PI )
                  IF ( GETCT ) THEN ! Obtain the CV discretised Mmat%CT eqations plus RHS on the boundary...


                      DO IDIM = 1, Mdims%ndim
                          CT_CON(IDIM,:) = LIMDT * suf_area * DIRECTION_NORM(IDIM) * INV_SIGMA_GI / DEN_ALL%val(1,:,JCV_NOD)
                      END DO
                      ! Put into Mmat%CT matrix...
                      COUNT2=0
                      DO COUNT = Mspars%CT%fin(JCV_NOD), Mspars%CT%fin(JCV_NOD+1)-1
                          IF ( Mspars%CT%col(COUNT)==JU_NOD ) COUNT2=COUNT
                      END DO
                      LOC_CT_RHS_U_ILOC = 0.0
                      DO IPHASE = wells_first_phase, final_phase*2
                          IF ( WIC_U_BC_ALL_NODS( IPHASE, JCV_NOD ) == WIC_U_BC_DIRICHLET ) THEN
                              DO IDIM = 1, Mdims%ndim
                                  LOC_CT_RHS_U_ILOC(IPHASE) = LOC_CT_RHS_U_ILOC(IPHASE) &
                                      - CT_CON(IDIM,IPHASE) * SUF_U_BC_ALL_NODS(IDIM,IPHASE,JCV_NOD)
                              END DO
                          ELSE
                              Mmat%CT( :, IPHASE, COUNT2 ) = Mmat%CT( :, IPHASE, COUNT2 ) + CT_CON( :, IPHASE )
                          END IF
                      END DO
                      call addto( Mmat%CT_RHS, Mdims%npres, jcv_nod, &
                          sum( LOC_CT_RHS_U_ILOC( 1+(Mdims%npres-1)*final_phase : Mdims%npres*final_phase ) ) )
                  END IF ! IF ( GETCT ) THEN

                  IF ( GETCV_DISC ) THEN ! this is on the boundary...
                      ! Put results into the RHS vector
                      LOC_CV_RHS_I = 0.0
                      LOC_MAT_II = 0.
                      DO IPHASE= 1, final_phase
                        compact_phase = iphase + (Mdims%npres - 1)*final_phase
                        LOC_CV_RHS_I( compact_phase ) =  LOC_CV_RHS_I( compact_phase ) &
                            ! subtract 1st order adv. soln.
                            + suf_area  * NDOTQ(compact_phase) * ( 1. - INCOME(compact_phase) ) * LIMD(compact_phase) * FVT(compact_phase)  &
                            - suf_area * NDOTQ(compact_phase) * LIMDT(compact_phase) ! hi order adv
                          if (.not.conservative_advection)  then
                            global_phase = iphase + (Mdims%npres - 1)*Mdims%n_in_pres
                            LOC_CV_RHS_I( compact_phase ) = LOC_CV_RHS_I( compact_phase ) &
                                - suf_area * NDOTQ(compact_phase) * LIMD(compact_phase) * T_ALL%val(1,global_phase,JCV_NOD) &
                                + suf_area * NDOTQ(compact_phase) * LIMD(compact_phase) * T_ALL%val(1,global_phase,JCV_NOD)
                          end if
                      end do
                      ! Put into matrix...
                      do iphase = wells_first_phase, final_phase*2
                        LOC_MAT_II(iphase) = LOC_MAT_II(iphase) + suf_area * NDOTQ(iphase) * ( 1. - INCOME(iphase) ) * LIMD(iphase)
                        if (.not.conservative_advection) LOC_MAT_II(iphase) = LOC_MAT_II(iphase) -suf_area * NDOTQ(iphase) * LIMD(iphase)
                      end do
                      do iphase = wells_first_phase, final_phase*2
                        assembly_phase = iphase
                        !For the RHS collapsing to assemble into phase 2 can be done just here
                        if (assemble_collapsed_to_one_phase) assembly_phase = 2
                        call addto(Mmat%CV_RHS,assembly_phase, JCV_NOD,LOC_CV_RHS_I(IPHASE))
                        !Introduce the information into the petsc_ACV matrix
                        call addto(Mmat%petsc_ACV,assembly_phase,assembly_phase,JCV_NOD,JCV_NOD, LOC_MAT_II(iphase) )
                      end do

                  ENDIF ! ENDOF IF ( GETCV_DISC ) THEN

                  !Finally store fluxes across all the boundaries either for mass conservation check or mass outflux
                  if (compute_outfluxes) then
                    sele = sele_from_cv_nod(Mdims, ndgln, JCV_NOD)!We need SELE for this, not ideal but this operation is not done much overall
                    call update_outfluxes(bcs_outfluxes, outfluxes, sele, JCV_NOD,  &
                        NDOTQ * suf_area * LIMT, NDOTQ * suf_area * LIMDT, &!Vol_flux and Mass_flux
                        T_ALL, temp_field, salt_field, wells_first_phase, final_phase*2 )
                  end if

              ENDIF ! ENDOF IF(JCV_NOD.NE.0) THEN
          END DO ! DO IPIPE2 = 1, NPIPES_IN_ELE
      END DO ! DO ELE = 1, Mdims%totele
      DO IPHASE = 1, final_phase
          INV_SIGMA(IPHASE,:) = INV_SIGMA(IPHASE,:) / MAX( MASS_PIPE, 1.E-15 )
      END DO
      IF ( GETCV_DISC ) THEN
          do iphase = wells_first_phase, final_phase*2
            assembly_phase = iphase
            if (assemble_collapsed_to_one_phase) assembly_phase = 2
              do cv_nodi = 1, Mdims%cv_nonods
                  if ( pipe_diameter%val(cv_nodi) <= 1e-8 ) then
                      cv_nodj = cv_nodi
                      i_indx = Mmat%petsc_ACV%row_numbering%gnn2unn( cv_nodi, assembly_phase )
                      j_indx = Mmat%petsc_ACV%column_numbering%gnn2unn( cv_nodj, assembly_phase )
                      call MatSetValue( Mmat%petsc_ACV, i_indx, j_indx, 1.0, INSERT_VALUES, ierr )
                  end if
              end do
          end do
       end if
  CONTAINS
    !>@brief: This sub calculates the limited face values TDADJ(1...SNGI) from the central
    !> difference face values TDCEN(1...SNGI) using a NVD shceme.
    !> INCOME(1...SNGI)=1 for incomming to element ELE  else =0.
    !> LIBETA is the flux limiting parameter.
    !> TDMAX(PELE)=maximum of the surrounding 6 element values of element PELE.
    !> TDMIN(PELE)=minimum of the surrounding 6 element values of element PELE.
    !> PELEOT=element at other side of current face.
    !> ELEOT2=element at other side of the element ELEOTH.
    !> ELESID=element next to oposing current face.
    !> DENOIN, CTILIN, DENOOU, CTILOU, FTILIN, FTILOU => memory
    !> The elements are arranged in this order: ELEOT2,ELE, PELEOT, ELESID.
    !> This sub finds the neighbouring elements. Suppose that this is the face IFACE.
    !>---------------------------------------------------
    !>|   ELEOT2   |   ELEOTH   |   ELE     |   ELESID   |
    !>---------------------------------------------------
    !> TAIN         THALF       TAOUT
    !>---------------------------------------------------
    !>TEXTIN
    !>TEXOUT<
    !>---------------------------------------------------
    PURE SUBROUTINE ONVDLIM_ANO_MANY( NFIELD, &
        TDLIM, TDCEN, INCOME, &
        ETDNEW_PELE, ETDNEW_PELEOT, XI_LIMIT,  &
        TUPWIN, TUPWI2, DENOIN, CTILIN, DENOOU, CTILOU, FTILIN, FTILOU )
        implicit none
        INTEGER, intent( in ) :: NFIELD
        REAL, DIMENSION( NFIELD ), intent( inout ) :: TDLIM
        REAL, DIMENSION( NFIELD ), intent( in ) :: TDCEN, INCOME, XI_LIMIT, TUPWIN, TUPWI2
        REAL, DIMENSION( NFIELD ), intent( in ) :: ETDNEW_PELE, ETDNEW_PELEOT
        real, dimension( NFIELD ), intent(inout) :: DENOIN, CTILIN, DENOOU, CTILOU, FTILIN, FTILOU
        ! Local variables
        REAL, PARAMETER :: TOLER=1.0E-10
        ! Calculate normalisation parameters for incomming velocities
        DENOIN = ( ETDNEW_PELE - TUPWIN )
        where( ABS( DENOIN ) < TOLER )
            DENOIN = SIGN( TOLER, DENOIN )
        end where
        CTILIN = ( ETDNEW_PELEOT - TUPWIN ) / DENOIN
        ! Calculate normalisation parameters for out going velocities
        DENOOU = ( ETDNEW_PELEOT - TUPWI2 )
        where( ABS( DENOOU ) < TOLER )
            DENOOU = SIGN( TOLER, DENOOU )
        end where
        CTILOU = ( ETDNEW_PELE - TUPWI2 ) / DENOOU
        FTILIN = ( TDCEN - TUPWIN ) / DENOIN
        FTILOU = ( TDCEN - TUPWI2 ) / DENOOU
        ! Velocity is going out of element
        TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW_MANY( FTILIN, CTILIN, XI_LIMIT ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW_MANY( FTILOU, CTILOU, XI_LIMIT ) * DENOOU )
        TDLIM = MAX( TDLIM, 0.0 )
        RETURN
    END SUBROUTINE ONVDLIM_ANO_MANY

    !!>@brief: Obtain sele from a cv_nod that is on the boundary
    !if not found then returns -1. Important to read BCs
    real function sele_from_cv_nod(Mdims, ndgln, cv_jnod)
        implicit none
        integer, intent(in) ::cv_jnod
        type(multi_ndgln), intent(in) :: ndgln
        type(multi_dimensions), intent(in) :: Mdims
        !Local variables
        integer :: sele, cv_siloc

        sele_from_cv_nod = -1
        do sele = 1, Mdims%stotel
            do cv_siloc = 1, Mdims%cv_snloc
                if (ndgln%suf_cv((sele-1)*Mdims%cv_snloc + cv_siloc) == cv_jnod) then
                    sele_from_cv_nod = sele
                    return
                end if
            end do
        end do

    end function sele_from_cv_nod

  END SUBROUTINE MOD_1D_CT_AND_ADV

  !>@brief: This sub modifies either Mmat%CT or the Advection-diffusion equation for 1D pipe modelling
  !> NOTE final_phase has to be define for the reservoir domain, i.e. for two phase flow it can be either 1 or 2, not 3 or 4.
  !> We define wells_first_phase as the first phase of the well domain
  subroutine ASSEMBLE_PIPE_TRANSPORT_AND_CTY( state, packed_state, tracer, den_all, denold_all, final_phase, Mdims, ndgln, DERIV, CV_P, &
                  SOURCT_ALL, ABSORBT_ALL, WIC_T_BC_ALL,WIC_D_BC_ALL, WIC_U_BC_ALL, SUF_T_BC_ALL,SUF_D_BC_ALL,SUF_U_BC_ALL,&
                  getcv_disc, getct, Mmat, Mspars, upwnd, GOT_T2, DT, pipes_aux, DIAG_SCALE_PRES_COUP, DIAG_SCALE_PRES, &
                  mean_pore_cv, eles_with_pipe, thermal, CV_BETA, MASS_CV, INV_B, MASS_ELE, bcs_outfluxes, outfluxes, porous_heat_coef, assemble_collapsed_to_one_phase )
      type(tensor_field), intent(inout) :: tracer
      type(state_type), intent(inout) :: packed_state
      type(state_type), dimension(:), intent(in) :: state
      type(multi_dimensions), intent(in) :: Mdims
      type(multi_ndgln), intent(in) :: ndgln
      type (multi_matrices), intent(inout) :: Mmat
      type (multi_sparsities), intent(in) :: Mspars
      type (porous_adv_coefs), intent(inout) :: upwnd
      type (multi_pipe_package), intent(in) :: pipes_aux
      REAL, DIMENSION( :, : ), intent( in) :: SOURCT_ALL
      REAL, DIMENSION( :, : ), intent( in ) :: DERIV ! (nphase,Mdims%cv_nonods)
      REAL, DIMENSION( :, :, : ), pointer, intent( in ) :: ABSORBT_ALL
      REAL, DIMENSION( :, :, : ), intent( in ) :: CV_P ! (1,Mdims%npres,Mdims%cv_nonods)
      integer, dimension(:,:,:), intent( in ) :: WIC_T_BC_ALL, WIC_D_BC_ALL, WIC_U_BC_ALL
      real, dimension(:,:,:), intent( inout ) :: SUF_T_BC_ALL,SUF_D_BC_ALL,SUF_U_BC_ALL, INV_B, DIAG_SCALE_PRES_COUP
      real, dimension(:,:),intent( inout ) :: den_all, denold_all, DIAG_SCALE_PRES
      REAL, DIMENSION( :, : ), intent( inout ) :: MEAN_PORE_CV ! (Mdims%npres,Mdims%cv_nonods)
      real, dimension(:),intent( inout ) :: MASS_ELE, MASS_CV
      logical, intent( in ) :: getcv_disc, getct, thermal, GOT_T2, assemble_collapsed_to_one_phase
      real, intent(in) :: DT, CV_BETA
      integer, intent(in) :: final_phase
      !Variables that are used to define the pipe pos.
      type(pipe_coords), dimension(:), intent(in):: eles_with_pipe
      REAL, DIMENSION( : ), optional, intent(in) :: porous_heat_coef
      !variables to store the pipe outfluxes if asked by the user
      real, dimension(:,:, :), allocatable, intent(inout):: bcs_outfluxes!<= if allocated then calculate outfluxes
      type (multi_outfluxes), intent(inout) :: outfluxes
      ! Local variables
      real :: auxR, cc, deltaP, rp, rp_nano, skin, h, h_nano, INV_SIGMA_ND, INV_SIGMA_NANO_N, w_sum_one1, w_sum_one2, one_m_cv_beta
      INTEGER :: u_lnloc, ele, i, jpres, u_iloc, x_iloc, idim, cv_lkloc, u_knod, u_lngi, &
          U_NOD, U_SILOC, MAT_NODI, ipres, k, CV_NODI, IPHASE, COUNT, cv_iloc, jphase, iphase_ipres, jphase_jpres,&
          compact_phase, global_phase, wells_first_phase

      !Variables for boundary conditions
      REAL, PARAMETER :: FEM_PIPE_CORRECTION = 0.035
      ! FEM_PIPE_CORRECTION is the FEM pipe correction factor used because the Peacement
      ! model is derived for a 7-point 3D finite difference stencil. This correction factor is obtained
      ! by correlating the IC-FERST production results with Eclipse results on a regular mesh
      ! of linear tetrahedra elements and a single well at steady state. =0.0 is no correction =0.35 recommended.

      type(tensor_field), pointer :: old_tracer
      real, dimension(:,:,:), pointer :: U_ALL
      real, dimension(:,:), pointer :: t2_all, T2OLD_ALL, T_ALL, TOLD_ALL
      real, dimension( :, :, : ), pointer :: fem_p
      real, dimension( Mdims%cv_nonods ) :: N
      !Variables for assemble_collapsed_to_one_phase; Note that diffusion parameters etc need
      integer :: assembly_phase, assembly_phase_2
      !Logical to check if we using a conservative method or not, to save cpu time
      logical :: conservative_advection, have_absorption
      !Parameters of the simulation
      real, parameter :: INFINY=1.0E+20
      integer, parameter :: WIC_B_BC_DIRICHLET = 1
      !
      ! !variables for pipes, allocatable because they are big and barely used
      type( scalar_field ), pointer :: pipe_Diameter, conductivity_pipes, well_thickness
      logical :: has_conductivity_pipes, Solve_all_phases
      REAL, DIMENSION( final_phase * 2 ) ::  RSUM_VEC, LOC_CV_RHS_I, LOC_CV_RHS_J, DIAG_SCALE_PRES_phase, &
                                    MEAN_PORE_CV_PHASE, ct_rhs_phase,r_pres, LOC_MAT_II
      real, dimension( final_phase * 2, final_phase * 2, Mdims%cv_nonods ) :: GAMMA_PRES_ABS2
      REAL , DIMENSION( final_phase * 2, Mdims%cv_nonods) :: opt_vel_upwind_coefs_new_cv
      real, dimension(final_phase * 2):: R_PHASE, Loc_DEN_I, LOC_PRES_I, LOC_SAT_I
      REAL , DIMENSION( final_phase * 2, Mdims%cv_nonods ) :: SIGMA_INV_APPROX!WE NEED FOR THIS SUBROUTINE, I THINK, TO BE DEFINED FOR RESERVOIR AND WELLS DOMAINS...
      real, dimension( Mdims%cv_nonods ) :: MASS_PIPE_FOR_COUP
      real, dimension( final_phase * 2, final_phase * 2, Mdims%cv_nonods ) :: A_GAMMA_PRES_ABS, PIPE_ABS
      !#################SET WORKING VARIABLES#################

      !Define phase where we start the assembly, this is the first phase of the well domains.
      wells_first_phase = 1 + final_phase

      have_absorption = associated( absorbt_all )
      one_m_cv_beta = 1.0 - cv_beta
      Solve_all_phases = .not. have_option("/numerical_methods/solve_nphases_minus_one")

      call get_var_from_packed_state(packed_state, NonlinearVelocity = U_ALL, FEPressure = FEM_P)
      !For every Field_selector value but 3 (saturation) we need U_ALL to be NU_ALL
      T_ALL =>tracer%val(1,:,:)
      old_tracer=>extract_tensor_field(packed_state,GetOldName(tracer))
      TOLD_ALL =>old_tracer%val(1,:,:)
      if (tracer%name == "PackedPhaseVolumeFraction") call get_var_from_packed_state(packed_state,Velocity = U_ALL)
      IF ( GOT_T2 .OR. THERMAL) call get_var_from_packed_state( packed_state, &
      PhaseVolumeFraction = T2_ALL, OldPhaseVolumeFraction = T2OLD_ALL )
      pipe_Diameter => extract_scalar_field( state(1), "DiameterPipe" )

      !For thermal retrieve, if present, the conductivity of the pipes to calculate the heat loss
      has_conductivity_pipes = .false.
      if (thermal .and. is_porous_media) then
          has_conductivity_pipes = have_option('/porous_media/wells_and_pipes/thermal_well_properties')
          if (has_conductivity_pipes) then
              conductivity_pipes => extract_scalar_field( state(1), "Conductivity" )
              well_thickness => extract_scalar_field( state(1), "well_thickness" )
          end if
      end if
    !################## END OF SET VARIABLES ##################
    conservative_advection = abs(cv_beta) > 0.99
    OPT_VEL_UPWIND_COEFS_NEW_CV=0.0 ; N=0.0
    DO ELE = 1, Mdims%totele
        DO CV_ILOC = 1, Mdims%cv_nloc
            CV_NODI = ndgln%cv( CV_ILOC + (ELE-1)*Mdims%cv_nloc )
            !Only go through the nodes that have a well
            if (PIPE_DIAMETER%val(cv_nodi) < 1e-8) cycle
            MAT_NODI = ndgln%mat( CV_ILOC + (ELE-1)*Mdims%cv_nloc )
            if (is_porous_media) then
                RSUM_VEC = 0.0
                do ipres = 1, Mdims%npres
                  do iphase = 1, final_phase
                    DO IDIM = 1, Mdims%ndim
                        RSUM_VEC(iphase + (ipres - 1)*final_phase) = RSUM_VEC(iphase + (ipres - 1)*final_phase) +&
                            upwnd%adv_coef( IDIM, IDIM, iphase + (ipres - 1)*Mdims%n_in_pres, MAT_NODI ) / REAL( Mdims%ndim )
                    END DO
                  end do
              end do
            else
                RSUM_VEC = 1.0
            end if
            OPT_VEL_UPWIND_COEFS_NEW_CV( :, CV_NODI ) = OPT_VEL_UPWIND_COEFS_NEW_CV( :, CV_NODI ) + RSUM_VEC * MASS_ELE( ELE )
            N( CV_NODI ) = N( CV_NODI ) + MASS_ELE( ELE )
        END DO
      END DO

      !Populate SIGMA_INV_APPROX depending on what type of problem we are trying to solve
      if (is_porous_media) then
          DO CV_NODI = 1, Mdims%cv_nonods
              !Only go through the nodes that have a well
              if (PIPE_DIAMETER%val(cv_nodi) < 1e-8) cycle
              SIGMA_INV_APPROX(:, CV_NODI) = 1.0 / ( OPT_VEL_UPWIND_COEFS_NEW_CV(:, CV_NODI) / N(CV_NODI) )
          end do
      else
          SIGMA_INV_APPROX = 1.0
      end if

      MASS_PIPE_FOR_COUP = 0.
      CALL MOD_1D_CT_AND_ADV( state, packed_state, final_phase, wells_first_phase, Mdims, ndgln, WIC_T_BC_ALL,WIC_D_BC_ALL, WIC_U_BC_ALL, SUF_T_BC_ALL,SUF_D_BC_ALL,SUF_U_BC_ALL, &
          getcv_disc, getct, Mmat, Mspars, DT, pipes_aux%MASS_CVFEM2PIPE, pipes_aux%MASS_PIPE2CVFEM, pipes_aux%MASS_CVFEM2PIPE_TRUE, pipes_aux%MASS_PIPE, MASS_PIPE_FOR_COUP, &
          SIGMA_INV_APPROX, upwnd%adv_coef, eles_with_pipe, THERMAL, cv_beta, bcs_outfluxes, outfluxes, assemble_collapsed_to_one_phase)

      GAMMA_PRES_ABS2 = 0.0
      !A_GAMMA_PRES_ABS only for compressible flow? sprint_to_do DO WE NEED TO DO THIS FOR Incompressible FLOW??
      A_GAMMA_PRES_ABS = 0.0
      DO CV_NODI = 1, Mdims%cv_nonods

          !Only go through the nodes that have a well
          if (PIPE_DIAMETER%val(cv_nodi) < 1e-8) cycle

          ! variables used in the edge approach
          h = pipes_aux%MASS_PIPE( cv_nodi )/( pi*(0.5*max(PIPE_DIAMETER%val(cv_nodi), 1.0e-10))**2 )
          h = max( h, 1.0e-10 )
          cc = 0.0
          if ( pipes_aux%MASS_PIPE( cv_nodi )>0.0 ) cc = 1.0 * (1.-FEM_PIPE_CORRECTION) ! to convert Peacement to FEM.
          h = (Mass_CV( cv_nodi )/h)**(1.0/(Mdims%ndim-1)) ! This is the lengthscale normal to the wells.
          rp = 0.14 * h
          Skin = 0.0
          !Create local memory to avoid slicing
          if (GOT_T2) then
            do ipres = 1, Mdims%npres
              DO IPHASE=1,final_phase
                LOC_SAT_I(iphase + (ipres - 1)*final_phase) = T2_ALL( iphase + (ipres - 1)*Mdims%n_in_pres, CV_NODI )
              end do
            end do
          else
            do ipres = 1, Mdims%npres
              DO IPHASE=1,final_phase
                LOC_SAT_I(iphase + (ipres - 1)*final_phase) = T_ALL( iphase + (ipres - 1)*Mdims%n_in_pres, CV_NODI )
              end do
            end do
          end if

          do ipres = 1, Mdims%npres
            DO IPHASE=1,final_phase
                Loc_DEN_I(iphase + (ipres - 1)*final_phase) =  DEN_ALL( iphase + (ipres - 1)*Mdims%n_in_pres, CV_NODI )
            end do
          end do
          DO IPRES = 1, Mdims%npres
             LOC_PRES_I(1+(ipres-1)*final_phase:ipres*final_phase) = FEM_P( 1, IPRES, CV_NODI )
          END DO


          DO IPHASE = 1, final_phase*2
              IPRES = 1 + INT( (IPHASE-1)/final_phase )
              DO JPHASE = 1, final_phase*2
                  JPRES = 1 + INT( (JPHASE-1)/final_phase )
                  ! This is the edge approach
                  ! We do NOT divide by r**2 here because we have not multiplied by r**2 in the pipes_aux%MASS_CVFEM2PIPE matrix (in MOD_1D_CT_AND_ADV)
                  IF ( IPRES /= JPRES ) THEN
                      !Peaceman correction
                      IF ( LOC_PRES_I(IPHASE) > LOC_PRES_I(JPHASE) ) THEN
                          GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI ) = pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * &
                              cc * LOC_SAT_I(IPHASE) * 2.0 * SIGMA_INV_APPROX( IPHASE, CV_NODI ) &
                              / ( 1.0*(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                      ELSE
                          GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI ) = pipes_aux%GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) * &
                              cc * LOC_SAT_I(JPHASE) * 2.0 * SIGMA_INV_APPROX( JPHASE, CV_NODI ) &
                              / ( 1.0*(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                      END IF
                  END IF ! IF ( IPRES /= JPRES ) THEN
              END DO
          END DO

         DO IPHASE = 1, final_phase*2
             DO JPHASE = 1, final_phase*2
                 IPRES = 1 + INT( (IPHASE-1)/final_phase )
                 JPRES = 1 + INT( (JPHASE-1)/final_phase )
                 IF( IPRES /= JPRES ) THEN
                     A_GAMMA_PRES_ABS( IPHASE, JPHASE, CV_NODI ) = - GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI )
                     A_GAMMA_PRES_ABS( IPHASE, IPHASE, CV_NODI ) = A_GAMMA_PRES_ABS( IPHASE, IPHASE, CV_NODI ) + GAMMA_PRES_ABS2( IPHASE, JPHASE, CV_NODI )
                 END IF
             END DO
         END DO
      END DO  ! DO CV_NODI = 1, Mdims%cv_nonods


      if(GETCV_DISC) then
          PIPE_ABS = 0.0
          DO CV_NODI = 1, Mdims%cv_nonods
              !Only go through the nodes that have a well
              if (PIPE_DIAMETER%val(cv_nodi) < 1e-8) cycle

              h = pipes_aux%MASS_PIPE( cv_nodi )/( pi*(0.5*max(PIPE_DIAMETER%val(cv_nodi),1.0e-10))**2)
              h = max( h, 1.0e-10 )
              h_nano = h
              cc = 0.0
              if ( pipes_aux%MASS_PIPE( cv_nodi )>0.0 ) cc = 1.0 * (1.-FEM_PIPE_CORRECTION) ! to convert Peacement to FEM.
              h = (Mass_CV( cv_nodi )/h)**(1.0/(Mdims%ndim-1))  ! This is the lengthscale normal to the wells.
              rp = 0.14 * h
              rp_NANO = 0.14 * h_NANO
              Skin = 0.0
              !Create local memory to reduce slicing
              do ipres = 1, Mdims%npres
                DO IPHASE=1,final_phase
                    Loc_DEN_I(iphase + (ipres - 1)*final_phase) =  DEN_ALL( iphase + (ipres - 1)*Mdims%n_in_pres, CV_NODI )
                end do
              end do
              DO IPRES = 1, Mdims%npres
                 LOC_PRES_I(1+(ipres-1)*final_phase:ipres*final_phase) = FEM_P( 1, IPRES, CV_NODI )
              END DO

              if (.not. conservative_advection) then
                !In the non-conservative method the concept is to only exchange the difference of the mass/Energy
                !instead of one side sending X and reciving Y and the other sending Y and receiving X
                do ipres = 1, Mdims%npres
                  do jpres = 1, Mdims%npres
                    IF ( IPRES /= JPRES ) THEN
                      DO iphase=1, final_phase
                        iphase_ipres = iphase + (ipres - 1)*final_phase
                        DO jphase=1, final_phase
                            jphase_jpres = jphase + (jpres - 1)*final_phase
                            IF ( LOC_PRES_I(iphase_ipres) < LOC_PRES_I(jphase_jpres)  ) THEN
                                DeltaP = LOC_PRES_I(iphase_ipres) - LOC_PRES_I(jphase_jpres)
                                !Add coupling in a non_conservative manner (mass conservative anyway is just the formulation)
                                PIPE_ABS( iphase_ipres, jphase_jpres, CV_NODI ) = PIPE_ABS( iphase_ipres, jphase_jpres, CV_NODI ) +&
                                    DeltaP * pipes_aux%GAMMA_PRES_ABS( iphase + (ipres - 1)*Mdims%n_in_pres, jphase + (jpres - 1)*Mdims%n_in_pres, CV_NODI ) * Loc_DEN_I( jphase_jpres ) * &
                                    cc * 2.0 * SIGMA_INV_APPROX( jphase_jpres, CV_NODI ) &
                                    / ( 1.0 *(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                               PIPE_ABS( iphase_ipres, iphase_ipres, CV_NODI ) = PIPE_ABS( iphase_ipres, iphase_ipres, CV_NODI ) - &!notice the negative sign
                                   DeltaP * pipes_aux%GAMMA_PRES_ABS( jphase + (jpres - 1)*Mdims%n_in_pres, iphase + (ipres - 1)*Mdims%n_in_pres, CV_NODI ) * Loc_DEN_I( jphase_jpres ) * &
                                   cc * 2.0 * SIGMA_INV_APPROX( jphase_jpres, CV_NODI ) &
                                   / ( 1.0 *(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                                 END IF
                           END DO
                         END DO
                       END IF
                     end do
                   end do
              else
                  do ipres = 1, Mdims%npres
                    do jpres = 1, Mdims%npres
                      IF ( IPRES /= JPRES ) THEN
                        DO iphase=1, final_phase
                          iphase_ipres = iphase + (ipres - 1)*final_phase
                          DO jphase=1, final_phase
                              jphase_jpres = jphase + (jpres - 1)*final_phase
                                DeltaP = LOC_PRES_I(iphase_ipres) - LOC_PRES_I(jphase + (jpres - 1)*final_phase)
                                ! We do NOT divide by r**2 here because we have not multiplied by r**2 in the pipes_aux%MASS_CVFEM2PIPE matrix (in MOD_1D_CT_AND_ADV)
                                IF ( LOC_PRES_I(iphase_ipres) >= LOC_PRES_I(jphase_jpres) ) THEN
                                    PIPE_ABS( iphase_ipres, iphase_ipres, CV_NODI ) = PIPE_ABS( iphase_ipres, iphase_ipres, CV_NODI ) + &
                                        DeltaP * pipes_aux%GAMMA_PRES_ABS( iphase + (ipres - 1)*Mdims%n_in_pres, jphase + (jpres - 1)*Mdims%n_in_pres, CV_NODI ) * Loc_DEN_I( iphase_ipres ) * &
                                        cc * 2.0 * SIGMA_INV_APPROX( iphase_ipres, CV_NODI ) &
                                        / (log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin)
                                ELSE
                                    PIPE_ABS( iphase_ipres, jphase_jpres, CV_NODI ) = PIPE_ABS( iphase_ipres, jphase_jpres, CV_NODI ) +&
                                        DeltaP * pipes_aux%GAMMA_PRES_ABS( iphase + (ipres - 1)*Mdims%n_in_pres, jphase + (jpres - 1)*Mdims%n_in_pres, CV_NODI ) * Loc_DEN_I( jphase_jpres ) * &
                                        cc * 2.0 * SIGMA_INV_APPROX( jphase_jpres, CV_NODI ) &
                                        / ( 1.0 *(log( rp / max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) ) + Skin) )
                                END IF
                          END DO
                        END DO
                      END IF
                    end do
                  end do
              end if

              !add the heat loss contribution due to diffusion to the nodes with pipes defined, even if it is closed
              !this might be true only for thermal and porous media
              if (has_conductivity_pipes) then
                  !Apply only where wells are closed, this is a good approximation
                  !Gamma should be the same for at least the well phases, so we check nphase
                  if (pipes_aux%GAMMA_PRES_ABS( Mdims%nphase, Mdims%nphase, CV_NODI )<1d-8) then
                      count = min(1,cv_nodi)
                      !Rp is the internal radius of the well
                      rp = max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) - well_thickness%val(count)
                      h = 1./(PI * (0.5*pipe_Diameter%val( cv_nodi )**2.))!height/volume pipe
                      !we apply Q = 1/WellVolume * (Tin-Tout) * 2 * PI * K * L/(ln(Rout/Rin))
                      auxR = conductivity_pipes%val(count) * 2.0 * PI * h / log( max( 0.5*pipe_Diameter%val( cv_nodi ), 1.0e-10 ) / rp )
                      DO IPHASE = 1, final_phase
                          jphase = iphase+final_phase
                          PIPE_ABS( IPHASE, IPHASE, CV_NODI ) = PIPE_ABS( IPHASE, IPHASE, CV_NODI ) + auxR
                          PIPE_ABS( iphase, jphase, CV_NODI ) = PIPE_ABS( IPHASE, JPHASE, CV_NODI ) - auxR
                      end do
                      !Need to apply variation of heat from the other domain as well
                      DO IPHASE = wells_first_phase, final_phase*2
                          jphase = iphase-final_phase
                          PIPE_ABS( IPHASE, IPHASE, CV_NODI ) = PIPE_ABS( IPHASE, IPHASE, CV_NODI ) + auxR
                          PIPE_ABS( iphase, jphase, CV_NODI ) = PIPE_ABS( IPHASE, JPHASE, CV_NODI ) - auxR
                      end do
                  end if
              end if
          END DO ! DO CV_NODI = 1, Mdims%cv_nonods
      endif ! if(GETCV_DISC) then
      IF ( GETCT ) THEN
          !SPRINT_TO_DO IS INV_B USED, it is used to calculate the CMC matrix ...
          !INV_B IS COMPLETELY REDUNDANT AND SHOULD BE REMOVED
            DO IPHASE = 1, Mdims%nphase
                INV_B( IPHASE, IPHASE, : ) = 1./DEN_ALL( IPHASE, : )
            END DO
      !Introduce terms into the matrix and the RHS
      end if
      Conditional_GETCV_DISC2: IF( GETCV_DISC ) THEN ! Obtain the CV discretised advection/diffusion equations
          Loop_CVNODI2: DO CV_NODI = 1, Mdims%cv_nonods ! Put onto the diagonal of the matrix
              !Only go through the nodes that have a well
              if (PIPE_DIAMETER%val(cv_nodi) < 1e-8) cycle

              LOC_CV_RHS_I=0.0; LOC_MAT_II = 0.
              R_PHASE(wells_first_phase:final_phase*2) = MEAN_PORE_CV( Mdims%npres, CV_NODI ) * pipes_aux%MASS_PIPE( cv_nodi ) / DT

              IF ( THERMAL .and. Mdims%npres == 1) THEN !gr never enters this if statement
                DO IPHASE = wells_first_phase, final_phase*2
                  LOC_CV_RHS_I(iphase)=LOC_CV_RHS_I(iphase) &
                      - CV_P(1, Mdims%npres, cv_nodi) * ( Mass_CV( CV_NODI ) / DT ) * ( T2_ALL( iphase, CV_NODI ) - T2OLD_ALL( iphase, CV_NODI ) )
                end do
              END IF

              IF ( GOT_T2 ) THEN
                  DO IPHASE = wells_first_phase, final_phase*2
                    LOC_CV_RHS_I(iphase)=LOC_CV_RHS_I(iphase)  &
                        + Mass_CV(CV_NODI) * SOURCT_ALL( iphase, CV_NODI )
                  end do
                  if (thermal .and. is_porous_media) then !gr  first term goes to 0?
                      !In this case for the time-integration term the effective rho Cp is a combination of the porous media
                      ! and the fluids. Here we add the porous media contribution
                      DO IPHASE = wells_first_phase, final_phase*2
                        LOC_MAT_II(iphase) = LOC_MAT_II(iphase) +  porous_heat_coef( CV_NODI ) * T2_ALL( IPHASE, CV_NODI ) &
                        * R_PHASE(IPHASE) * (1-MEAN_PORE_CV( 1, CV_NODI ))/MEAN_PORE_CV( 1, CV_NODI )
                              !R_PHASE includes the porosity. Since in this case we are interested in what is NOT porous
                                  !we divide to remove that term and multiply by the correct term (1-porosity)
                          LOC_CV_RHS_I(iphase)=LOC_CV_RHS_I(iphase)  &
                              + (CV_BETA * porous_heat_coef( CV_NODI ) * T2OLD_ALL( iphase, CV_NODI ) &
                              + (ONE_M_CV_BETA) * porous_heat_coef( CV_NODI ) * T2_ALL( iphase, CV_NODI ) ) &
                              * R_PHASE(iphase) * TOLD_ALL( iphase, CV_NODI )* (1-MEAN_PORE_CV( 1, CV_NODI ))/MEAN_PORE_CV( 1, CV_NODI )
                      END DO
                  end if

                  DO IPHASE=1 , final_phase !gr time derivative
                    global_phase = iphase + Mdims%n_in_pres
                    compact_phase = iphase + final_phase
                    LOC_MAT_II(compact_phase) = LOC_MAT_II(compact_phase) + DEN_ALL( global_phase, CV_NODI ) * T2_ALL( global_phase, CV_NODI ) &
                    * R_PHASE(compact_phase)
                  END DO
                 DO IPHASE = 1 , final_phase
                   global_phase = iphase + Mdims%n_in_pres
                   compact_phase = iphase + final_phase
                   LOC_CV_RHS_I(compact_phase)=LOC_CV_RHS_I(compact_phase)  &
                        + (CV_BETA * DENOLD_ALL( global_phase, CV_NODI ) * T2OLD_ALL( global_phase, CV_NODI ) &
                        + (ONE_M_CV_BETA) * DEN_ALL( global_phase, CV_NODI ) * T2_ALL( global_phase, CV_NODI ) ) &
                        * R_PHASE(compact_phase) * TOLD_ALL( global_phase, CV_NODI )
                 end do
              ELSE
                DO IPHASE=1 , final_phase
                  global_phase = iphase + Mdims%n_in_pres
                  compact_phase = iphase + final_phase
                  LOC_MAT_II(compact_phase) = LOC_MAT_II(compact_phase) + DEN_ALL( global_phase, CV_NODI )  &
                      * R_PHASE(compact_phase)

                  LOC_CV_RHS_I(compact_phase)=LOC_CV_RHS_I(compact_phase)  &
                    + Mass_CV( CV_NODI ) * SOURCT_ALL( global_phase, CV_NODI )&
                    + ( CV_BETA * DENOLD_ALL( global_phase, CV_NODI ) &
                    + (ONE_M_CV_BETA) * DEN_ALL( global_phase, CV_NODI ) ) &
                    * R_PHASE(compact_phase) * TOLD_ALL( global_phase, CV_NODI )
                  END DO
              END IF

              !Here we add the coupling term between phases for the wells, so we have to iterate over all the phases
              do ipres = 1, Mdims%npres
                DO jphase=1 , final_phase
                  ! global_phase = jphase + (ipres - 1)*Mdims%n_in_pres
                  compact_phase = jphase + (ipres - 1)*final_phase
                  assembly_phase = compact_phase
                  if (assemble_collapsed_to_one_phase) assembly_phase = ipres!1 + (jphase-1)/Mdims%n_in_pres
                    do iphase=1 , final_phase*2
                      !Implicit method in all the cases
                      assembly_phase_2 = iphase
                      if (assemble_collapsed_to_one_phase) assembly_phase_2 = 1 + (iphase-1)/Mdims%n_in_pres
                      call addto(Mmat%petsc_ACV,assembly_phase_2,assembly_phase, &
                          cv_nodi, cv_nodi, &
                          MASS_PIPE_FOR_COUP( CV_NODI ) * PIPE_ABS( iphase, compact_phase, CV_NODI ))
                          ! if(.not.conservative_advection) then ! original method - all implicit (may be unstable in some cases 12/07/2017)
                          !     call addto(Mmat%petsc_ACV,iphase,compact_phase, &
                          !         cv_nodi, cv_nodi, &
                          !         MASS_PIPE_FOR_COUP( CV_NODI ) * PIPE_ABS( iphase, compact_phase, CV_NODI ))
                          ! else
                          !   !Some sort of predictor corrector. Implicit may be unstable, therefore the explicit acts as a Dumping term
                          !   call addto(Mmat%petsc_ACV,iphase,compact_phase, &
                          !       cv_nodi, cv_nodi, &
                          !         (1. + dumping_well_factor) * MASS_PIPE_FOR_COUP( CV_NODI ) * PIPE_ABS( iphase, compact_phase, CV_NODI ))
                          !
                          !       ! ! well mass exchange is introduced in the RHS so it is treated explictly
                          !       LOC_CV_RHS_I(IPHASE)=LOC_CV_RHS_I(IPHASE)  &
                          !       + dumping_well_factor * MASS_PIPE_FOR_COUP( CV_NODI ) * PIPE_ABS( iphase, compact_phase, CV_NODI )* T_ALL(global_phase, CV_NODI)
                          !
                          !     ! ! ! well mass exchange is introduced in the RHS so it is treated explictly <+ all explicit is too slow
                          !     ! LOC_CV_RHS_I(IPHASE)=LOC_CV_RHS_I(IPHASE)  &
                          !     ! - MASS_PIPE_FOR_COUP( CV_NODI ) * PIPE_ABS( iphase, jphase, CV_NODI )* T_ALL(global_phase, CV_NODI)
                          ! endif
                    end do
                end do
              end do

              if ( have_absorption ) then
                DO jphase=1 , final_phase
                  compact_phase = jphase + final_phase
                  assembly_phase = compact_phase
                  if (assemble_collapsed_to_one_phase) assembly_phase = 1 + (jphase-1)/Mdims%n_in_pres
                    do iphase=wells_first_phase, final_phase*2
                      assembly_phase_2 = iphase
                      if (assemble_collapsed_to_one_phase) assembly_phase_2 = 1 + (iphase-1)/Mdims%n_in_pres
                       call addto(Mmat%petsc_ACV,assembly_phase_2,assembly_phase, cv_nodi, cv_nodi, &
                           Mass_CV( CV_NODI ) * ABSORBT_ALL( iphase, compact_phase, CV_NODI ))
                    end do
                end do
            end if

            !Assemble into the matrix/RHS
            DO IPHASE= 1, final_phase*2
                assembly_phase = iphase
                !For the RHS collapsing to assemble into phase 1 can be done just here
                if (assemble_collapsed_to_one_phase) assembly_phase = 2
                call addto(Mmat%CV_RHS,assembly_phase, CV_NODI,LOC_CV_RHS_I(IPHASE))
                !Introduce the information into the petsc_ACV matrix
                call addto(Mmat%petsc_ACV,assembly_phase,assembly_phase,cv_nodi,cv_nodi, LOC_MAT_II(iphase) )
            end do

            ! call addto(Mmat%CV_RHS,CV_NODI,LOC_CV_RHS_I)
          END DO Loop_CVNODI2
      END IF Conditional_GETCV_DISC2

      IF ( GETCT ) THEN

          W_SUM_ONE1 = 0.0 !If == 1.0 applies constraint to T
          if (Solve_all_phases) W_SUM_ONE1 = 1.0
          W_SUM_ONE2 = 0.0 !If == 1.0 applies constraint to TOLD !sprint_to_do Unnecessary, should be removed
          DIAG_SCALE_PRES = 0.0
          DIAG_SCALE_PRES_COUP=0.0
          DO CV_NODI = 1, Mdims%cv_nonods
            !Only go through the nodes that have a well
            if (PIPE_DIAMETER%val(cv_nodi) < 1e-8) cycle
              ct_rhs_phase=0.0 ; DIAG_SCALE_PRES_phase=0.0
              IPRES= Mdims%npres
              R_PRES(IPRES) = pipes_aux%MASS_PIPE( cv_nodi ) * MEAN_PORE_CV( IPRES, CV_NODI ) / DT
              ! Add constraint to force sum of volume fracts to be unity...
                 ! W_SUM_ONE==1 applies the constraint
                 ! W_SUM_ONE==0 does NOT apply the constraint
              if ( Solve_all_phases) call addto(Mmat%CT_RHS,IPRES,cv_nodi,&
                  - ( W_SUM_ONE1 - W_SUM_ONE2 ) * R_PRES(IPRES))
              R_PHASE(1+(ipres-1)*final_phase:ipres*final_phase)=R_PRES(IPRES)!SPRINT_TO_DO WHY COPIYING MEMORY HERE???
              MEAN_PORE_CV_PHASE(1+(ipres-1)*final_phase:ipres*final_phase) = MEAN_PORE_CV( IPRES, CV_NODI )!SPRINT_TO_DO WHY COPIYING MEMORY HERE???
              do iphase = wells_first_phase, final_phase*2
                ct_rhs_phase(iphase) = ct_rhs_phase(iphase) &
                    - R_PHASE(iphase) * ( &
                    + (1.0-W_SUM_ONE1) * T_ALL( iphase, CV_NODI ) - (1.0-W_SUM_ONE2) * TOLD_ALL( iphase, CV_NODI ) &
                    + ( TOLD_ALL( iphase, CV_NODI ) * ( DEN_ALL( iphase, CV_NODI ) - DENOLD_ALL( iphase, CV_NODI ) ) &
                    - DERIV( iphase, CV_NODI ) * CV_P(1, Mdims%npres, cv_nodi ) * T_ALL( iphase, CV_NODI ) ) / DEN_ALL(iphase, CV_NODI ) )
                DIAG_SCALE_PRES_phase( iphase ) = DIAG_SCALE_PRES_phase( iphase ) &
                    +  MEAN_PORE_CV_PHASE(iphase) * T_ALL( iphase, CV_NODI ) * DERIV( iphase, CV_NODI ) &
                    / ( DT * DEN_ALL(iphase, CV_NODI) )
                ct_rhs_phase(iphase)=ct_rhs_phase(iphase) &
                    + Mass_CV(CV_NODI ) * SOURCT_ALL( iphase, CV_NODI ) / DEN_ALL(iphase, CV_NODI)
                IF ( HAVE_ABSORPTION ) THEN!No absorption in the wells for the time being
                   DO JPHASE = wells_first_phase, final_phase*2
                      ct_rhs_phase(iphase)=ct_rhs_phase(iphase)  &
                         - Mass_CV( CV_NODI ) * ABSORBT_ALL( iphase, JPHASE, CV_NODI ) * T_ALL( JPHASE, CV_NODI ) / DEN_ALL(iphase, CV_NODI)
                   END DO
                END IF
            end do
            ! scaling coefficient...
            DO IPRES=1,Mdims%npres
                DO JPRES=1,Mdims%npres
                    DO jphase=1+(jpres-1)*final_phase, jpres*final_phase
                        ! dont divid the pipe to reservoir mass exchange term by density...
                        DIAG_SCALE_PRES_COUP(IPRES,JPRES, cv_nodi) = DIAG_SCALE_PRES_COUP(IPRES,JPRES, cv_nodi)  &
                            + sum( A_GAMMA_PRES_ABS(1+(ipres-1)*final_phase:ipres*final_phase,JPHASE, CV_NODI )  )
                    end do
                END DO
            END DO
            IPRES= Mdims%npres
            !Introduce into RHS and scaling coefficient
            call addto(Mmat%CT_RHS, IPRES, cv_nodi, SUM( ct_rhs_phase(1+(ipres-1)*final_phase:ipres*final_phase)))
            DIAG_SCALE_PRES( IPRES,CV_NODI ) = DIAG_SCALE_PRES( IPRES,CV_NODI ) &
                + sum( DIAG_SCALE_PRES_phase(1+(ipres-1)*final_phase:ipres*final_phase))
          END DO  ! endof DO CV_NODI = 1, Mdims%cv_nonods
      ENDIF ! ENDOF IF ( GETCT ) THEN

  end subroutine ASSEMBLE_PIPE_TRANSPORT_AND_CTY


  !>@brief: This sub modifies Mmat%C for 1D pipe modelling
    SUBROUTINE MOD_1D_FORCE_BAL_C( STATE, packed_state, Mdims, Mspars, Mmat, ndgln, eles_with_pipe, GET_PIVIT_MAT, &
        &                         WIC_P_BC_ALL,SUF_P_BC_ALL, SIGMA, NU_ALL, &
        &                         U_SOURCE, U_SOURCE_CV )
        IMPLICIT NONE
        TYPE(STATE_TYPE), DIMENSION( : ), INTENT( IN ) :: STATE
        TYPE(STATE_TYPE), INTENT( IN ) :: packed_STATE
        type(multi_dimensions), intent(in) :: Mdims
        type (multi_sparsities), intent(in) :: Mspars
        type(multi_ndgln), intent(in) :: ndgln
        type (multi_matrices), intent(inout) :: Mmat
        !REAL, DIMENSION( :, :, : ), INTENT( IN ) :: SUF_P_BC_ALL, U_SOURCE, U_SOURCE_CV
        REAL, DIMENSION( :, :, : ), INTENT( IN ) :: SUF_P_BC_ALL, U_SOURCE_CV
        type( multi_field ), INTENT( IN ) :: U_SOURCE
        REAL, DIMENSION( :, : ), INTENT( IN ) :: SIGMA
        REAL, DIMENSION( :, :, : ), INTENT( IN ) :: NU_ALL
        INTEGER, DIMENSION( :,:,: ), INTENT( IN ) :: WIC_P_BC_ALL
        LOGICAL, INTENT( IN ) :: GET_PIVIT_MAT
        !Variables that are used to define the pipe pos.
        type(pipe_coords), dimension(:), intent(in):: eles_with_pipe

        !Local variables
        LOGICAL :: CV_QUADRATIC, U_QUADRATIC, ELE_HAS_PIPE, PIPE_MIN_DIAM, U_P0DG
        LOGICAL :: CALC_SIGMA_PIPE
        INTEGER :: ELE, PIPE_NOD_COUNT, ICORNER, &
            &     CV_ILOC, U_ILOC, CV_NODI, IPIPE, CV_LILOC, U_LILOC, CV_LNLOC, U_LNLOC, CV_KNOD, MAT_KNOD, IDIM, &
            &     IU_NOD, P_LJLOC, JCV_NOD, COUNT, COUNT2, IPHASE
        integer, dimension(Mdims%ndim + 1) :: CV_LOC_CORNER, U_LOC_CORNER
        INTEGER, DIMENSION(:), ALLOCATABLE :: CV_GL_LOC, CV_GL_GL, X_GL_GL, MAT_GL_GL, U_GL_LOC, U_GL_GL
        INTEGER, DIMENSION(Mdims%npres, Mdims%cv_nonods) :: WIC_P_BC_ALL_NODS
        INTEGER, DIMENSION(Mdims%ndim + 1, Mdims%ndim + 1) :: CV_MID_SIDE, U_MID_SIDE
        TYPE(SCALAR_FIELD), POINTER :: PIPE_DIAMETER, WD
        TYPE(VECTOR_FIELD), POINTER :: X
        TYPE(TENSOR_FIELD), POINTER :: WM, CV_VOL_FRAC
        REAL, DIMENSION( :, : ), ALLOCATABLE :: scvfen, scvfenslx, scvfensly, &
            &                                  scvfenlx, scvfenly, scvfenlz, &
            &                                  sufen, sufenslx, sufensly, &
            &                                  sufenlx, sufenly, sufenlz
        REAL, DIMENSION( : ), ALLOCATABLE :: scvfeweigh
        REAL, DIMENSION( :, :, : ), ALLOCATABLE :: L_CVFENX_ALL_REVERSED
        REAL, DIMENSION( :, : ), ALLOCATABLE :: L_CVFENX_ALL, L_UFENX_ALL, L_UFEN_REVERSED
        REAL, DIMENSION( :, : ), ALLOCATABLE :: X_ALL_CORN
        REAL, DIMENSION(Mdims%ndim, Mdims%nphase) :: LOC_U_RHS_U_ILOC
        REAL, DIMENSION(Mdims%npres, Mdims%cv_nonods) ::SUF_P_BC_ALL_NODS, RVEC_SUM
        REAL, DIMENSION( : ), ALLOCATABLE :: DETWEI, PIPE_DIAM_GI, NMX_ALL, WELL_DENSITY, WELL_VISCOSITY
        REAL, DIMENSION( :, : ), ALLOCATABLE :: SIGMA_GI, SIGMA_ON_OFF_GI
        TYPE( SCALAR_FIELD ), POINTER :: PHASE_EXCLUDE_PIPE_SAT_MIN, PHASE_EXCLUDE_PIPE_SAT_MAX, SIGMA_SWITCH_ON_OFF_PIPE
        INTEGER :: PHASE_EXCLUDE
        LOGICAL, DIMENSION( Mdims%ndim + 1 ) :: PIPE_INDEX_LOGICAL
        REAL, dimension(Mdims%ndim) :: DIRECTION, DIRECTION_NORM
        REAL :: DX, ELE_ANGLE, NN, NM, suf_area, PIPE_DIAM_END, MIN_DIAM, U_GI, E_ROUGHNESS
        REAL :: S_WATER, S_WATER_MIN, S_WATER_MAX, SIGMA_SWITCH_ON_OFF_PIPE_GI, PIPE_SWITCH
        INTEGER :: ncorner, scvngi, k, &
            &     i_indx, j_indx, jdim, jphase, u_ljloc, u_jloc, ICORNER1, ICORNER2, ICORNER3, ICORNER4
        INTEGER :: SELE, CV_SILOC, JCV_NOD1, JCV_NOD2, IPRES, JU_NOD, CV_NOD, CV_LOC1, CV_LOC2, IPHASE_IN_PIPE, GI, IWATER

        ncorner = Mdims%ndim + 1
        PIPE_MIN_DIAM=.TRUE. ! Take the min diamter of the pipe as the real diameter.
        CALC_SIGMA_PIPE = have_option("/porous_media/wells_and_pipes/well_options/calculate_sigma_pipe")
        call get_option("/porous_media/wells_and_pipes/well_options/calculate_sigma_pipe/pipe_roughness", E_ROUGHNESS, default=1.0E-6)
        ! Add the sigma associated with the switch to switch the pipe flow on and off...
        if ( CALC_SIGMA_PIPE ) then
            allocate( well_density(Mdims%nphase), well_viscosity(Mdims%nphase) )
            do iphase = Mdims%n_in_pres+1, Mdims%nphase
                wd => extract_scalar_field( state(iphase), "Density" )
                wm => extract_tensor_field( state(iphase), "Viscosity" )
                well_density( iphase ) = wd%val(1)
                well_viscosity( iphase ) = wm%val(1,1,1)
            end do
        end if
        ! Set rhs of the force balce equation to zero just for the pipes...
        Mmat%U_RHS( :, Mdims%n_in_pres+1:Mdims%nphase, : ) = 0.0
        IF ( GET_PIVIT_MAT ) THEN
            DO U_ILOC = 1, Mdims%u_nloc
                DO U_JLOC = 1, Mdims%u_nloc
                    DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                        JPHASE = IPHASE
                        DO IDIM = 1, Mdims%ndim
                            JDIM = IDIM
                            i_indx = IDIM + (IPHASE-1)*Mdims%ndim + (U_ILOC-1)*Mdims%ndim*Mdims%nphase
                            j_indx = JDIM + (JPHASE-1)*Mdims%ndim + (U_JLOC-1)*Mdims%ndim*Mdims%nphase
                            Mmat%PIVIT_MAT(i_indx, j_indx, :) = 0.0
                        END DO
                    END DO
                END DO
            END DO
        END IF
        CV_QUADRATIC = (Mdims%cv_nloc==6.AND.Mdims%ndim==2) .OR. (Mdims%cv_nloc==10.AND.Mdims%ndim==3)
        U_QUADRATIC = (Mdims%u_nloc==6.AND.Mdims%ndim==2) .OR. (Mdims%u_nloc==10.AND.Mdims%ndim==3)
        U_P0DG = Mdims%u_nloc == 1
        ! The following is for 1d integration...
        scvngi = 2
        if ( CV_QUADRATIC ) then
            cv_lnloc = 3
            scvngi = 3
        else
            cv_lnloc = 2
        end if
        if ( U_QUADRATIC ) then
            u_lnloc = 3
            scvngi = 3
        else
            u_lnloc = 2
            if (U_P0DG) u_lnloc = 1
        end if

        PIPE_Diameter => EXTRACT_SCALAR_FIELD(state(1), "DiameterPipe")
        X => EXTRACT_VECTOR_FIELD( PACKED_STATE, "PressureCoordinate" )
        allocate( PIPE_DIAM_GI(scvngi) )
        allocate( SIGMA_GI(Mdims%nphase,scvngi), SIGMA_ON_OFF_GI(Mdims%nphase,scvngi) )
        ! Get the 1D shape functions...
        allocate( scvfeweigh(scvngi), &
            scvfen(cv_lnloc, scvngi), scvfenslx(cv_lnloc, scvngi), scvfensly(cv_lnloc, scvngi), &
            scvfenlx(cv_lnloc, scvngi), scvfenly(cv_lnloc, scvngi), scvfenlz(cv_lnloc, scvngi), &
            sufen(u_lnloc, scvngi), sufenslx(u_lnloc, scvngi), sufensly(u_lnloc, scvngi), &
            sufenlx(u_lnloc, scvngi), sufenly(u_lnloc, scvngi), sufenlz(u_lnloc, scvngi) )
        allocate( detwei(scvngi), &
            l_cvfenx_all(cv_lnloc, scvngi), l_ufenx_all(u_lnloc, scvngi), &
            l_cvfenx_all_reversed(Mdims%ndim, scvngi, cv_lnloc) , &
            l_ufen_reversed(scvngi, u_lnloc), &
            nmx_all( Mdims%ndim ), X_ALL_CORN(Mdims%ndim, ncorner) )
        call fv_1d_quad( scvngi, cv_lnloc, scvfen, scvfenslx, scvfensly, scvfeweigh, &
            scvfenlx, scvfenly, scvfenlz ) ! For scalar fields
        call fv_1d_quad( scvngi, u_lnloc, sufen, sufenslx, sufensly, scvfeweigh, &
            sufenlx, sufenly, sufenlz ) ! For U fields

        ! Calculate the local corner nodes...
        CALL CALC_CORNER_NODS(CV_LOC_CORNER, Mdims%ndim, Mdims%cv_nloc, CV_QUADRATIC, CV_MID_SIDE)
        CALL CALC_CORNER_NODS(U_LOC_CORNER, Mdims%ndim, Mdims%u_nloc, U_QUADRATIC, U_MID_SIDE)

        allocate( CV_GL_LOC(cv_lnloc), CV_GL_GL(cv_lnloc), X_GL_GL(cv_lnloc), MAT_GL_GL(cv_lnloc) )
        allocate( U_GL_LOC(U_lnloc), U_GL_GL(U_lnloc) )
        ! set up the surface B.Mmat%C.'s
        WIC_P_BC_ALL_NODS=0 ; RVEC_SUM=0.0 ; SUF_P_BC_ALL_NODS=0.0
        DO SELE = 1, Mdims%stotel
            DO IPRES = 2, Mdims%npres
                IF ( WIC_P_BC_ALL( 1, IPRES, SELE ) == WIC_P_BC_DIRICHLET ) THEN
                    DO CV_SILOC = 1, Mdims%cv_snloc
                        CV_NOD = ndgln%suf_p((SELE-1)*Mdims%cv_snloc + CV_SILOC )
                        WIC_P_BC_ALL_NODS( IPRES, CV_NOD ) = WIC_P_BC_DIRICHLET
                        SUF_P_BC_ALL_NODS( IPRES, CV_NOD ) = SUF_P_BC_ALL_NODS( IPRES, CV_NOD ) + SUF_P_BC_ALL( 1,IPRES, (SELE-1)*Mdims%cv_snloc + CV_SILOC )
                        RVEC_SUM( IPRES, CV_NOD ) = RVEC_SUM( IPRES, CV_NOD ) + 1.0
                    END DO
                END IF
            END DO
        END DO
        DO CV_NOD = 1, Mdims%cv_nonods
            DO IPRES = 2, Mdims%npres
                IF ( WIC_P_BC_ALL_NODS(IPRES, CV_NOD) /= 0 ) SUF_P_BC_ALL_NODS(IPRES, CV_NOD) = SUF_P_BC_ALL_NODS(IPRES, CV_NOD) / RVEC_SUM(IPRES, CV_NOD)
            END DO
        END DO
        Mmat%C( :, Mdims%n_in_pres+1:Mdims%nphase, : ) = 0.0

        DO k = 1, size(eles_with_pipe)
            ELE = eles_with_pipe(k)%ele!Element with pipe

            X_ALL_CORN(:,1:NCORNER) = x%val(:, ndgln%x( ( ELE - 1 ) * Mdims%cv_nloc + CV_LOC_CORNER(1:NCORNER)) )
            DO IPIPE = 1, eles_with_pipe(k)%npipes
                ! DEFINE CV_LILOC:
                CV_LILOC = 1
                ICORNER = eles_with_pipe(k)%pipe_corner_nds1( IPIPE )
                ICORNER1 = ICORNER
                CV_ILOC = CV_LOC_CORNER( ICORNER )
                CV_GL_LOC( CV_LILOC ) = CV_ILOC
                CV_LILOC = CV_LNLOC
                ICORNER = eles_with_pipe(k)%pipe_corner_nds2( IPIPE )
                ICORNER2 = ICORNER
                CV_ILOC = CV_LOC_CORNER( ICORNER )
                CV_GL_LOC( CV_LILOC ) = CV_ILOC
                IF ( CV_QUADRATIC ) CV_GL_LOC(2) = CV_MID_SIDE( ICORNER1, ICORNER2 )
                CV_GL_GL( : ) = ndgln%cv( (ELE-1)*Mdims%cv_nloc + CV_GL_LOC )
                MAT_GL_GL( : ) = ndgln%mat( (ELE-1)*Mdims%cv_nloc + CV_GL_LOC )
                ! DEFINE U_LILOC:
                if (U_P0DG) then
                    U_GL_LOC = 1
                else
                    U_LILOC = 1
                    ICORNER = eles_with_pipe(k)%pipe_corner_nds1( IPIPE )
                    U_ILOC = U_LOC_CORNER( ICORNER )
                    U_GL_LOC( U_LILOC ) = U_ILOC
                    U_LILOC = U_LNLOC
                    ICORNER = eles_with_pipe(k)%pipe_corner_nds2( IPIPE )
                    U_ILOC = U_LOC_CORNER( ICORNER )
                    U_GL_LOC( U_LILOC ) = U_ILOC
                    IF ( U_QUADRATIC ) U_GL_LOC( 2 ) = U_MID_SIDE( ICORNER1, ICORNER2 )
                end if
                DO U_LILOC = 1, U_LNLOC
                    U_ILOC = U_GL_LOC( U_LILOC )
                    U_GL_GL( U_LILOC ) = ndgln%u( (ELE-1)*Mdims%u_nloc + U_ILOC )
                END DO

                ! Calculate element angle sweeped out by element and pipe
                IF ( Mdims%ndim == 2 ) THEN
                    ELE_ANGLE = PI
                ELSE
                    ! find the nodes other than the pipe end corner nodes...
                    ICORNER3 = 0
                    DO ICORNER = 1, NCORNER
                        IF ( ICORNER /= ICORNER1 ) THEN
                            IF ( ICORNER /= ICORNER2 ) THEN
                                IF ( ICORNER3 == 0 ) THEN
                                    ICORNER3 = ICORNER
                                ELSE
                                    ICORNER4 = ICORNER
                                END IF
                            END IF
                        END IF
                    END DO
                    ELE_ANGLE = CALC_ELE_ANGLE_3D( X_ALL_CORN(:, ICORNER1), X_ALL_CORN(:, ICORNER2) , &
                        &                         X_ALL_CORN(:, ICORNER3), X_ALL_CORN(:, ICORNER4) )
                END IF
                DIRECTION = X_ALL_CORN( :, ICORNER2 ) - X_ALL_CORN( :, ICORNER1 )
                DX = SQRT( SUM( DIRECTION**2 ) )
                DIRECTION = DIRECTION / DX

                IF ( PIPE_MIN_DIAM ) THEN
                    MIN_DIAM = MINVAL( PIPE_diameter%val( CV_GL_GL( : ) ) )
                    PIPE_DIAM_GI = MIN_DIAM
                    DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                        SIGMA_GI(IPHASE,:) = MINVAL( SIGMA(IPHASE, MAT_GL_GL( : ) ) )
                    END DO
                ELSE
                    PIPE_DIAM_GI = 0.0
                    DO CV_LILOC = 1, CV_LNLOC
                        CV_KNOD = CV_GL_GL( CV_LILOC )
                        PIPE_DIAM_GI = PIPE_DIAM_GI + PIPE_DIAMETER%val( CV_KNOD ) * SCVFEN( CV_LILOC, : )
                    END DO
                    PIPE_DIAM_GI = MAX(PIPE_DIAM_GI, 0.0)
                    SIGMA_GI(:,:)=0.0
                    DO IPHASE=Mdims%n_in_pres+1,Mdims%nphase
                        DO CV_LILOC = 1, CV_LNLOC
                            MAT_KNOD = MAT_GL_GL( CV_LILOC )
                            SIGMA_GI(IPHASE,:) = SIGMA_GI(IPHASE,:) + SIGMA( IPHASE, MAT_KNOD ) * SCVFEN( CV_LILOC, : )
                        END DO
                    END DO
                END IF
                ! Recalculate SIGMA if we need to...
                IF ( CALC_SIGMA_PIPE ) THEN
                    MIN_DIAM = MINVAL( PIPE_diameter%val( CV_GL_GL( : ) ) )
                    DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                        IPHASE_IN_PIPE=IPHASE-Mdims%n_in_pres
                        DO GI = 1, scvngi
                            U_GI = 0.0
                            DO IDIM = 1, Mdims%ndim
                                U_GI = U_GI + SUM( SUFEN( : , GI ) * NU_ALL( IDIM, IPHASE, U_GL_GL( : ) ) ) * DIRECTION( IDIM )
                            END DO
                            CALL SIGMA_PIPE_FRICTION( SIGMA_GI( IPHASE, GI ), U_GI, MIN_DIAM, WELL_DENSITY( IPHASE ), WELL_VISCOSITY( IPHASE ), E_ROUGHNESS )
                        END DO
                    END DO
                END IF
                ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
                ! Adjust according to the volume of the pipe...
                DETWEI = SCVFEWEIGH * 0.5* DX * PI * ( (0.5*PIPE_DIAM_GI)**2 ) * ELE_ANGLE / ( 2.0 * PI )
                L_CVFENX_ALL = 2.0 * SCVFENLX / DX
                L_UFENX_ALL = 2.0 * SUFENLX / DX

                DO CV_LILOC = 1, CV_LNLOC
                    DO IDIM = 1, Mdims%ndim
                        L_CVFENX_ALL_REVERSED( IDIM, :, CV_LILOC ) = L_CVFENX_ALL( CV_LILOC, : ) * DIRECTION( IDIM )
                    END DO
                END DO

                DO U_LILOC = 1, U_LNLOC
                    L_UFEN_REVERSED( :, U_LILOC ) = SUFEN( U_LILOC, : )
                END DO

                ! Add in Mmat%C matrix contribution: (DG velocities)
                DO U_LILOC = 1, U_LNLOC
                    IU_NOD = U_GL_GL( U_LILOC ) ! ndgln%u( ( ELE - 1 ) * Mdims%u_nloc + U_ILOC )
                    DO P_LJLOC = 1, CV_LNLOC
                        JCV_NOD = CV_GL_GL(P_LJLOC) ! ndgln%p( ( ELE - 1 ) * Mdims%p_nloc + P_JLOC )
                        ! In this section we multiply the shape functions over the GI points. i.e: we perform the integration
                        ! over the element of the pressure like source term.
                        ! Put into matrix
                        DO COUNT2 = Mspars%C%fin(IU_NOD), Mspars%C%fin(IU_NOD+1)-1
                            IF ( Mspars%C%col(COUNT2) == JCV_NOD ) COUNT = COUNT2
                        END DO
                        ! Prepare aid variable NMX_ALL to improve the speed of the calculations
                        NMX_ALL( : ) = matmul( L_CVFENX_ALL_REVERSED( :, :, P_LJLOC ), DETWEI( : ) * L_UFEN_REVERSED( :, U_LILOC ) )
                        NM = sum( L_UFEN_REVERSED( :, U_LILOC ) * SCVFEN( P_LJLOC, : ) * DETWEI( : ) )
                        ! Put into matrix
                        DO IDIM = 1, Mdims%ndim
                            Mmat%C( IDIM, Mdims%n_in_pres+1:Mdims%nphase, COUNT ) = Mmat%C( IDIM, Mdims%n_in_pres+1:Mdims%nphase, COUNT ) - NMX_ALL( IDIM )
                        END DO
                        DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                            DO IDIM = 1, Mdims%ndim
                                ! This is ( S \dot n ) n
                                Mmat%U_RHS( IDIM, IPHASE, IU_NOD ) =  Mmat%U_RHS( IDIM, IPHASE, IU_NOD ) + &
                                    NM * dot_product(U_SOURCE_CV( :, IPHASE, JCV_NOD ), DIRECTION( : ) ) * DIRECTION( IDIM )
                            END DO
                        END DO
                    END DO ! DO P_LJLOC = 1, CV_LNLOC
                    U_ILOC = U_GL_LOC( U_LILOC )
                    DO U_LJLOC = 1, U_LNLOC
                        U_JLOC = U_GL_LOC( U_LJLOC )
                        JU_NOD = U_GL_GL( U_LJLOC )
                        IF ( GET_PIVIT_MAT ) THEN
                            DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                                NN = SUM( L_UFEN_REVERSED( :, U_LILOC ) * L_UFEN_REVERSED( :, U_LJLOC ) * DETWEI( : ) * SIGMA_GI( IPHASE , : ) )
                                JPHASE = IPHASE
                                DO IDIM = 1, Mdims%ndim
                                    DO JDIM = 1, Mdims%ndim
                                        i_indx = IDIM + (IPHASE-1)*Mdims%ndim + (U_ILOC-1)*Mdims%ndim*Mdims%nphase
                                        j_indx = JDIM + (JPHASE-1)*Mdims%ndim + (U_JLOC-1)*Mdims%ndim*Mdims%nphase
                                        Mmat%PIVIT_MAT( i_indx, j_indx, ele ) = Mmat%PIVIT_MAT( i_indx, j_indx, ele ) + NN * DIRECTION( IDIM ) * DIRECTION( JDIM )
                                    END DO
                                END DO
                            END DO
                        END IF ! GET_PIVIT_MAT
                        if ( u_source%have_field .and. .false.) then!DISABLED SOURCES
                            NN = sum( L_UFEN_REVERSED( :, U_LILOC ) * L_UFEN_REVERSED( :, U_LJLOC ) * DETWEI( : ) )
                            DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                                DO IDIM = 1, Mdims%ndim
                                    ! This is ( S \dot n ) n
                                    Mmat%U_RHS( IDIM, IPHASE, IU_NOD ) =  Mmat%U_RHS( IDIM, IPHASE, IU_NOD ) + &
                                        NN * SUM(U_SOURCE%VAL( :, IPHASE, 1, JU_NOD ) * DIRECTION( : ) ) * DIRECTION( IDIM )
                                END DO
                            END DO
                        end if
                    END DO ! DO U_LJLOC = 1, U_LNLOC
                END DO ! DO U_LILOC = 1, U_LNLOC

                ! Add pressure b.c to matrix and rhs...
                CV_LOC1 = CV_LOC_CORNER( ICORNER1 )
                JCV_NOD1 = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_LOC1 )
                CV_LOC2 = CV_LOC_CORNER( ICORNER2 )
                JCV_NOD2 = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_LOC2 )
                DO IPRES = 2, Mdims%npres
                    JCV_NOD = 0
                    IF ( WIC_P_BC_ALL_NODS( IPRES, JCV_NOD1 ) == WIC_P_BC_DIRICHLET ) THEN
                        CV_LILOC = 1
                        JCV_NOD = JCV_NOD1
                        U_LILOC = 1
                        JU_NOD = U_GL_GL( U_LILOC )
                        direction_norm = - direction ! for the b.c it must be negative at the bottom of element
                    END IF
                    IF ( WIC_P_BC_ALL_NODS( IPRES, JCV_NOD2 ) == WIC_P_BC_DIRICHLET ) THEN
                        CV_LILOC = CV_LNLOC
                        JCV_NOD = JCV_NOD2
                        U_LILOC = U_LNLOC
                        JU_NOD = U_GL_GL( U_LILOC )
                        direction_norm = +direction ! for the b.c it must be positive at the top of element
                    END IF
                    IF ( JCV_NOD /= 0 ) THEN
                        ! Add in Mmat%C matrix contribution: (DG velocities)
                        ! In this section we multiply the shape functions over the GI points. i.e: we perform the integration
                        ! over the element of the pressure like source term.
                        ! Put into matrix
                        PIPE_DIAM_END = PIPE_diameter%val( JCV_NOD )
                        DO COUNT2 = Mspars%C%fin(JU_NOD), Mspars%C%fin(JU_NOD+1)-1
                            IF ( Mspars%C%col(COUNT2) == JCV_NOD ) COUNT = COUNT2
                        END DO
                        ! Prepare aid variable NMX_ALL to improve the speed of the calculations
                        suf_area = PI * ( (0.5*PIPE_DIAM_END)**2. ) * ELE_ANGLE/ ( 2.0 * PI )
                        NMX_ALL( : ) = direction_norm* suf_area
                        LOC_U_RHS_U_ILOC = 0.0
                        DO IPHASE = 1+(IPRES-1)*Mdims%n_in_pres, IPRES*Mdims%n_in_pres
                            DO IDIM = 1, Mdims%ndim
                                Mmat%C( IDIM, IPHASE, COUNT ) = Mmat%C( IDIM, IPHASE, COUNT ) + NMX_ALL( IDIM )
                                LOC_U_RHS_U_ILOC( IDIM, IPHASE) =  LOC_U_RHS_U_ILOC( IDIM, IPHASE ) &
                                    - NMX_ALL( IDIM ) * SUF_P_BC_ALL_NODS( IPRES,JCV_NOD )
                            END DO
                        END DO
                        Mmat%U_RHS( :, Mdims%n_in_pres+1:Mdims%nphase, JU_NOD ) = Mmat%U_RHS( :, Mdims%n_in_pres+1:Mdims%nphase, JU_NOD ) + LOC_U_RHS_U_ILOC( :, Mdims%n_in_pres+1:Mdims%nphase)
                    END IF
                END DO ! DO IPRES = 2, Mdims%npres
            END DO ! DO IPIPE = 1, NPIPES
        END DO ! DO ELE = 1, Mdims%totele
        if ( GET_PIVIT_MAT ) then
            DO U_ILOC = 1, Mdims%u_nloc
                U_JLOC = U_ILOC
                DO IPHASE = Mdims%n_in_pres+1, Mdims%nphase
                    JPHASE = IPHASE
                    DO IDIM = 1, Mdims%ndim
                        JDIM = IDIM
                        i_indx = IDIM + (IPHASE-1)*Mdims%ndim + (U_ILOC-1)*Mdims%ndim*Mdims%nphase
                        j_indx = JDIM + (JPHASE-1)*Mdims%ndim + (U_JLOC-1)*Mdims%ndim*Mdims%nphase
                        WHERE ( abs( Mmat%PIVIT_MAT(i_indx, j_indx, :) ) < 1.0e-10 )
                            Mmat%PIVIT_MAT(i_indx, j_indx, :) = 1.0
                        END WHERE
                    END DO
                END DO
            END DO
        end if
        RETURN
    END SUBROUTINE MOD_1D_FORCE_BAL_C

    !>@brief: Introduces friction within the pipes based on a Moody diagram approach and based on a given roughness on the pipes
    SUBROUTINE SIGMA_PIPE_FRICTION( SIGMA, U, DIAM, DEN, VISC, E_ROUGHNESS )
        IMPLICIT NONE
        REAL, INTENT( IN ) :: U,DIAM,DEN,VISC,E_ROUGHNESS
        REAL, INTENT( INOUT ) :: SIGMA

        REAL, PARAMETER :: TOLER=1.E-10
        REAL :: RE, A, MAX_RE, F, F1, F2, W, MAX_U

        MAX_U = MAX( TOLER, ABS( U ) )
        RE = DEN * MAX_U * DIAM / VISC

        MAX_RE = RE
        A = LOG10( 6.9/MAX_RE + (E_ROUGHNESS/(3.7*DIAM))**(10.0/9.0) )

        F1 = 16.0 / MAX_RE
        F2 = 1.0 / MAX( -3.6 * A, TOLER )

        W = MAX( 0.0,  MIN( 1.0, (RE-2000.0)/2000.0 ) ) ! Relaxation factor between 2 friction factor expressions

        F = ( 1 - W ) * F1 + W * F2

        SIGMA = ( F / DIAM ) * 2.0 * DEN * ABS( U ) ! Based on Fanning friction factor
        !SIGMA = ( F / DIAM ) * 2.0 * DEN * 1.0 ! Based on Fanning friction factor
        !SIGMA = ( F / DIAM ) * 2.0 * DEN * MAX_U ! Based on Fanning friction factor

        RETURN
    END SUBROUTINE SIGMA_PIPE_FRICTION





    !>@brief: Calculate the local corner nodes...
    !> CV_MID_SIDE(ICORN,JCORN)= CV_ILOC local node number for node between these two corner nodes
    SUBROUTINE CALC_CORNER_NODS( CV_LOC_CORNER, NDIM, CV_NLOC, CV_QUADRATIC, CV_MID_SIDE )
        IMPLICIT NONE
        INTEGER, INTENT( IN ) :: CV_NLOC, NDIM
        INTEGER, DIMENSION( : ), INTENT( INOUT ) :: CV_LOC_CORNER
        INTEGER, DIMENSION( :, : ), INTENT( INOUT ), optional :: CV_MID_SIDE
        LOGICAL, INTENT( INOUT ), optional :: CV_QUADRATIC
        INTEGER :: CV_ILOC, ICORN, JCORN, CV_NCORNER

        CV_NCORNER = NDIM + 1
        if ( CV_NLOC==1) then!For P0DG the simplest case
            if (present(CV_QUADRATIC) ) CV_QUADRATIC=.FALSE.
            CV_LOC_CORNER(1)=1
        else IF( ((CV_NLOC==3).AND.(NDIM==2)) .OR. ((CV_NLOC==4).AND.(NDIM==3))) THEN ! assume linear...
            if (present(CV_QUADRATIC) ) CV_QUADRATIC=.FALSE.
            DO CV_ILOC=1,CV_NLOC
                CV_LOC_CORNER(CV_ILOC)=CV_ILOC
            END DO
        ELSE ! assume quadratic...
            if (present(CV_QUADRATIC) ) CV_QUADRATIC=.TRUE.
            CV_LOC_CORNER(1)=1
            CV_LOC_CORNER(2)=3
            CV_LOC_CORNER(3)=6

            CV_MID_SIDE=0
            CV_MID_SIDE(1,2)=2
            CV_MID_SIDE(1,3)=4
            CV_MID_SIDE(2,3)=5

            IF(NDIM==3) THEN
                CV_MID_SIDE(1,4)=7
                CV_MID_SIDE(2,4)=8
                CV_MID_SIDE(3,4)=9
                CV_LOC_CORNER(4)=10
            END IF
            if (present(CV_MID_SIDE)) then
                DO ICORN=1,CV_NCORNER
                    DO JCORN=ICORN+1,CV_NCORNER
                        CV_MID_SIDE(JCORN,ICORN)=CV_MID_SIDE(ICORN,JCORN)
                    END DO
                END DO
            end if
        END IF

        RETURN
    END SUBROUTINE CALC_CORNER_NODS

    !>@brief: Calculate element angle sweeped out by element and pipe
    !> X_ALL_CORN_PIPE1, X_ALL_CORN_PIPE2 are the coordinates of the ends of the pipe within an element.
    !> X_ALL_CORN_PIPE3, X_ALL_CORN_PIPE4 are the other corner 2 nodes of an element.
    REAL FUNCTION CALC_ELE_ANGLE_3D( X_ALL_CORN_PIPE1, X_ALL_CORN_PIPE2, X_ALL_CORN_PIPE3, X_ALL_CORN_PIPE4 )
        IMPLICIT NONE
        REAL, intent( in ) :: X_ALL_CORN_PIPE1(3), X_ALL_CORN_PIPE2(3),  X_ALL_CORN_PIPE3(3), X_ALL_CORN_PIPE4(3)
        REAL :: X_PIPE1(3), X_PIPE2(3), X_PIPE3(3), X_PIPE4(3)
        REAL :: X_PIPE3_N(3), X_PIPE4_N(3)
        REAL :: X_PIPE3_2D(2), X_PIPE4_2D(2)
        REAL :: NORM(3), T1(3), T2(3), A(3,3)
        INTEGER :: IDIM

        ! Translate
        X_PIPE1 = 0.0
        X_PIPE2 = X_ALL_CORN_PIPE2 - X_ALL_CORN_PIPE1
        X_PIPE3 = X_ALL_CORN_PIPE3 - X_ALL_CORN_PIPE1
        X_PIPE4 = X_ALL_CORN_PIPE4 - X_ALL_CORN_PIPE1

        ! Rotate to be aligned with x-axis
        Norm = X_PIPE2 / SQRT( SUM( X_PIPE2**2 ) )

        ! Determine the tangent and bi-normal vectors from the normal NormX, NormY, NormZ
        CALL GET_TANG_BINORM( Norm(1), Norm(2), Norm(3), T1(1), T1(2), T1(3), T2(1), T2(2), T2(3), 1 )

        A(1,:) = NORM
        A(2,:) = T1
        A(3,:) = T2

        DO IDIM = 1, 3
            X_PIPE3_N(IDIM) = SUM( A(IDIM,:) * X_PIPE3 )
            X_PIPE4_N(IDIM) = SUM( A(IDIM,:) * X_PIPE4 )
        END DO

        X_PIPE3_2D(1:2) = X_PIPE3_N(2:3)
        X_PIPE4_2D(1:2) = X_PIPE4_N(2:3)

        ! Use cosine rule....
        CALC_ELE_ANGLE_3D = ACOS( SUM( X_PIPE3_2D * X_PIPE4_2D ) / ( SQRT( SUM( X_PIPE3_2D**2 ) ) * SQRT( SUM( X_PIPE4_2D**2 ) ) ) )

        RETURN
    END FUNCTION CALC_ELE_ANGLE_3D


    !>@brief: In this subroutine the elements that contain pipes are identified
    !> The pipes can either be defined using python (DEPRECATED)
    !> or a nastran file defining the trajectory of each well by points
    subroutine retrieve_pipes_coords(state, packed_state, Mdims, ndgln, eles_with_pipe)
        implicit none
        type(state_type), dimension(:), intent(inout) :: state
        type(state_type), intent(in) :: packed_state
        type(pipe_coords), dimension(:), allocatable, intent(inout) :: eles_with_pipe!>allocated inside
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_ndgln), intent(in) :: ndgln
        !Local variables
        type(vector_field), pointer :: X
        integer :: ele, k, ICORNER, CV_ILOC, CV_NODI, PIPE_NOD_COUNT, NCORNER, x_iloc, x_iloc2
        real, dimension(Mdims%cv_nonods) :: diameter_of_the_pipe_aux
!        integer, parameter :: NCORNER = Mdims%ndim + 1!Is this because we only consider P1 tets and triangles?
        type(pipe_coords), dimension(Mdims%totele):: AUX_eles_with_pipe
        type(scalar_field), pointer :: pipe_diameter
        integer, dimension(Mdims%ndim + 1) :: CV_LOC_CORNER!NCORNER
        logical, dimension(Mdims%ndim + 1) :: PIPE_INDEX_LOGICAL!NCORNER
        !Variables to read from input files
        integer :: number_well_files, edge, stat, nb
        real :: diam
        character( len = option_path_len ):: file_path
        real, dimension(:,:), allocatable :: nodes
        integer, dimension(:,:), allocatable :: edges
        integer, dimension(:), allocatable :: pipe_seeds
        type( tensor_field ), pointer:: tfield
        type (scalar_field), pointer :: well_domains
        logical :: skip
        logical, save :: dump_vtu_zero = .true.
        integer, dimension(:), pointer :: neighbours
        !Initialise
        AUX_eles_with_pipe%ele = -1
        k = 1; NCORNER = Mdims%ndim + 1
        PIPE_DIAMETER => extract_scalar_field( state(1), "DiameterPipe" )
        X => EXTRACT_VECTOR_FIELD( packed_state, "PressureCoordinate" )
        !Create list of corners
        call CALC_CORNER_NODS( CV_LOC_CORNER, Mdims%NDIM, Mdims%CV_NLOC)
        !Retrieve, if there are any number of input .bdf files
        number_well_files = option_count("/porous_media/wells_and_pipes/well_from_file")
        if (number_well_files > 0) then
            !Need the mesh to get neighbouring elements
            tfield => extract_tensor_field( packed_state, "PackedFEPressure" )

            !for the time being we re-use from diamond
            diam = maxval(PIPE_DIAMETER%val)
            !Clean eles_with_pipes before reading the files
            if (allocated(eles_with_pipe)) deallocate(eles_with_pipe)
            diameter_of_the_pipe_aux = 0.

            well_domains =>  extract_scalar_field( state(1), "Well_domains" , stat)
            if (stat /= 0) then
                if (getprocno() == 1)then
                    FLExit( "ERROR: Wells defined by a file requires the well_volumes_ids to be specified")
                end if
            end if
            do k = 1, number_well_files
                !First identify the well trajectory
                call get_option("/porous_media/wells_and_pipes/well_from_file["// int2str(k-1) //"]/file_path", file_path)
                call read_nastran_file(file_path, nodes, edges)
                call find_pipe_seeds(well_domains, X%val, nodes, edges, pipe_seeds)
                !Only if a seed is found then the well is constructed
                if ( size(pipe_seeds)>0 ) call find_nodes_of_well(X%val, nodes, edges, pipe_seeds, eles_with_pipe, diameter_of_the_pipe_aux)
                deallocate(nodes, edges)!because nodes and edges are allocated inside read_nastran_file
                deallocate(pipe_seeds)
            end do

            if (.not.allocated(eles_with_pipe)) allocate(eles_with_pipe(0)) !This if is important for parallel so it exists and the loops are skipped
            !Re-populate properly PIPE_DIAMETER
            PIPE_DIAMETER%val = 0.
            !Copy values back to PIPE_DIAMETER from diameter_of_the_pipe_aux. This should go over less than 1% of the nodes

            do ele = 1, size(eles_with_pipe)
              !For parallel we only assign a diameter value if the node is owned by the processor
              if (IsParallel()) then
                  if (.not. assemble_ele(PIPE_DIAMETER,ele)) then
                      skip=.true.
                      neighbours=>ele_neigh(PIPE_DIAMETER,ele)
                      do nb=1,size(neighbours)
                          if (neighbours(nb)<=0) cycle
                          if (assemble_ele(PIPE_DIAMETER,neighbours(nb))) then
                              skip=.false.
                              exit
                          end if
                      end do
                      if (skip) cycle
                  end if
              end if

                do k = 1, eles_with_pipe(ele)%npipes
                    x_iloc = eles_with_pipe(ele)%pipe_corner_nds1(k)
                    PIPE_DIAMETER%val(ndgln%cv( ( eles_with_pipe(ele)%ele - 1 ) * Mdims%cv_nloc + x_iloc )) = &
                        diameter_of_the_pipe_aux(ndgln%cv( ( eles_with_pipe(ele)%ele - 1 ) * Mdims%cv_nloc + x_iloc ))
                    x_iloc = eles_with_pipe(ele)%pipe_corner_nds2(k)
                    PIPE_DIAMETER%val(ndgln%cv( ( eles_with_pipe(ele)%ele - 1 ) * Mdims%cv_nloc + x_iloc )) = &
                     diameter_of_the_pipe_aux(ndgln%cv( ( eles_with_pipe(ele)%ele - 1 ) * Mdims%cv_nloc + x_iloc ))
                end do
            end do

        else
            do ele = 1, Mdims%totele

              !For parallel we only assign a diameter value if the node is owned by the processor
              if (IsParallel()) then
                  if (.not. assemble_ele(PIPE_DIAMETER,ele)) then
                      skip=.true.
                      neighbours=>ele_neigh(PIPE_DIAMETER,ele)
                      do nb=1,size(neighbours)
                          if (neighbours(nb)<=0) cycle
                          if (assemble_ele(PIPE_DIAMETER,neighbours(nb))) then
                              skip=.false.
                              exit
                          end if
                      end do
                      if (skip) cycle
                  end if
              end if


                ! Look for pipe indicator in element:
                PIPE_INDEX_LOGICAL=.FALSE.
                PIPE_NOD_COUNT=0
                DO ICORNER = 1, NCORNER ! Number of corner nodes in element
                    CV_ILOC = CV_LOC_CORNER( ICORNER )
                    CV_NODI = ndgln%cv( ( ELE - 1 ) * Mdims%cv_nloc + CV_ILOC )
                    IF ( PIPE_diameter%val( CV_NODI ) > 0.0 ) THEN
                        PIPE_NOD_COUNT = PIPE_NOD_COUNT + 1
                        PIPE_INDEX_LOGICAL(ICORNER) = .TRUE.
                    END IF
                END DO
                !If it has pipe, then store the element
                if (PIPE_NOD_COUNT-1 >= 1) then
                    AUX_eles_with_pipe(k)%ele = ele
                    AUX_eles_with_pipe(k)%npipes = 1!just one pipe !PIPE_NOD_COUNT - 1
                    allocate(AUX_eles_with_pipe(k)%pipe_index(NCORNER))
                    AUX_eles_with_pipe(k)%pipe_index = PIPE_INDEX_LOGICAL
                    k = k + 1
                end if
            end do
            !Store only the required elements
            if(allocated(eles_with_pipe)) deallocate(eles_with_pipe)!re-adjust if required
            allocate(eles_with_pipe(k-1))

            do k = 1, size(eles_with_pipe)
                allocate(eles_with_pipe(k)%pipe_index(NCORNER))
                eles_with_pipe(k) = AUX_eles_with_pipe(k)
                allocate(eles_with_pipe(k)%pipe_corner_nds1(eles_with_pipe(k)%npipes))
                allocate(eles_with_pipe(k)%pipe_corner_nds2(eles_with_pipe(k)%npipes))
                ! Calculate the pipes within an element...
                ! we return the pipe corner nodes for each pipe in
                ! If we have more than one pipe then choose the 2 edges with the shortest sides
                ! and have a maximum of 2 pipes per element...
                call CALC_PIPES_IN_ELE( x%val(:, ndgln%x( ( eles_with_pipe(k)%ele - 1 ) * Mdims%cv_nloc + CV_LOC_CORNER(1:NCORNER)) ),&
                                    eles_with_pipe(k)%pipe_index, Mdims%ndim, eles_with_pipe(k)%pipe_corner_nds1,   &
                                    eles_with_pipe(k)%pipe_corner_nds2, eles_with_pipe(k)%npipes )
            end do
        end if
        if (dump_vtu_zero) then
            !Overwrite the first dump with the values of the diameter of the pipe, don't remove this section
            ewrite(1,*) "Overwritting the initial vtu file with the updated values of the diameter of the well"
            k = 0
            call write_state( k, state )
            dump_vtu_zero = .false.
        end if
    contains

        !!>@brief: Check whether a node is within a pipe
        !> We define a virtual cylinder and we check whether the node falls within it.
        logical function is_within_pipe(P, v1, v2, tol)
            implicit none
            real, dimension(:), intent(in) :: P, v1, v2
            real, intent(in) :: tol
            !local variables
            real, dimension(Mdims%ndim) :: vec1, vec2, Vaux
            real :: c1, c2, distance, Saux1, Saux2, diam
            !To compare if a mesh node falls within the well section we need
            !to use a tighter tolerance than when checking if the node is within a cirtual well cylinder
            !this is because we are not affected by the lack of precision of the nastran .bdf files
            real, parameter :: single_precision_tol = 1e-8
            !Initialiase variables
            is_within_pipe = .false.

            vec1 = v2(1:Mdims%ndim) - v1(1:Mdims%ndim)
            vec2 = P(1:Mdims%ndim) - v1(1:Mdims%ndim)

            c1 = dot_product(Vec2,Vec1)
            c2 = dot_product(Vec1,Vec1)!<=lenght of the section**2
            !First we check that the point is between the two vertexes
            if (c1 <= single_precision_tol)then!Before v1
                distance = sqrt(dot_product(Vec2,Vec2))
            else if ( c2 - c1 <= single_precision_tol)then!after v2
                Vec2 = P(1:Mdims%ndim)-v2(1:Mdims%ndim)
                distance = sqrt(dot_product(Vec2,Vec2))
            else !Calculate distance to a line
                Vaux = abs(cross_product(P(1:Mdims%ndim)-v1(1:Mdims%ndim),P(1:Mdims%ndim)-v2(1:Mdims%ndim)))
                Saux1 = sqrt(dot_product(Vaux,Vaux))

                Vaux = abs(v2(1:Mdims%ndim)-v1(1:Mdims%ndim))
                Saux2 = sqrt(dot_product(Vaux,Vaux))
                distance = Saux1/max(Saux2,1d-16)

           end if

!           is_within_pipe = (distance <= diameter)
            !Use a relative tolerance to the lenght of the section, as mesh adaptivity allows a bit of movement
            diam = max(0.01 * sqrt(c2), tol)
           is_within_pipe = (distance <= diam)
        end function is_within_pipe


        !>@brief: For a given NASTRAN file find a node that is within the pipe
        !> It uses brute force but only within the regions that contain a well, the sleeves defined in diamond
        subroutine find_pipe_seeds(well_domains, X, nodes, edges, pipe_seeds)
            implicit none
            type (scalar_field), pointer :: well_domains
            real, dimension(:,:), intent(in) :: X
            real, dimension(:,:), allocatable, intent(in) :: nodes
            integer, dimension(:,:), allocatable, intent(in) :: edges
            integer, dimension(:), allocatable, intent(inout) :: pipe_seeds
            !Local variables
            logical, save :: first_time =.true.
            logical :: found
            integer :: i, j, l, count
            integer :: sele, siloc, sinod, ele, iloc, inod
            real, dimension(Mdims%stotel) :: aux_pipe_seeds
            aux_pipe_seeds = -1
            !Initialise tolerancePipe just once per simulation
            if (first_time) then
                first_time = .false.
                if (have_option('/porous_media/wells_and_pipes/well_options/wells_bdf_tolerance')) then
                    call get_option('/porous_media/wells_and_pipes/well_options/wells_bdf_tolerance', tolerancePipe)
                end if
            end if

            !Use brute force through the whole mesh but selecting only elements that live within the wells region ids to find a seed
            !This should definetively work in parallel and be cheap
            l = 0; aux_pipe_seeds = -1
            element_loop: do ele = 1, size(well_domains%val)
                if (well_domains%val(ele) <= 1e-8) cycle!Only look at elements within the well regions, the rest is set to 0
                do iloc  = 1, Mdims%cv_nloc
                    inod = ndgln%p( ( ele - 1 ) * Mdims%cv_nloc + iloc )
                    do edge = 1, size(edges,2)
                        if (is_within_pipe(X(:,inod), nodes(:,edges(1,edge)), nodes(:,edges(2,edge)), tolerancePipe)) then
                            found = .false.
                            do j = 1, size(aux_pipe_seeds)!Make sure that we do not store the same position many times
                                if (aux_pipe_seeds(j)==inod) found = .true.
                                if (aux_pipe_seeds(j) < 0 .or. found) exit
                            end do
                            if (.not.found) then
                                l = l + 1
                                aux_pipe_seeds(l) = inod!Store the global node number
                                exit element_loop !just one seed per well
                            end if
                        end if
                    end do
                end do
            end do element_loop
            !OLD METHOD, June/2018. REMOVEME IF THIS DATE IS TOO OLD...
!            !Use brute force through the surface of the domain. This only works in serial...
!            !The cost is minimised by looping over the boundary of the domain only
!            l = 0; aux_pipe_seeds = -1
!            do sele = 1, Mdims%stotel
!                do siloc = 1, Mdims%p_snloc
!                    sinod = ndgln%suf_p( ( sele - 1 ) * Mdims%p_snloc + siloc )
!                    do edge = 1, size(edges,2)
!                        if (is_within_pipe(X(:,sinod), nodes(:,edges(1,edge)), nodes(:,edges(2,edge)), tolerancePipe)) then
!
!                            found = .false.
!                            do j = 1, size(aux_pipe_seeds)!Make sure that we do not store the same position many times
!                                if (aux_pipe_seeds(j)==sinod) found = .true.
!                                if (aux_pipe_seeds(j) < 0 .or. found) exit
!                            end do
!                            if (.not.found) then
!                                l = l + 1
!                                aux_pipe_seeds(l) = sinod!Store the global node number
!                            end if
!                        end if
!                    end do
!                end do
!            end do

            allocate(pipe_seeds(l))
            if (size(pipe_seeds)>0) pipe_seeds = aux_pipe_seeds(1:l)
        end subroutine find_pipe_seeds

        !>@brief: Once a seed for a node is found, this subroutine finds all the nodes that form a well by searching neirhbouring elements only
        subroutine find_nodes_of_well(X, nodes, edges, pipe_seeds, eles_with_pipe, diameter_of_the_pipe_aux)
            implicit none
            real, dimension(:,:), intent(in) :: X
            real, dimension(:,:), allocatable, intent(in) :: nodes
            integer, dimension(:,:), allocatable, intent(in) :: edges
            integer, dimension(:), intent(in) :: pipe_seeds
            type(pipe_coords), dimension(:), allocatable, intent(inout) :: eles_with_pipe
            real, dimension(:), intent(inout) :: diameter_of_the_pipe_aux
            !Local variables
            logical, save :: first_time = .true.
            integer, dimension(2, Mdims%totele) :: visited_eles!Number of element visited and neigbours used
            integer :: starting_ele, edge, neig, ele_bak, neig_bak, visit_counter, seed
            integer :: ele, ele2, inode, k, i, j, x_iloc, x_iloc2, x_inod, ipipe, first_node, first_loc, neighbours_left
            real :: aux
            real, dimension(Mdims%ndim) :: Vaux
            logical :: got_new_ele, touching_well, continue_looking, found, resize
            type(pipe_coords), dimension(:), allocatable :: AUX_eles_with_pipe, BAK_eles_with_pipe
            !Initialise AUX_eles_with_pipe
            allocate(AUX_eles_with_pipe(Mdims%totele))
            AUX_eles_with_pipe%ele = -1
            do j = 1, Mdims%totele!We can study this, but it is VERY unlikely that this goes beyond a 10%
                AUX_eles_with_pipe(j)%npipes = 0
                allocate(AUX_eles_with_pipe(j)%pipe_index(Mdims%ndim + 1))
                allocate(AUX_eles_with_pipe(j)%pipe_corner_nds1(Mdims%ndim))!Maximum of number of dimension pipes per element
                allocate(AUX_eles_with_pipe(j)%pipe_corner_nds2(Mdims%ndim))
            end do

            !First retrieve the first seed of the well
            seeds_loop: do seed = 1, size(pipe_seeds)
                do ele = 1, Mdims%totele
                    do x_iloc = 1, Mdims%x_nloc
                        x_inod = ndgln%x( ( ele - 1 ) * Mdims%x_nloc + x_iloc )
                        if (x_inod == pipe_seeds(seed)) then!element found
                            starting_ele = ele
                        end if
                    end do
                end do
                ele = starting_ele
                visited_eles(1,:) = -1; visited_eles(2,:) = 0
                visited_eles(1,1) = ele; visited_eles(2,1) = 1

                !Once we have the starting node we use that to go through the neighbouring nodes to build up the well and the connections
                k = 1; ele2 = ele
                ele_loop: do while (.true.)
                    neig = 0
                    neig_loop: do while (neig < Mdims%ndim + 1)!A tet/triangle element can have one neighbour more than dimensions it has
                        neig = neig + 1
                        touching_well = .false.
                        i = 1
                        loc_loop: do x_iloc = 1, Mdims%x_nloc
                            x_inod = ndgln%x( ( ele2 - 1 ) * Mdims%x_nloc + x_iloc )
                            do edge = 1, size(edges,2)!<= this can be optimised if we know that there is one well only defined per edges array
                                if (is_within_pipe(X(:,x_inod), nodes(:,edges(1,edge)), nodes(:,edges(2,edge)),  tolerancePipe)) then
                                    select case (i)
                                        case (1)!First true
                                            first_node = x_inod
                                            first_loc = x_iloc
                                            touching_well = .true.
                                            i = i + 1
                                            !backup just in case this element is not an in-between-element
                                            ele_bak = ele; neig_bak = neig!<= has to be BEFORE updating ele
                                            !Update position
                                            ele = ele2; neig = max(1, visited_eles(2,get_pos(ele2, visited_eles)))
                                            !And visited list
                                            visit_counter = get_pos(ele, visited_eles)
                                            visited_eles(1, visit_counter) = ele
                                            visited_eles(2, visit_counter) = max(visited_eles(2,visit_counter),neig)
                                            exit!to ensure that we are not in a joint and X(:,x_inod) is not considered twice
                                        case (2)!Second true, we got a well in the element!
                                            !Need to test if the nodes are already stored
                                            found = .false.
                                            j = 1
                                            do while (AUX_eles_with_pipe(j)%ele > 0)
                                                do ipipe = 1, AUX_eles_with_pipe(j)%npipes!Test all the available pipes
                                                    if ((first_node == AUX_eles_with_pipe(j)%pipe_corner_nds1(ipipe) .or.&
                                                        first_node == AUX_eles_with_pipe(j)%pipe_corner_nds2(ipipe)) .and. &
                                                        (x_inod == AUX_eles_with_pipe(j)%pipe_corner_nds1(ipipe) .or.&
                                                        x_inod == AUX_eles_with_pipe(j)%pipe_corner_nds2(ipipe))) then
                                                        found = .true.
                                                        touching_well = .false.
                                                        !Do not move to this new element to pivot around it
                                                        ele = ele_bak; neig = neig_bak
                                                        !Don't consider it again
                                                        visit_counter = get_pos(ele2, visited_eles)
                                                        visited_eles(1, visit_counter) = ele2
                                                        visited_eles(2, visit_counter) = 100
                                                        exit loc_loop
                                                    end if
                                                    j = j + 1
                                                end do
                                            end do
                                            if (.not.found) then
!                                                ipipe = AUX_eles_with_pipe(j)%npipes + 1 !Add one pipe
                                                ipipe = 1 !One pipe per element for the time being
                                                AUX_eles_with_pipe(j)%npipes = ipipe
                                                AUX_eles_with_pipe(k)%ele = ele2
                                                AUX_eles_with_pipe(k)%pipe_index(first_loc) = .true.!Don't know if necessary now...
                                                AUX_eles_with_pipe(k)%pipe_corner_nds1(ipipe) = first_loc!first_node
                                                AUX_eles_with_pipe(k)%pipe_index(x_iloc) = .true.!Don't know if necessary now...
                                                AUX_eles_with_pipe(k)%pipe_corner_nds2(ipipe) = x_iloc!x_inod
                                                k = k + 1
                                                exit loc_loop!Change reference to this element
                                            end if
                                    end select
                                end if
                            end do
                        end do loc_loop
                        if (.not.touching_well) then
                            !Don't consider it again
                            visit_counter = get_pos(ele2, visited_eles)
                            visited_eles(1, visit_counter) = ele2
                            visited_eles(2, visit_counter) = 100
                        end if
                        !Look for new proposed element
                        do while (neig <= Mdims%ndim + 1)
                            ele2 = max(ele_neigh(tfield%mesh, ele, neig),0)!This is the cv mesh; test next neighbour
                            !Update neighbour used in the list
                            i = get_pos(ele, visited_eles)
                            visited_eles(2,i) = max(visited_eles(2,i),neig+1)
                            got_new_ele = .true.
                            !Store element about to be inspected
                            i = get_pos(ele2, visited_eles)
                            visited_eles(1,i) = ele2
                            !If proposed element does not have available elements to study then find another
                            if (visited_eles(2,i) > Mdims%ndim + 1 .or. ele2 == 0) then!Ignore boundary
                                got_new_ele = .false.
                                !Advance to the next possible neighbour
                                neig = neig + 1
                                if (neig > Mdims%ndim + 1) then!If not more to look at then
                                    !Start to go back along the elements that still have neigbours to look at
                                    j = 1
                                    do while (visited_eles(2,j) > Mdims%ndim + 1)
                                        if (visited_eles(1,j) <= 0) then
                                            !Impossible to continue the search
                                            print *, "WARNING: Exit due to visited_list full"
                                            exit ele_loop
                                        end if
                                        j = j + 1
                                    end do
                                    visit_counter = max(j,1)
                                    if (visited_eles(1,j) > 0) then
                                        ele = visited_eles(1,j)
                                        neig = max(1, visited_eles(2,j))
                                    end if
                                end if
                            !if proposed element already visited (and was touching a well)
                            !then we can know if it was touching the well,
                            !in which case we move to it and continue the search
                            else if (visited_eles(2,i) >= 1 .and. visited_eles(2,i) <= Mdims%ndim + 1) then
                                !element already studied for wells, now we just want to find a neighbour
                                !Update position
                                ele = ele2
                                neig = visited_eles(2,get_pos(ele, visited_eles))
                                got_new_ele = .false.
                            end if
                            if (got_new_ele) exit


                        end do
                        if (neig > Mdims%ndim + 1) then
                            exit ele_loop!Can't continue searching so well finished
                        end if
                    end do neig_loop

                end do ele_loop
            end do seeds_loop


!print *, "ELEMENT COORDINATES"
! j = 1
! do while (visited_eles(1,j) > 0)
! ele = visited_eles(1,j)
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 1 ))
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 2 ))
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 1 ))
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 3 ))
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 1 ))
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 4 ))
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 2 ))
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 3 ))
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 2 ))
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 4 ))
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 3 ))
! print *, X(:,ndgln%x( ( ele - 1 ) * Mdims%x_nloc + 4 ))
! j = j + 1
! end do
! read*
            !Count useful values
            j = 0
            do while (AUX_eles_with_pipe(j+1)%ele > 0)
                j = j + 1
            end do
            if (j==0) return
            !Create copy if required of the original list
            resize = allocated(eles_with_pipe)
            k = 0
            if (resize) then
                k = size(eles_with_pipe)
                call copy_from_pipe_coords(eles_with_pipe, BAK_eles_with_pipe, 1, k, siz = k)
            end if
            !Now copy the new elements into the beginning of eles_with_pipes
            call copy_from_pipe_coords(Aux_eles_with_pipe, eles_with_pipe, 1, j, siz = k + j)
            !Finally if required, copy bak_eles back into eles_with_pipes
            if (resize) call copy_from_pipe_coords(BAK_eles_with_pipe, eles_with_pipe, j+1, k + j)




            !#######################################################################
            !####THIS NEEDS TO BE REVISITED ONCE THE MEMORY IS CORRECTLY CREATED####
            !Now, introduce the value of the diameter only in the correct regions
            !This should be temporary until it is being read from the well file as well.

            !Populate here diameter_of_the_pipe_aux
            !This method enables to use different diameters
            if (size(PIPE_DIAMETER%val) > 1) then
                !Find nodes that have wells and tag them. This first part should go over less than 1% of the nodes
                do ele = 1, size(eles_with_pipe)
                    do ipipe = 1, eles_with_pipe(ele)%npipes
                        !get the two ends of the pipe
                        x_iloc  = ndgln%cv( ( eles_with_pipe(ele)%ele - 1 ) * Mdims%cv_nloc + eles_with_pipe(ele)%pipe_corner_nds1(ipipe) )
                        x_iloc2 = ndgln%cv( ( eles_with_pipe(ele)%ele - 1 ) * Mdims%cv_nloc + eles_with_pipe(ele)%pipe_corner_nds2(ipipe) )
                        !ensure consistency within the pipe section
                        aux = max(PIPE_DIAMETER%val(x_iloc), PIPE_DIAMETER%val(x_iloc2 ) )
                        !Store value
                        diameter_of_the_pipe_aux(x_iloc)  = aux
                        diameter_of_the_pipe_aux(x_iloc2) = aux
                    end do
                end do
            end if

            ! !We have to ensure that the python prescribed field is not recalculated
            ! !this needs to be removed once the memory is properly allocated
            ! if (first_time .and. getprocno() == 1) then
            !     first_time = .false.
            !     if (have_option("porous_media/wells_and_pipes/scalar_field::DiameterPipe/prescribed")) &
            !         call add_option("porous_media/wells_and_pipes/scalar_field::DiameterPipe/prescribed/do_not_recalculate", stat = k)
            ! end if
            !#######################################################################
!    !To test the results gnuplot and the run spl'test' w linesp
!to compare with well plotted from multi_tools: spl'test'  using 1:2:3 with lines palette title "Eles", 'well_coords' with lines
!do j = 1, size(eles_with_pipe)
!print *, X(:,eles_with_pipe(j)%pipe_corner_nds1(1))
!print *, X(:,eles_with_pipe(j)%pipe_corner_nds2(1))
!end do
!read*
        end subroutine find_nodes_of_well

        subroutine copy_from_pipe_coords(original, copy, start, end, siz)
            !This subroutine copies from an input pipe_coords structure into another
            !from an starting point to a final point, allocating all the internal
            !variables and deallocating the original
            !siz is the size of the copy file if it is required ot allocate the file
            Implicit none
            integer, intent(in) :: start, end
            type(pipe_coords), dimension(:), allocatable, intent(inout) :: original, copy
            integer, optional, intent(in) :: siz
            !Local variables
            integer :: k, k_orig
            if (present(siz)) allocate(copy(siz))

            k_orig = 1
            do k = start, end
                allocate(copy(k)%pipe_index(Mdims%ndim + 1))
                allocate(copy(k)%pipe_corner_nds1(Mdims%ndim))
                allocate(copy(k)%pipe_corner_nds2(Mdims%ndim))
                copy(k)%ele = original(k_orig)%ele
                copy(k)%npipes = original(k_orig)%npipes
                copy(k)%pipe_index(1:size(copy(k)%pipe_index)) =&
                     original(k_orig)%pipe_index(1:size(copy(k)%pipe_index))
                copy(k)%pipe_corner_nds1(1:size(copy(k)%pipe_corner_nds1)) =&
                     original(k_orig)%pipe_corner_nds1(1:size(copy(k)%pipe_corner_nds1))
                copy(k)%pipe_corner_nds2(1:size(copy(k)%pipe_corner_nds2)) = &
                    original(k_orig)%pipe_corner_nds2(1:size(copy(k)%pipe_corner_nds2))
                deallocate(original(k_orig)%pipe_index)
                deallocate(original(k_orig)%pipe_corner_nds1, original(k_orig)%pipe_corner_nds2)
                k_orig = k_orig + 1
            end do
            deallocate(original)
        end subroutine copy_from_pipe_coords


        integer function get_pos(ele, visited_eles)
            !This either gives the position in the list where element is,
            ! or the new position to store the information in
            integer, intent(in) :: ele
            integer, dimension(2, Mdims%totele), intent(in) :: visited_eles
            !Local variables
            integer :: i
            do i = 1, size(visited_eles,2)
                if (visited_eles(1,i)==ele .or. visited_eles(1,i) <=0) then
                    get_pos = i
                    return
                end if
            end do
        end function get_pos


        SUBROUTINE CALC_PIPES_IN_ELE( X_ALL_CORN, PIPE_INDEX_LOGICAL, NDIM, &
            pipe_corner_nds1, pipe_corner_nds2, npipes )
            ! Calculate the pipes within an element...
            ! Return the pipe corner nodes for each pipe in element.
            IMPLICIT NONE
            INTEGER, intent( in ) :: NDIM, npipes
            REAL, intent( in ) :: X_ALL_CORN(NDIM,NDIM+1)
            LOGICAL, intent( in ) :: PIPE_INDEX_LOGICAL(NDIM+1)
            integer, intent( inout ) :: pipe_corner_nds1(npipes), pipe_corner_nds2(npipes)
            ! logal variables...
            REAL :: cp(ndim)
            REAL :: area_sqr(ndim+1)
            INTEGER :: I,J,K, II(ndim+1),JJ(ndim+1),KK(ndim+1), ipipe, i_nd,j_nd, iSTORE(4), jSTORE(4), kSTORE(4)
            INTEGER :: icorn, iface, iii, ik, il
            logical :: found_i, found_j
            real :: max_area

            IF ( NPIPES==1 ) THEN

                i_nd = 0
                DO ICORN = 1, NDIM+1
                    IF ( PIPE_INDEX_LOGICAL(ICORN) ) THEN
                        IF ( I_ND == 0 ) THEN
                            I_ND = ICORN
                        ELSE
                            J_ND = ICORN
                        END IF
                    END IF
                END DO
                ipipe = 1
                pipe_corner_nds1(ipipe) = i_nd
                pipe_corner_nds2(ipipe) = j_nd

            ELSE  ! IF ( NPIPES==1 ) THEN

                IF ( NDIM==2 ) THEN

                    i=1 ; j=2
                    area_sqr(1) = sum( ( X_ALL_CORN(:,1) - X_ALL_CORN(:,2) )**2 ) ! THIS IS THE LENGTH SQUARED
                    ii(1)=i ; jj(1)=j

                    i=1 ; j=3
                    area_sqr(2) = sum( ( X_ALL_CORN(:,1) - X_ALL_CORN(:,3) )**2 )
                    ii(2)=i ; jj(2)=j

                    i=2 ; j=3
                    area_sqr(3) = sum( ( X_ALL_CORN(:,2) - X_ALL_CORN(:,3) )**2 )
                    ii(3)=i ; jj(3)=j

                    kk=0

                ELSE ! ENDOF IF(NDIM==2) THEN...(THIS IS FOR 3D)

                    ! NB. The cross product is directly related to the area of the triangle...
                    iSTORE(1)=1 ; jSTORE(1)=2 ; kSTORE(1)=3
                    iSTORE(2)=1 ; jSTORE(2)=2 ; kSTORE(2)=4
                    iSTORE(3)=2 ; jSTORE(3)=3 ; kSTORE(3)=4
                    iSTORE(4)=1 ; jSTORE(4)=3 ; kSTORE(4)=4

                    IFACE = 0
                    DO III = 1, 4
                        i=iSTORE(III) ; j=jSTORE(III) ; k=kSTORE(III)
                        IF ( PIPE_INDEX_LOGICAL(I) .AND. PIPE_INDEX_LOGICAL(J) .AND. PIPE_INDEX_LOGICAL(K) ) THEN
                            IFACE = IFACE+1
                            call CrossProduct( ndim, cp,  ( X_ALL_CORN(:,j) - X_ALL_CORN(:,i) ), ( X_ALL_CORN(:,k) - X_ALL_CORN(:,i) ) )
                            area_sqr(IFACE) =  sum(cp**2 )
                            ii(IFACE)=i ; jj(IFACE)=j ; kk(IFACE)=k
                        END IF
                    END DO ! ENDOF DO III=1,4

                    IF(IFACE==1) THEN ! A 3d adjustment for when we have only 2 pipes...
                        i=ii(IFACE) ; j=jj(IFACE) ; k=kk(IFACE) ! nodes on the face

                        area_sqr(1) = sum( ( X_ALL_CORN(:,i) - X_ALL_CORN(:,j) )**2 ) ! THIS IS THE LENGTH SQUARED

                        area_sqr(2) = sum( ( X_ALL_CORN(:,i) - X_ALL_CORN(:,k) )**2 ) ! THIS IS THE LENGTH SQUARED

                        area_sqr(3) = sum( ( X_ALL_CORN(:,j) - X_ALL_CORN(:,k) )**2 ) ! THIS IS THE LENGTH SQUARED

                        ik = 1
                        max_area = area_sqr(ik)
                        do il = 2, 3
                            if ( max_area < area_sqr(il) ) then
                                max_area = area_sqr(il)
                                ik = il
                            end if
                        end do

                        if(ik==1) then
                            ipipe = 1
                            pipe_corner_nds1(ipipe) = i
                            pipe_corner_nds2(ipipe) = k
                            ipipe = 2
                            pipe_corner_nds1(ipipe) = j
                            pipe_corner_nds2(ipipe) = k
                        else if(ik==2) then
                            ipipe = 1
                            pipe_corner_nds1(ipipe) = i
                            pipe_corner_nds2(ipipe) = j
                            ipipe = 2
                            pipe_corner_nds1(ipipe) = j
                            pipe_corner_nds2(ipipe) = k
                        else
                            ipipe = 1
                            pipe_corner_nds1(ipipe) = i
                            pipe_corner_nds2(ipipe) = j
                            ipipe = 2
                            pipe_corner_nds1(ipipe) = i
                            pipe_corner_nds2(ipipe) = k
                        endif
                        RETURN

                    ENDIF  ! IF(IFACE==1) THEN

                END IF ! ENDOF IF(NDIM==2) THEN

                ! Find the maximum surface area..
                ik = 1
                max_area = area_sqr(ik)
                do i = 2, NDIM+1
                    if ( max_area < area_sqr(i) ) then
                        max_area = area_sqr(i)
                        ik = i
                    end if
                end do

                ! Now work out the pipes that do not contain the edges on this maximum surface
                ! and contain
                ipipe = 0
                do i_nd = 1, ndim+1
                    do j_nd = i_nd+1, ndim+1
                        ! If edge nodes i_nd with j_nd are both in the ii,jj or kk
                        IF ( PIPE_INDEX_LOGICAL(i_nd) .AND. PIPE_INDEX_LOGICAL(J_ND) ) THEN
                            ! Make sure both nodes do not belong to the face with the maximum area...
                            found_i = .false.
                            found_j = .false.

                            if ( i_nd==ii(ik) ) found_i = .true.
                            if ( i_nd==jj(ik) ) found_i = .true.
                            if ( i_nd==kk(ik) ) found_i = .true.

                            if ( j_nd==ii(ik) ) found_j = .true.
                            if ( j_nd==jj(ik) ) found_j = .true.
                            if ( j_nd==kk(ik) ) found_j = .true.

                            if ( (.not.found_i) .and. (.not.found_j) ) then
                                ipipe = ipipe + 1
                                pipe_corner_nds1(ipipe) = i_nd
                                pipe_corner_nds2(ipipe) = j_nd
                            end if

                        end if
                    end do
                end do

            END IF ! IF ( NPIPES==1 ) THEN ELSE

            RETURN
        END SUBROUTINE CALC_PIPES_IN_ELE

    end subroutine retrieve_pipes_coords

    subroutine initialize_pipes_package_and_gamma(state, pipes_aux, Mdims, Mspars)
        implicit none
        type(state_type), dimension(:), intent(in) :: state
        type (multi_pipe_package), intent(inout) :: pipes_aux
        type (multi_dimensions), intent(in)  ::Mdims
        type (multi_sparsities), intent(in) :: Mspars
         !Local variables
        type( scalar_field ), pointer :: sfield
        !Variables to initialize pipes_aux
        integer :: ipres, iphase, jpres, jphase, iphase_real, jphase_real
        !Initialize memory, despite we call it, if it is already allocated no memory is re-allocated
        call allocate_multi_pipe_package(pipes_aux, Mdims, Mspars)

        !Initialize gamma
        sfield=>extract_scalar_field(state(1),"Gamma",ipres)
        pipes_aux%GAMMA_PRES_ABS = 0.0
        do ipres = 1, Mdims%npres
           do iphase = 1+(ipres-1)*Mdims%n_in_pres, ipres*Mdims%n_in_pres
              do jpres = 1, Mdims%npres
                 if ( ipres /= jpres ) then
                    do jphase = 1+(jpres-1)*Mdims%n_in_pres, jpres*Mdims%n_in_pres
                       iphase_real = iphase-(ipres-1)*Mdims%n_in_pres
                       jphase_real = jphase-(jpres-1)*Mdims%n_in_pres
                       if ( iphase_real == jphase_real ) then
                          call assign_val(pipes_aux%GAMMA_PRES_ABS(IPHASE,JPHASE,:),sfield%val)
                       end if
                    end do
                 end if
              end do
           end do
        end do
        pipes_aux%GAMMA_PRES_ABS_NANO = pipes_aux%GAMMA_PRES_ABS
    end subroutine initialize_pipes_package_and_gamma

end module multi_pipes
