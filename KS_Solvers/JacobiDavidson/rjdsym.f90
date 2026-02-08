!
! Copyright (C) 2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE rjdsym( h_psi_ptr, s_psi_ptr, uspp, &
                   npw, npwx, nvec, nvecx, evc, ethr, &
                   g2kin, e, btype, notcnv, jd_iter, nhpsi )
  !----------------------------------------------------------------------------
  !
  ! ... Jacobi-Davidson iterative diagonalization of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is a real symmetric operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector
  ! ... (real wavefunctions with only half plane waves stored).
  !
  ! ... Gamma-point version: uses real BLAS with 2*npw trick.
  ! ... Processes one eigenvalue at a time with explicit deflation.
  ! ... Uses TPA preconditioner with default shift: t = r / (g2kin + 1).
  !
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, &
                            nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp_bands_util, ONLY : gstart
  USE mp,            ONLY : mp_sum, mp_bcast
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx,nvec)
    !  evc contains the refined estimates of the eigenvectors
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : do not calculate S|psi>
  REAL(DP), INTENT(IN) :: g2kin(npwx)
    ! kinetic energy of each G-vector, used for preconditioning
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: jd_iter, notcnv
    ! integer number of iterations performed
    ! number of unconverged roots
  INTEGER, INTENT(OUT) :: nhpsi
    ! total number of individual hpsi
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 200
    ! maximum number of iterations
  REAL(DP), PARAMETER :: default_shift = 1.0_DP
    ! TPA default shift for preconditioner
  !
  INTEGER :: j, nconv, iter, npw2, npwx2, ierr, jmin
    ! current subspace dimension
    ! number of converged eigenvalues
    ! iteration counter
    ! double-length dimensions for real BLAS
  INTEGER :: i, ig
    ! loop counters
  REAL(DP) :: nr, nt, tol, empty_ethr
    ! residual norm
    ! correction norm
    ! convergence tolerance
    ! threshold for empty bands
  REAL(DP) :: coeff_u_r, rdot
    ! projection coefficient (real)
    ! dot product
  !
  COMPLEX(DP), ALLOCATABLE :: V(:,:), W(:,:), SW(:,:)
    ! search space basis vectors
    ! H * V
    ! S * V (if uspp)
  REAL(DP), ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:)
    ! Hamiltonian on the reduced basis (real for gamma)
    ! S matrix on the reduced basis (real for gamma)
    ! eigenvectors of the reduced problem (real for gamma)
  REAL(DP), ALLOCATABLE :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: u(:), Au(:), Su(:), r(:), t(:)
    ! Ritz vector
    ! H * Ritz vector
    ! S * Ritz vector
    ! residual
    ! correction vector
  REAL(DP), ALLOCATABLE :: rwork(:)
    ! real work array for projections
  COMPLEX(DP), ALLOCATABLE :: Vtmp(:,:), Wtmp(:,:), SWtmp(:,:)
    ! temporary arrays for restart/deflation
  !
  EXTERNAL  h_psi_ptr, s_psi_ptr
    ! h_psi_ptr(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi_ptr(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
  !
  nhpsi = 0
  CALL start_clock( 'rjdsym' )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'rjdsym', 'nvecx is too small', 1 )
  IF ( gstart == -1 ) CALL errore( 'rjdsym', 'gstart variable not initialized', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  ! ... convergence tolerance (scaled as in Julia reference)
  !
  tol = ethr / SQRT( DBLE(nvec) )
  !
  npw2  = 2*npw
  npwx2 = 2*npwx
  !
  ! ... jmin for restart
  !
  jmin = MAX( 5, nvec / 2 )
  IF ( jmin > nvecx / 2 ) jmin = nvecx / 2
  !
  ! ... Allocate search space and projected matrices
  !
  ALLOCATE( V( npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'rjdsym', 'cannot allocate V', ABS(ierr) )
  ALLOCATE( W( npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'rjdsym', 'cannot allocate W', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( SW( npwx, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( 'rjdsym', 'cannot allocate SW', ABS(ierr) )
  END IF
  !
  ! ... Reduced matrices are REAL for gamma-point
  !
  ALLOCATE( hr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'rjdsym', 'cannot allocate hr', ABS(ierr) )
  ALLOCATE( sr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'rjdsym', 'cannot allocate sr', ABS(ierr) )
  ALLOCATE( vr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'rjdsym', 'cannot allocate vr', ABS(ierr) )
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'rjdsym', 'cannot allocate ew', ABS(ierr) )
  !
  ALLOCATE( u( npwx ), Au( npwx ), Su( npwx ) )
  ALLOCATE( r( npwx ), t( npwx ) )
  ALLOCATE( rwork( nvecx ) )
  ALLOCATE( Vtmp( npwx, nvecx ) )
  ALLOCATE( Wtmp( npwx, nvecx ) )
  IF ( uspp ) ALLOCATE( SWtmp( npwx, nvecx ) )
  !
  ! ... Initialize
  !
  nconv = 0
  j = nvec
  hr = 0.0_DP
  sr = 0.0_DP
  vr = 0.0_DP
  ew = 0.0_DP
  !
  ! ... Copy initial vectors to search space
  !
  V(:,:) = ZERO
  W(:,:) = ZERO
  V(1:npwx, 1:nvec) = evc(1:npwx, 1:nvec)
  !
  ! ... Set Im[psi(G=0)] = 0 for numerical stability
  !
  IF ( gstart == 2 ) THEN
     DO i = 1, nvec
        V(1,i) = CMPLX( DBLE( V(1,i) ), 0.0_DP, kind=DP )
     END DO
  END IF
  !
  ! ... Orthogonalize initial search space (modified Gram-Schmidt)
  ! ... This ensures the projected overlap matrix is positive definite,
  ! ... which is critical for the first SCF iteration when initial
  ! ... wavefunctions may not be orthonormal.
  !
  DO i = 1, nvec
     !
     ! ... Orthogonalize V(:,i) against V(:,1:i-1)
     !
     IF ( i > 1 ) THEN
        CALL DGEMM( 'T', 'N', i-1, 1, npw2, 2.0_DP, V, npwx2, &
                    V(1,i), npwx2, 0.0_DP, rwork, i-1 )
        IF ( gstart == 2 ) CALL MYDGER( i-1, 1, -1.0_DP, V, npwx2, &
                                         V(1,i), npwx2, rwork, i-1 )
        CALL mp_sum( rwork(1:i-1), intra_bgrp_comm )
        CALL DGEMM( 'N', 'N', npw2, 1, i-1, -1.0_DP, V, npwx2, &
                    rwork, i-1, 1.0_DP, V(1,i), npwx2 )
        ! ... Repeat for numerical stability
        CALL DGEMM( 'T', 'N', i-1, 1, npw2, 2.0_DP, V, npwx2, &
                    V(1,i), npwx2, 0.0_DP, rwork, i-1 )
        IF ( gstart == 2 ) CALL MYDGER( i-1, 1, -1.0_DP, V, npwx2, &
                                         V(1,i), npwx2, rwork, i-1 )
        CALL mp_sum( rwork(1:i-1), intra_bgrp_comm )
        CALL DGEMM( 'N', 'N', npw2, 1, i-1, -1.0_DP, V, npwx2, &
                    rwork, i-1, 1.0_DP, V(1,i), npwx2 )
     END IF
     !
     ! ... Normalize
     !
     nr = 0.0_DP
     CALL DGEMM( 'T', 'N', 1, 1, npw2, 2.0_DP, V(1,i), npwx2, &
                 V(1,i), npwx2, 0.0_DP, nr, 1 )
     IF ( gstart == 2 ) nr = nr - DBLE( CONJG(V(1,i)) * V(1,i) )
     CALL mp_sum( nr, intra_bgrp_comm )
     nr = SQRT( nr )
     IF ( nr > 1.0D-14 ) THEN
        V(1:npw,i) = V(1:npw,i) / nr
     END IF
     !
     ! ... Zero padding and fix G=0
     !
     IF ( npw < npwx ) V(npw+1:npwx,i) = ZERO
     IF ( gstart == 2 ) V(1,i) = CMPLX( DBLE(V(1,i)), 0.0_DP, kind=DP )
     !
  END DO
  !
  ! ... Compute H*V and S*V for initial vectors
  !
  CALL h_psi_ptr( npwx, npw, nvec, V, W )
  nhpsi = nhpsi + nvec
  !
  IF ( uspp ) THEN
     SW(:,:) = ZERO
     CALL s_psi_ptr( npwx, npw, nvec, V, SW )
  END IF
  !
  ! ... Build initial projected Hamiltonian hr = V^T * W (real, using gamma trick)
  ! ... <a|b> = 2 * Re(a^H * b) - Re(a(G=0)^* * b(G=0))
  !
  CALL DGEMM( 'T', 'N', j, j, npw2, 2.0_DP, V, npwx2, W, npwx2, &
              0.0_DP, hr, nvecx )
  IF ( gstart == 2 ) CALL MYDGER( j, j, -1.0_DP, V, npwx2, W, npwx2, hr, nvecx )
  CALL mp_sum( hr(1:j, 1:j), intra_bgrp_comm )
  !
  ! ... Build initial projected overlap sr
  !
  IF ( uspp ) THEN
     CALL DGEMM( 'T', 'N', j, j, npw2, 2.0_DP, V, npwx2, SW, npwx2, &
                 0.0_DP, sr, nvecx )
     IF ( gstart == 2 ) CALL MYDGER( j, j, -1.0_DP, V, npwx2, SW, npwx2, sr, nvecx )
  ELSE
     CALL DGEMM( 'T', 'N', j, j, npw2, 2.0_DP, V, npwx2, V, npwx2, &
                 0.0_DP, sr, nvecx )
     IF ( gstart == 2 ) CALL MYDGER( j, j, -1.0_DP, V, npwx2, V, npwx2, sr, nvecx )
  END IF
  CALL mp_sum( sr(1:j, 1:j), intra_bgrp_comm )
  !
  ! ... Initialize eigenvalue estimates
  !
  e = 0.0_DP
  !
  ! ====================================================================
  ! ... Main Jacobi-Davidson loop
  ! ====================================================================
  !
  iterate: DO iter = 1, maxter
     !
     jd_iter = iter
     !
     ! ... Diagonalize the projected Hamiltonian
     !
     CALL start_clock( 'rjdsym:diag' )
     IF ( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( j, MIN(nvec - nconv, j), hr, sr, nvecx, ew, vr, &
                     me_bgrp, root_bgrp, intra_bgrp_comm )
     END IF
     IF ( nbgrp > 1 ) THEN
        CALL mp_bcast( vr, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     END IF
     CALL stop_clock( 'rjdsym:diag' )
     !
     ! ... Extract best Ritz pair (smallest eigenvalue)
     ! ... u = V * vr(:,1), Au = W * vr(:,1), Su = SW * vr(:,1)
     !
     u = ZERO
     CALL DGEMM( 'N', 'N', npw2, 1, j, 1.0_DP, V, npwx2, &
                 vr(1,1), nvecx, 0.0_DP, u, npwx2 )
     !
     Au = ZERO
     CALL DGEMM( 'N', 'N', npw2, 1, j, 1.0_DP, W, npwx2, &
                 vr(1,1), nvecx, 0.0_DP, Au, npwx2 )
     !
     IF ( uspp ) THEN
        Su = ZERO
        CALL DGEMM( 'N', 'N', npw2, 1, j, 1.0_DP, SW, npwx2, &
                    vr(1,1), nvecx, 0.0_DP, Su, npwx2 )
     ELSE
        Su = u
     END IF
     !
     ! ... Compute residual: r = Au - theta * Su
     !
     r(1:npw) = Au(1:npw) - ew(1) * Su(1:npw)
     !
     ! ... Set Im[r(G=0)] = 0 for numerical stability
     !
     IF ( gstart == 2 ) r(1) = CMPLX( DBLE(r(1)), 0.0_DP, kind=DP )
     !
     ! ... Orthogonalize residual against converged eigenvectors (double for stability)
     !
     IF ( nconv > 0 ) THEN
        DO i = 1, 2
           CALL DGEMM( 'T', 'N', nconv, 1, npw2, 2.0_DP, evc, npwx2, &
                       r, npwx2, 0.0_DP, rwork, nconv )
           IF ( gstart == 2 ) CALL MYDGER( nconv, 1, -1.0_DP, evc, npwx2, &
                                            r, npwx2, rwork, nconv )
           CALL mp_sum( rwork(1:nconv), intra_bgrp_comm )
           CALL DGEMM( 'N', 'N', npw2, 1, nconv, -1.0_DP, evc, npwx2, &
                       rwork, nconv, 1.0_DP, r, npwx2 )
        END DO
     END IF
     !
     ! ... Compute residual norm: ||r||^2 = 2*sum|r_G|^2 - |r_0|^2
     !
     nr = 0.0_DP
     CALL DGEMM( 'T', 'N', 1, 1, npw2, 2.0_DP, r, npwx2, &
                 r, npwx2, 0.0_DP, nr, 1 )
     IF ( gstart == 2 ) nr = nr - DBLE( CONJG(r(1)) * r(1) )
     CALL mp_sum( nr, intra_bgrp_comm )
     nr = SQRT( nr )
     !
     ! ... Check convergence
     !
     IF ( btype(nconv+1) == 1 ) THEN
        IF ( nr < tol ) THEN
           !
           ! ... Converged! Store eigenpair
           !
           nconv = nconv + 1
           evc(1:npwx, nconv) = u(1:npwx)
           e(nconv) = ew(1)
           !
           IF ( nconv >= nvec ) EXIT iterate
           !
           ! ... Deflate: remove converged component from subspace
           !
           IF ( j > 1 ) THEN
              !
              CALL DGEMM( 'N', 'N', npw2, j-1, j, 1.0_DP, V, npwx2, &
                          vr(1,2), nvecx, 0.0_DP, Vtmp, npwx2 )
              V(1:npwx, 1:j-1) = Vtmp(1:npwx, 1:j-1)
              !
              CALL DGEMM( 'N', 'N', npw2, j-1, j, 1.0_DP, W, npwx2, &
                          vr(1,2), nvecx, 0.0_DP, Wtmp, npwx2 )
              W(1:npwx, 1:j-1) = Wtmp(1:npwx, 1:j-1)
              !
              IF ( uspp ) THEN
                 CALL DGEMM( 'N', 'N', npw2, j-1, j, 1.0_DP, SW, npwx2, &
                             vr(1,2), nvecx, 0.0_DP, SWtmp, npwx2 )
                 SW(1:npwx, 1:j-1) = SWtmp(1:npwx, 1:j-1)
              END IF
              !
              j = j - 1
              !
              hr = 0.0_DP
              sr = 0.0_DP
              DO i = 1, j
                 hr(i,i) = ew(i+1)
                 sr(i,i) = 1.0_DP
              END DO
              !
           ELSE
              CALL rjd_reinit_subspace()
           END IF
           !
           CYCLE iterate
           !
        END IF
     ELSE
        ! ... empty band: relaxed threshold
        IF ( nr < empty_ethr / SQRT( DBLE(nvec) ) ) THEN
           !
           nconv = nconv + 1
           evc(1:npwx, nconv) = u(1:npwx)
           e(nconv) = ew(1)
           !
           IF ( nconv >= nvec ) EXIT iterate
           !
           IF ( j > 1 ) THEN
              CALL DGEMM( 'N', 'N', npw2, j-1, j, 1.0_DP, V, npwx2, &
                          vr(1,2), nvecx, 0.0_DP, Vtmp, npwx2 )
              V(1:npwx, 1:j-1) = Vtmp(1:npwx, 1:j-1)
              !
              CALL DGEMM( 'N', 'N', npw2, j-1, j, 1.0_DP, W, npwx2, &
                          vr(1,2), nvecx, 0.0_DP, Wtmp, npwx2 )
              W(1:npwx, 1:j-1) = Wtmp(1:npwx, 1:j-1)
              !
              IF ( uspp ) THEN
                 CALL DGEMM( 'N', 'N', npw2, j-1, j, 1.0_DP, SW, npwx2, &
                             vr(1,2), nvecx, 0.0_DP, SWtmp, npwx2 )
                 SW(1:npwx, 1:j-1) = SWtmp(1:npwx, 1:j-1)
              END IF
              !
              j = j - 1
              hr = 0.0_DP
              sr = 0.0_DP
              DO i = 1, j
                 hr(i,i) = ew(i+1)
                 sr(i,i) = 1.0_DP
              END DO
           ELSE
              CALL rjd_reinit_subspace()
           END IF
           !
           CYCLE iterate
        END IF
     END IF
     !
     ! ... Restart if subspace is too large
     !
     IF ( j >= nvecx ) THEN
        !
        CALL start_clock( 'rjdsym:restart' )
        !
        CALL DGEMM( 'N', 'N', npw2, jmin, j, 1.0_DP, V, npwx2, &
                    vr(1,1), nvecx, 0.0_DP, Vtmp, npwx2 )
        V(1:npwx, 1:jmin) = Vtmp(1:npwx, 1:jmin)
        !
        CALL DGEMM( 'N', 'N', npw2, jmin, j, 1.0_DP, W, npwx2, &
                    vr(1,1), nvecx, 0.0_DP, Wtmp, npwx2 )
        W(1:npwx, 1:jmin) = Wtmp(1:npwx, 1:jmin)
        !
        IF ( uspp ) THEN
           CALL DGEMM( 'N', 'N', npw2, jmin, j, 1.0_DP, SW, npwx2, &
                       vr(1,1), nvecx, 0.0_DP, SWtmp, npwx2 )
           SW(1:npwx, 1:jmin) = SWtmp(1:npwx, 1:jmin)
        END IF
        !
        j = jmin
        !
        hr = 0.0_DP
        sr = 0.0_DP
        DO i = 1, j
           hr(i,i) = ew(i)
           sr(i,i) = 1.0_DP
        END DO
        !
        CALL stop_clock( 'rjdsym:restart' )
        !
     END IF
     !
     ! ================================================================
     ! ... Solve correction equation (simplified Jacobi-Davidson)
     ! ... t = M^{-1} * r, then project out [converged, u]
     ! ================================================================
     !
     CALL start_clock( 'rjdsym:correction' )
     !
     ! ... Apply TPA preconditioner with default shift: t = r / (g2kin + 1)
     !
     t = ZERO
     DO ig = 1, npw
        t(ig) = r(ig) / ( g2kin(ig) + default_shift )
     END DO
     !
     ! ... Set Im[t(G=0)] = 0 for numerical stability
     !
     IF ( gstart == 2 ) t(1) = CMPLX( DBLE(t(1)), 0.0_DP, kind=DP )
     !
     ! ... Project out converged eigenvectors (double for stability)
     !
     IF ( nconv > 0 ) THEN
        DO i = 1, 2
           CALL DGEMM( 'T', 'N', nconv, 1, npw2, 2.0_DP, evc, npwx2, &
                       t, npwx2, 0.0_DP, rwork, nconv )
           IF ( gstart == 2 ) CALL MYDGER( nconv, 1, -1.0_DP, evc, npwx2, &
                                            t, npwx2, rwork, nconv )
           CALL mp_sum( rwork(1:nconv), intra_bgrp_comm )
           CALL DGEMM( 'N', 'N', npw2, 1, nconv, -1.0_DP, evc, npwx2, &
                       rwork, nconv, 1.0_DP, t, npwx2 )
        END DO
     END IF
     !
     ! ... Project out current Ritz vector u
     !
     ! ... dot_ut = <u|t>, dot_uu = <u|u>
     !
     rdot = 0.0_DP
     CALL DGEMM( 'T', 'N', 1, 1, npw2, 2.0_DP, u, npwx2, &
                 t, npwx2, 0.0_DP, rdot, 1 )
     IF ( gstart == 2 ) rdot = rdot - DBLE( CONJG(u(1)) * t(1) )
     CALL mp_sum( rdot, intra_bgrp_comm )
     !
     coeff_u_r = 0.0_DP
     CALL DGEMM( 'T', 'N', 1, 1, npw2, 2.0_DP, u, npwx2, &
                 u, npwx2, 0.0_DP, coeff_u_r, 1 )
     IF ( gstart == 2 ) coeff_u_r = coeff_u_r - DBLE( CONJG(u(1)) * u(1) )
     CALL mp_sum( coeff_u_r, intra_bgrp_comm )
     !
     IF ( coeff_u_r > 1.0D-30 ) THEN
        coeff_u_r = rdot / coeff_u_r
        t(1:npw) = t(1:npw) - coeff_u_r * u(1:npw)
     END IF
     !
     CALL stop_clock( 'rjdsym:correction' )
     !
     ! ... Orthogonalize t against V (double Gram-Schmidt for stability)
     !
     CALL start_clock( 'rjdsym:ortho' )
     !
     DO i = 1, 2
        CALL DGEMM( 'T', 'N', j, 1, npw2, 2.0_DP, V, npwx2, &
                    t, npwx2, 0.0_DP, rwork, j )
        IF ( gstart == 2 ) CALL MYDGER( j, 1, -1.0_DP, V, npwx2, &
                                         t, npwx2, rwork, j )
        CALL mp_sum( rwork(1:j), intra_bgrp_comm )
        CALL DGEMM( 'N', 'N', npw2, 1, j, -1.0_DP, V, npwx2, &
                    rwork, j, 1.0_DP, t, npwx2 )
     END DO
     !
     ! ... Normalize t
     !
     nt = 0.0_DP
     CALL DGEMM( 'T', 'N', 1, 1, npw2, 2.0_DP, t, npwx2, &
                 t, npwx2, 0.0_DP, nt, 1 )
     IF ( gstart == 2 ) nt = nt - DBLE( CONJG(t(1)) * t(1) )
     CALL mp_sum( nt, intra_bgrp_comm )
     nt = SQRT( nt )
     !
     CALL stop_clock( 'rjdsym:ortho' )
     !
     IF ( nt < 1.0D-14 ) CYCLE iterate    ! correction too small, skip
     !
     t(1:npw) = t(1:npw) / nt
     ! ... zero out padding
     IF ( npw < npwx ) t(npw+1:npwx) = ZERO
     ! ... Set Im[t(G=0)] = 0
     IF ( gstart == 2 ) t(1) = CMPLX( DBLE(t(1)), 0.0_DP, kind=DP )
     !
     ! ================================================================
     ! ... Expand subspace
     ! ================================================================
     !
     j = j + 1
     V(1:npwx, j) = t(1:npwx)
     !
     ! ... Compute H*t and S*t
     !
     CALL h_psi_ptr( npwx, npw, 1, V(1,j), W(1,j) )
     nhpsi = nhpsi + 1
     !
     IF ( uspp ) CALL s_psi_ptr( npwx, npw, 1, V(1,j), SW(1,j) )
     !
     ! ... Update projected Hamiltonian: new column hr(1:j, j)
     !
     CALL start_clock( 'rjdsym:overlap' )
     !
     CALL DGEMM( 'T', 'N', j, 1, npw2, 2.0_DP, V, npwx2, &
                 W(1,j), npwx2, 0.0_DP, hr(1,j), nvecx )
     IF ( gstart == 2 ) CALL MYDGER( j, 1, -1.0_DP, V, npwx2, &
                                      W(1,j), npwx2, hr(1,j), nvecx )
     CALL mp_sum( hr(1:j, j), intra_bgrp_comm )
     !
     ! ... Symmetry: hr(j, 1:j-1) = hr(1:j-1, j)
     !
     DO i = 1, j - 1
        hr(j,i) = hr(i,j)
     END DO
     !
     ! ... Update projected overlap: new column sr(1:j, j)
     !
     IF ( uspp ) THEN
        CALL DGEMM( 'T', 'N', j, 1, npw2, 2.0_DP, V, npwx2, &
                    SW(1,j), npwx2, 0.0_DP, sr(1,j), nvecx )
        IF ( gstart == 2 ) CALL MYDGER( j, 1, -1.0_DP, V, npwx2, &
                                         SW(1,j), npwx2, sr(1,j), nvecx )
     ELSE
        CALL DGEMM( 'T', 'N', j, 1, npw2, 2.0_DP, V, npwx2, &
                    V(1,j), npwx2, 0.0_DP, sr(1,j), nvecx )
        IF ( gstart == 2 ) CALL MYDGER( j, 1, -1.0_DP, V, npwx2, &
                                         V(1,j), npwx2, sr(1,j), nvecx )
     END IF
     CALL mp_sum( sr(1:j, j), intra_bgrp_comm )
     !
     DO i = 1, j - 1
        sr(j,i) = sr(i,j)
     END DO
     !
     CALL stop_clock( 'rjdsym:overlap' )
     !
  END DO iterate
  !
  ! ... Set output
  !
  notcnv = nvec - nconv
  !
  ! ... For any remaining unconverged eigenvalues, use best Ritz approximation
  !
  IF ( notcnv > 0 .AND. j > 0 ) THEN
     DO i = 1, MIN(notcnv, j)
        CALL DGEMM( 'N', 'N', npw2, 1, j, 1.0_DP, V, npwx2, &
                    vr(1,i), nvecx, 0.0_DP, evc(1,nconv+i), npwx2 )
        e(nconv+i) = ew(i)
     END DO
     DO i = j + 1, notcnv
        e(nconv+i) = ew(MIN(i,j))
     END DO
  END IF
  !
  ! ... Deallocate
  !
  IF ( uspp ) DEALLOCATE( SWtmp )
  DEALLOCATE( Wtmp )
  DEALLOCATE( Vtmp )
  DEALLOCATE( rwork )
  DEALLOCATE( t, r, Su, Au, u )
  DEALLOCATE( ew )
  DEALLOCATE( vr, sr, hr )
  IF ( uspp ) DEALLOCATE( SW )
  DEALLOCATE( W, V )
  !
  CALL stop_clock( 'rjdsym' )
  !
  RETURN
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE rjd_reinit_subspace()
    !-----------------------------------------------------------------------
    !
    ! ... Reinitialize the search space when it becomes empty after deflation.
    ! ... Creates a random vector orthogonal to converged eigenvectors.
    !
    IMPLICIT NONE
    INTEGER :: ig2
    REAL(DP) :: rr, ri
    !
    j = 1
    V(1:npwx, 1) = ZERO
    DO ig2 = 1, npw
       rr = DBLE(MOD(ig2 + nconv*137, 1000)) / 1000.0_DP
       ri = DBLE(MOD(ig2 + nconv*251, 1000)) / 1000.0_DP
       V(ig2, 1) = CMPLX( rr, ri, kind=DP )
    END DO
    IF ( gstart == 2 ) V(1,1) = CMPLX( DBLE(V(1,1)), 0.0_DP, kind=DP )
    !
    ! ... Orthogonalize against converged eigenvectors (double)
    !
    DO i = 1, 2
       CALL DGEMM( 'T', 'N', nconv, 1, npw2, 2.0_DP, evc, npwx2, &
                   V(1,1), npwx2, 0.0_DP, rwork, nconv )
       IF ( gstart == 2 ) CALL MYDGER( nconv, 1, -1.0_DP, evc, npwx2, &
                                        V(1,1), npwx2, rwork, nconv )
       CALL mp_sum( rwork(1:nconv), intra_bgrp_comm )
       CALL DGEMM( 'N', 'N', npw2, 1, nconv, -1.0_DP, evc, npwx2, &
                   rwork, nconv, 1.0_DP, V(1,1), npwx2 )
    END DO
    !
    ! ... Normalize
    !
    nt = 0.0_DP
    CALL DGEMM( 'T', 'N', 1, 1, npw2, 2.0_DP, V(1,1), npwx2, &
                V(1,1), npwx2, 0.0_DP, nt, 1 )
    IF ( gstart == 2 ) nt = nt - DBLE( CONJG(V(1,1)) * V(1,1) )
    CALL mp_sum( nt, intra_bgrp_comm )
    nt = SQRT( nt )
    IF ( nt > 1.0D-14 ) THEN
       V(1:npw, 1) = V(1:npw, 1) / nt
    END IF
    IF ( npw < npwx ) V(npw+1:npwx, 1) = ZERO
    IF ( gstart == 2 ) V(1,1) = CMPLX( DBLE(V(1,1)), 0.0_DP, kind=DP )
    !
    ! ... Compute H*V and S*V
    !
    CALL h_psi_ptr( npwx, npw, 1, V(1,1), W(1,1) )
    nhpsi = nhpsi + 1
    IF ( uspp ) CALL s_psi_ptr( npwx, npw, 1, V(1,1), SW(1,1) )
    !
    ! ... Build projected matrices (1x1)
    !
    hr = 0.0_DP
    sr = 0.0_DP
    CALL DGEMM( 'T', 'N', 1, 1, npw2, 2.0_DP, V(1,1), npwx2, &
                W(1,1), npwx2, 0.0_DP, hr(1,1), nvecx )
    IF ( gstart == 2 ) hr(1,1) = hr(1,1) - DBLE( CONJG(V(1,1)) * W(1,1) )
    CALL mp_sum( hr(1:1, 1:1), intra_bgrp_comm )
    !
    IF ( uspp ) THEN
       CALL DGEMM( 'T', 'N', 1, 1, npw2, 2.0_DP, V(1,1), npwx2, &
                   SW(1,1), npwx2, 0.0_DP, sr(1,1), nvecx )
       IF ( gstart == 2 ) sr(1,1) = sr(1,1) - DBLE( CONJG(V(1,1)) * SW(1,1) )
    ELSE
       CALL DGEMM( 'T', 'N', 1, 1, npw2, 2.0_DP, V(1,1), npwx2, &
                   V(1,1), npwx2, 0.0_DP, sr(1,1), nvecx )
       IF ( gstart == 2 ) sr(1,1) = sr(1,1) - DBLE( CONJG(V(1,1)) * V(1,1) )
    END IF
    CALL mp_sum( sr(1:1, 1:1), intra_bgrp_comm )
    !
  END SUBROUTINE rjd_reinit_subspace
  !
END SUBROUTINE rjdsym
