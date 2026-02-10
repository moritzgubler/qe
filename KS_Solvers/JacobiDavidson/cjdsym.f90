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
SUBROUTINE cjdsym( h_psi_ptr, s_psi_ptr, uspp, g_psi_ptr, &
                   npw, npwx, nvec, nvecx, npol, evc, ethr, &
                   g2kin, e, btype, notcnv, jd_iter, nhpsi )
  !----------------------------------------------------------------------------
  !
  ! ... Blocked Jacobi-Davidson iterative diagonalization:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... In each iteration, nblock correction vectors are computed and
  ! ... added to the search space for better BLAS3 utilization.
  ! ... Preconditioner: g_psi_ptr (use_g_psi=T) or TPA (use_g_psi=F).
  !
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, &
                            nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp,            ONLY : mp_sum, mp_bcast
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx*npol,nvec)
  REAL(DP), INTENT(IN) :: ethr
  LOGICAL, INTENT(IN) :: uspp
  REAL(DP), INTENT(IN) :: g2kin(npwx)
  INTEGER, INTENT(IN) :: btype(nvec)
  REAL(DP), INTENT(OUT) :: e(nvec)
  INTEGER, INTENT(OUT) :: jd_iter, notcnv
  INTEGER, INTENT(OUT) :: nhpsi
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
  INTEGER :: nblock
    ! number of correction vectors per iteration
  REAL(DP), PARAMETER :: default_shift = 1.0_DP
  LOGICAL, PARAMETER :: use_g_psi = .TRUE.
  !
  INTEGER :: j, nconv, iter, kdim, kdmx, ierr, jmin
  INTEGER :: i, ig, ipol, ib, nb, nact, nconv_new, k
  REAL(DP) :: tol, empty_ethr, norm_t
  LOGICAL :: lprint
  REAL(DP) :: rnorms(nvec)
  !
  COMPLEX(DP), ALLOCATABLE :: V(:,:), W(:,:), SW(:,:)
    ! search space basis vectors / H * V / S * V
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
    ! projected Hamiltonian / overlap / eigenvectors
  REAL(DP), ALLOCATABLE :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: ub(:,:), rb(:,:), tb(:,:)
    ! Ritz vectors / residuals (then corrections) / scratch
  COMPLEX(DP), ALLOCATABLE :: work2d(:,:)
    ! workspace for projections
  COMPLEX(DP), ALLOCATABLE :: Vtmp(:,:), Wtmp(:,:), SWtmp(:,:)
    ! temporary arrays for restart/deflation
  !
  EXTERNAL  h_psi_ptr, s_psi_ptr, g_psi_ptr
    ! h_psi_ptr(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi_ptr(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    ! g_psi_ptr(npwx,npw,notcnv,npol,psi,e)
    !     calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
  !
  nblock = nvec
  nhpsi = 0
  lprint = .FALSE.
  CALL start_clock( 'cjdsym' )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'cjdsym', 'nvecx is too small', 1 )

  print*, npw, npwx, nvec, nvecx
  !
  empty_ethr = sqrt(MAX( ( ethr * 5.D0 ), 1.D-5 ))
  tol = sqrt(ethr)
  !
  IF ( npol == 1 ) THEN
     kdim = npw
     kdmx = npwx
  ELSE
     kdim = npwx*npol
     kdmx = npwx*npol
  END IF
  !
  jmin = MAX( npw, nvec + 5 )
  IF ( jmin > nvecx / 2 ) jmin = nvecx / 2
  !
  ! ... Allocate workspace
  !
  ALLOCATE( V( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'cjdsym', 'cannot allocate V', ABS(ierr) )
  ALLOCATE( W( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'cjdsym', 'cannot allocate W', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( SW( npwx*npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( 'cjdsym', 'cannot allocate SW', ABS(ierr) )
  END IF
  !
  ALLOCATE( hc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'cjdsym', 'cannot allocate hc', ABS(ierr) )
  ALLOCATE( sc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'cjdsym', 'cannot allocate sc', ABS(ierr) )
  ALLOCATE( vc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'cjdsym', 'cannot allocate vc', ABS(ierr) )
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( 'cjdsym', 'cannot allocate ew', ABS(ierr) )
  !
  ALLOCATE( ub( npwx*npol, nblock ) )
  ALLOCATE( rb( npwx*npol, nblock ) )
  ALLOCATE( tb( npwx*npol, nblock ) )
  ALLOCATE( work2d( nvecx, nblock ) )
  ALLOCATE( Vtmp( npwx*npol, nvecx ) )
  ALLOCATE( Wtmp( npwx*npol, nvecx ) )
  IF ( uspp ) ALLOCATE( SWtmp( npwx*npol, nvecx ) )
  !
  ! ... Initialize subspace with input eigenvectors
  !
  CALL cjd_init_subspace()
  !
  e = 0.0_DP
  !
  ! ====================================================================
  ! ... Main blocked Jacobi-Davidson loop
  ! ====================================================================
  !
  iterate: DO iter = 1, maxter
     !
     jd_iter = iter
     !
     ! ... Diagonalize projected problem (must be done before restart
     ! ... so that vc is fresh and matches current j)
     !
     CALL cjd_diag_projected()
     !
     ! ... Block size for this iteration
     !
     nb = MIN( nblock, nvec - nconv, j )
     IF ( nb < 1 ) nb = 1
     !
     ! ... Restart if not enough room for nb new vectors
     !
     IF ( j + nb > nvecx ) THEN
        CALL cjd_restart()
        nb = MIN( nblock, nvec - nconv, j, nvecx - j )
        IF ( nb < 1 ) nb = 1
        ! ... Re-diag so vc matches the restarted (smaller) j
        CALL cjd_diag_projected()
     END IF
     !
     ! ... Compute Ritz pairs and residuals for nb lowest
     !
     CALL cjd_compute_residuals_block( nb )
     !
     ! ... Check convergence (consecutive from pair 1)
     !
     nconv_new = 0
     DO ib = 1, nb
        IF ( btype(nconv+ib) == 1 ) THEN
           IF ( rnorms(ib) < tol ) THEN
              nconv_new = nconv_new + 1
           ELSE
              EXIT
           END IF
        ELSE
           IF ( rnorms(ib) < empty_ethr ) THEN
              nconv_new = nconv_new + 1
           ELSE
              EXIT
           END IF
        END IF
     END DO
     !
     IF ( lprint ) THEN
        WRITE(6, '(5X,"cjdsym it=",I4," nconv=",I3," j=",I3,' // &
             '" nb=",I2," rnorms=",4ES10.3)') &
             iter, nconv, j, nb, (rnorms(ib), ib=1, MIN(nb,4))
        FLUSH(6)
     END IF
     !
     ! ... Store converged eigenpairs and deflate
     !
     IF ( nconv_new > 0 ) THEN
        !
        DO ib = 1, nconv_new
           nconv = nconv + 1
           evc(1:npwx*npol, nconv) = ub(1:npwx*npol, ib)
           e(nconv) = ew(ib)
           IF ( lprint ) THEN
              WRITE(6, '(5X,"cjdsym >>> band ",I4," CONVERGED: e=",F14.8,' // &
                   '" |r|=",ES10.3," at iter ",I4)') nconv, ew(ib), rnorms(ib), iter
              FLUSH(6)
           END IF
        END DO
        !
        IF ( nconv >= nvec ) EXIT iterate
        !
        CALL cjd_deflate_block( nconv_new )
        !
        ! ... Use remaining unconverged residuals instead of wasting an iteration
        !
        nb = nb - nconv_new
        IF ( nb > 0 ) THEN
           DO ib = 1, nb
              rb(1:npwx*npol, ib) = rb(1:npwx*npol, nconv_new + ib)
              ew(ib) = ew(nconv_new + ib)
           END DO
        ELSE
           CYCLE iterate
        END IF
        !
     END IF
     !
     ! ... Solve nb correction equations
     !
     CALL cjd_solve_corrections_block( nb )
     !
     ! ... Expand subspace with correction vectors
     !
     CALL cjd_expand_subspace_block( nb )
     !
  END DO iterate
  !
  ! ... Set output
  !
  notcnv = nvec - nconv
  !
  WRITE(6, '(5X,"cjdsym: finished. nconv=",I4," notcnv=",I4,' // &
       '" iter=",I4," nhpsi=",I6)') nconv, notcnv, jd_iter, nhpsi
  FLUSH(6)
  !
  ! ... For any remaining unconverged eigenvalues, use best Ritz approximation
  !
  IF ( notcnv > 0 .AND. j > 0 ) THEN
     DO i = 1, MIN(notcnv, j)
        CALL ZGEMV( 'N', kdim, j, ONE, V, kdmx, vc(1,i), 1, ZERO, &
                    evc(1,nconv+i), 1 )
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
  DEALLOCATE( Wtmp, Vtmp )
  DEALLOCATE( work2d )
  DEALLOCATE( tb, rb, ub )
  DEALLOCATE( ew )
  DEALLOCATE( vc, sc, hc )
  IF ( uspp ) DEALLOCATE( SW )
  DEALLOCATE( W, V )
  !
  CALL stop_clock( 'cjdsym' )
  !
  RETURN
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_init_subspace()
    !-----------------------------------------------------------------------
    !
    ! ... Initialize the search subspace: copy input eigenvectors,
    ! ... orthogonalize via modified Gram-Schmidt, compute H*V and S*V,
    ! ... and build the projected Hamiltonian and overlap matrices.
    !
    IMPLICIT NONE
    !
    nconv = 0
    j = nvec
    hc = ZERO
    sc = ZERO
    vc = ZERO
    ew = 0.0_DP
    !
    V(:,:) = ZERO
    W(:,:) = ZERO
    V(1:npwx*npol, 1:nvec) = evc(1:npwx*npol, 1:nvec)
    !
    ! ... Orthogonalize initial search space (modified Gram-Schmidt)
    !
    DO i = 1, nvec
       !
       IF ( i > 1 ) THEN
          CALL ZGEMV( 'C', kdim, i-1, ONE, V, kdmx, V(1,i), 1, ZERO, work2d(1,1), 1 )
          CALL mp_sum( work2d(1:i-1,1), intra_bgrp_comm )
          CALL ZGEMV( 'N', kdim, i-1, -ONE, V, kdmx, work2d(1,1), 1, ONE, V(1,i), 1 )
          ! ... Repeat for numerical stability
          CALL ZGEMV( 'C', kdim, i-1, ONE, V, kdmx, V(1,i), 1, ZERO, work2d(1,1), 1 )
          CALL mp_sum( work2d(1:i-1,1), intra_bgrp_comm )
          CALL ZGEMV( 'N', kdim, i-1, -ONE, V, kdmx, work2d(1,1), 1, ONE, V(1,i), 1 )
       END IF
       !
       norm_t = 0.0_DP
       DO ig = 1, kdim
          norm_t = norm_t + DBLE( CONJG(V(ig,i)) * V(ig,i) )
       END DO
       CALL mp_sum( norm_t, intra_bgrp_comm )
       norm_t = SQRT( norm_t )
       IF ( norm_t > 1.0D-14 ) THEN
          V(1:kdim,i) = V(1:kdim,i) / norm_t
       END IF
       !
       IF ( npol == 1 .AND. npw < npwx ) V(npw+1:npwx,i) = ZERO
       IF ( npol == 2 .AND. npw < npwx ) THEN
          V(npw+1:npwx,i) = ZERO
          V(npwx+npw+1:2*npwx,i) = ZERO
       END IF
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
    ! ... Build projected Hamiltonian hc = V^H * W
    !
    CALL ZGEMM( 'C', 'N', j, j, kdim, ONE, V, kdmx, W, kdmx, ZERO, hc, nvecx )
    CALL mp_sum( hc(1:j, 1:j), intra_bgrp_comm )
    !
    ! ... Build projected overlap sc = V^H * S*V (or V^H * V)
    !
    IF ( uspp ) THEN
       CALL ZGEMM( 'C', 'N', j, j, kdim, ONE, V, kdmx, SW, kdmx, ZERO, sc, nvecx )
    ELSE
       CALL ZGEMM( 'C', 'N', j, j, kdim, ONE, V, kdmx, V, kdmx, ZERO, sc, nvecx )
    END IF
    CALL mp_sum( sc(1:j, 1:j), intra_bgrp_comm )
    !
    ! ... Symmetrize (ensure Hermitian)
    !
    DO i = 1, j
       hc(i,i) = CMPLX( REAL( hc(i,i) ), 0.0_DP, kind=DP )
       sc(i,i) = CMPLX( REAL( sc(i,i) ), 0.0_DP, kind=DP )
       DO k = i + 1, j
          hc(i,k) = CONJG( hc(k,i) )
          sc(i,k) = CONJG( sc(k,i) )
       END DO
    END DO
    !
  END SUBROUTINE cjd_init_subspace
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_diag_projected()
    !-----------------------------------------------------------------------
    !
    ! ... Diagonalize the projected eigenproblem hc * vc = sc * vc * diag(ew).
    !
    IMPLICIT NONE
    !
    CALL start_clock( 'cjdsym:diag' )
    IF ( my_bgrp_id == root_bgrp_id ) THEN
       CALL diaghg( j, j, hc, sc, nvecx, ew, vc, &
                    me_bgrp, root_bgrp, intra_bgrp_comm )
    END IF
    IF ( nbgrp > 1 ) THEN
       CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
       CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
    END IF
    CALL stop_clock( 'cjdsym:diag' )
    !
  END SUBROUTINE cjd_diag_projected
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_compute_residuals_block( nb_in )
    !-----------------------------------------------------------------------
    !
    ! ... Compute nb_in Ritz vectors, their residuals, and residual norms.
    ! ... Uses ZGEMM for blocked matrix-vector products.
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nb_in
    !
    ! ... Compute Ritz vectors: ub(:,1:nb) = V * vc(:,1:nb)
    !
    CALL ZGEMM( 'N', 'N', kdim, nb_in, j, ONE, V, kdmx, &
                vc(1,1), nvecx, ZERO, ub, kdmx )
    !
    ! ... Compute H * Ritz vectors: rb = W * vc(:,1:nb) (temporary)
    !
    CALL ZGEMM( 'N', 'N', kdim, nb_in, j, ONE, W, kdmx, &
                vc(1,1), nvecx, ZERO, rb, kdmx )
    !
    ! ... Compute S * Ritz vectors (into tb) and form residuals
    !
    IF ( uspp ) THEN
       CALL ZGEMM( 'N', 'N', kdim, nb_in, j, ONE, SW, kdmx, &
                   vc(1,1), nvecx, ZERO, tb, kdmx )
       DO ib = 1, nb_in
          rb(1:kdim,ib) = rb(1:kdim,ib) - ew(ib) * tb(1:kdim,ib)
       END DO
    ELSE
       DO ib = 1, nb_in
          rb(1:kdim,ib) = rb(1:kdim,ib) - ew(ib) * ub(1:kdim,ib)
       END DO
    END IF
    !
    ! ... Orthogonalize residuals against converged eigenvectors
    !
    IF ( nconv > 0 ) THEN
       DO i = 1, 2
          CALL ZGEMM( 'C', 'N', nconv, nb_in, kdim, ONE, evc, kdmx, &
                      rb, kdmx, ZERO, work2d, nvecx )
          CALL mp_sum( work2d(1:nconv, 1:nb_in), intra_bgrp_comm )
          CALL ZGEMM( 'N', 'N', kdim, nb_in, nconv, -ONE, evc, kdmx, &
                      work2d, nvecx, ONE, rb, kdmx )
       END DO
    END IF
    !
    ! ... Compute residual norms
    !
    DO ib = 1, nb_in
       rnorms(ib) = 0.0_DP
       DO ig = 1, kdim
          rnorms(ib) = rnorms(ib) + DBLE( CONJG(rb(ig,ib)) * rb(ig,ib) )
       END DO
    END DO
    CALL mp_sum( rnorms(1:nb_in), intra_bgrp_comm )
    rnorms(1:nb_in) = SQRT( rnorms(1:nb_in) )
    !
  END SUBROUTINE cjd_compute_residuals_block
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_deflate_block( nrem )
    !-----------------------------------------------------------------------
    !
    ! ... Remove nrem converged Ritz vectors from the search subspace.
    ! ... Rotates V, W, SW to keep Ritz vectors nrem+1:j.
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nrem
    INTEGER :: jnew
    !
    IF ( j > nrem ) THEN
       !
       jnew = j - nrem
       !
       CALL ZGEMM( 'N', 'N', kdim, jnew, j, ONE, V, kdmx, &
                   vc(1,nrem+1), nvecx, ZERO, Vtmp, kdmx )
       V(1:npwx*npol, 1:jnew) = Vtmp(1:npwx*npol, 1:jnew)
       !
       CALL ZGEMM( 'N', 'N', kdim, jnew, j, ONE, W, kdmx, &
                   vc(1,nrem+1), nvecx, ZERO, Wtmp, kdmx )
       W(1:npwx*npol, 1:jnew) = Wtmp(1:npwx*npol, 1:jnew)
       !
       IF ( uspp ) THEN
          CALL ZGEMM( 'N', 'N', kdim, jnew, j, ONE, SW, kdmx, &
                      vc(1,nrem+1), nvecx, ZERO, SWtmp, kdmx )
          SW(1:npwx*npol, 1:jnew) = SWtmp(1:npwx*npol, 1:jnew)
       END IF
       !
       j = jnew
       !
       ! ... After rotation by eigenvectors: hc = diag(ew(nrem+1:)), sc = I
       !
       hc = ZERO
       sc = ZERO
       DO i = 1, j
          hc(i,i) = CMPLX( ew(nrem+i), 0.0_DP, kind=DP )
          sc(i,i) = ONE
       END DO
       !
    ELSE
       !
       CALL cjd_reinit_subspace()
       !
    END IF
    !
  END SUBROUTINE cjd_deflate_block
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_restart()
    !-----------------------------------------------------------------------
    !
    ! ... Restart the subspace by keeping the best jmin Ritz vectors.
    !
    IMPLICIT NONE
    !
    CALL start_clock( 'cjdsym:restart' )
    !
    IF ( lprint ) THEN
       WRITE(6, '(5X,"cjdsym: RESTART j=",I4," -> jmin=",I4)') j, jmin
       FLUSH(6)
    END IF
    !
    CALL ZGEMM( 'N', 'N', kdim, jmin, j, ONE, V, kdmx, &
                vc(1,1), nvecx, ZERO, Vtmp, kdmx )
    V(1:npwx*npol, 1:jmin) = Vtmp(1:npwx*npol, 1:jmin)
    !
    CALL ZGEMM( 'N', 'N', kdim, jmin, j, ONE, W, kdmx, &
                vc(1,1), nvecx, ZERO, Wtmp, kdmx )
    W(1:npwx*npol, 1:jmin) = Wtmp(1:npwx*npol, 1:jmin)
    !
    IF ( uspp ) THEN
       CALL ZGEMM( 'N', 'N', kdim, jmin, j, ONE, SW, kdmx, &
                   vc(1,1), nvecx, ZERO, SWtmp, kdmx )
       SW(1:npwx*npol, 1:jmin) = SWtmp(1:npwx*npol, 1:jmin)
    END IF
    !
    j = jmin
    !
    hc = ZERO
    sc = ZERO
    DO i = 1, j
       hc(i,i) = CMPLX( ew(i), 0.0_DP, kind=DP )
       sc(i,i) = ONE
    END DO
    !
    CALL stop_clock( 'cjdsym:restart' )
    !
  END SUBROUTINE cjd_restart
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_solve_corrections_block( nb_in )
    !-----------------------------------------------------------------------
    !
    ! ... Apply preconditioner to nb_in residuals (in-place in rb),
    ! ... then project out converged eigenvectors.
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nb_in
    !
    CALL start_clock( 'cjdsym:correction' )
    !
    ! ... Apply preconditioner
    !
    IF ( use_g_psi ) THEN
       ! ... g_psi_ptr handles all nb_in vectors at once
       CALL g_psi_ptr( npwx, npw, nb_in, npol, rb, ew )
    ELSE
       ! ... TPA preconditioner: rb = rb / (g2kin + shift - ew)
       DO ib = 1, nb_in
          DO ipol = 1, npol
             DO ig = 1, npw
                rb(ig+(ipol-1)*npwx,ib) = rb(ig+(ipol-1)*npwx,ib) / &
                   ( g2kin(ig) + default_shift - ew(ib) )
             END DO
          END DO
       END DO
    END IF
    !
    ! ... Project out converged eigenvectors (double for stability)
    !
    IF ( nconv > 0 ) THEN
       DO i = 1, 2
          CALL ZGEMM( 'C', 'N', nconv, nb_in, kdim, ONE, evc, kdmx, &
                      rb, kdmx, ZERO, work2d, nvecx )
          CALL mp_sum( work2d(1:nconv, 1:nb_in), intra_bgrp_comm )
          CALL ZGEMM( 'N', 'N', kdim, nb_in, nconv, -ONE, evc, kdmx, &
                      work2d, nvecx, ONE, rb, kdmx )
       END DO
    END IF
    !
    CALL stop_clock( 'cjdsym:correction' )
    !
  END SUBROUTINE cjd_solve_corrections_block
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_expand_subspace_block( nb_in )
    !-----------------------------------------------------------------------
    !
    ! ... Orthogonalize nb_in correction vectors (in rb) against V and
    ! ... each other, normalize, add to V, compute H*V and S*V for the
    ! ... new vectors, and update the projected matrices hc and sc.
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nb_in
    INTEGER :: jj, ii
    !
    CALL start_clock( 'cjdsym:ortho' )
    !
    nact = 0
    !
    DO ib = 1, nb_in
       !
       ! ... Orthogonalize rb(:,ib) against V(:,1:j+nact) (double Gram-Schmidt)
       ! ... This includes previously accepted corrections in V(:,j+1:j+nact)
       !
       DO i = 1, 2
          CALL ZGEMV( 'C', kdim, j+nact, ONE, V, kdmx, rb(1,ib), 1, &
                      ZERO, work2d(1,1), 1 )
          CALL mp_sum( work2d(1:j+nact,1), intra_bgrp_comm )
          CALL ZGEMV( 'N', kdim, j+nact, -ONE, V, kdmx, work2d(1,1), 1, &
                      ONE, rb(1,ib), 1 )
       END DO
       !
       ! ... Normalize
       !
       norm_t = 0.0_DP
       DO ig = 1, kdim
          norm_t = norm_t + DBLE( CONJG(rb(ig,ib)) * rb(ig,ib) )
       END DO
       CALL mp_sum( norm_t, intra_bgrp_comm )
       norm_t = SQRT( norm_t )
       !
       IF ( norm_t < 1.0D-14 ) CYCLE  ! skip this correction
       !
       rb(1:kdim,ib) = rb(1:kdim,ib) / norm_t
       !
       ! ... Zero padding
       !
       IF ( npol == 1 .AND. npw < npwx ) rb(npw+1:npwx,ib) = ZERO
       IF ( npol == 2 .AND. npw < npwx ) THEN
          rb(npw+1:npwx,ib) = ZERO
          rb(npwx+npw+1:2*npwx,ib) = ZERO
       END IF
       !
       ! ... Accept: place into search space
       !
       nact = nact + 1
       V(1:npwx*npol, j+nact) = rb(1:npwx*npol, ib)
       !
    END DO
    !
    CALL stop_clock( 'cjdsym:ortho' )
    !
    IF ( nact == 0 ) RETURN
    !
    ! ... Compute H*V and S*V for new columns (single blocked call)
    !
    CALL h_psi_ptr( npwx, npw, nact, V(1,j+1), W(1,j+1) )
    nhpsi = nhpsi + nact
    !
    IF ( uspp ) CALL s_psi_ptr( npwx, npw, nact, V(1,j+1), SW(1,j+1) )
    !
    ! ... Update projected Hamiltonian (blocked ZGEMM)
    !
    CALL start_clock( 'cjdsym:overlap' )
    !
    CALL ZGEMM( 'C', 'N', j+nact, nact, kdim, ONE, V, kdmx, &
                W(1,j+1), kdmx, ZERO, hc(1,j+1), nvecx )
    CALL mp_sum( hc(1:j+nact, j+1:j+nact), intra_bgrp_comm )
    !
    DO ib = 1, nact
       jj = j + ib
       DO ii = 1, jj - 1
          hc(jj,ii) = CONJG( hc(ii,jj) )
       END DO
       hc(jj,jj) = CMPLX( REAL( hc(jj,jj) ), 0.0_DP, kind=DP )
    END DO
    !
    ! ... Update projected overlap (blocked ZGEMM)
    !
    IF ( uspp ) THEN
       CALL ZGEMM( 'C', 'N', j+nact, nact, kdim, ONE, V, kdmx, &
                   SW(1,j+1), kdmx, ZERO, sc(1,j+1), nvecx )
    ELSE
       CALL ZGEMM( 'C', 'N', j+nact, nact, kdim, ONE, V, kdmx, &
                   V(1,j+1), kdmx, ZERO, sc(1,j+1), nvecx )
    END IF
    CALL mp_sum( sc(1:j+nact, j+1:j+nact), intra_bgrp_comm )
    !
    DO ib = 1, nact
       jj = j + ib
       DO ii = 1, jj - 1
          sc(jj,ii) = CONJG( sc(ii,jj) )
       END DO
       sc(jj,jj) = CMPLX( REAL( sc(jj,jj) ), 0.0_DP, kind=DP )
    END DO
    !
    j = j + nact
    !
    CALL stop_clock( 'cjdsym:overlap' )
    !
  END SUBROUTINE cjd_expand_subspace_block
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_reinit_subspace()
    !-----------------------------------------------------------------------
    !
    ! ... Reinitialize the search space when it becomes empty after deflation.
    !
    IMPLICIT NONE
    INTEGER :: ig2, ipol2
    REAL(DP) :: rr, ri
    !
    j = 1
    V(1:npwx*npol, 1) = ZERO
    DO ipol2 = 1, npol
       DO ig2 = 1, npw
          rr = DBLE(MOD(ig2 + nconv*137, 1000)) / 1000.0_DP
          ri = DBLE(MOD(ig2 + nconv*251, 1000)) / 1000.0_DP
          V(ig2 + (ipol2-1)*npwx, 1) = CMPLX( rr, ri, kind=DP )
       END DO
    END DO
    !
    DO i = 1, 2
       CALL ZGEMV( 'C', kdim, nconv, ONE, evc, kdmx, V(1,1), 1, ZERO, work2d(1,1), 1 )
       CALL mp_sum( work2d(1:nconv,1), intra_bgrp_comm )
       CALL ZGEMV( 'N', kdim, nconv, -ONE, evc, kdmx, work2d(1,1), 1, ONE, V(1,1), 1 )
    END DO
    !
    norm_t = 0.0_DP
    DO ig2 = 1, kdim
       norm_t = norm_t + DBLE( CONJG(V(ig2,1)) * V(ig2,1) )
    END DO
    CALL mp_sum( norm_t, intra_bgrp_comm )
    norm_t = SQRT( norm_t )
    IF ( norm_t > 1.0D-14 ) THEN
       V(1:kdim, 1) = V(1:kdim, 1) / norm_t
    END IF
    !
    CALL h_psi_ptr( npwx, npw, 1, V(1,1), W(1,1) )
    nhpsi = nhpsi + 1
    IF ( uspp ) CALL s_psi_ptr( npwx, npw, 1, V(1,1), SW(1,1) )
    !
    hc = ZERO
    sc = ZERO
    CALL ZGEMV( 'C', kdim, 1, ONE, V, kdmx, W(1,1), 1, ZERO, hc(1,1), 1 )
    CALL mp_sum( hc(1:1, 1:1), intra_bgrp_comm )
    hc(1,1) = CMPLX( REAL( hc(1,1) ), 0.0_DP, kind=DP )
    !
    IF ( uspp ) THEN
       CALL ZGEMV( 'C', kdim, 1, ONE, V, kdmx, SW(1,1), 1, ZERO, sc(1,1), 1 )
    ELSE
       CALL ZGEMV( 'C', kdim, 1, ONE, V, kdmx, V(1,1), 1, ZERO, sc(1,1), 1 )
    END IF
    CALL mp_sum( sc(1:1, 1:1), intra_bgrp_comm )
    sc(1,1) = CMPLX( REAL( sc(1,1) ), 0.0_DP, kind=DP )
    !
  END SUBROUTINE cjd_reinit_subspace
  !
END SUBROUTINE cjdsym
