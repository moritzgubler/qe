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
SUBROUTINE cjdsym( h_psi_ptr, s_psi_ptr, uspp, &
                   npw, npwx, nvec, nvecx, npol, evc, ethr, &
                   g2kin, e, btype, notcnv, jd_iter, nhpsi )
  !----------------------------------------------------------------------------
  !
  ! ... Jacobi-Davidson iterative diagonalization of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is a Hermitian operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector.
  !
  ! ... Processes one eigenvalue at a time with explicit deflation.
  ! ... Uses TPA preconditioner: t = r / (g2kin + default_shift).
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
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    ! number of spin polarizations
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx*npol,nvec)
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
  INTEGER, PARAMETER :: maxter = 400
    ! maximum number of iterations
  REAL(DP), PARAMETER :: default_shift = 1.0_DP
    ! minimum denominator for preconditioner
  !
  INTEGER :: j, nconv, iter, kdim, kdmx, ierr, jmin
  INTEGER :: i, ig, ipol
  REAL(DP) :: nr, nt, tol, empty_ethr
  COMPLEX(DP) :: coeff_u, cdot
  REAL(DP) :: rdot
  LOGICAL :: lprint, skip
  !
  COMPLEX(DP), ALLOCATABLE :: V(:,:), W(:,:), SW(:,:)
    ! search space basis vectors / H * V / S * V
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
    ! projected Hamiltonian / overlap / eigenvectors
  REAL(DP), ALLOCATABLE :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: u(:), Au(:), Su(:), r(:), t(:)
    ! Ritz vector / H*u / S*u / residual / correction
  COMPLEX(DP), ALLOCATABLE :: work(:)
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
  lprint = .FALSE.
  CALL start_clock( 'cjdsym' )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'cjdsym', 'nvecx is too small', 1 )
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  tol = ethr
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
  IF ( lprint ) THEN
     WRITE(6, '(5X,"cjdsym: npw=",I8," npwx=",I8," nvec=",I4,' // &
          '" nvecx=",I4," npol=",I2)') npw, npwx, nvec, nvecx, npol
     WRITE(6, '(5X,"cjdsym: ethr=",ES10.3," tol=",ES10.3,' // &
          '" empty_ethr=",ES10.3," uspp=",L2)') ethr, tol, empty_ethr, uspp
     WRITE(6, '(5X,"cjdsym: jmin=",I4," maxter=",I6,' // &
          '" default_shift=",F6.2)') jmin, maxter, default_shift
     WRITE(6, '(5X,"cjdsym: btype(1:nvec)=",20I2)') btype(1:nvec)
     FLUSH(6)
  END IF
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
  ALLOCATE( u( npwx*npol ), Au( npwx*npol ), Su( npwx*npol ) )
  ALLOCATE( r( npwx*npol ), t( npwx*npol ) )
  ALLOCATE( work( nvecx ) )
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
  ! ... Main Jacobi-Davidson loop
  ! ====================================================================
  !
  IF ( lprint ) THEN
     WRITE(6, '(5X,"cjdsym: entering main loop, initial subspace dim j=",I4)') j
     FLUSH(6)
  END IF
  !
  iterate: DO iter = 1, maxter
     !
     jd_iter = iter
     !
     ! ... Diagonalize projected problem and extract best Ritz pair
     !
     CALL cjd_extract_ritz_pair( nr )
     !
     ! ... Check convergence (relaxed threshold for empty bands)
     !
     IF ( ( btype(nconv+1) == 1 .AND. nr < tol ) .OR. &
          ( btype(nconv+1) /= 1 .AND. nr < empty_ethr ) ) THEN
        !
        nconv = nconv + 1
        evc(1:npwx*npol, nconv) = u(1:npwx*npol)
        e(nconv) = ew(1)
        !
        IF ( lprint ) THEN
           WRITE(6, '(5X,"cjdsym >>> band ",I4," CONVERGED: e=",F14.8,' // &
                '" |r|=",ES10.3," at iter ",I4)') nconv, ew(1), nr, iter
           FLUSH(6)
        END IF
        !
        IF ( nconv >= nvec ) EXIT iterate
        !
        CALL cjd_deflate()
        CYCLE iterate
        !
     END IF
     !
     ! ... Restart if subspace is full
     !
     IF ( j >= nvecx ) CALL cjd_restart()
     !
     ! ... Solve correction equation
     !
     CALL cjd_solve_correction()
     !
     ! ... Expand subspace with correction vector
     !
     CALL cjd_expand_subspace( skip )
     IF ( skip ) CYCLE iterate
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
     IF ( lprint ) THEN
        WRITE(6, '(5X,"cjdsym: unconverged Ritz values: ",8F12.6)') &
             (e(nconv+i), i=1, MIN(notcnv, 8))
     END IF
  END IF
  !
  ! ... Deallocate
  !
  IF ( uspp ) DEALLOCATE( SWtmp )
  DEALLOCATE( Wtmp )
  DEALLOCATE( Vtmp )
  DEALLOCATE( work )
  DEALLOCATE( t, r, Su, Au, u )
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
    INTEGER :: k
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
          CALL ZGEMV( 'C', kdim, i-1, ONE, V, kdmx, V(1,i), 1, ZERO, work, 1 )
          CALL mp_sum( work(1:i-1), intra_bgrp_comm )
          CALL ZGEMV( 'N', kdim, i-1, -ONE, V, kdmx, work, 1, ONE, V(1,i), 1 )
          ! ... Repeat for numerical stability
          CALL ZGEMV( 'C', kdim, i-1, ONE, V, kdmx, V(1,i), 1, ZERO, work, 1 )
          CALL mp_sum( work(1:i-1), intra_bgrp_comm )
          CALL ZGEMV( 'N', kdim, i-1, -ONE, V, kdmx, work, 1, ONE, V(1,i), 1 )
       END IF
       !
       nr = 0.0_DP
       DO ig = 1, kdim
          nr = nr + DBLE( CONJG(V(ig,i)) * V(ig,i) )
       END DO
       CALL mp_sum( nr, intra_bgrp_comm )
       nr = SQRT( nr )
       IF ( nr > 1.0D-14 ) THEN
          V(1:kdim,i) = V(1:kdim,i) / nr
       ELSE IF ( lprint ) THEN
          WRITE(6, '(5X,"cjdsym WARNING: initial vector ",I4,' // &
               '" has near-zero norm after orthogonalization: ",ES10.3)') i, nr
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
    IF ( lprint ) THEN
       WRITE(6, '(5X,"cjdsym: initial hc diag(1:min(5,j)): ",5F12.6)') &
            (REAL(hc(i,i)), i=1, MIN(5,j))
       FLUSH(6)
    END IF
    !
  END SUBROUTINE cjd_init_subspace
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_extract_ritz_pair( rnorm )
    !-----------------------------------------------------------------------
    !
    ! ... Diagonalize the projected eigenproblem, extract the best Ritz
    ! ... pair (u, theta), and compute the residual norm.
    !
    IMPLICIT NONE
    REAL(DP), INTENT(OUT) :: rnorm
    !
    ! ... Diagonalize hc * vc = sc * vc * diag(ew)
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
    ! ... Compute Ritz vector u = V * vc(:,1)
    !
    u = ZERO
    CALL ZGEMV( 'N', kdim, j, ONE, V, kdmx, vc(1,1), 1, ZERO, u, 1 )
    !
    ! ... Compute Au = W * vc(:,1) = H*u
    !
    Au = ZERO
    CALL ZGEMV( 'N', kdim, j, ONE, W, kdmx, vc(1,1), 1, ZERO, Au, 1 )
    !
    ! ... Compute Su = SW * vc(:,1) = S*u (or u if no uspp)
    !
    IF ( uspp ) THEN
       Su = ZERO
       CALL ZGEMV( 'N', kdim, j, ONE, SW, kdmx, vc(1,1), 1, ZERO, Su, 1 )
    ELSE
       Su = u
    END IF
    !
    ! ... Compute residual: r = Au - theta * Su
    !
    r(1:kdim) = Au(1:kdim) - ew(1) * Su(1:kdim)
    !
    ! ... Orthogonalize residual against converged eigenvectors
    !
    IF ( nconv > 0 ) THEN
       DO i = 1, 2
          CALL ZGEMV( 'C', kdim, nconv, ONE, evc, kdmx, r, 1, ZERO, work, 1 )
          CALL mp_sum( work(1:nconv), intra_bgrp_comm )
          CALL ZGEMV( 'N', kdim, nconv, -ONE, evc, kdmx, work, 1, ONE, r, 1 )
       END DO
    END IF
    !
    ! ... Compute residual norm
    !
    rnorm = 0.0_DP
    DO i = 1, kdim
       rnorm = rnorm + DBLE( CONJG(r(i)) * r(i) )
    END DO
    CALL mp_sum( rnorm, intra_bgrp_comm )
    rnorm = SQRT( rnorm )
    !
    IF ( lprint .AND. ( MOD(iter, 1) == 0 .OR. iter <= 10 ) ) THEN
       WRITE(6, '(5X,"cjdsym it=",I4," nconv=",I3," j=",I3,' // &
            '" theta=",F14.8," |r|=",ES10.3," tol=",ES10.3)') &
            iter, nconv, j, ew(1), rnorm, tol
       IF ( j >= 2 ) THEN
          WRITE(6, '(5X,"  ew(1:",I2,")=",8F12.6)') &
               MIN(j, 8), (ew(i), i=1, MIN(j, 8))
       END IF
       FLUSH(6)
    END IF
    !
  END SUBROUTINE cjd_extract_ritz_pair
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_deflate()
    !-----------------------------------------------------------------------
    !
    ! ... Remove the converged component from the search subspace.
    ! ... Rotates V, W, SW to keep Ritz vectors 2:j, or reinitializes
    ! ... if the subspace becomes empty.
    !
    IMPLICIT NONE
    !
    IF ( j > 1 ) THEN
       !
       ! ... Rotate subspace: V_new = V * vc(:,2:j), etc.
       !
       CALL ZGEMM( 'N', 'N', kdim, j-1, j, ONE, V, kdmx, &
                   vc(1,2), nvecx, ZERO, Vtmp, kdmx )
       V(1:npwx*npol, 1:j-1) = Vtmp(1:npwx*npol, 1:j-1)
       !
       CALL ZGEMM( 'N', 'N', kdim, j-1, j, ONE, W, kdmx, &
                   vc(1,2), nvecx, ZERO, Wtmp, kdmx )
       W(1:npwx*npol, 1:j-1) = Wtmp(1:npwx*npol, 1:j-1)
       !
       IF ( uspp ) THEN
          CALL ZGEMM( 'N', 'N', kdim, j-1, j, ONE, SW, kdmx, &
                      vc(1,2), nvecx, ZERO, SWtmp, kdmx )
          SW(1:npwx*npol, 1:j-1) = SWtmp(1:npwx*npol, 1:j-1)
       END IF
       !
       j = j - 1
       !
       ! ... After rotation by eigenvectors: hc = diag(ew(2:)), sc = I
       !
       hc = ZERO
       sc = ZERO
       DO i = 1, j
          hc(i,i) = CMPLX( ew(i+1), 0.0_DP, kind=DP )
          sc(i,i) = ONE
       END DO
       !
       IF ( lprint ) THEN
          WRITE(6, '(5X,"cjdsym: deflated, new subspace dim j=",I4)') j
          FLUSH(6)
       END IF
       !
    ELSE
       !
       ! ... Subspace empty after deflation, reinitialize
       !
       IF ( lprint ) THEN
          WRITE(6, '(5X,"cjdsym: subspace empty, reinitializing")')
          FLUSH(6)
       END IF
       CALL cjd_reinit_subspace()
       !
    END IF
    !
  END SUBROUTINE cjd_deflate
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_restart()
    !-----------------------------------------------------------------------
    !
    ! ... Restart the subspace by keeping the best jmin Ritz vectors
    ! ... when the subspace dimension reaches nvecx.
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
    ! ... After rotation by eigenvectors: hc = diag(ew), sc = I
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
  SUBROUTINE cjd_solve_correction()
    !-----------------------------------------------------------------------
    !
    ! ... Solve the correction equation (simplified Jacobi-Davidson):
    ! ... Apply TPA preconditioner t = r / (g2kin + shift), then
    ! ... project out converged eigenvectors and the current Ritz vector.
    !
    IMPLICIT NONE
    !
    CALL start_clock( 'cjdsym:correction' )
    !
    ! ... Apply TPA preconditioner
    !
    t = ZERO
    DO ipol = 1, npol
       DO ig = 1, npw
        t(ig + (ipol-1)*npwx) = r(ig + (ipol-1)*npwx) / (g2kin(ig) + default_shift)
       END DO
    END DO
    !
    ! ... Project out converged eigenvectors (double for stability)
    !
    IF ( nconv > 0 ) THEN
       DO i = 1, 2
          CALL ZGEMV( 'C', kdim, nconv, ONE, evc, kdmx, t, 1, ZERO, work, 1 )
          CALL mp_sum( work(1:nconv), intra_bgrp_comm )
          CALL ZGEMV( 'N', kdim, nconv, -ONE, evc, kdmx, work, 1, ONE, t, 1 )
       END DO
    END IF
    !
    ! ... Project out current Ritz vector u
    !
    cdot = ZERO
    DO i = 1, kdim
       cdot = cdot + CONJG(u(i)) * t(i)
    END DO
    CALL mp_sum( cdot, intra_bgrp_comm )
    !
    rdot = 0.0_DP
    DO i = 1, kdim
       rdot = rdot + DBLE( CONJG(u(i)) * u(i) )
    END DO
    CALL mp_sum( rdot, intra_bgrp_comm )
    !
    IF ( rdot > 1.0D-30 ) THEN
       coeff_u = cdot / CMPLX( rdot, 0.0_DP, kind=DP )
       t(1:kdim) = t(1:kdim) - coeff_u * u(1:kdim)
    END IF
    !
    CALL stop_clock( 'cjdsym:correction' )
    !
  END SUBROUTINE cjd_solve_correction
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_expand_subspace( skip_expand )
    !-----------------------------------------------------------------------
    !
    ! ... Orthogonalize the correction vector t against the search space V,
    ! ... normalize it, add it to V, compute H*t and S*t, and update the
    ! ... projected matrices hc and sc.
    !
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: skip_expand
    !
    skip_expand = .FALSE.
    !
    ! ... Orthogonalize t against V (double Gram-Schmidt for stability)
    !
    CALL start_clock( 'cjdsym:ortho' )
    !
    DO i = 1, 2
       CALL ZGEMV( 'C', kdim, j, ONE, V, kdmx, t, 1, ZERO, work, 1 )
       CALL mp_sum( work(1:j), intra_bgrp_comm )
       CALL ZGEMV( 'N', kdim, j, -ONE, V, kdmx, work, 1, ONE, t, 1 )
    END DO
    !
    ! ... Normalize t
    !
    nt = 0.0_DP
    DO i = 1, kdim
       nt = nt + DBLE( CONJG(t(i)) * t(i) )
    END DO
    CALL mp_sum( nt, intra_bgrp_comm )
    nt = SQRT( nt )
    !
    CALL stop_clock( 'cjdsym:ortho' )
    !
    IF ( nt < 1.0D-14 ) THEN
       IF ( lprint ) THEN
          WRITE(6, '(5X,"cjdsym: correction too small, |t|=",ES10.3,' // &
               '" skipping")') nt
          FLUSH(6)
       END IF
       skip_expand = .TRUE.
       RETURN
    END IF
    !
    t(1:kdim) = t(1:kdim) / nt
    IF ( npol == 1 .AND. npw < npwx ) t(npw+1:npwx) = ZERO
    IF ( npol == 2 .AND. npw < npwx ) THEN
       t(npw+1:npwx) = ZERO
       t(npwx+npw+1:2*npwx) = ZERO
    END IF
    !
    ! ... Add correction to search space
    !
    j = j + 1
    V(1:npwx*npol, j) = t(1:npwx*npol)
    !
    ! ... Compute H*t and S*t
    !
    CALL h_psi_ptr( npwx, npw, 1, V(1,j), W(1,j) )
    nhpsi = nhpsi + 1
    !
    IF ( uspp ) CALL s_psi_ptr( npwx, npw, 1, V(1,j), SW(1,j) )
    !
    ! ... Update projected Hamiltonian: new column hc(1:j, j)
    !
    CALL start_clock( 'cjdsym:overlap' )
    !
    CALL ZGEMV( 'C', kdim, j, ONE, V, kdmx, W(1,j), 1, ZERO, hc(1,j), 1 )
    CALL mp_sum( hc(1:j, j), intra_bgrp_comm )
    !
    DO i = 1, j - 1
       hc(j,i) = CONJG( hc(i,j) )
    END DO
    hc(j,j) = CMPLX( REAL( hc(j,j) ), 0.0_DP, kind=DP )
    !
    ! ... Update projected overlap: new column sc(1:j, j)
    !
    IF ( uspp ) THEN
       CALL ZGEMV( 'C', kdim, j, ONE, V, kdmx, SW(1,j), 1, ZERO, sc(1,j), 1 )
    ELSE
       CALL ZGEMV( 'C', kdim, j, ONE, V, kdmx, V(1,j), 1, ZERO, sc(1,j), 1 )
    END IF
    CALL mp_sum( sc(1:j, j), intra_bgrp_comm )
    !
    DO i = 1, j - 1
       sc(j,i) = CONJG( sc(i,j) )
    END DO
    sc(j,j) = CMPLX( REAL( sc(j,j) ), 0.0_DP, kind=DP )
    !
    CALL stop_clock( 'cjdsym:overlap' )
    !
  END SUBROUTINE cjd_expand_subspace
  !
  !-----------------------------------------------------------------------
  SUBROUTINE cjd_reinit_subspace()
    !-----------------------------------------------------------------------
    !
    ! ... Reinitialize the search space when it becomes empty after deflation.
    ! ... Creates a random vector orthogonal to converged eigenvectors.
    !
    IMPLICIT NONE
    INTEGER :: ig2, ipol2
    REAL(DP) :: rr, ri
    !
    ! ... Create a pseudo-random initial vector
    !
    j = 1
    V(1:npwx*npol, 1) = ZERO
    DO ipol2 = 1, npol
       DO ig2 = 1, npw
          ! ... Simple deterministic initialization
          rr = DBLE(MOD(ig2 + nconv*137, 1000)) / 1000.0_DP
          ri = DBLE(MOD(ig2 + nconv*251, 1000)) / 1000.0_DP
          V(ig2 + (ipol2-1)*npwx, 1) = CMPLX( rr, ri, kind=DP )
       END DO
    END DO
    !
    ! ... Orthogonalize against converged eigenvectors (double)
    !
    DO i = 1, 2
       CALL ZGEMV( 'C', kdim, nconv, ONE, evc, kdmx, V(1,1), 1, ZERO, work, 1 )
       CALL mp_sum( work(1:nconv), intra_bgrp_comm )
       CALL ZGEMV( 'N', kdim, nconv, -ONE, evc, kdmx, work, 1, ONE, V(1,1), 1 )
    END DO
    !
    ! ... Normalize
    !
    nt = 0.0_DP
    DO ig2 = 1, kdim
       nt = nt + DBLE( CONJG(V(ig2,1)) * V(ig2,1) )
    END DO
    CALL mp_sum( nt, intra_bgrp_comm )
    nt = SQRT( nt )
    IF ( nt > 1.0D-14 ) THEN
       V(1:kdim, 1) = V(1:kdim, 1) / nt
    END IF
    !
    ! ... Compute H*V and S*V
    !
    CALL h_psi_ptr( npwx, npw, 1, V(1,1), W(1,1) )
    nhpsi = nhpsi + 1
    IF ( uspp ) CALL s_psi_ptr( npwx, npw, 1, V(1,1), SW(1,1) )
    !
    ! ... Build projected matrices (1x1)
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
    IF ( lprint ) THEN
       WRITE(6, '(5X,"cjdsym: reinit done, hc(1,1)=",F14.8," sc(1,1)=",F14.8)') &
            REAL(hc(1,1)), REAL(sc(1,1))
       FLUSH(6)
    END IF
    !
  END SUBROUTINE cjd_reinit_subspace
  !
END SUBROUTINE cjdsym
