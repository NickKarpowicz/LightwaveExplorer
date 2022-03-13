!===============================================================================
! Copyright 1999-2019 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!  Content:
!      Intel(R) Math Kernel Library (Intel(R) MKL) Fortran 90 interface for
!	   Cluster Sparse Solver
!*******************************************************************************

!DEC$ IF .NOT. DEFINED( __MKL_CLUSTER_SPARSE_SOLVER_F90 )

!DEC$ DEFINE __MKL_CLUSTER_SPARSE_SOLVER_F90

      MODULE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
        TYPE MKL_CLUSTER_SPARSE_SOLVER_HANDLE; INTEGER(KIND=8) DUMMY; END TYPE
      END MODULE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE

      MODULE MKL_CLUSTER_SPARSE_SOLVER
        USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE

!
! Subroutine prototype for CLUSTER_SPARSE_SOLVER
!

      INTERFACE CLUSTER_SPARSE_SOLVER
        SUBROUTINE CLUSTER_SPARSE_SOLVER_D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          REAL(KIND=8),     INTENT(IN)    :: A(*)
          REAL(KIND=8),     INTENT(INOUT) :: B(*)
          REAL(KIND=8),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_S(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          REAL(KIND=4),     INTENT(IN)    :: A(*)
          REAL(KIND=4),     INTENT(INOUT) :: B(*)
          REAL(KIND=4),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_S

        SUBROUTINE CLUSTER_SPARSE_SOLVER_DC(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),     INTENT(INOUT) :: B(*)
          COMPLEX(KIND=8),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_DC

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SC(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=4),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),     INTENT(INOUT) :: B(*)
          COMPLEX(KIND=4),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SC

        SUBROUTINE CLUSTER_SPARSE_SOLVER_D_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          REAL(KIND=8),     INTENT(IN)    :: A(*)
          REAL(KIND=8),     INTENT(INOUT) :: B(N,*)
          REAL(KIND=8),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_D_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_S_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          REAL(KIND=4),     INTENT(IN)    :: A(*)
          REAL(KIND=4),     INTENT(INOUT) :: B(N,*)
          REAL(KIND=4),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_S_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),     INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=8),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=4),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),     INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=4),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_2D
      END INTERFACE

!
! Subroutine prototype for CLUSTER_SPARSE_SOLVER_64
!

      INTERFACE CLUSTER_SPARSE_SOLVER_64
        SUBROUTINE CLUSTER_SPARSE_SOLVER_D_64(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=8),     INTENT(IN)    :: A(*)
          REAL(KIND=8),     INTENT(INOUT) :: B(*)
          REAL(KIND=8),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_D_64

        SUBROUTINE CLUSTER_SPARSE_SOLVER_S_64(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=4),     INTENT(IN)    :: A(*)
          REAL(KIND=4),     INTENT(INOUT) :: B(*)
          REAL(KIND=4),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_S_64

        SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_64(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),     INTENT(INOUT) :: B(*)
          COMPLEX(KIND=8),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_64

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_64(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=4),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),     INTENT(INOUT) :: B(*)
          COMPLEX(KIND=4),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_64

        SUBROUTINE CLUSTER_SPARSE_SOLVER_D_64_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=8),     INTENT(IN)    :: A(*)
          REAL(KIND=8),     INTENT(INOUT) :: B(N,*)
          REAL(KIND=8),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_D_64_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_S_64_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=4),     INTENT(IN)    :: A(*)
          REAL(KIND=4),     INTENT(INOUT) :: B(N,*)
          REAL(KIND=4),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_S_64_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_64_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),     INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=8),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_64_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_64_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=4),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),     INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=4),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_64_2D
      END INTERFACE

!   operation names for the export routines
      INTEGER*4, PARAMETER :: SPARSE_PTLUQT  = 0
      INTEGER*4, PARAMETER :: SPARSE_DPTLUQT = 1

!   data names for the export routines
      INTEGER*4, PARAMETER :: SPARSE_PTLUQT_L  = 0
      INTEGER*4, PARAMETER :: SPARSE_PTLUQT_U  = 1
      INTEGER*4, PARAMETER :: SPARSE_PTLUQT_P  = 2
      INTEGER*4, PARAMETER :: SPARSE_PTLUQT_Q  = 3
      INTEGER*4, PARAMETER :: SPARSE_DPTLUQT_L = 4
      INTEGER*4, PARAMETER :: SPARSE_DPTLUQT_U = 5
      INTEGER*4, PARAMETER :: SPARSE_DPTLUQT_P = 6
      INTEGER*4, PARAMETER :: SPARSE_DPTLUQT_Q = 7
      INTEGER*4, PARAMETER :: SPARSE_DPTLUQT_D = 8

      INTERFACE CLUSTER_SPARSE_SOLVER_GET_CSR_SIZE
        SUBROUTINE CLUSTER_SPARSE_SOLVER_GET_CSR_SIZE_F(PT,DATANAME,LOCAL_NROWS,LOCAL_NNZ,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER*4,        INTENT(IN)    :: DATANAME
          INTEGER,          INTENT(OUT)   :: LOCAL_NROWS
          INTEGER,          INTENT(OUT)   :: LOCAL_NNZ
          INTEGER*4,        INTENT(IN)    :: COMM
          INTEGER,          INTENT(OUT)   :: ERROR
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_GET_CSR_SIZE_F
      END INTERFACE

      INTERFACE CLUSTER_SPARSE_SOLVER_SET_CSR_PTRS
        SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_CSR_PTRS_D(PT,DATANAME,ROWPTR,COLINDX,VALS,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER*4,        INTENT(IN)    :: DATANAME
          INTEGER,          INTENT(IN)    :: ROWPTR(*)
          INTEGER,          INTENT(IN)    :: COLINDX(*)
          REAL(KIND=8),     INTENT(IN)    :: VALS(*)
          INTEGER*4,        INTENT(IN)    :: COMM
          INTEGER,          INTENT(OUT)   :: ERROR
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_CSR_PTRS_D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_CSR_PTRS_S(PT,DATANAME,ROWPTR,COLINDX,VALS,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER*4,        INTENT(IN)    :: DATANAME
          INTEGER,          INTENT(IN)    :: ROWPTR(*)
          INTEGER,          INTENT(IN)    :: COLINDX(*)
          REAL(KIND=4),     INTENT(IN)    :: VALS(*)
          INTEGER*4,        INTENT(IN)    :: COMM
          INTEGER,          INTENT(OUT)   :: ERROR
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_CSR_PTRS_S

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_CSR_PTRS_DC(PT,DATANAME,ROWPTR,COLINDX,VALS,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER*4,        INTENT(IN)    :: DATANAME
          INTEGER,          INTENT(IN)    :: ROWPTR(*)
          INTEGER,          INTENT(IN)    :: COLINDX(*)
          COMPLEX(KIND=8),  INTENT(IN)    :: VALS(*)
          INTEGER*4,        INTENT(IN)    :: COMM
          INTEGER,          INTENT(OUT)   :: ERROR
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_CSR_PTRS_DC

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_CSR_PTRS_SC(PT,DATANAME,ROWPTR,COLINDX,VALS,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER*4,        INTENT(IN)    :: DATANAME
          INTEGER,          INTENT(IN)    :: ROWPTR(*)
          INTEGER,          INTENT(IN)    :: COLINDX(*)
          COMPLEX(KIND=4),  INTENT(IN)    :: VALS(*)
          INTEGER*4,        INTENT(IN)    :: COMM
          INTEGER,          INTENT(OUT)   :: ERROR
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_CSR_PTRS_SC

      END INTERFACE

      INTERFACE CLUSTER_SPARSE_SOLVER_SET_PTR
        SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_PTR_D(PT,DATANAME,PTR,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER*4,        INTENT(IN)    :: DATANAME
          REAL(KIND=8),     INTENT(IN)    :: PTR(*)
          INTEGER*4,        INTENT(IN)    :: COMM
          INTEGER,          INTENT(OUT)   :: ERROR
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_PTR_D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_PTR_S(PT,DATANAME,PTR,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER*4,        INTENT(IN)    :: DATANAME
          REAL(KIND=4),     INTENT(IN)    :: PTR(*)
          INTEGER*4,        INTENT(IN)    :: COMM
          INTEGER,          INTENT(OUT)   :: ERROR
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_PTR_S

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_PTR_DC(PT,DATANAME,PTR,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER*4,        INTENT(IN)    :: DATANAME
          COMPLEX(KIND=8),  INTENT(IN)    :: PTR(*)
          INTEGER*4,        INTENT(IN)    :: COMM
          INTEGER,          INTENT(OUT)   :: ERROR
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_PTR_DC

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_PTR_SC(PT,DATANAME,PTR,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER*4,        INTENT(IN)    :: DATANAME
          COMPLEX(KIND=4),  INTENT(IN)    :: PTR(*)
          INTEGER*4,        INTENT(IN)    :: COMM
          INTEGER,          INTENT(OUT)   :: ERROR
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_PTR_SC

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_PTR_I(PT,DATANAME,PTR,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER*4,        INTENT(IN)    :: DATANAME
          INTEGER,          INTENT(IN)    :: PTR(*)
          INTEGER*4,        INTENT(IN)    :: COMM
          INTEGER,          INTENT(OUT)   :: ERROR
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SET_PTR_I

      END INTERFACE

      INTERFACE CLUSTER_SPARSE_SOLVER_EXPORT
        SUBROUTINE CLUSTER_SPARSE_SOLVER_EXPORT_F(PT,OPERATION,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER*4,        INTENT(IN)    :: OPERATION
          INTEGER*4,        INTENT(IN)    :: COMM
          INTEGER,          INTENT(OUT)   :: ERROR
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_EXPORT_F
      END INTERFACE

      END MODULE MKL_CLUSTER_SPARSE_SOLVER

!DEC$ ENDIF
