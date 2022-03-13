!===============================================================================
! Copyright 2014-2019 Intel Corporation.
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

!
!   This file contains preliminary version of new SpBLAS API which supports
!   two-step execution (inspector-executor) model.
!
!*******************************************************************************

!*******************************************************************************
!*********************************** Basic types and constants *****************
!*******************************************************************************

    MODULE MKL_SOLVERS_EE

    USE MKL_SPBLAS
    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INTPTR_T, C_INT

!*************************************************************************************************
!*** Opaque structure for sparse matrix in internal format, further D - means double precision ***
!*************************************************************************************************

!    struct  sparse_matrix;
!    typedef struct sparse_matrix *sparse_matrix_t;

!    TYPE, BIND(C) :: SPARSE_MATRIX_T
!        INTEGER(C_INTPTR_T) :: PTR
!    END TYPE SPARSE_MATRIX_T

!   descriptor of main sparse matrix properties
!    TYPE, BIND(C) :: MATRIX_DESCR
!        INTEGER(C_INT) :: TYPE
!        INTEGER(C_INT) :: MODE
!        INTEGER(C_INT) :: DIAG
!    END TYPE MATRIX_DESCR

    INTERFACE

!****************************************************************************************
!****************************** Computational routines **********************************
!****************************************************************************************
        FUNCTION MKL_SPARSE_EE_INIT(pm) &
                 BIND(C, name='MKL_SPARSE_EE_INIT')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            INTEGER, INTENT(IN), DIMENSION(*) :: pm
            INTEGER(C_INT) MKL_SPARSE_EE_INIT
        END FUNCTION

        FUNCTION MKL_SPARSE_S_EV(which, pm, A, descrA, k0, k, E, X, res) &
                 BIND(C, name='MKL_SPARSE_S_EV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            CHARACTER           , INTENT(IN) :: which
            INTEGER, INTENT(IN), DIMENSION(*) :: pm
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descrA
            INTEGER, INTENT(IN) :: k0
            INTEGER, INTENT(INOUT) :: k
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: E
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: X
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: res
            INTEGER(C_INT) MKL_SPARSE_S_EV
        END FUNCTION

        FUNCTION MKL_SPARSE_D_EV(which, pm, A, descrA, k0, k, E, X, res) &
                 BIND(C, name='MKL_SPARSE_D_EV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            CHARACTER           , INTENT(IN) :: which
            INTEGER, INTENT(IN), DIMENSION(*) :: pm
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descrA
            INTEGER, INTENT(IN) :: k0
            INTEGER, INTENT(INOUT) :: k
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: E
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: X
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: res
            INTEGER(C_INT) MKL_SPARSE_D_EV
        END FUNCTION

        FUNCTION MKL_SPARSE_D_GV(which, pm, A, descrA, B, descrB, k0, k, E, X, res) &
                 BIND(C, name='MKL_SPARSE_D_GV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            CHARACTER           , INTENT(IN) :: which
            INTEGER, INTENT(IN), DIMENSION(*) :: pm
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descrA
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descrB
            INTEGER, INTENT(IN) :: k0
            INTEGER, INTENT(INOUT) :: k
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: E
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: X
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: res
            INTEGER(C_INT) MKL_SPARSE_D_GV
        END FUNCTION

        FUNCTION MKL_SPARSE_S_GV(which, pm, A, descrA, B, descrB, k0, k, E, X, res) &
                 BIND(C, name='MKL_SPARSE_S_GV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            CHARACTER           , INTENT(IN) :: which
            INTEGER, INTENT(IN), DIMENSION(*) :: pm
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descrA
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descrB
            INTEGER, INTENT(IN) :: k0
            INTEGER, INTENT(INOUT) :: k
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: E
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: X
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: res
            INTEGER(C_INT) MKL_SPARSE_S_GV
        END FUNCTION

        FUNCTION MKL_SPARSE_S_SVD(whichE, whichV, pm, A, descrA, k0, k, E, XL, XR, res) &
                 BIND(C, name='MKL_SPARSE_S_SVD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            CHARACTER           , INTENT(IN) :: whichE
            CHARACTER           , INTENT(IN) :: whichV            
            INTEGER, INTENT(IN), DIMENSION(*) :: pm
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descrA
            INTEGER, INTENT(IN) :: k0
            INTEGER, INTENT(INOUT) :: k
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: E
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: XL
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: XR
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: res
            INTEGER(C_INT) MKL_SPARSE_S_SVD
        END FUNCTION

        FUNCTION MKL_SPARSE_D_SVD(whichE, whichV, pm, A, descrA, k0, k, E, XL, XR, res) &
                 BIND(C, name='MKL_SPARSE_D_SVD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            CHARACTER           , INTENT(IN) :: whichE
            CHARACTER           , INTENT(IN) :: whichV
            INTEGER, INTENT(IN), DIMENSION(*) :: pm
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descrA
            INTEGER, INTENT(IN) :: k0
            INTEGER, INTENT(INOUT) :: k
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: E
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: XL
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: XR
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: res
            INTEGER(C_INT) MKL_SPARSE_D_SVD
        END FUNCTION


    END INTERFACE

    END MODULE MKL_SOLVERS_EE
