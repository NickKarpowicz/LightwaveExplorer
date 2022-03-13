!===============================================================================
! Copyright 2018-2019 Intel Corporation.
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
!   Content:
!           Intel(R) Math Kernel Library (Intel(R) MKL) Sparse QR FORTRAN header file
!
!           Contains interface to: MKL_SPARSE_X_QR
!                                  MKL_SPARSE_QR_REORDER
!                                  MKL_SPARSE_X_QR_FACTORIZE
!                                  MKL_SPARSE_X_QR_SOLVE
!                                  MKL_SPARSE_X_QR_QMULT
!                                  MKL_SPARSE_X_QR_RSOLVE
!                                  MKL_SPARSE_SET_QR_HINT
!
!===============================================================================

MODULE MKL_SPARSE_QR

    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INTPTR_T, C_INT
    USE MKL_SPBLAS, ONLY : SPARSE_MATRIX_T, MATRIX_DESCR

    ENUM, BIND(C)
        ENUMERATOR :: SPARSE_QR_WITH_PIVOTS
    END ENUM

    INTERFACE

    FUNCTION MKL_SPARSE_SET_QR_HINT(A, HINT) &
        BIND(C, name='MKL_SPARSE_SET_QR_HINT')
        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
        IMPORT SPARSE_MATRIX_T
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
        INTEGER(C_INT)       , INTENT(IN)    :: HINT
        INTEGER(C_INT) MKL_SPARSE_SET_QR_HINT
    END FUNCTION

    FUNCTION MKL_SPARSE_D_QR(operation, A, descr, layout, columns, x, ldx, b, ldb) &
        BIND(C, name='MKL_SPARSE_D_QR')

        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
        IMPORT SPARSE_MATRIX_T
        IMPORT MATRIX_DESCR
        INTEGER(C_INT)       , INTENT(IN)               :: operation
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT)            :: A
        TYPE(MATRIX_DESCR)   , INTENT(IN)               :: descr
        INTEGER(C_INT)       , INTENT(IN)               :: layout
        INTEGER              , INTENT(IN)               :: columns
        REAL(C_DOUBLE)       , INTENT(OUT), DIMENSION(*):: x
        INTEGER              , INTENT(IN)               :: ldx
        REAL(C_DOUBLE)       , INTENT(IN) , DIMENSION(*):: b
        INTEGER              , INTENT(IN)               :: ldb

        INTEGER(C_INT) MKL_SPARSE_D_QR
    END FUNCTION

    FUNCTION MKL_SPARSE_S_QR(operation,A,descr, layout, columns, x, ldx, b, ldb) &
        BIND(C, name='MKL_SPARSE_S_QR')

        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
        IMPORT SPARSE_MATRIX_T
        IMPORT MATRIX_DESCR
        INTEGER(C_INT)       , INTENT(IN)               :: operation
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT)            :: A
        TYPE(MATRIX_DESCR)   , INTENT(IN)               :: descr
        INTEGER(C_INT)       , INTENT(IN)               :: layout
        INTEGER              , INTENT(IN)               :: columns
        REAL(C_FLOAT)        , INTENT(OUT), DIMENSION(*):: x
        INTEGER              , INTENT(IN)               :: ldx
        REAL(C_FLOAT)        , INTENT(IN) , DIMENSION(*):: b
        INTEGER              , INTENT(IN)               :: ldb

        INTEGER(C_INT) MKL_SPARSE_S_QR
    END FUNCTION

    FUNCTION MKL_SPARSE_QR_REORDER(A, descr) &
        BIND(C, name='MKL_SPARSE_QR_REORDER')

        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
        IMPORT SPARSE_MATRIX_T
        IMPORT MATRIX_DESCR
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT)            :: A
        TYPE(MATRIX_DESCR)   , INTENT(IN)               :: descr

        INTEGER(C_INT) MKL_SPARSE_QR_REORDER
    END FUNCTION

    FUNCTION MKL_SPARSE_D_QR_FACTORIZE(A, alt_values) &
        BIND(C, name='MKL_SPARSE_D_QR_FACTORIZE')
        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
        IMPORT SPARSE_MATRIX_T
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT)            :: A
        REAL(C_DOUBLE)       , INTENT(IN) , DIMENSION(*):: alt_values

        INTEGER(C_INT) MKL_SPARSE_D_QR_FACTORIZE
    END FUNCTION

    FUNCTION MKL_SPARSE_S_QR_FACTORIZE(A, alt_values) &
        BIND(C, name='MKL_SPARSE_S_QR_FACTORIZE')
        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
        IMPORT SPARSE_MATRIX_T
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT)            :: A
        REAL(C_FLOAT)        , INTENT(IN) , DIMENSION(*):: alt_values

        INTEGER(C_INT) MKL_SPARSE_S_QR_FACTORIZE
    END FUNCTION

    FUNCTION MKL_SPARSE_D_QR_SOLVE(operation, A, alt_values, layout, columns, x, ldx, b, ldb) &
        BIND(C, name='MKL_SPARSE_D_QR_SOLVE')

        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
        IMPORT SPARSE_MATRIX_T
        IMPORT MATRIX_DESCR
        INTEGER(C_INT)       , INTENT(IN)               :: operation
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT)            :: A
        REAL(C_DOUBLE)       , INTENT(OUT), DIMENSION(*):: alt_values
        INTEGER(C_INT)       , INTENT(IN)               :: layout
        INTEGER              , INTENT(IN)               :: columns
        REAL(C_DOUBLE)       , INTENT(OUT), DIMENSION(*):: x
        INTEGER              , INTENT(IN)               :: ldx
        REAL(C_DOUBLE)       , INTENT(IN) , DIMENSION(*):: b
        INTEGER              , INTENT(IN)               :: ldb

        INTEGER(C_INT) MKL_SPARSE_D_QR_SOLVE
    END FUNCTION

    FUNCTION MKL_SPARSE_S_QR_SOLVE(operation, A, alt_values, layout, columns, x, ldx, b, ldb) &
        BIND(C, name='MKL_SPARSE_S_QR_SOLVE')

        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
        IMPORT SPARSE_MATRIX_T
        IMPORT MATRIX_DESCR
        INTEGER(C_INT)       , INTENT(IN)               :: operation
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT)            :: A
        REAL(C_FLOAT)        , INTENT(OUT), DIMENSION(*):: alt_values
        INTEGER(C_INT)       , INTENT(IN)               :: layout
        INTEGER              , INTENT(IN)               :: columns
        REAL(C_FLOAT)        , INTENT(OUT), DIMENSION(*):: x
        INTEGER              , INTENT(IN)               :: ldx
        REAL(C_FLOAT)        , INTENT(IN) , DIMENSION(*):: b
        INTEGER              , INTENT(IN)               :: ldb

        INTEGER(C_INT) MKL_SPARSE_S_QR_SOLVE
    END FUNCTION

    FUNCTION MKL_SPARSE_D_QR_QMULT(operation, A, layout, columns, x, ldx, b, ldb) &
        BIND(C, name='MKL_SPARSE_D_QR_QMULT')

        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
        IMPORT SPARSE_MATRIX_T
        IMPORT MATRIX_DESCR
        INTEGER(C_INT)       , INTENT(IN)               :: operation
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT)            :: A
        INTEGER(C_INT)       , INTENT(IN)               :: layout
        INTEGER              , INTENT(IN)               :: columns
        REAL(C_DOUBLE)       , INTENT(OUT), DIMENSION(*):: x
        INTEGER              , INTENT(IN)               :: ldx
        REAL(C_DOUBLE)       , INTENT(IN) , DIMENSION(*):: b
        INTEGER              , INTENT(IN)               :: ldb

        INTEGER(C_INT) MKL_SPARSE_D_QR_QMULT
    END FUNCTION

    FUNCTION MKL_SPARSE_S_QR_QMULT(operation, A, layout, columns, x, ldx, b, ldb) &
        BIND(C, name='MKL_SPARSE_S_QR_QMULT')

        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
        IMPORT SPARSE_MATRIX_T
        IMPORT MATRIX_DESCR
        INTEGER(C_INT)       , INTENT(IN)               :: operation
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT)            :: A
        INTEGER(C_INT)       , INTENT(IN)               :: layout
        INTEGER              , INTENT(IN)               :: columns
        REAL(C_FLOAT)        , INTENT(OUT), DIMENSION(*):: x
        INTEGER              , INTENT(IN)               :: ldx
        REAL(C_FLOAT)        , INTENT(IN) , DIMENSION(*):: b
        INTEGER              , INTENT(IN)               :: ldb

        INTEGER(C_INT) MKL_SPARSE_S_QR_QMULT
    END FUNCTION

    FUNCTION MKL_SPARSE_D_QR_RSOLVE(operation, A, layout, columns, x, ldx, b, ldb) &
        BIND(C, name='MKL_SPARSE_D_QR_RSOLVE')

        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
        IMPORT SPARSE_MATRIX_T
        IMPORT MATRIX_DESCR
        INTEGER(C_INT)       , INTENT(IN)               :: operation
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT)            :: A
        INTEGER(C_INT)       , INTENT(IN)               :: layout
        INTEGER              , INTENT(IN)               :: columns
        REAL(C_DOUBLE)       , INTENT(OUT), DIMENSION(*):: x
        INTEGER              , INTENT(IN)               :: ldx
        REAL(C_DOUBLE)       , INTENT(IN) , DIMENSION(*):: b
        INTEGER              , INTENT(IN)               :: ldb

        INTEGER(C_INT) MKL_SPARSE_D_QR_RSOLVE
    END FUNCTION

    FUNCTION MKL_SPARSE_S_QR_RSOLVE(operation, A, layout, columns, x, ldx, b, ldb) &
        BIND(C, name='MKL_SPARSE_S_QR_RSOLVE')

        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
        IMPORT SPARSE_MATRIX_T
        IMPORT MATRIX_DESCR
        INTEGER(C_INT)       , INTENT(IN)               :: operation
        TYPE(SPARSE_MATRIX_T), INTENT(INOUT)            :: A
        INTEGER(C_INT)       , INTENT(IN)               :: layout
        INTEGER              , INTENT(IN)               :: columns
        REAL(C_FLOAT)        , INTENT(OUT), DIMENSION(*):: x
        INTEGER              , INTENT(IN)               :: ldx
        REAL(C_FLOAT)        , INTENT(IN) , DIMENSION(*):: b
        INTEGER              , INTENT(IN)               :: ldb

        INTEGER(C_INT) MKL_SPARSE_S_QR_RSOLVE
    END FUNCTION

    END INTERFACE

END MODULE
