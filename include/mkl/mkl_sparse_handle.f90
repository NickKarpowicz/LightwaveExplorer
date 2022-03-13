!===============================================================================
! Copyright 2004-2019 Intel Corporation.
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

!   Content:
!           Intel(R) Math Kernel Library (Intel(R) MKL) DSS Fortran header file
!
!           Contains more detailed information on internal datatypes and
!           constants used by DSS interface to PARDISO.
!
!*******************************************************************************

      MODULE MKL_SPARSE_HANDLE

      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INTPTR_T

      ENUM, BIND(C)
          ENUMERATOR :: MKL_ZERO_BASED, MKL_ONE_BASED
      END ENUM

      ENUM, BIND(C)
          ENUMERATOR :: MKL_C_STYLE, MKL_FORTRAN_STYLE
      END ENUM

      ENUM, BIND(C)
          ENUMERATOR :: MKL_NO_PRINT, MKL_PRINT
      END ENUM

      ENUM, BIND(C)
          ENUMERATOR :: MKL_GENERAL_STRUCTURE, MKL_UPPER_TRIANGULAR, MKL_LOWER_TRIANGULAR, MKL_STRUCTURAL_SYMMETRIC
      END ENUM

      ENUM, BIND(C)
          ENUMERATOR :: MKL_CSR
      END ENUM

      TYPE, BIND(C) :: SPARSE_STRUCT
          INTEGER N
          INTEGER (C_INTPTR_T) :: CSR_IA
          INTEGER (C_INTPTR_T) :: CSR_JA
          INTEGER CHECK_RESULT(3)
          INTEGER(KIND=4) INDEXING
          INTEGER(KIND=4) MATRIX_STRUCTURE
          INTEGER(KIND=4) MATRIX_FORMAT
          INTEGER(KIND=4) MESSAGE_LEVEL
          INTEGER(KIND=4) PRINT_STYLE
      END TYPE SPARSE_STRUCT

      INTERFACE

          FUNCTION sparse_matrix_checker(PT)
              IMPORT SPARSE_STRUCT
              TYPE(SPARSE_STRUCT), INTENT(INOUT) :: PT
              INTEGER sparse_matrix_checker
          END

          SUBROUTINE sparse_matrix_checker_init(PT)
              IMPORT SPARSE_STRUCT
              TYPE(SPARSE_STRUCT), INTENT(INOUT) :: PT
          END

      END INTERFACE

      INTEGER, PARAMETER :: MKL_SPARSE_CHECKER_SUCCESS = 0
      INTEGER, PARAMETER :: MKL_SPARSE_CHECKER_NON_MONOTONIC = 21
      INTEGER, PARAMETER :: MKL_SPARSE_CHECKER_OUT_OF_RANGE = 22
      INTEGER, PARAMETER :: MKL_SPARSE_CHECKER_NONTRIANGULAR = 23
      INTEGER, PARAMETER :: MKL_SPARSE_CHECKER_NONORDERED = 24


      END MODULE MKL_SPARSE_HANDLE
