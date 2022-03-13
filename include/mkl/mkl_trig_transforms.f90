!===============================================================================
! Copyright 2006-2019 Intel Corporation.
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
!      Intel(R) Math Kernel Library (Intel(R) MKL) interface for TT routines
!*******************************************************************************

      MODULE  MKL_TT_TYPE

! Parameters definitions for the kind of the Trigonometric Transform
      INTEGER, PARAMETER :: MKL_SINE_TRANSFORM               = 0
      INTEGER, PARAMETER :: MKL_COSINE_TRANSFORM             = 1
      INTEGER, PARAMETER :: MKL_STAGGERED_COSINE_TRANSFORM   = 2
      INTEGER, PARAMETER :: MKL_STAGGERED_SINE_TRANSFORM     = 3
      INTEGER, PARAMETER :: MKL_STAGGERED2_COSINE_TRANSFORM  = 4
      INTEGER, PARAMETER :: MKL_STAGGERED2_SINE_TRANSFORM    = 5

      END MODULE MKL_TT_TYPE

      MODULE  MKL_TRIG_TRANSFORMS

      USE MKL_TT_TYPE
      USE MKL_DFTI

      INTERFACE

        SUBROUTINE D_INIT_TRIG_TRANSFORM(n, tt_type, ipar,dpar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_d_init_trig_transform' :: D_INIT_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: n
            !MS$ATTRIBUTES REFERENCE :: tt_type
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: dpar
            !MS$ATTRIBUTES REFERENCE :: stat
            INTEGER, INTENT(IN) :: n, tt_type
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(8), INTENT(INOUT) :: dpar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE D_INIT_TRIG_TRANSFORM

        SUBROUTINE D_COMMIT_TRIG_TRANSFORM(f, handle, ipar,dpar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_d_commit_trig_transform' :: D_COMMIT_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: dpar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(8), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(8), INTENT(OUT) :: dpar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE D_COMMIT_TRIG_TRANSFORM

        SUBROUTINE D_FORWARD_TRIG_TRANSFORM(f, handle, ipar,dpar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_d_forward_trig_transform' :: D_FORWARD_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: dpar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(8), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(8), INTENT(IN) :: dpar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE D_FORWARD_TRIG_TRANSFORM

        SUBROUTINE D_BACKWARD_TRIG_TRANSFORM(f, handle, ipar,dpar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_d_backward_trig_transform' :: D_BACKWARD_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: dpar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(8), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(8), INTENT(IN) :: dpar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE D_BACKWARD_TRIG_TRANSFORM

        SUBROUTINE S_INIT_TRIG_TRANSFORM(n, tt_type, ipar,spar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_s_init_trig_transform' :: S_INIT_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: n
            !MS$ATTRIBUTES REFERENCE :: tt_type
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: spar
            !MS$ATTRIBUTES REFERENCE :: stat
            INTEGER, INTENT(IN) :: n, tt_type
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(4), INTENT(INOUT) :: spar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE S_INIT_TRIG_TRANSFORM

        SUBROUTINE S_COMMIT_TRIG_TRANSFORM(f, handle, ipar,spar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_s_commit_trig_transform' :: S_COMMIT_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: spar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(4), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(4), INTENT(OUT) :: spar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE S_COMMIT_TRIG_TRANSFORM

        SUBROUTINE S_FORWARD_TRIG_TRANSFORM(f, handle, ipar,spar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_s_forward_trig_transform' :: S_FORWARD_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: spar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(4), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(4), INTENT(IN) :: spar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE S_FORWARD_TRIG_TRANSFORM

        SUBROUTINE S_BACKWARD_TRIG_TRANSFORM(f, handle, ipar,spar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_s_backward_trig_transform' :: S_BACKWARD_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: spar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(4), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(4), INTENT(IN) :: spar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE S_BACKWARD_TRIG_TRANSFORM

        SUBROUTINE FREE_TRIG_TRANSFORM(handle, ipar,stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_free_trig_transform' :: FREE_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: stat
            INTEGER, INTENT(INOUT) :: ipar(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE FREE_TRIG_TRANSFORM

      END INTERFACE

      END MODULE MKL_TRIG_TRANSFORMS
