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

MODULE MKL_POISSON

USE MKL_DFTI

INTERFACE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!INTERFACES FOR 3D CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE D_INIT_HELMHOLTZ_3D(AX,BX,AY,BY,AZ,BZ,NX,NY,NZ,BCTYPE,Q,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, NZ, STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION AX,BX,AY,BY,AZ,BZ,Q
   CHARACTER(6) BCTYPE
   DOUBLE PRECISION DPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_COMMIT_HELMHOLTZ_3D(F,BD_AX,BD_BX,BD_AY,BD_BY,BD_AZ,BD_BZ,XHANDLE,YHANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION DPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,IPAR(12)+1,*)
   DOUBLE PRECISION BD_AX(IPAR(12)+1,*),BD_BX(IPAR(12)+1,*),BD_AY(IPAR(11)+1,*),BD_BY(IPAR(11)+1,*)
   DOUBLE PRECISION BD_AZ(IPAR(11)+1,*),BD_BZ(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: XHANDLE, YHANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_HELMHOLTZ_3D(F,BD_AX,BD_BX,BD_AY,BD_BY,BD_AZ,BD_BZ,XHANDLE,YHANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,IPAR(12)+1,*)
   DOUBLE PRECISION BD_AX(IPAR(12)+1,*),BD_BX(IPAR(12)+1,*),BD_AY(IPAR(11)+1,*),BD_BY(IPAR(11)+1,*)
   DOUBLE PRECISION BD_AZ(IPAR(11)+1,*),BD_BZ(IPAR(11)+1,*)
   DOUBLE PRECISION DPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: XHANDLE, YHANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE FREE_HELMHOLTZ_3D(XHANDLE,YHANDLE,IPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: XHANDLE, YHANDLE
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!INTERFACES FOR 2D CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE D_INIT_HELMHOLTZ_2D(AX,BX,AY,BY,NX,NY,BCTYPE,Q,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION AX,BX,AY,BY,Q
   CHARACTER(4) BCTYPE
   DOUBLE PRECISION DPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_COMMIT_HELMHOLTZ_2D(F,BD_AX,BD_BX,BD_AY,BD_BY,HANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   DOUBLE PRECISION BD_AX(*),BD_BX(*),BD_AY(*),BD_BY(*)
   DOUBLE PRECISION DPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_HELMHOLTZ_2D(F,BD_AX,BD_BX,BD_AY,BD_BY,HANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   DOUBLE PRECISION BD_AX(*),BD_BX(*),BD_AY(*),BD_BY(*)
   DOUBLE PRECISION DPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE FREE_HELMHOLTZ_2D(HANDLE,IPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE
!**********************************************************************
!***********************SINGLE*****************************************
!**********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!INTERFACES FOR 3D CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE S_INIT_HELMHOLTZ_3D(AX,BX,AY,BY,AZ,BZ,NX,NY,NZ,BCTYPE,Q,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, NZ, STAT
   INTEGER IPAR(*)
   REAL AX,BX,AY,BY,AZ,BZ,Q
   CHARACTER(6) BCTYPE
   REAL SPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_COMMIT_HELMHOLTZ_3D(F,BD_AX,BD_BX,BD_AY,BD_BY,BD_AZ,BD_BZ,XHANDLE,YHANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL SPAR(*)
   REAL F(IPAR(11)+1,IPAR(12)+1,*)
   REAL BD_AX(IPAR(12)+1,*),BD_BX(IPAR(12)+1,*),BD_AY(IPAR(11)+1,*),BD_BY(IPAR(11)+1,*)
   REAL BD_AZ(IPAR(11)+1,*),BD_BZ(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: XHANDLE, YHANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_HELMHOLTZ_3D(F,BD_AX,BD_BX,BD_AY,BD_BY,BD_AZ,BD_BZ,XHANDLE,YHANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL F(IPAR(11)+1,IPAR(12)+1,*)
   REAL BD_AX(IPAR(12)+1,*),BD_BX(IPAR(12)+1,*),BD_AY(IPAR(11)+1,*),BD_BY(IPAR(11)+1,*)
   REAL BD_AZ(IPAR(11)+1,*),BD_BZ(IPAR(11)+1,*)
   REAL SPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: XHANDLE, YHANDLE
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!INTERFACES FOR 2D CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE S_INIT_HELMHOLTZ_2D(AX,BX,AY,BY,NX,NY,BCTYPE,Q,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   REAL AX,BX,AY,BY,Q
   CHARACTER(4) BCTYPE
   REAL SPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_COMMIT_HELMHOLTZ_2D(F,BD_AX,BD_BX,BD_AY,BD_BY,HANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL F(IPAR(11)+1,*)
   REAL BD_AX(*),BD_BX(*),BD_AY(*),BD_BY(*)
   REAL SPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_HELMHOLTZ_2D(F,BD_AX,BD_BX,BD_AY,BD_BY,HANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL F(IPAR(11)+1,*)
   REAL BD_AX(*),BD_BX(*),BD_AY(*),BD_BY(*)
   REAL SPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

SUBROUTINE D_INIT_SPH_P(AX,BX,AY,BY,NX,NY,Q,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION AX,BX,AY,BY,Q
   DOUBLE PRECISION DPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_COMMIT_SPH_P(F,HANDLE_S,HANDLE_C,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION DPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE_C, HANDLE_S
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_SPH_P(F,HANDLE_S,HANDLE_C,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION DPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE_C, HANDLE_S
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE FREE_SPH_P(HANDLE_S,HANDLE_C,IPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE_S, HANDLE_C
END SUBROUTINE

!---------------------------------------------------------------------
SUBROUTINE D_INIT_SPH_NP(AX,BX,AY,BY,NX,NY,Q,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION AX,BX,AY,BY,Q
   DOUBLE PRECISION DPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_COMMIT_SPH_NP(F,HANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION DPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_SPH_NP(F,HANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION DPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE FREE_SPH_NP(HANDLE,IPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!======================================================================

SUBROUTINE S_INIT_SPH_P(AX,BX,AY,BY,NX,NY,Q,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   REAL AX,BX,AY,BY,Q
   REAL SPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_COMMIT_SPH_P(F,HANDLE_S,HANDLE_C,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL SPAR(*)
   REAL F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE_C, HANDLE_S
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_SPH_P(F,HANDLE_S,HANDLE_C,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL SPAR(*)
   REAL F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE_C, HANDLE_S
END SUBROUTINE

!--------------------------------------------------------------------

SUBROUTINE S_INIT_SPH_NP(AX,BX,AY,BY,NX,NY,Q,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   REAL AX,BX,AY,BY,Q
   REAL SPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_COMMIT_SPH_NP(F,HANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL SPAR(*)
   REAL F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_SPH_NP(F,HANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL SPAR(*)
   REAL F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE


END INTERFACE
END MODULE
