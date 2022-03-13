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
!      Intel(R) Math Kernel Library (Intel(R) MKL) FORTRAN interface for JIT
!      BLAS routines
!*******************************************************************************

  MODULE MKL_JIT_BLAS

    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INTPTR_T, C_PTR, C_INT, C_FUNPTR, &
                                             C_FLOAT, C_DOUBLE, C_FLOAT_COMPLEX, C_DOUBLE_COMPLEX

!   return status of the routines
    ENUM, BIND(C)
       ENUMERATOR :: MKL_JIT_SUCCESS             = 0,  &  ! jitter was created and kernel jitted
                     MKL_NO_JIT                  = 1,  &  ! jitter was created but no kernel jitted, will use standard GEMM
                     MKL_JIT_ERROR               = 2      ! jitter was not created
          
    END ENUM

!     define corresponding fortran type of jit_get_?gemm_ptr returned function pointer 
    
    ABSTRACT INTERFACE
       subroutine sgemm_jit_kernel_t ( jitter, a, b, c ) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_float, c_intptr_t, c_ptr
         TYPE(C_PTR), INTENT(IN), VALUE  :: jitter
         REAL :: a(*), b(*), c(*)
       end subroutine sgemm_jit_kernel_t

       subroutine dgemm_jit_kernel_t ( jitter, a, b, c ) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_double, c_intptr_t, c_ptr
         TYPE(C_PTR), INTENT(IN), VALUE  :: jitter
         DOUBLE PRECISION :: a(*), b(*), c(*)
       end subroutine dgemm_jit_kernel_t

       subroutine cgemm_jit_kernel_t ( jitter, a, b, c ) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_float_complex, c_intptr_t, c_ptr
         TYPE(C_PTR), INTENT(IN), VALUE :: jitter
         COMPLEX :: a(*), b(*), c(*)
       end subroutine cgemm_jit_kernel_t
       
       subroutine zgemm_jit_kernel_t ( jitter, a, b, c ) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_double_complex, c_intptr_t, c_ptr
         TYPE(C_PTR), INTENT(IN), VALUE :: jitter
         DOUBLE COMPLEX :: a(*), b(*), c(*)
       end subroutine zgemm_jit_kernel_t
    END INTERFACE

!   JIT API interface

    INTERFACE

!      create a jitter, store it in first argument, generate the corresponding GEMM kernel (can be a call to standard GEMM), return status is either MKL_JIT_ERROR, MKL_JIT_SUCCESS, MKL_NO_JIT
       function mkl_jit_create_dgemm ( jitter, transa, transb, m, n, k, alpha, lda, ldb, beta, ldc ) RESULT (status) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_intptr_t, c_ptr, c_int, c_double
         integer :: status
         TYPE(C_PTR) :: jitter
         character*1 :: transa, transb
         integer :: m, n, k, lda, ldb, ldc
         double precision :: alpha, beta
       END function mkl_jit_create_dgemm

       function mkl_jit_create_sgemm ( jitter, transa, transb, m, n, k, alpha, lda, ldb, beta, ldc ) RESULT (status) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_intptr_t, c_ptr, c_int, c_float
         integer :: status
         TYPE(C_PTR) :: jitter
         character*1 :: transa, transb
         integer :: m, n, k, lda, ldb, ldc
         real :: alpha, beta
       END function mkl_jit_create_sgemm

       function mkl_jit_create_cgemm ( jitter, transa, transb, m, n, k, alpha, lda, ldb, beta, ldc ) RESULT (status) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_intptr_t, c_ptr, c_int, c_float_complex
         integer :: status
         TYPE(C_PTR) :: jitter
         character*1 :: transa, transb
         integer :: m, n, k, lda, ldb, ldc
         complex :: alpha, beta
       END function mkl_jit_create_cgemm

       function mkl_jit_create_zgemm ( jitter, transa, transb, m, n, k, alpha, lda, ldb, beta, ldc ) RESULT (status) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_intptr_t, c_ptr, c_int, c_double_complex
         integer :: status
         TYPE(C_PTR) :: jitter
         character*1 :: transa, transb
         integer :: m, n, k, lda, ldb, ldc
         double complex :: alpha, beta
       END function mkl_jit_create_zgemm

       ! destroy jitter and free memory, return status is either MKL_JIT_SUCCESS or MKL_JIT_ERROR (if given pointer is not a handle on a jitter)
       function mkl_jit_destroy ( jitter ) RESULT (status) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_intptr_t, c_ptr, c_int
         integer :: status
         TYPE(C_PTR), INTENT(IN), VALUE :: jitter
       END function mkl_jit_destroy
       
       ! return a C procedure pointer to the generated GEMM kernel
       ! this pointer needs to be converted to a Fortran procedure pointer using ?gemm_kernel_t interfaces above 

       function mkl_jit_get_dgemm_ptr ( jitter ) RESULT (ptr) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_funptr, c_intptr_t, c_ptr
         TYPE(C_FUNPTR) :: ptr
         TYPE(C_PTR), INTENT(IN), VALUE :: jitter
       END function mkl_jit_get_dgemm_ptr

       function mkl_jit_get_sgemm_ptr ( jitter ) RESULT (ptr) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_funptr, c_intptr_t, c_ptr
         TYPE(C_FUNPTR) :: ptr
         TYPE(C_PTR), INTENT(IN), VALUE :: jitter
       END function mkl_jit_get_sgemm_ptr

       function mkl_jit_get_cgemm_ptr ( jitter ) RESULT (ptr) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_funptr, c_intptr_t, c_ptr
         TYPE(C_FUNPTR) :: ptr
         TYPE(C_PTR), INTENT(IN), VALUE :: jitter
       END function mkl_jit_get_cgemm_ptr

       function mkl_jit_get_zgemm_ptr ( jitter ) RESULT (ptr) BIND(C)
         use, intrinsic :: ISO_C_BINDING, only : c_funptr, c_intptr_t, c_ptr
         TYPE(C_FUNPTR) :: ptr
         TYPE(C_PTR), INTENT(IN), VALUE :: jitter
       END function mkl_jit_get_zgemm_ptr

    END INTERFACE

  END MODULE MKL_JIT_BLAS
