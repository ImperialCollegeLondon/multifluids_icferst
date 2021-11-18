!----------------------------------------------------------------------------------------
!> @author Vinicius L S Silva
!> @brief Interface to call XGBoost library C API from fortran
!> This code is a inspired in https://github.com/christophkeller/fortran2xgb 
!----------------------------------------------------------------------------------------

module xgb_interface

   use iso_c_binding

   implicit none

   interface
   ! Create matrix content from dense matrix
   ! Returns 0 when success, -1 when failure happens
      integer(c_int) function fortran_XGDMatrixCreateFromMat(data, nrow, ncol, missing, out) bind(c, name="XGDMatrixCreateFromMat")
         use iso_c_binding, only: c_int, c_ptr, c_float, c_int64_t
         ! Parameters
         real(c_float)             :: data(*) ! Pointer to the data space
         integer(c_int64_t), value :: nrow    ! Number of rows
         integer(c_int64_t), value :: ncol    ! Number columns
         real(c_float), value      :: missing ! Which value to represent missing value
         type(c_ptr)               :: out     ! Created dmatrix
         ! End parameters 
      end function
   end interface    

   interface
   ! Create xgboost learner
   ! Returns 0 when success, -1 when failure happens
      integer(c_int) function fortran_XGBoosterCreate(dmats, len, out) bind(c, name="XGBoosterCreate")
         use iso_c_binding, only: c_int, c_ptr, c_int64_t
         ! Parameters
         type(c_ptr), value        :: dmats ! Matrices that are set to be cached
         integer(c_int64_t), value :: len   ! Length of dmats
         type(c_ptr)               :: out   ! Handle to the result booster
         ! End parameters 
      end function
   end interface   
   
   interface
   ! Load model from existing file.
   ! Returns 0 when success, -1 when failure happens
      integer(c_int) function fortran_XGBoosterLoadModel(handle, fname) bind(C, name="XGBoosterLoadModel")
         use iso_c_binding, only: c_int, c_char, c_ptr
         ! Parameters
         type(c_ptr), value                          :: handle ! Handle
         character(len=1, kind=c_char), dimension(*) :: fname  ! File URI or file name (*needs to be in the right c format - use trim(name_in_fortran)//c_null_char)
         ! End parameters 
      end function
   end interface
   
   interface
   ! Make prediction based on dmatrix
   ! Returns 0 when success, -1 when failure happens
      integer(c_int) function fortran_XGBoosterPredict(handle, dmat, option_mask, ntree_limit, out_len, out_result) bind(C, name="XGBoosterPredict")
         use iso_c_binding, only: c_int, c_ptr, c_int64_t
         ! Parameters
         type(c_ptr), value        :: handle       ! Handle
         type(c_ptr), value        :: dmat         ! Data matrix
         integer(c_int), value     :: option_mask  ! Bit-mask of options taken in prediction, possible values 0:normal prediction 1:output margin instead of transformed value 2:output leaf index of trees instead of leaf value, note leaf index is unique per tree 4:output feature contributions to individual predictions
         integer(c_int), value     :: ntree_limit  ! Limit number of trees used for prediction, this is only valid for boosted trees when the parameter is set to 0, we will use all the trees
         integer(c_int64_t)        :: out_len      ! Used to store length of returning result
         type(c_ptr)               :: out_result   ! Used to set a pointer to array
         ! End parameters 
      end function
   end interface
   
   interface
   ! Free obj in handle
   ! Returns 0 when success, -1 when failure happens
      integer(c_int) function fortran_XGBoosterFree(handle) bind(c, name="XGBoosterFree")
         use iso_c_binding, only: c_int, c_ptr
         ! Parameters
         type(c_ptr), value :: handle ! Handle to be freed
         ! End parameters 
      end function
   end interface
   
end module