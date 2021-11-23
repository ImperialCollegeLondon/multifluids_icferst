!    Copyright (C) 2020 Imperial College London and others.
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Affero General Public License
!    as published by the Free Software Foundation,
!    version 3.0 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without seven the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

!----------------------------------------------------------------------------------------
!> @author Vinicius L S Silva
!> @brief Module to load and call a XGBoost model 
!----------------------------------------------------------------------------------------
module multi_machine_learning

  use iso_c_binding
  use xgb_interface
  use spud

  implicit none

  ! Module private variables
  type(c_ptr), private, save   :: xgb_model ! XGBoost model 
  character(len=255), private :: name_xgb_model ! Path for XGBoost model (xgb_model.bin if not defined in Diamond)
  integer(c_int64_t), parameter, private :: nrow = 1 ! Number of rows in the model input (usually 1)
  integer(c_int64_t), parameter, private :: ncol = 17 ! Number of colunms (features) in the model input

  contains
    
    !----------------------------------------------------------------------------------------
    !> @author Vinicius L S Silva
    !> @brief Load the XGBoost model as -> xgb_model (private module variable)
    !----------------------------------------------------------------------------------------
    subroutine xgboost_load_model()

      implicit none
      
      ! Local variables
      type(c_ptr)                :: dmatrix         
      real(c_float), allocatable :: xgb_input(:)    
      integer                    :: error  
      logical                    :: there = .false.
      
      ! Constants
      real(c_float), parameter      :: missing = -999.0
      integer(c_int64_t), parameter :: dmatrix_len = 0

      call get_option("/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/ML_model_path", name_xgb_model, default='xgb_model.bin')
      inquire( FILE=name_xgb_model, EXIST=there)
      if ( .not. there ) write(*,*) 'Machine learning model do not exist! filepath: ./', name_xgb_model

      !!! Load XGB model !!!
      ! Create  a dummy  XGDMatrix 
      allocate(xgb_input(ncol))
      xgb_input(:) = 0.0
      error = fortran_XGDMatrixCreateFromMat(xgb_input, nrow, ncol, missing, dmatrix)
      ! Create XGBooster object
      error = fortran_XGBoosterCreate(dmatrix, dmatrix_len, xgb_model)
      ! Load XGBooster model from binary file
      write(*,*) 'Reading '//trim(name_xgb_model)
      ! Always use "trim(name)//c_null_char" to pass the file name 
      error = fortran_XGBoosterLoadModel(xgb_model, trim(name_xgb_model)//c_null_char)
      ! Forces the XGBoost model to use only 1 thread for prediction (faster than using all of them)
      error = fortran_XGBoosterSetParam(xgb_model, trim('nthread')//c_null_char, trim('1')//c_null_char)

      !--- Cleanup
      deallocate(xgb_input)

    end subroutine xgboost_load_model

    !----------------------------------------------------------------------------------------
    !> @author Vinicius L S Silva
    !> @brief Predict using the loaded XGBoost model 
    !> xgboost_load_model() needs to be run first
    !----------------------------------------------------------------------------------------
    subroutine xgboost_predict(raw_input, out_result)

      implicit none

      real(c_float), dimension(17), intent( in )  :: raw_input
      real(c_float), pointer, intent( inout )     :: out_result(:)
      
      ! Local variables
      type(c_ptr)                   :: dmatrix         
      real(c_float), allocatable    :: xgb_input(:)  
      integer(c_int64_t)            :: out_len 
      type(c_ptr)                   :: out_result_c  
      integer                       :: error  
      real(c_float), dimension(17)  :: mean, std, norm_input 
      
      ! Constants
      integer(c_int), parameter  :: option_mask = 0
      integer(c_int), parameter  :: ntree_limit = 0
      real(c_float), parameter   :: missing = -999.0 
      

      ! Initialize variables
      mean = (/4.606498E+00, 4.134381E+01, 1.312609E+01, 1.869254E-02, 7.138944E+02, 3.206514E-06, 1.012754E+00, 3.133900E-03, 8.245164E-02, 9.482336E-02, 1.752161E-01, 2.600834E+00, -2.709370E+09, 2.655941E-06, 1.865938E-05, 8.406276E-01, 1.148470E+01 /)
      std = (/1.685133E+00, 2.110226E+03, 6.861566E+02, 2.188420E-02, 1.111934E+03, 9.776763E-06, 2.756072E+00, 1.341650E-02, 4.978913E-01, 3.173449E-01, 1.437997E+00, 9.703361E+00, 3.664636E+11, 6.481185E-05, 1.378213E-03, 9.822488E+01, 1.795359E+01 /)
      norm_input = (raw_input - mean)/std 

      !!! Make predictions !!!
      ! Create XGDMatrix using the input values 
      allocate(xgb_input(ncol))
      xgb_input = norm_input
      error = fortran_XGDMatrixCreateFromMat(xgb_input, nrow, ncol, missing, dmatrix)
      ! Make prediction. The result will be stored in c pointer out_result_c 
      error = fortran_XGBoosterPredict(xgb_model, dmatrix, option_mask, ntree_limit, out_len, out_result_c)
      ! Link to fortran pointer out_result 
      call c_f_pointer(out_result_c, out_result, [out_len])
      !write(*,*) 'XGB model Prediction: ',out_result
      
      ! Cleanup
      deallocate(xgb_input)
    
    end subroutine xgboost_predict  

    !----------------------------------------------------------------------------------------
    !> @author Vinicius L S Silva
    !> @brief Free the loaded XGBoost model from memory
    !----------------------------------------------------------------------------------------
    subroutine xgboost_free_model()

      implicit none

      integer :: error 
      
      ! Free XGB model
      error = fortran_XGBoosterFree(xgb_model)
    
    end subroutine xgboost_free_model

    !----------------------------------------------------------------------------------------
    !> @author Vinicius L S Silva
    !> @brief Teste the XGBoost coupling
    !----------------------------------------------------------------------------------------
    subroutine test_xgboost()

      implicit none
      
      real(c_float), dimension(17)  :: raw_input 
      real(c_float), pointer        :: out_result(:)

      !--- Initialize variables
      raw_input = (/5.000000E+00, 4.890481E+00, 3.937589E-01, 2.464229E-02, 2.388089E+02, 3.333715E-07, 9.755054E-01, 2.556481E-03, 6.391201E-02, 1.238833E-01, 3.023031E-01, 3.097082E+00, -1.180094E+06, 5.637608E-08, 1.044071E-07, 5.399639E-01, 3.000000E+00 /)

      call xgboost_load_model()

      call xgboost_predict(raw_input, out_result)
      write(*,*) 'XGB model prediction: ',out_result           

      nullify(out_result)
      call xgboost_free_model()

    end subroutine test_xgboost

end module multi_machine_learning