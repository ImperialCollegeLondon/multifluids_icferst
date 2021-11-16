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

module multi_machine_learning

  use iso_c_binding
  use xgb_interface

  implicit none

  contains

  subroutine init_xgboost()

    implicit none
    integer                    :: error   
    type(c_ptr)                :: xgb_model
    type(c_ptr)                :: dmatrix         
    real(c_float), allocatable :: xgb_input(:)    
    type(c_ptr)                :: out_result_c
    integer(c_int64_t)         :: out_len 
    real(c_float), pointer     :: out_result(:)
    
    ! Constants
    integer(c_int)             :: option_mask = 0
    integer(c_int)             :: ntree_limit = 0
    real(c_float), parameter   :: missing = -999.0
    integer(c_int64_t)         :: dmatrix_len = 0
    integer(c_int64_t)         :: nrow = 1
    integer(c_int64_t)         :: ncol = 17 
    character(len=255)         :: name_xgb_model = 'XGBRFmodel_orig_lessfeatures_xgbformat.bin'
    
    ! Local variables
    logical :: there = .false.
    real(c_float), dimension(17) :: mean, std, raw_input, norm_input 

    !--- Initialize variables
    mean = (/4.606498E+00, 4.134381E+01, 1.312609E+01, 1.869254E-02, 7.138944E+02, 3.206514E-06, 1.012754E+00, 3.133900E-03, 8.245164E-02, 9.482336E-02, 1.752161E-01, 2.600834E+00, -2.709370E+09, 2.655941E-06, 1.865938E-05, 8.406276E-01, 1.148470E+01 /)
    std = (/1.685133E+00, 2.110226E+03, 6.861566E+02, 2.188420E-02, 1.111934E+03, 9.776763E-06, 2.756072E+00, 1.341650E-02, 4.978913E-01, 3.173449E-01, 1.437997E+00, 9.703361E+00, 3.664636E+11, 6.481185E-05, 1.378213E-03, 9.822488E+01, 1.795359E+01 /)
    raw_input = (/5.000000E+00, 4.890481E+00, 3.937589E-01, 2.464229E-02, 2.388089E+02, 3.333715E-07, 9.755054E-01, 2.556481E-03, 6.391201E-02, 1.238833E-01, 3.023031E-01, 3.097082E+00, -1.180094E+06, 5.637608E-08, 1.044071E-07, 5.399639E-01, 3.000000E+00 /)
    norm_input = (raw_input - mean)/std 

    !--- Starts here
    write(*,*) 'Starting XGBoost program'
    inquire( FILE=name_xgb_model, EXIST=there) 
    if ( .not. there ) write(*,*) 'Machine learning model do not exist!'
    
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

    
    !!! Make predictions !!!
    ! Create XGDMatrix using the input values 
    xgb_input = norm_input
    error = fortran_XGDMatrixCreateFromMat(xgb_input, nrow, ncol, missing, dmatrix)
    ! Make prediction. The result will be stored in c pointer out_result_c 
    error = fortran_XGBoosterPredict(xgb_model, dmatrix, option_mask, ntree_limit, out_len, out_result_c)
    ! Link to fortran pointer out_result 
    call c_f_pointer(out_result_c, out_result, [out_len])
    write(*,*) 'XGB model Prediction: ',out_result
    
    !--- Cleanup
    error = fortran_XGBoosterFree(xgb_model)
    deallocate(xgb_input)
    nullify(out_result)
  
  end subroutine init_xgboost
end module multi_machine_learning