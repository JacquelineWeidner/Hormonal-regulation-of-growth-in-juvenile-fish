!*******************************************************************************
! PURPOSE: Write arrays/matrices to a binary file
! include use BINARY_WRITE before using the module!!
! $Id: binary_write.f90 2996 2017-03-10 10:03:21Z cje012 $
!*******************************************************************************

! NOTE: form='binary' is not portable and should **never be used** in any
!       new codes. The outdated 'binary' functionality is emulated in
!       F2003 standard combination of 'unformatted' and access='stream';
!       status='replace' makes sure the file is opened strictly for writing
!       and not for appending (default for stream is system dependent).

module binary_write

  use parameter_values, only : PREC_REAL

  implicit none

  ! Global constants.
  integer, parameter, private :: unit_number = 110
  character(len=*), parameter, private :: file_form = 'unformatted'
  ! Folder path from file
  character(100), parameter, private :: foldername = 'D:\' 

  ! Generic interface, can use matrix_to_binary in the calling program.
  interface matrix_to_binary
    module procedure matrix_to_binary_real8_3d
    module procedure matrix_to_binary_int_4d
  end interface matrix_to_binary

  contains

  ! Write to binary for 3D array with real values
  subroutine matrix_to_binary_real8_3d(matrix, file_name, funit)
    integer, parameter :: matrix_rank = 3                                          
    real(kind=PREC_REAL), dimension(:,:,:), intent(in) :: matrix                   
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    ! Local variables
    integer, dimension(matrix_rank) :: dim_length
    integer :: i, funit_loc
    
    ! Check if optional unit name is provided, if yes, use it. If not use default unit_number
    if (present(funit)) then;funit_loc=funit;else;funit_loc=unit_number;end if
    ! Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(foldername) // file_name, form=file_form, access='stream',    &
          status='replace')
    dim_length(1:matrix_rank) =  shape(matrix)
    write(funit_loc) matrix_rank, (dim_length(i), i = 1, matrix_rank), matrix
    close(funit_loc)

  end subroutine matrix_to_binary_real8_3d


  ! Write to binary for 4D array with integer values
  subroutine matrix_to_binary_int_4d(matrix, file_name, funit)
    integer, parameter :: matrix_rank = 4                                         
    integer, dimension(:,:,:,:), intent(in) :: matrix                             
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    ! Local variables
    integer, dimension(matrix_rank) :: dim_length
    integer :: i, funit_loc
    
    ! Check if optional unit name is provided, if yes, use it. If not use default unit_number.
    if (present(funit)) then;funit_loc=funit;else;funit_loc=unit_number;end if
    ! Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(foldername) // file_name, form=file_form, access='stream',    &
          status='replace')
    dim_length(1:matrix_rank) =  shape(matrix)
    write(funit_loc) matrix_rank, (dim_length(i), i = 1, matrix_rank), matrix
    close(funit_loc)

  end subroutine matrix_to_binary_int_4d

end module binary_write

