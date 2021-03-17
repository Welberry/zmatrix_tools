program zmat_rotate

  use mol2_class
  use zmatrix_class
  use file_functions
  use string_functions, only: real, int, operator(.isnumber.)
  use iso_varying_string
  use cmdline_arguments
  use fundamental_constants, only: radian
  use vector_class
  use quaternion_class

  implicit none

  type (mol2_object)      :: mol2
  type (zmatrix_object)   :: zmat
  type (quaternion)       :: q
  type (rotmatrix)        :: R

  real, allocatable, dimension(:,:) :: new_coords, rot_coords
  integer, dimension(:), pointer :: ordering

  integer :: i, j, zmat_entry, param_id, num_steps, second_atom(1)
  real    :: original_value(1), delta, inc, param(1), rotation_vector(3)
  
  type (varying_string) ::  fname
  character(len=6), allocatable, dimension(:) :: atomlabels

  if ( .not. have_args() ) then
     write(stderr,*) "Must provide a zmatrix filename as command line argument"
     stop
  end if

  ! Grab the file name from the command line
  fname = next_arg()

  zmat = char(fname)

  ! call print(zmat)
  
  if ( .not. have_args() ) then
     write(stderr,*) "Must provide a z-matrix entry as first command line argument"
     stop
  end if
  zmat_entry = int(next_arg())
  
  if ( .not. have_args() ) then
     write(stderr,*) "Must provide a rotation angle the parameter as second command line argument"
     stop
  end if
  delta = real(next_arg())

  if (param_id /= 1) delta = delta / radian

  if ( .not. have_args() ) then
     write(stderr,*) "Must provide a value for the number of steps as third command line argument"
     stop
  end if
  num_steps = int(next_arg())

  allocate(new_coords(3,num(zmat)))
  allocate(rot_coords(3,num(zmat)))
  allocate(atomlabels(num(zmat)))

  ! Generate cartesian coordinates from the z-matrix
  new_coords = as_xyz(zmat)

  second_atom = connectivity(zmat, (/ zmat_entry /), (/  1 /))

  rotation_vector = new_coords(:,zmat_entry) - new_coords(:, second_atom(1)) 

  ! print *,zmat_entry, second_atom, unit(rotation_vector)

  atomlabels = labels(zmat)

  do i = 1, size(atomlabels)
     ! No element has a label of more than two letters 
     atomlabels(i)(3:) = ''
     if (.isnumber. atomlabels(i)(2:2)) atomlabels(i)(2:2) = ''
  end do

  inc = abs(delta)/real(num_steps)

  write(0,*) original_value(1), delta, inc
  
  ! Change the paramater to parameter + delta ...
  do j = 0,num_steps 
     call write_frame(j)
  end do

  ! .. and now down to parameter - delta  ...
  do j = num_steps-1,-num_steps,-1
     call write_frame(j)
  end do

  ! .. and back to the original value
  do j = -num_steps,0
     call write_frame(j)
  end do

  contains
    
    subroutine write_frame (j)

      integer, intent(in) :: j

      param = parameter(zmat, (/ zmat_entry /), (/ param_id /), (/ original_value(1) + real(j)*inc /))
      q = as_quaternion(rotation_vector,real(j)*inc)

      R = q

      rot_coords = new_coords

      call rotate(R,rot_coords)

      print *,num(zmat)
      print *,'Frame number ',j,' angle ',real(j)*inc*radian

      do i = 1,size(rot_coords,2)
         print *,atomlabels(i),rot_coords(:,i)
      end do

    end subroutine write_frame

end program zmat_rotate
