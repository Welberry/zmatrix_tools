program zmatchk

  use simpose
  use fundamental_constants
  use cmdline_arguments, only: have_args, next_arg
  use mol2_class
  use quaternion_class
  use rotmatrix_class
  use zmatrix_class
  use file_functions
  use string_functions, only: real
  use iso_varying_string

  implicit none

  type (mol2_object)      :: mol2
  type (zmatrix_object)   :: zmat
  type (rotmatrix)        :: rmat

  type (quaternion)  :: q, q1, q2
  real, dimension(3) :: trans, trans1, trans2

  real, pointer, dimension(:,:) :: mol2_coords, new_coords
  integer, dimension(:), pointer :: ordering

  integer :: i, j
  real    :: rmsd, old_rmsd, array(4)

  character(len=100) ::  fname
  character(len=6), allocatable, dimension(:) :: atomlabels

  if ( .not. have_args() ) then
     write(stderr,*) "Must provide a zmatrix filename as command line argument"
     stop
  end if

  ! Grab the file name from the command line
  fname = next_arg()

  zmat = fname

  ! call print(zmat)
  
  do i = 1,4
     if ( .not. have_args() ) then
        write(stderr,*) "Must provide a quaternion as command line argument"
        stop
     end if
     array(i) = real(next_arg())
  end do

  q = array

  do i = 1,3
     if ( .not. have_args() ) then
        write(stderr,*) "Must provide an xyz as command line argument"
        stop
     end if
     trans(i) = real(next_arg())
  end do

  allocate(new_coords(3,num(zmat)))
  allocate(atomlabels(num(zmat)))

  atomlabels = labels(zmat)

  ! Regenerate cartesian coordinates from the z-matrix
  new_coords = as_xyz(zmat)

  rmat = as_rotmatrix(q)

  call rotate(rmat,new_coords)
 
  print '(I0)',size(new_coords,2)
  print '(A)','output of zmat_chk'
  do i = 1,size(new_coords,2)
    new_coords(:,i) = new_coords(:,i) + trans
    print '(A,3F9.3)',atomlabels(i)(1:1),new_coords(:,i)
  end do

end program zmatchk
