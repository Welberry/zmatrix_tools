program mol2fixer

  use fundamental_constants
  use cmdline_arguments
  use mol2_class
  use file_functions
  use string_functions, only: real
  use crystallography_class

  implicit none

  type (mol2_object)      :: mol2

  type(crystallography_object) :: xtal

  real, allocatable, dimension(:,:) :: mol2_coords, new_coords

  integer :: i, j

  character(len=100) ::  fname
  character(len=6), allocatable, dimension(:) :: atomlabels

  if ( .not. have_args() ) then
     write(stderr,*) "Must provide a mol2 filename as command line argument"
     stop
  end if

  ! Grab the file name from the command line
  fname = next_arg()

  mol2 = fname

  xtal = mol2

  allocate(mol2_coords(3,atom_num(mol2,1)), new_coords(3,atom_num(mol2,1)))

  mol2_coords = coords(mol2,1) 

  new_coords = as_fractional(xtal,mol2_coords)

  do i = 1,size(new_coords,2)
    write(*,'(6F8.4)') mol2_coords(:,i),new_coords(:,i)
  end do

end program mol2fixer

