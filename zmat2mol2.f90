program zmat2mol2

  use simpose
  use fundamental_constants
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use mol2_class
  use quaternion_class
  use rotmatrix_class
  use zmatrix_class
  use file_functions
  use string_functions
  use iso_varying_string
  use variable_array

  implicit none

  type (mol2_object)      :: mol2
  type (zmatrix_object)   :: zmat
  type (rotmatrix)        :: rmat

  type (quaternion)  :: q
  real, dimension(3) :: trans, trans1, trans2

  real, pointer, dimension(:,:) :: mol2_coords, new_coords
  integer, dimension(:), pointer :: ordering

  integer :: i, j, error
  logical :: inversion

  real, pointer :: array(:)

  character(len=100) ::  fname
  character(len=6), allocatable, dimension(:) :: atomlabels

  type (varying_string) :: arraystring, myoptions(5)

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  myoptions(1) = 'help'
  myoptions(2) = 'q'
  myoptions(3) = 'quiet'
  myoptions(4) = 'com'
  myoptions(5) = 'improper'

  ! This call parses the command line arguments for command line options
  call get_options(myoptions, error)

  ! Check we weren't passed duff options -- spit the dummy if we were
  if (error > 0) then
     write(stderr,*) 'ERROR! Unknown options: ',join(bad_options()," ")
     call usage
     STOP
  end if

  ! Check if we just want to print the usage
  if (option_exists('help')) then
     call usage
     STOP
  end if

  if ( .not. have_args() ) then
     write(stderr,*) "Must provide a zmatrix filename as command line argument"
     call usage
     stop
  end if

  inversion = option_exists('improper')

  if (option_exists('q')) then
     ! Make sure we have a value
     if (.NOT. has_value('q')) then
        write(stderr,*) 'Option q must have 4 comma separated values (no spaces)!'
        call usage
        stop
     end if
     arraystring = get_value('q')
     i = push(array, real(split(arraystring,",")))
  else
     i = push(array, (/ 1., 0., 0., 0. /))
  end if

  call new(q, array, unit=.true., improper=inversion)

  i = splice(array,0)

  if (option_exists('com')) then
     ! Make sure we have a value
     if (.NOT. has_value('com')) then
        write(stderr,*) 'Option com must have 3 comma separated values (no spaces)!'
        call usage
        stop
     end if
     arraystring = get_value('com')
     i = push(array, real(split(arraystring,",")))
  else
     i = push(array, (/ 0., 0., 0. /))
  end if

  trans = array

  ! Grab the file name from the command line
  fname = next_arg()

  zmat = fname

  ! call print(zmat)
  
  allocate(new_coords(3,num(zmat)))
  allocate(atomlabels(num(zmat)))

  atomlabels = labels(zmat)

  ! Regenerate cartesian coordinates from the z-matrix
  new_coords = as_xyz(zmat)

  rmat = as_rotmatrix(q)

  call rotate(rmat,new_coords)
 
  print '(I0)',size(new_coords,2)
  print '(3A,2(F0.4,","),F0.4,A,3(F0.4,","),F0.4,A,L1)','Converted by zmat2xyz from ',&
       trim(fname),' com=',trans,' q=',as_array(q),' improper=',inversion
  do i = 1,size(new_coords,2)
    new_coords(:,i) = new_coords(:,i) + trans
    print '(A2,3F10.4)',label_to_atomtype(atomlabels(i)),new_coords(:,i)
  end do

contains

  function label_to_atomtype (label) result(atomtype)

    ! Interface variable
    character(len=*), intent(in) :: label

    ! Return value
    character(len=2) atomtype

    integer ::  i

    atomtype = ' '

    i = scan(label,'0123456789')

    if (i == 0) then 
       ! No numbers in label, assume both characters form a valid atom type
       i = 2
    else
       ! Want to make i the first place where numbers ARE NOT, so decrement
       i = i - 1
       ! Make sure we have a maximum of two characters
       i = min(2,i)
    end if

    atomtype(1:i) = label(1:i)

  end function label_to_atomtype

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'zmat2mol2 takes a zmatrix file and outputs a mol2 file'
    write(stderr,*)
    write(stderr,*) 'Usage: zmat2xyz [--q=q1,q2,q3,q4] [--improper] [--com=x,y,z] [--help]'
    write(stderr,*) '                [--qxyz=file] zmatrixfile'
    write(stderr,*)
    write(stderr,*) '  --q         - quaternion rotation (default q=1,0,0,0)'
    write(stderr,*) '  --improper  - rotation is improper'
    write(stderr,*) '  --com       - centre of mass translation (default com=0,0,0)'
    write(stderr,*) '  --qxyz=file - file containing quaternions and centres of mass'
    write(stderr,*) '                (as output from zmat_maker)'
    write(stderr,*) '  --help      - print this message'
    write(stderr,*)

  end subroutine usage

end program zmat2mol2
