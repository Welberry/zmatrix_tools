program zmat_anim

  use mol2_class
  use zmatrix_class
  use file_functions
  use string_functions, only: real, int, operator(.isnumber.), join
  use iso_varying_string
  use cmdline_arguments
  use fundamental_constants, only: radian

  implicit none

  !! $Log: zmat_anim.f90,v $

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: zmat_anim.f90,v 1.1 2006/05/10 05:06:50 aidan Exp aidan $"

  type (mol2_object)      :: mol2
  type (zmatrix_object)   :: zmat

  real, dimension(3) :: trans, trans1, trans2

  real, allocatable, dimension(:,:) :: new_coords
  integer, dimension(:), pointer :: ordering

  integer :: i, j, zmat_entry, param_id, num_steps, error
  real    :: original_value(1), delta, inc, param(1)
  
  type (varying_string) ::  fname, myoptions(4)
  character(len=6), allocatable, dimension(:) :: atomlabels

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'delta'
  myoptions(3) = 'param'
  myoptions(4) = 'nstep'

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

  ! See if we have specified a number of steps in our animation
  if (option_exists('nstep')) then
     ! Make sure we have a value
     if (.NOT. has_value('nstep')) then
        write(stderr,*) 'Option nstep must have a value!'
        call usage
        stop
     end if
     ! Get the number of steps
     num_steps = get_value('nstep')
  else
     num_steps = 20
  end if

  ! Must have at least two steps ...
  num_steps = max(num_steps,2)

  ! See if we have specified the parameter we will change in our animation
  if (option_exists('param')) then
     ! Make sure we have a value
     if (.NOT. has_value('param')) then
        write(stderr,*) 'Option param must have a value!'
        call usage
        stop
     end if
     ! Get the number of steps
     param_id = get_value('param')
  else
     param_id = 3
  end if

  ! See if we have specified the amount we will change our parameter
  if (option_exists('delta')) then
     ! Make sure we have a value
     if (.NOT. has_value('delta')) then
        write(stderr,*) 'Option delta must have a value!'
        call usage
        stop
     end if
     ! Get the number of steps
     delta = get_value('delta')
  else
     delta = 0.1
  end if

  ! if (param_id /= 1) delta = delta / radian

  ! Grab the file name from the command line
  fname = next_arg()

  zmat = char(fname)

  ! call print(zmat)
  
  if ( .not. have_args() ) then
     write(stderr,*) "ERROR! Must provide a zmatrix entry (line number) as a command line argument"
     call usage()
     stop
  end if
  zmat_entry = int(next_arg())
  
  allocate(new_coords(3,num(zmat)))
  allocate(atomlabels(num(zmat)))

  atomlabels = labels(zmat)

  do i = 1, size(atomlabels)
     ! No element has a label of more than two letters 
     atomlabels(i)(3:) = ''
     if (.isnumber. atomlabels(i)(2:2)) atomlabels(i)(2:2) = ''
  end do

  original_value = parameter(zmat, (/ zmat_entry /), (/ param_id /))
  
  if (param_id == 1) then
     ! If we are dealing with a bond length then delta is a fraction of that length
     delta = delta * original_value(1)
  else
     ! If we are dealing with an angle then delta is a fraction of 180 degrees
     delta = (delta * 180. / radian)
  end if

  inc = abs(delta)/real(num_steps)

!!$  ! Change the paramater to parameter + delta ...
!!$  do j = 0,num_steps 
!!$     call write_frame(j)
!!$  end do
!!$
!!$  ! .. and now down to parameter - delta  ...
!!$  do j = num_steps-1,-num_steps,-1
!!$     call write_frame(j)
!!$  end do
!!$
!!$  ! .. and back to the original value
!!$  do j = -num_steps,0
!!$     call write_frame(j)
!!$  end do

  do j = -num_steps/2,num_steps/2
     call write_frame(j)
  end do

contains
    
  subroutine write_frame (j)
    
    integer, intent(in) :: j
      
    param = parameter(zmat, (/ zmat_entry /), (/ param_id /), (/ original_value(1) + real(j)*inc /))
    
    ! Regenerate cartesian coordinates from the z-matrix
    new_coords = as_xyz(zmat)
    
    print *,num(zmat)
    print '(A,I0,A,F)','Frame number ',j,', parameter: ',param(1)
    
    do i = 1,size(new_coords,2)
       new_coords(:,i) = new_coords(:,i) + trans
       print '(A6,3F12.6)',atomlabels(i),new_coords(:,i)
    end do
    
  end subroutine write_frame
  
  subroutine usage
    
    write(stderr,*)
    write(stderr,*) 'zmat_anim creats a series of frames by changing a specified parameter'
    write(stderr,*) 'in a z-matrix input file, to check that the z-matrix is properly constructed'
    write(stderr,*) 'The result is series of "xyz" structures printed to standard output '
    write(stderr,*)
    write(stderr,*) 'Usage: zmat_anim [OPTIONS] <zmatrix file> <atom_number>'
    write(stderr,*)
    write(stderr,*) '  where atom_number is the entry in the zmatrix file that will change'
    write(stderr,*)
    write(stderr,*) '  Options:'
    write(stderr,*)
    write(stderr,*) '  --help          - print this message'
    write(stderr,*) '  --param=<1,2,3> - specify the parameter: 1=bond, 2=bond angle,'
    write(stderr,*) '                    3=torsion angle (defaults to 3)'
    write(stderr,*) '  --delta=<value> - the fraction the parameter should change (default=0.1 (of'
    write(stderr,*) '                    bond length for param=1, otherwise a fraction of 180 degrees))'
    write(stderr,*) '  --nstep=<value> - sets the number of frames (steps) in the animation (default=20)'
    write(stderr,*)

  end subroutine usage

end program zmat_anim
