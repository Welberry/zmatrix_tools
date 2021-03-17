program zmat_anim

  use simpose
  use zmatrix_class
  use quaternion_class
  use rotmatrix_class
  use file_functions
  use string_functions, only: real, int, operator(.isnumber.), join
  use iso_varying_string
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use fundamental_constants, only: radian

  implicit none

  !! $Log: zmat_anim.f90,v $
  !! Revision 1.1  2006/05/10 05:07:28  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: zmat_anim.f90,v 1.1 2006/05/10 05:07:28 aidan Exp aidan $"

  type (zmatrix_object)   :: zmat

  real, dimension(3) :: trans, trans1, trans2

  real, allocatable, dimension(:,:) :: new_coords
  integer, dimension(:), pointer :: ordering

  integer :: i, j, zmat_entry, param_id, num_steps, error, num_atoms
  real    :: original_value(1), delta, inc, param(1)
  
  type (varying_string) ::  fname, myoptions(5)
  character(len=6), allocatable, dimension(:) :: atomlabels

  logical :: palindrome

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'delta'
  myoptions(3) = 'param'
  myoptions(4) = 'nstep'
  myoptions(5) = 'palindrome'

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

  palindrome = option_exists('palindrome')

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
  
  num_atoms = num(zmat)
  allocate(new_coords(3,num_atoms))
  allocate(atomlabels(num_atoms))

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

  if (palindrome) then
     write(stderr,*) 'Yes palindrome'
     do j = num_steps/2,-num_steps/2,-1
        call write_frame(j)
     end do
  end if

contains
    
  subroutine write_frame (j)
    
    integer, intent(in) :: j
    
    type (quaternion) :: q
    real :: trans(3), rsd

    type (zmatrix_object)   :: zmat_tmp
    type (rotmatrix)        :: rmat

    ! zmat_tmp = zmat
    call copy(zmat_tmp, zmat)

    ! print *,parameter(zmat, (/ zmat_entry /), (/ param_id /) )
    param = parameter(zmat_tmp, (/ zmat_entry /), (/ param_id /), (/ original_value(1) + real(j)*inc /))
    ! print *,parameter(zmat, (/ zmat_entry /), (/ param_id /) )
    ! print *,parameter(zmat_tmp, (/ zmat_entry /), (/ param_id /) )

    if (zmat == zmat_tmp) print*,'zmatrices are equal!'
    
    ! Regenerate cartesian coordinates from the z-matrix
    new_coords = as_xyz(zmat_tmp)

    ! new_coords(1,:) = new_coords(1,:) + 1.3
    ! print *,'coods: ',new_coords(:,zmat_entry)
    
    rsd = fit(zmat, new_coords, (/ (i,i=1,num_atoms) /), q, trans)
        
    print *,num_atoms
    write(*,'(A,I0,A,F6.3,A,F4.3,A,3F6.3,A,4F6.3,A,L)') 'Frame number ',j,', parameter: ',param(1),' RSD = ',rsd,' trans = ',trans,' q = ',as_array(q),' improper = ',improper(q)

    rmat = q

    ! call rotate(rmat,new_coords)
    
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
    write(stderr,*) '  --palindrome    - make the loop go backwards at the end'
    write(stderr,*) '                    which results in smoother animation'
    write(stderr,*)

  end subroutine usage

  real function fit(zmat, coords, order, q, trans)

    ! Calculate the best fit quaternion and translation required
    ! to map a z-matrix on to a set of real space coordinates (NB: the 
    ! coordinates dimensions should be (3,n), n=number of atoms). 
    ! Returns the root mean square deviation between the real space 
    ! coordinates and the fitted z-matrix

    ! Interface variables
    type (zmatrix_object), intent(in) :: zmat
    real, intent(in)                  :: coords(:,:)
    integer, intent(in)               :: order(:)
    type (quaternion), intent(out)    :: q
    real, dimension(3), intent(out)   :: trans

    ! Local variables
    real             :: new_coords(size(coords,1),size(coords,2)) 
    real             :: deviation, deviation_invert
    type(quaternion) :: q1, q2
    real             :: trans1(3), trans2(3)

    ! Regenerate cartesian coordinates from the z-matrix
    new_coords = as_xyz(zmat)

    ! Find quaternion and translation relating z-matrix to cartesian coordinates
    deviation = superimpose(coords(:,order), new_coords, q1, trans1)
    ! print *,as_array(q1),deviation
        
    ! Check to see we don't have a roto-inversion ... invert the z-matrix coordinates
    deviation_invert = superimpose(coords(:,order), -1*new_coords, q2, trans2)
    ! print *,as_array(q2),deviation_invert
        
    if (deviation > deviation_invert) then
       ! Accept the roto inverted solution
       q = q2
       trans = trans2
       fit = deviation_invert 
       ! Set the improper flag on the quaternion (we use a "faux" if 
       ! statement to satisfy the need for a return value)
       if(improper(q,.true.)) continue
    else
       q = q1
       trans = trans1
       fit = deviation
    end if

    ! Make it a mean square deviation
    ! fit = fit / real(size(new_coords))
    
  end function fit
    
end program zmat_anim
