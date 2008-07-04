program zmat_maker

  use simpose
  use fundamental_constants
  use cmdline_arguments
  use mol2_class
  use quaternion_class
  use rotmatrix_class
  use zmatrix_class
  use file_functions, only: stderr, stdout, log, delete
  use string_functions, only: join
  use iso_varying_string
  use variable_array

  implicit none

  !! $Log: zmat_maker.f90,v $
  !! Revision 1.3  2007/02/12 03:48:25  aidan
  !! Big code clean up. Took the fitting routine out into a separate function.
  !! Removed unused variables. Added command line arguments.
  !!
  !! Revision 1.2  2006/05/18 06:05:17  aidan
  !! Working version that can automatically output multiple z-matrices.
  !!
  !! Revision 1.1  2006/05/18 01:31:34  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: zmat_maker.f90,v 1.3 2007/02/12 03:48:25 aidan Exp $"

  type (mol2_object)    :: mol2
  type (zmatrix_object) :: zmat(100), zmat_template
  type (rotmatrix)      :: rmat

  type (quaternion)  :: q

  real    :: rsd, trans(3), rsd_limit, com(3), loc_coincident
  integer :: i, j, k, iz, numzmat, numloc, error, ind, loc, moltype, dontcare
  logical :: verbose, have_template, new_location

  character(len=1000)   :: fname, outname, buffer
  type (varying_string) :: myoptions(4)

  type location
     real :: coords(3)
     integer, pointer :: zmats(:), number(:)
  end type location

  type(location), allocatable :: locations(:)

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  myoptions(1) = 'help'
  myoptions(2) = 'rsd'
  myoptions(3) = 'quiet'
  myoptions(4) = 'locationtol'

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
     write(stderr,*) "Must provide a mol2 filename as command line argument"
     call usage
     stop
  end if

  verbose = .NOT. option_exists('quiet')

  ! See if we have a maximum rsd
  if (option_exists('rsd')) then
     ! Make sure we have a value
     if (.NOT. has_value('rsd')) then
        write(stderr,*) 'Option rsd must have a value!'
        call usage
        stop
     end if
     rsd_limit = get_value('rsd')
  else
     rsd_limit = 0.01
  end if

  ! See if we have a tolerance for deciding if molecules are on the same location
  if (option_exists('locationtol')) then
     ! Make sure we have a value
     if (.NOT. has_value('locationtol')) then
        write(stderr,*) 'Option locationtol must have a value!'
        call usage
        stop
     end if
     loc_coincident = get_value('locationtol')
  else
     loc_coincident = 4.
  end if

!!$  ! See if we have specified a template z-matrix file
!!$  if (option_exists('template')) then
!!$     ! Make sure we have a value
!!$     if (.NOT. has_value('template')) then
!!$        write(stderr,*) 'Must specify a filename with the option template!'
!!$        call usage
!!$        stop
!!$     end if
!!$     fname = get_value('template')
!!$     zmat_template = fname
!!$     have_template = .TRUE.
!!$  else
!!$     have_template = .FALSE.
!!$  end if

  ! Grab the mol2 file name from the command line
  fname = next_arg()

  ! Make a mol2 object from it
  mol2 = fname

  ! Use the input filename as a template for the output filename.
  ! Search for a ".mol2" suffix 
  ind = index(fname,".mol2")
  if (ind == 0) then
     ! Couldn't find a matching suffix, so look for the end of the
     ! filename string in the character variable
     ind = verify(fname," ",back=.TRUE.) + 1
  end if
  outname = fname
  ! Append new suffix
  outname(ind:) = ".qxyz"

  ! Delete an existing quaternion output file (ignoring an error 
  ! if one did not already exist)
  call delete(trim(outname),error)
     
  numzmat = 0
  numloc = 0

  allocate(locations(sub_num(mol2)))

  ! For each substructure we attempt to find a quaternion
  ! and associated centre of mass translations that will
  ! map the cartesian coordinates converted from a z-matrix 
  ! onto those in the other substructure positions in the mol2 file
  do i = 1, sub_num(mol2)

     rsd = rsd_limit + 1e6
     iz = 0

     ! Loop over various options of z-matrix until we find a match or
     ! make a new one!
     do while (rsd > rsd_limit)
        iz = iz + 1
        if (iz > sub_num(mol2)) then
           write(stdout,'(A,I0)') 'WOAH! Something nasty has happened .. have failed to fit substructure ',i
           write(stdout,'(A)') 'below the rsd limit. There is probably something wrong with the z-matrix'
           write(stdout,'(A)') 'this program will produce. Examine it for any coincident atoms, reorder '
           write(stdout,'(A)') 'the definition for that atom and re-run.'
           ! write(stdout,'(A)') 'the definition for that atom and use this z-matrix as a template and re-run.'
           exit
        end if
        if (iz > numzmat) then
           ! Make a new z-matrix, as none of the other ones are a good fit
           numzmat = numzmat + 1
           if (have_template) then
              ! zmat(numzmat) = as_zmatrix(zmat_template,mol2,i)
           else
              zmat(numzmat) = as_zmatrix(mol2,i)
           end if
        end if
        rsd = fit(zmat(iz), coords(mol2,i), order(mol2,i), q, trans)
     end do

     com = sum(coords(mol2,i),dim=2)/real(atom_num(mol2,i))
     new_location = .true.

     NUMLOOP: do j = 1, numloc
        if (sum(abs(com - locations(j)%coords)) < loc_coincident) then
           ! These molecules are in the same place.
           loc = j
           new_location = .false.
           do k = 1, size(locations(j)%zmats)
              if (locations(j)%zmats(k) == iz) then
                 locations(j)%number(k) = locations(j)%number(k) + 1
                 moltype = locations(j)%number(k) 
                 exit NUMLOOP
              end if
           end do
           ! If we got here then there were no zmatrices of this type
           ! at this location. So we put a new zmat on this location, 
           ! with population 1
           dontcare = push(locations(j)%zmats,iz)
           dontcare = push(locations(j)%number,1)
           moltype = 1
        end if
     end do NUMLOOP

     if (new_location) then
        numloc = numloc + 1
        if (push(locations(numloc)%zmats,iz) /= 1) stop 'Error initialising location!'
        moltype = push(locations(numloc)%number,1)
        locations(numloc)%coords = com
        loc = numloc
     end if

     write(buffer,'(4I6, 4F10.6,L,3F11.6)') i, loc, iz, moltype, as_array(q), improper(q), trans
     call log(outname, trim(buffer))

     if (verbose) then
        write(stdout,'("  Sub structure: ",I0)') i
        write(stdout,'("       Z-matrix: ",I0)') iz
        write(stdout,'("     Quaternion: ",4F11.6)') as_array(q)
        write(stdout,'("       Improper: ",L4)') improper(q)
        write(stdout,'("COM Translation: ",3F11.6)') trans
        write(stdout,'("            RSD: ", F11.6)') rsd
     end if
     
  end do

  if (verbose) print *,'Saved quaternion data to ',trim(outname)

  do i = 1, numzmat
     outname = fname
     if (numzmat > 1) then
        write(outname(ind:),'(A,I0,A)') "_",i,".zmat"
     else
        outname(ind:) = ".zmat"
     end if
     call print(trim(outname),zmat(i))
     if (verbose) print *,'Writing ',trim(outname)
  end do
     
contains

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
        
    ! Check to see we don't have a roto-inversion ... invert the z-matrix coordinates
    deviation_invert = superimpose(coords(:,order), -1*new_coords, q2, trans2)
        
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
    
  subroutine usage

    write(stderr,*)
    write(stderr,*) 'zmat_maker takes a mol2 file and fits a z-matrix to each'
    write(stderr,*) 'residue (i.e. molecule). The minimum number of z-matrices'
    write(stderr,*) 'are generated and a quaternion and centre of mass translation'
    write(stderr,*) 'is output for each residue. A new z-matrix is required when' 
    write(stderr,*) 'a fit exceeds the root square deviation (RSD) limit.'
    write(stderr,*)
    write(stderr,*) 'Usage: zmat_maker [--help] [--quiet] [--rsd=value] mol2file'
    ! write(stderr,*) 'Usage: zmat_maker [--help] [--quiet] [--rsd=value] [--template=zmatrix_file] mol2file'
    write(stderr,*)
    write(stderr,*) '  --help     - print this message'
    write(stderr,*) '  --quiet    - less verbose output'
    write(stderr,*) '  --rsd      - the value at which a new z-matrix is created'
    ! write(stderr,*) '  --template - specify a template z-matrix file'
    write(stderr,*)

  end subroutine usage

end program zmat_maker
