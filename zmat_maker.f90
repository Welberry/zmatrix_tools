program zmat_maker

  use simpose
  use fundamental_constants
  use cmdline_arguments, only: have_args, next_arg
  use mol2_class
  use quaternion_class
  use rotmatrix_class
  use zmatrix_class
  use file_functions
  use iso_varying_string

  implicit none

  !! $Log: zmat_maker.f90,v $
  !! Revision 1.1  2006/05/18 01:31:34  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: zmat_maker.f90,v 1.1 2006/05/18 01:31:34 aidan Exp aidan $"

  type (mol2_object)             :: mol2
  type (zmatrix_object)          :: template
  type (zmatrix_object), pointer :: zmat
  type (zmatrix_object), target  :: zmat_library(100)
  type (zmatrix_object)          :: zmat2
  type (rotmatrix)               :: rmat

  type (quaternion)  :: q, q1, q2
  real, dimension(3) :: trans, trans1, trans2

  real, dimension(:,:), allocatable  :: mol2_coords, new_coords
  integer, dimension(:), allocatable :: ordering

  integer :: i, j, iz, numzmat
  real    :: rmsd, old_rmsd

  character(len=100) ::  fname

  if ( .not. have_args() ) then
     write(stderr,*) "Must provide a mol2 filename as command line argument"
     stop
  end if

  ! Grab the file name from the command line
  mol2 = char(next_arg())

  if (have_args()) then
     fname = ''
     fname = next_arg()
     write(stdout,*) "Using ",trim(fname)," as template"
     ! Make a z-matrix from these coordinates
     template = trim(fname)
     ! zmat = as_zmatrix(template,mol2)
  else
     ! Make a z-matrix from these coordinates
     write(*,'()',advance='no') ! Need this bogus write otherwise I get an address error!

     ! zmat = as_zmatrix(mol2)
     ! By default we grab the first substructure for our zmat 
     numzmat = 1
     zmat_library(numzmat) = mol2
  end if

  call print(zmat_library(numzmat))

  ! For each substructure we attempt to find a quaternion
  ! and associated centre of mass translations that will
  ! map the cartesian coordinates converted from a z-matrix 
  ! onto those in the other substructure positions in the mol2 file
  do i = 1, sub_num(mol2)

     if (allocated(new_coords)) deallocate(new_coords)
     if (allocated(mol2_coords)) deallocate(mol2_coords)
     if (allocated(ordering)) deallocate(ordering)

     allocate(new_coords(3,atom_num(mol2,i)))
     allocate(mol2_coords(3,atom_num(mol2,i)))
     allocate(ordering(atom_num(mol2,i)))

     rmsd = 10000.

     iz = 0

     ! Loop over various options of z-matrix until we find a match or
     ! make a new one!
     do while (rmsd > 0.01)

        iz = iz + 1

        print *,'Trying z-matrix ',iz

        if (iz > numzmat) then
           write (stderr,*) "ERROR! We should never be here .. that is consoling n'cest pas?"
        end if

        zmat => zmat_library(iz)

        ! Regenerate cartesian coordinates from the z-matrix
        new_coords = as_xyz(zmat)

        ! Grab the coordinates for substructure i
        mol2_coords = coords(mol2,i)

        ! And the ordering array for substructure i
        ordering = order(mol2,i)

        ! Find quaternion and translation relating z-matrix to cartesian
        ! coordinates
        rmsd = superimpose(mol2_coords(:,ordering), new_coords, q1, trans1)
        
        old_rmsd = rmsd
        
        ! Check to see we don't have a roto-inversion ...
        
        ! Invert the z-matrix coordinates
        rmsd = superimpose(mol2_coords(:,ordering), -1*new_coords, q2, trans2)
        
        if (old_rmsd > rmsd) then
           if (rmsd < 0.01) then
              ! Accept this -- set the improper flag on the quaternion
              ! (we use a "faux" if statement to satisfy the need for a 
              ! return value)
              q = q2
              trans = trans2
              if(improper(q,.true.)) continue
           else
              ! write(stderr,*) 'Something bad has happened -- cannot fit z-matrix to coordinates'
              ! write(stderr,*) '           RMSD: ',old_rmsd
              ! write(stderr,*) 'RMSD (inverted): ',rmsd
              ! stop
              if (iz == numzmat) then
                 numzmat = numzmat + 1
                 zmat_library(numzmat) = as_zmatrix(mol2,i)
                 call print(zmat_library(numzmat))
                 ! print *,'New z-matrix: ',numzmat
              end if
              cycle
           end if
        else
           q = q1
           trans = trans1
           rmsd = old_rmsd
        end if
        
        exit
     
     end do

     ! Make a rotation matrix out of this quaternion
     rmat = as_rotmatrix(q)
        
     ! Rotate the coordinates using this rotation matrix
     call rotate(rmat,new_coords)
        
     ! Add the translation vector (we use the spread function to
     ! create a temporary matrix of the translation vector 
     ! replicated the same number of times as the number of coordinates)
     new_coords = new_coords + spread(trans,2,size(new_coords,2))
     
     ! if (sum(abs(new_coords - mol2_coords)) > 1e-3) then
     !    print *,sum(abs(new_coords - mol2_coords))
     !    stop
     ! end if
     
     write(stdout,'("  Sub structure: ",I0)') i
     write(stdout,'("       Z-matrix: ",I0)') iz
     write(stdout,'("     Quaternion: ",4F11.6)') as_array(q)
     write(stdout,'("       Improper: ",L4)') improper(q)
     write(stdout,'("COM Translation: ",3F11.6)') trans
     write(stdout,'("           RMSD: ", F11.6)') rmsd
     
  end do
     
end program zmat_maker
