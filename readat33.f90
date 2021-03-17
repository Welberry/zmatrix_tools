!---------------------------------------------------------------------------------------------!
!August 12, 2003: Darren Goossens                                                             !
! readat for 33dimethoxybenzil cdeveloped from MC program as at 1:00pm                        !
!---------------------------------------------------------------------------------------------!
!August 12, 2003: Darren Goossens                                                             !
! 33dimethoxybenzil MC program as at 1:00 pm.  Needs to be tested for run-time errors         !
!---------------------------------------------------------------------------------------------!
!August 11, 2003: Darren Goossens                                                             !
! 33dimethoxybenzil MC program as at 5:00 pm.  Unfinished                                     !
!---------------------------------------------------------------------------------------------!
!August 7, 2003: Darren Goossens                                                              !
!  Start on turning the contact plotting program into an MC program.                          !
! In the getcarts subroutine,                                                                 !
!          NOTE that shiftvec is commented out, UNLIKE in the contact vector calculation      !
!---------------------------------------------------------------------------------------------!
!August 5, 2003: Darren Goossens                                                              !
!Start on contact plotting program using contacts33benz.f90 as a basis.  Also use TRW's       !
!ps_routines, since I have some familiarity. And some code from my own plotcarts.f            !
!  --------------------------------------------------------------------------------------     !
!August 4, 2003: Darren Goossens                                                              !
!Write contacts33benz.f90 to contacts22benz_bak.f90 and start modifying contacts33benz.f90    !
!to incorporate some sorting of contact vectors.  Sort by which mol. they start from, which   !
!they go to and which rotating fragment of the molecule they go from/to.                      !
!  ----------------------------------------------------------------------------------------   !
!July 30, 2003.  Darren Goossens                                                              !
!This is a butchered version of Aidan's z-matrix-handling test program benzil.test.           !
!Whereas that program was to read in a mol2 and generate the quaternions, origin translations !
!and z-matrix for 33dimethoxy benzil, this program is to read in the z-matrix and various     !
!information and replicate the molecules through space and work out a plausuble set of contact!
!vectors.  It will also therefore form the basis of an MC program.                            !
!---------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!          Module declared by Darren - change all fixed size arrays here       !
!------------------------------------------------------------------------------!

module myparams
  !Integer Parameters!
  integer, parameter:: zm_max = 1      !max number of different z-mats
  integer, parameter:: mol_max = 4     !max number of molecules of each z-mat
  integer, parameter:: amax = 32       !max number of unit cells along a dirn
  integer, parameter:: bmax = 32       !max number of unit cells along b dirn
  integer, parameter:: cmax = 32       !max number of unit cells along c dirn
  integer, parameter:: atom_max = 21   !max number of atoms in a z-matrix
  integer, parameter:: con_max = 200   !max number of contacts per molecule
  integer, parameter:: ty_max = 70     !max number of types of contacts
  integer, parameter:: in_max = 7      !max number of internal variables
  integer, parameter:: cell_atoms_max = 90      !max number of atoms for whole unit cell
end module myparams

!------------------------------------------------------------------------------!
!          Program proper starts here                                          !
!------------------------------------------------------------------------------!

subroutine readat(LLL,isite,cat,inum,xxx,yyy,zzz)

  !------------------------------------------------------------------------------!
  !          Use these modules -- requires compile flag -lmodules                !
  !------------------------------------------------------------------------------!

  use cmdline_arguments, only: have_args, next_arg
  Use file_functions
  use quaternion_class
  use rotmatrix_class
  use zmatrix_class
  use crystallography_class
  use myparams

  implicit none

  !------------------------------------------------------------------------------!
  !          These variables are for dovetailing with diffuse.f                  !
  !------------------------------------------------------------------------------!

  real, dimension(amax,bmax,cmax,cell_atoms_max,3):: xyz
  real, intent(out)                               :: xxx,yyy,zzz
  CHARACTER(len=2)                                :: ty1, ty2
  CHARACTER(len=2),dimension(cell_atoms_max)      :: occc
  INTEGER,intent(in)                              :: isite
  INTEGER,intent(out)                             :: inum
  INTEGER                                         :: itest,ia,ib,ic
  INTEGER,dimension(3),intent(in)                 :: LLL
  CHARACTER(len=4),intent(in)                     :: cat

  !------------------------------------------------------------------------------!
  !                      End of declarations                                     !
  !------------------------------------------------------------------------------!

  !------------------------------------------------------------------------------!
  !                      Need to save a couple of arrays which we set up on      !
  !       first pass through the readat.f
  !------------------------------------------------------------------------------!

  save occc, xyz,itest

  !------------------------------------------------------------------------------!
  !      If positions already established, skip over this part                   !
  !------------------------------------------------------------------------------!

  IF (itest.eq.1234) goto 657

  call getcoords(xyz,occc)

657 itest = 1234

  ia=LLL(1)
  ib=LLL(2)
  ic=LLL(3)
  xxx = xyz(ia,ib,ic,isite,1)
  yyy = xyz(ia,ib,ic,isite,2)
  zzz = xyz(ia,ib,ic,isite,3)
  inum = 0
  ty1 = occc(isite)
  ty2 = cat(1:2)
  if(ty1.eq.ty2) inum = 1
  return

end subroutine readat

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
subroutine getcoords(xyz,occc)

  !------------------------------------------------------------------------------!
  !          Use these modules -- requires compile flag -lmodules                !
  !------------------------------------------------------------------------------!

  use cmdline_arguments, only: have_args, next_arg
  Use file_functions
  use quaternion_class
  use rotmatrix_class
  use zmatrix_class
  use crystallography_class
  use myparams

  implicit none

  !------------------------------------------------------------------------------!
  !          Variables declared by Darren                                        !
  !------------------------------------------------------------------------------!

  !Characters!
  character(len=100)                   :: outname, fname,oname,vname, ename, header
  ! 	fname is master file, vname is control 
  !       variables file, and oname is output filename

  !Integers!
  integer                                  :: dunit, mcell, i1,i2
  integer                                  :: typ, ntyp,ncy,mtyp, icy
  integer                                  :: num_zm ,k,m,iz
  integer, dimension(zm_max)               :: num_mol,numat,nintern 
  integer, dimension(3)                    :: simsize
  integer                                  :: iat
  integer,dimension(mol_max,zm_max)        :: conmol
  integer,dimension(zm_max,mol_max,con_max)::da,db,dc,dz,dm,ty,oat,dat
  integer,dimension(in_max,2,zm_max)       :: intern
  integer,dimension(-amax:amax+amax,3)     :: wrap
  integer                                  :: oa, ob, oc, oz, om, nin, irej, iacc, iup

  !Reals!
  real, dimension(zm_max)                               :: rnum_mol
  real, dimension(6)                                    :: cellpar
  real                                                  ::v_min,v_max,temp
  real                                                  :: edif,new_e,old_e
  real, dimension(mol_max,3)                            :: translate
  real, dimension(amax,bmax,cmax,zm_max,mol_max,3)      :: xyzcrystal
  real, dimension(amax,bmax,cmax,zm_max,mol_max,in_max) :: ncrystal
  real, dimension(in_max)                               :: tempin
  real, dimension(3)                                    :: rsimsize, tempxyz
  real, dimension(3,3)                                  :: celltrans
  real,dimension(3,amax,bmax,cmax,mol_max,atom_max,zm_max):: carts
  real,dimension(3,atom_max)                            :: tempcarts
  real,dimension(zm_max,mol_max,con_max)                ::length
  real                                                  :: rvar,rnum_zm
  real, dimension(ty_max)                               :: sprcon
  real,dimension(7)                                     :: widths
  real,dimension(7+in_max)                              :: tweak,initw
  real,dimension(in_max,zm_max)                         :: inwidths
  real,  dimension(amax,bmax,cmax,zm_max,mol_max,in_max):: incrystal

  !Custom Types!
  type (quaternion), dimension(mol_max)                          :: quatern
  type (quaternion), dimension(amax,bmax,cmax,zm_max,mol_max)    :: qcrystal
  type (crystallography_object)                                  :: xtal
  type (quaternion)                                              :: tempq

  !Logicals!

  !------------------------------------------------------------------------------!
  !          These variables are for dovetailing with diffuse.f                  !
  !------------------------------------------------------------------------------!

  real, dimension(amax,bmax,cmax,cell_atoms_max,3),intent(out):: xyz
  real xxx,yyy,zzz
  CHARACTER(len=2),dimension(cell_atoms_max),intent(out)      :: occc
  CHARACTER(len=2)                                :: ty1, ty2
  INTEGER                                         :: itest, isite,inum
  INTEGER,dimension(3)                            :: LLL
  CHARACTER(len=4)                                :: cat

  !------------------------------------------------------------------------------!
  !          Variables declared by Aidan in the original program                 !
  !------------------------------------------------------------------------------!

  integer                                 :: i

  type (rotmatrix)                        :: rmat
  type (zmatrix_object),dimension(zm_max) :: zmat
  type (zmatrix_object)                   :: z1

  real, dimension(3,atom_max) :: new_coords !new_coords is a 3xN array, effectively.

  !------------------------------------------------------------------------------!
  !                      End of declarations                                     !
  !------------------------------------------------------------------------------!

  !------------------------------------------------------------------------------!
  ! To convert the MC to a readat, we remove most of it and read in the molecule !
  ! variables instead of replicating them across the simulation.                 !
  ! Can use the same input files                                                 !
  !------------------------------------------------------------------------------!

  if ( .not. have_args() ) then
     write(stdout,*)'---------------------------------------------------------------------'
     write(stderr,*) "Must provide a master file filename as command line argument"
     write(stderr,*) 'Program exiting.'
     write(stdout,*)'---------------------------------------------------------------------'
     stop
  end if   ! arg test
  fname = next_arg()  

  !----------------------------------------------------------------------------------------!
  !Now we read in the masterfile (fname), which at present contains                        !
  !name of the file from which we will read the contact vectors, and the name of the       ! 
  !control file which hold the size of the simulation, the unit cell, and also gives       !
  !some parameters for contact vector ploitting and the zmatrices and CofM/quaternion sets.!
  !----------------------------------------------------------------------------------------!
  ! Note; masterfile is probably benzil33MC.master,                                        !
  !which is contact33.master with with an extra filename in it for the MC controls.        !
  ! But for readat we'll only use the controlfile and the name of the outputted variable   !
  ! file, which we can get from the MC file, but which is the ONLY thing we need to get    !
  ! from that file.                                                                        !
  !----------------------------------------------------------------------------------------!

  write(stdout,*)'--------------------------------------------------------------------------'
  dunit = open(fname, status='old') 
  read(dunit,*)vname
  write(stdout,'(" Crystal variables will be read from           : ",a34)')vname
  read(dunit,*)oname
  write(stdout,'(" Contact vectors iare in but will not be read from: ",a34)')oname
  read(dunit,*)ename
  write(stdout,'(" MC params will be read from                   : ",a34)')ename
  close(dunit)
  write(stdout,'(" Master file closed                            : ",a34)')fname
  write(stdout,*)'--------------------------------------------------------------------------'

  !-------------------------------------------------------------------------------------------!
  !                      Now we read in the infomation                                        !
  !-------------------------------------------------------------------------------------------!

  call readfiles(vname,num_zm,num_mol,simsize,cellpar,v_min,v_max,translate,quatern,zmat)

  !-------------------------------------------------------------------------------------------!
  !                           Now we read in the MC infomation                                !
  ! this is MCreadat, not MC read since we only need it for the filename, nintern and intern  !
  !-------------------------------------------------------------------------------------------!

  call MCreadat(ename,outname, intern, nintern,num_zm)

  !-------------------------------------------------------------------------------------------!
  !        Instead of replicator routines, we read in the outname file and fill the           !
  ! qcrystal, xyzcrystal, incrystal arrays.                                                   !
  !-------------------------------------------------------------------------------------------!

  call readcrystal(outname, simsize, num_zm, num_mol, qcrystal, xyzcrystal, incrystal, nintern)

  xtal = cellpar(:)

  !-------------------------------------------------------------------------------------------!
  !        Now we have to get all the cartesian positions.                                    !
  !-------------------------------------------------------------------------------------------!

  call getcarts(simsize,xtal,num_zm,num_mol, zmat, qcrystal, xyzcrystal, carts, &
       numat,incrystal,nintern,intern)

  !-------------------------------------------------------------------------------------------!
  !        Now we need to turn the cartesians into fractional and give each atom a number     !
  ! and type across the whole unit cell, since this is what diffuse.f wants                   !
  !-------------------------------------------------------------------------------------------!

  call cart2frac(carts, simsize,num_zm,num_mol,xtal,zmat,occc,numat,xyz)

  write(stdout,*)'--------------------------------------------------------------------------'
  write(stdout,*)'Initial cartesian positions calculated.'
  write(stdout,*)'----------------------------------------------------------------------------'
  write(stdout,*)'Intialisation completed.'
  write(stdout,*)'----------------------------------------------------------------------------'
  write(stdout,*)'Getcoords exiting.'
  write(stdout,*)'----------------------------------------------------------------------------'

  !-------------------------------------------------------------------------------------------!
  !   End of position calculation.                                                            !
  !-------------------------------------------------------------------------------------------!

  return

end subroutine getcoords

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
subroutine readfiles(vname,num_zm,num_mol,simsize,cellpar,v_min,v_max,translate,quatern,zmat)

  use file_functions
  use zmatrix_class
  use myparams

  implicit none

  !------------------------------------------------------------------------------!
  !          Variables declared by Darren                                        !
  !------------------------------------------------------------------------------!

  !Characters!
  character(len=100),dimension(zm_max)    ::  zname,qname
  character(len=100),intent(in)           ::  vname
  ! 	zname is z-matrix file, qname contains quaternions and trans
  ! 	vname is control variable file

  !Integers!
  integer                                           :: dunit     !unit nos for files
  integer, intent(out)                              :: num_zm    !number of z matrices
  integer, intent(out), dimension(zm_max)           :: num_mol   !no. of mol./cell per z-matrix
  integer, intent(out), dimension(3)                :: simsize
  integer                                           :: i,j,ii,k

  !Reals!
  real, intent(out), dimension(6)                        :: cellpar
  real, dimension(4)                                     :: tempquat
  real, intent(out)                                      :: v_min, v_max
  real, intent(out), dimension(mol_max,3)                :: translate

  !Custom Types!
  type (quaternion),  dimension(mol_max)                   :: quatern
  type (zmatrix_object), dimension(zm_max)                 :: zmat

  !Logicals!
  logical                                                  :: improp

  !------------------------------------------------------------------------------!
  !          Read in vname and print to stdout                                   !
  !------------------------------------------------------------------------------!

  dunit = open(vname, status='old') 
  read(dunit,*)num_zm
  if(num_zm.gt.zm_max) then
     write(stdout,*)'-------------------------------------------------------------------------'
     write(stdout,*)'Too many z-matrices, change limits in program.'
     write(stdout,*)'Program exiting.'
     write(stdout,*)'-------------------------------------------------------------------------'
     stop
  end if  ! num_zm test
  do i = 1,num_zm
     read(dunit,*) num_mol(i)
     if(num_mol(i).gt.mol_max) then
        write(stdout,*)'----------------------------------------------------------------------'
        write(stdout,*)'Too many molecules, change limits in program.'
        write(stdout,*)'Program exiting.'
        write(stdout,*)'----------------------------------------------------------------------'
        stop
     end if  ! num_mol test
     read(dunit,*) zname(i)
     write(stdout,'(" Z-matrix ",i3," of ",i3," will be read from                : ",a34)')  &
          i,num_zm,zname(i)
     read(dunit,*)qname(i)
     write(stdout,'(" Q and trans for z-matrix ",i3," of ",i3," will be read from: ",a34)')  &
          i,num_zm,qname(i)
     write(stdout,'(" There are ",i3," molecules of type ",i3," in the unit cell.")') &
          num_mol(i),i
  end do   !i = 1 , num_zm

  read(dunit,*)(cellpar(i),i=1,6)
  write(stdout,'(" Cell parameters are              : ",6f8.4)')(cellpar(i),i=1,6)
  read(dunit,*)(simsize(i),i=1,3)
  if((simsize(1).gt.amax).or.(simsize(2).gt.bmax).or.(simsize(3).gt.cmax)) then
     write(stdout,*)'-------------------------------------------------------------------------'
     write(stdout,*)'Model crystal too big, change limits in program or use smaller crystal.'
     write(stdout,*)'Program exiting.'
     write(stdout,*)'-------------------------------------------------------------------------'
     stop
  end if   !if simsize...
  write(stdout,'(" Simulation dimensions are        : ",3i4)')(simsize(i),i=1,3)
  read(dunit,*)v_min,v_max
  close(dunit)   !vname
  write(stdout,'(" Min. contact vector length is : ",f8.4,", max. is ",f8.4," Angstrom.")') &
       v_min,v_max
  write(stdout,'(" Control variable file closed     : ",a34)')vname
  write(stdout,*)'----------------------------------------------------------------------------'

  !------------------------------------------------------------------------------!
  !          Read in the z-matrix(ces) from zname and print to stdout            !
  !------------------------------------------------------------------------------!

  do i = 1,num_zm
     zmat(i) = load_zmatrix(zname(i))
     write(stdout,'(" Reading z-matrix ",i3," of ",i3," from : ",a34)')i,num_zm,zname(i)
     call print(stdout,zmat(i))
  end do !i = 1 , num_zm

  if (num_zm.lt.2)then
     write(stdout,*)'Z-matrix read in.'
  else
     write(stdout,*)'Z-matrices read in.'
  end if  !(num_zm.lt.2)

  !------------------------------------------------------------------------------!
  !          Read in the quaternions and xyz's from qname                        !
  !         The writing out part here is a direct steal from Aidan               !
  !------------------------------------------------------------------------------!

  write(stdout,*)'----------------------------------------------------------------------------'
  do i = 1,num_zm
     dunit = open(qname(i), status='old') 
     do ii = 1,num_mol(i)
        Read(dunit,*)j
        if(j.ne.ii)then
           write(stdout,*)'-------------------------------------------------------------------'
           write(stdout,*)'Something funny is going on with quaternion order.'
           write(stdout,*)'Program exiting.'
           write(stdout,*)'-------------------------------------------------------------------'
           stop
        end if    ! if(j.ne.ii)
        read(dunit,*)(tempquat(k),k=1,4)
        read(dunit,*)improp
        quatern(ii) = tempquat
        improp=improper(quatern(ii),improp)
        read(dunit,*)(translate(ii,k),k=1,3)
        write(stdout,'("  Sub structure: ",I4)') ii
        write(stdout,'("     Quaternion: ",4F11.6)') as_array(quatern(ii))
        write(stdout,'("       Improper: ",L4)') improp
        write(stdout,'("COM Translation: ",3F11.6)') translate(ii,:)
     end do !  ii = 1,num_mol(i)
     close(dunit)   !qname(i)
  end do   !i = 1 , num_zm
  return

end subroutine readfiles

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------

subroutine getcarts(simsize,xtal,num_zm,num_mol, zmat, qcrystal, xyzcrystal, carts, &
     numat,incrystal,nintern,intern)

  use quaternion_class
  use rotmatrix_class
  use zmatrix_class
  use crystallography_class
  use myparams

  implicit none

  !Characters!

  !Integers!
  integer,intent(in)                             :: num_zm !number of z matrices
  integer,intent(in), dimension(zm_max)          :: num_mol! no. of mol. in cell per z-matrix
  integer,intent(out), dimension(zm_max)         :: numat  ! number of atoms per z-matrix
  integer,intent(in), dimension(3)               :: simsize
  integer,intent(in), dimension(zm_max)          :: nintern 
  integer,intent(in),dimension(in_max,2,zm_max)  :: intern
  integer                                        :: ia,ib,ic,i,j,iz,m,in

  !Reals!
  real, intent(in),dimension(amax,bmax,cmax,zm_max,mol_max,3)         :: xyzcrystal
  real,intent(out),dimension(3,amax,bmax,cmax,mol_max,atom_max,zm_max):: carts
  real, dimension(3)                                                  :: shiftvec
  real,  dimension(3,atom_max)                                        ::  new_coords 
  real,intent(inout),dimension(amax,bmax,cmax,zm_max,mol_max,in_max)  :: incrystal

  !Custom Types!
  type (quaternion),intent(in),dimension(amax,bmax,cmax,zm_max,mol_max) :: qcrystal
  type (crystallography_object), intent(in)                             :: xtal
  type (zmatrix_object),intent(in),dimension(zm_max)                    :: zmat
  type (zmatrix_object)                                                 :: z1
  type (rotmatrix)                                                      :: rmat
  type (quaternion)::q

  !Logicals!
  logical                                                                    :: improp

  !------------------------------------------------------------------------------!
  !          First, fill the numat array: we will not using tha actual values    !
  !          this puts into new_coords.                                          !
  !------------------------------------------------------------------------------!

  do iz = 1,num_zm
     new_coords = as_xyz(zmat(1))
     numat(iz)=size(new_coords,2)
  end do ! iz = 1,num_zm

  !------------------------------------------------------------------------------!
  !          Now just fill up the carts array.  This needs updating the z-matrix !
  !          for each molecule, then getting its rmat and so on                  !
  !------------------------------------------------------------------------------!

  do ia = 1,simsize(1)
     do ib = 1,simsize(2)
        do ic = 1,simsize(3)
           ! shiftvec = some function of ia, ib, ic
           ! NOTE that shiftvec is commented out, UNLIKE in the contact vector calculation
           !shiftvec=as_cartesian(xtal,real((/ia,ib,ic/)))
           do iz = 1,num_zm
              z1 =zmat(iz)
              do m = 1,num_mol(iz)
                 in = nintern(iz)
                 incrystal(ia,ib,ic,iz,m,1:in)= &
                      parameter(z1,intern(1:in,1,iz),intern(1:in,2,iz), & 
                      incrystal(ia,ib,ic,iz,m,1:in))
                 new_coords = as_xyz(z1)
                 ! Make a rotation matrix out of this quaternion
                 rmat = as_rotmatrix(qcrystal(ia,ib,ic,iz,m))
                 ! Pre-rotate the coordinates using this rotation matrix
                 call rotate(rmat,new_coords)
                 ! Add on the translation. and the shiftvec
                 do j = 1,size(new_coords,2)
                    do i = 1,3
                       carts(i,ia,ib,ic,m,j,iz)=new_coords(i,j)+  &  
                            xyzcrystal(ia,ib,ic,iz,m,i) ! + shiftvec(i)
                    end do ! i = 1,3
                 end do ! j = 1,size
              end do !  m=1,num_mol
           end do ! iz = 1,num_zm
        end do ! ic
     end do ! ib
  end do ! ia
  return

end subroutine getcarts

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------

subroutine MCreadat(ename,outname, intern, nintern,num_zm)

  use myparams
  use file_functions

  implicit none

  !Characters!
  character(len=100),intent(in)                               ::  ename
  character(len=100),intent(out)                              :: outname
  character(len=100)                                          :: header
  character(len = 12)                                         :: fred

  !Reals!
  real                                                        :: temp
  real,            dimension(7)                               :: widths
  real,            dimension(in_max,zm_max)                   :: inwidths
  real,            dimension(ty_max)                          :: sprcon
  real,            dimension(7+in_max)                        :: initw

  !Integers!
  integer                                                     :: num_zm
  integer                                                     :: ntyp,ncy
  integer, intent(out),dimension(zm_max)                      :: nintern
  integer, intent(out), dimension(in_max,2,zm_max)            :: intern
  integer                                                     :: dunit, i,iz,iiz,miz

  !------------------------------------------------------------------------------!
  ! Open the file with MC input in it                                            !
  ! It's header will go into the output file to keep track of simulations        !
  !------------------------------------------------------------------------------!
  dunit = open(ename, status = 'old') 
  read(dunit,'(a100)')header
  !  write(stdout,*)'----------------------------------------------------------------------------'
  !  write(stdout,'("Header of MC file ",a20, "is ")')ename
  !  write(stdout,'(a60)')header

  !------------------------------------------------------------------------------!
  ! outname is where the outputted variables will go.                            !
  !------------------------------------------------------------------------------!

  read(dunit,*)outname
  write(stdout,*)'Will read variables from ',outname
  !------------------------------------------------------------------------------!
  !   Here we read in a bunch of stuff relevant to the MC part of the program    !
  !------------------------------------------------------------------------------!

  read(dunit,*)temp
  !  write(stdout,'("Simulation temperature is ",f4.2)')temp
  read(dunit,*)ncy
  !  write(stdout,'("Will do ",i4," cycles per molecule")')ncy
  read(dunit,*)ntyp
  !  write(stdout,'("There are ",i4," types of contact vectors of stengths...")')ntyp
  do i = 1,ntyp
     read(dunit,*)sprcon(i)
     !     write(stdout,'(f7.2)',advance = 'no')sprcon(i)
  end do
  !  write(stdout,*)
  do i = 1,7
     read(dunit,*)widths(i),initw(i)
  end do
  !  write(stdout,'("Widths of ditributions for increments on x,y,z are       ",3f7.5)') &
  !       (widths(i),i=1,3)
  !  write(stdout,'("Widths of ditributions for increments on q1,q2,q3,q4 are ",4f7.5)') &
  !       (widths(i),i=4,7)
  !  write(stdout,'("Widths of initial distributions for x,y,z are       ",3f7.5)') &
  !       (initw(i),i=1,3)
  !  write(stdout,'("Widths of initial distributions for q1,q2,q3,q4 are ",4f7.5)') &
  !       (initw(i),i=4,7)
  miz = 0
  do iiz = 1,num_zm
     read(dunit,*)iz, nintern(iz)
     if(iz.gt.miz)miz=iz
     write(stdout,'("There are ",i4," internal d.f. for z-matrix ",i2," which are:")')  &
          nintern(iz),iz
     do i = 1,nintern(iz)
        read(dunit,*)intern(i,1,iz),intern(i,2,iz),inwidths(i,iz),initw(7+i)
        if(intern(i,2,iz).eq.1) fred = 'bond length'
        if(intern(i,2,iz).eq.2) fred = 'bond angle'
        if(intern(i,2,iz).eq.3) fred = 'dihedral'
        if((intern(i,2,iz).gt.3).or.(intern(i,2,iz).lt.1)) then
           write(stdout,*)'--------------------------------------------------------------'
           write(stdout,*)'A non-existant internal degree of freedom has been specified.'
           write(stdout,*)'Program exiting.'
           write(stdout,*)'--------------------------------------------------------------'
           stop
        end if
        write(stdout,'(a12," of atom ",i3," of z-mat ",i2)') &
             fred,intern(i,1,iz),iz
        ! Note that initial widths for internals are not echoed to stdout.
        ! If intern(,2,) is type 2 or 3, convert the increment into one in radians
        if(intern(i,2,iz).gt.1) inwidths(i,iz)=inwidths(i,iz)*(3.14159/180.0)
     end do
  end do
  close(dunit)
  if(miz.ne.num_zm) then
     write(stdout,*)'--------------------------------------------------------------'
     write(stdout,*)'Number of internal df not/wrongly specified.'
     write(stdout,*)'Program exiting.'
     write(stdout,*)'--------------------------------------------------------------'
     stop
  end if
  return

End subroutine MCreadat


!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------

subroutine readcrystal(outname, simsize, num_zm, num_mol, qcrystal, xyzcrystal, incrystal, nintern)

  Use file_functions
  use quaternion_class
  use myparams

  implicit none

  !Characters!
  character(len=100), intent(in)                                             :: outname
  character(len=100)                                                         :: header

  !Integers!
  integer :: ia,ib,ic,im,iz,i,dunit,mcell,irow,trow
  integer, intent(in),dimension(zm_max)                                      :: num_mol
  integer, dimension(3)                                                      :: simsize
  integer, intent(in)                                                        :: num_zm
  integer, dimension(zm_max),intent(in)                                      :: nintern 

  !Reals!
  real, intent(out),                dimension(amax,bmax,cmax,zm_max,mol_max,3):: xyzcrystal
  real,  intent(out),          dimension(amax,bmax,cmax,zm_max,mol_max,in_max):: incrystal
  real, dimension(4)::ddd
  !Custom Types!
  type (quaternion),  dimension(amax,bmax,cmax,zm_max,mol_max)               :: qcrystal

  !Logicals!
  logical                                                                    :: improp

  ! First get the number of rows in the file.
  ! This should equal the number of molecules, which is a*b*c*mol per cell
  !

  mcell = 0
  do i = 1,num_zm
     mcell = mcell+num_mol(i)
  end do ! i = 1,num_zm
  trow = simsize(1)*simsize(2)*simsize(3)*mcell
  !------------------------------------------------------------------------------!
  ! Loop over crystal reading out the xyz, q and internals for each molecule     !
  !------------------------------------------------------------------------------!
  dunit = open(outname)
  read(dunit,'(a70)')header
  write(stdout,*)'Header of crystal variable file is....'
  write(stdout,*)header
  do irow = 1,trow
     read(dunit,'(5i3," ",7(f8.5," "),L3)', advance='no')ia,ib,ic,iz,im,&
          (xyzcrystal(ia,ib,ic,iz,im,i),i = 1,3), (ddd(i),i=1,4),improp
     read(dunit,*)(incrystal(ia,ib,ic,iz,im,i),i=1,nintern(iz))
     qcrystal(ia,ib,ic,iz,im)=as_quaternion(ddd)
     improp = improper(qcrystal(ia,ib,ic,iz,im),improp)
!!$  !    Test the reading in by writing out to some other unit.
!!$          improp = improper(qcrystal(ia,ib,ic,iz,im))
!!$          write(8,'(5i3," ",7(f8.5," "),L3," ")', advance = 'no')ia,ib,ic,iz,im,&
!!$               xyzcrystal(ia,ib,ic,iz,im,:),&
!!$               as_array(qcrystal(ia,ib,ic,iz,im)),improp
!!$          write(8,'(5f8.5)')incrystal(ia,ib,ic,iz,im,1:nintern(iz))
  end do ! irow
  close(dunit)
  return

end subroutine readcrystal
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------

subroutine cart2frac(carts, simsize,num_zm,num_mol,xtal,zmat,occc,numat,xyz)

  use crystallography_class
  use zmatrix_class
  use myparams
  use file_functions

  implicit none

  !Integers!
  integer,intent(in)                                                  :: num_zm 
  integer, dimension(zm_max),intent(in)                               :: num_mol,numat
  integer, dimension(3),intent(in)                                    :: simsize
  integer                                         :: ia,ib,ic,im,iz, dummycount
  integer                                         :: sitecount,iat,atoms_cell,ind,in2

  !Reals!
  real, dimension(amax,bmax,cmax,cell_atoms_max,3),intent(out)        :: xyz
  real, dimension(cell_atoms_max,3)                                   :: cvec, fvec 
  real,dimension(3,amax,bmax,cmax,mol_max,atom_max,zm_max) ,intent(in):: carts

  !Custom types!
  type (crystallography_object),intent(in)                            :: xtal
  type (zmatrix_object),dimension(zm_max),intent(in)                  :: zmat
  type (zmatrix_object)                                               :: testzmat

  !Characters!
  CHARACTER(len=2), intent(out) ,dimension(cell_atoms_max)            :: occc
  character(len=6),dimension(zm_max,atom_max)                         :: lable
  character(len=2),dimension(zm_max,atom_max)                         :: lab_short

  atoms_cell = 0
  do iz = 1,num_zm
     atoms_cell = atoms_cell+(num_mol(iz)*numat(iz))
  end do

  do iz = 1,num_zm
     testzmat = zmat(iz)
     ! Convert to len=2
     do iat = 1,numat(iz)
        lable(iz,iat)=labels(testzmat,iat)
        in2 = max((scan(lable(iz,iat),'0123456789')-1),1)
        ind = min(2,in2)
        lab_short(iz,iat)=lable(iz,iat)(1:ind)
     end do
  end do
  dummycount = 0
  sitecount=0
  do iz = 1,num_zm  
     do im = 1,num_mol(iz)
        do iat = 1,numat(iz)
           if((lab_short(iz,iat)(1:1)).eq.'x') then
              dummycount = dummycount +1
              ! do nothing
           else
              sitecount = sitecount+1
              occc(sitecount)=lab_short(iz,iat)
           end if
        end do
     end do
  end do
  if(sitecount.ne.atoms_cell-dummycount) then
     write(stdout,*)sitecount,dummycount,atoms_cell
     write(stdout,*)'--------------------------------------------------------------'
     write(stdout,*)'Something funny is going on with the number of sites.'
     write(stdout,*)'Program exiting.'
     write(stdout,*)'--------------------------------------------------------------'
     stop
  end if

  do ia = 1,simsize(1)
     do ib = 1,simsize(2)
        do ic = 1,simsize(3)
           sitecount = 0
           dummycount = 0
           do iz =1,num_zm
              do im = 1,num_mol(iz)
                 do iat = 1,numat(iz)
                    if((lab_short(iz,iat)(1:1)).eq.'x') then
                       dummycount = dummycount +1
                       !do nothing
                    else
                       sitecount = sitecount+1
                       cvec(sitecount,:)=carts(:,ia,ib,ic,im,iat,iz)
                    end if
                 end do  ! iat
              end do ! im
           end do ! iz
           ! Now cvec is full of site coords for one unit cell, convert to fractional
           if(sitecount.ne.atoms_cell-dummycount) then
              write(stdout,*)sitecount,dummycount,atoms_cell
              write(stdout,*)'--------------------------------------------------------------'
              write(stdout,*)'Something funny is going on with the number of sites.'
              write(stdout,*)'Program exiting.'
              write(stdout,*)'--------------------------------------------------------------'
              stop
           end if
           do iat = 1,sitecount
              xyz(ia,ib,ic,iat,:)=as_fractional(xtal,cvec(iat,:))
           end do ! iat
        end do ! ic
     end do ! ib
  end do ! ia

  return  

end subroutine cart2frac


