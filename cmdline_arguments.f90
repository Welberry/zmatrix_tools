module cmdline_arguments

  implicit none

  private

  ! This module provides an easy to use interface to the iargc and getarg
  ! compatibilty functions provided by the Intel Fortran Compiler (and
  ! supported by most Fortran compilers). These are used to get access to
  ! command line arguments.
  !
  ! Global variables are used to keep track of the argument object for each
  ! particular invocation. This is justifiable as no program has more than 
  ! one set of command line arguments (well I can't think of an instance
  ! where this would be the case, so I'm not supporting it). This makes the
  ! module much simpler and easier to use.

  !! $Log: cmdline_arguments.f90,v $
  !! Revision 1.1  2003/09/11 04:30:13  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: cmdline_arguments.f90,v 1.1 2003/09/11 04:30:13 aidan Exp $"

  integer, parameter :: arglength = 50

  ! We define an argument object which contains an array for the argument
  ! list and two integers, the number of arguments (the length of the 
  ! argument list) and the current argument, which is used to keep track
  ! of our current position when traversing the list
  type argument_object
     integer :: number_of_arguments, current_argument
     character(len=arglength), pointer, dimension(:) :: argument_list
  end type argument_object

  ! We define a number of interfaces to internal routines -- provides 
  ! easy to type and remember versions of longer names
  interface nargs
     module procedure number_of_arguments
  end interface

  interface current
     module procedure current_argument
  end interface

  interface have_args
     module procedure arguments_remaining
  end interface

  interface next_arg
     module procedure next_argument
  end interface

  ! Need to provide an explicit interface to the iargc function otherwise
  ! the compiler will complain due to our use of implicit none
  interface 
     integer(kind=4) function iargc()
     end function iargc
  end interface

  ! Do not *need* an interface to getarg, but it is good practice to 
  ! provide one
  interface getarg
     subroutine getarg(argindex, argument)
       integer(4), intent(in) :: argindex
       character(len=*), intent(out) :: argument
     end subroutine getarg
  end interface

  !!!!!!!!!!!!!!!!
  ! Public stuff !
  !!!!!!!!!!!!!!!!

  ! Only provide two public routines -- these are a pidgeon pair of 
  ! routines which traverse the argument array one by one returning 
  ! command line arguments
  public :: have_args, next_arg

  !!!!!!!!!!!!!!!!!!!!
  ! Global variables !
  !!!!!!!!!!!!!!!!!!!!

  ! We have a command line argument object which we save so that we can 
  ! initialise it once and then change the status of the initialised flag
  type (argument_object), save :: this_argument_object
  logical, save :: initialised = .false.

contains

  logical function arguments_remaining()

    ! This is the "main" function in this module. It is the first 
    ! function that is called, so it initialises our global argument
    ! object on the first pass, and then tells the user if there are
    ! any arguments. Subsequent invocations are usually after
    ! the 'next_arg' function, which returns the current argument and
    ! then increments the current argument counter.

    ! This function is used in loops like so
    !
    !         while (have_args()) do
    !            argument = next_arg()
    !         end do
    ! 
    ! or in single statements like
    !
    !         if (have_args()) argument = next_arg()
    ! 

    ! On the first pass through this function we will initialise the
    ! the global argument object (this_argument_object)
    if (.not. initialised) then
       call initialise_argument_object(this_argument_object)
       initialised = .true.
    end if

    ! Set default for arguments_remaining
    arguments_remaining = .true.

    ! Set to false if the number of the current object is greater than
    ! or equal to the number of arguments
    if (current(this_argument_object) >= nargs(this_argument_object)) &
         arguments_remaining = .false.

  end function arguments_remaining


  function next_argument()

    ! Return the next argument from the argument array of the globally defined 
    ! argument object 'this_argument_object'

    character(len=arglength) :: next_argument

    ! This next line (split over several lines to aid readability) does the 
    ! following tasks: 
    !
    !     o gets the current argument number and increments it by one
    !     o passes this incremented value back into the argument object
    !     o uses the return value of setting the incremented value as
    !       an index of the argument array
    !     o sets the next argument to this member of the argument array
    !
    next_argument = this_argument_object%argument_list(                      &
                                 current(this_argument_object,               &
                                         current(this_argument_object) + 1   &
                                        )                                    &
                                 )

  end function next_argument


  subroutine initialise_argument_object(arg_object)

    ! This is the workhorse routine which grabs the command line
    ! arguments and stores them internally in the argument object
    ! passed to it from 'arguments_remaining'
    
    type (argument_object), intent(out) :: arg_object

    integer           :: i, error, current_argument
    character(len=arglength) :: arg

    ! Initialise the argument list to a null pointer
    nullify(arg_object%argument_list)
    
    ! Set the current argument index to zero
    current_argument = 0
    current_argument = current(arg_object,current_argument)

    ! Grab the number of arguments and assign this internally in
    ! one step (nargs returns the number we put in as a matter of course)
    if (nargs(arg_object,iargc()) > 0) then

       ! Allocate some memory for the argument array
       allocate(arg_object%argument_list(nargs(arg_object)), stat=error)
       if (error /= 0) stop 'Error allocating memory for argument list'

       ! Cycle through the arguments and save in the argument array
       do i=1,nargs(arg_object) 
          call getarg(i,arg_object%argument_list(i))
       enddo

    endif
    
  end subroutine initialise_argument_object


  integer function number_of_arguments(arg_object, number)

    ! Simple function to set and/or return the number of arguments

    type (argument_object), intent(out) :: arg_object
    integer, intent(in), optional     :: number

    if (present(number)) then
       arg_object%number_of_arguments = number
    end if
    number_of_arguments = arg_object%number_of_arguments 

  end function number_of_arguments

  integer function current_argument(arg_object, number)

    ! Allows us to return (and optionally) set the number of the
    ! current argument -- this is for traversing the argument array

    type (argument_object), intent(out) :: arg_object
    integer, intent(in), optional     :: number

    ! If we specify a number then set this to tbe the current argument
    if (present(number)) then
       arg_object%current_argument = number
    end if

    ! Always return the value in the current argument
    current_argument = arg_object%current_argument 

  end function current_argument

end module cmdline_arguments
