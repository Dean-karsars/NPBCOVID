module mod_io
    use iso_fortran_env, only:stdin => input_unit, stdout => output_unit, &
                              real32, real64, int32
    
    use mod_events, only: event
    use csv_module

    implicit none
    private
    public :: read_event_data
contains
    subroutine read_event_data(filename, time, outcome, group)
    ! Read processed weather data from a CSV file.
    
      character(len=*), intent(in) :: filename
      real(real32), allocatable, intent(out) :: time(:)
      integer(int32), allocatable, intent(out) :: outcome(:)
      integer(int32), allocatable, intent(out) :: group(:)

      type(csv_file) :: f
      logical,dimension(:),allocatable :: t
      character(len=30),dimension(:),allocatable :: header
      logical :: status_ok
      integer,dimension(:),allocatable :: itypes

      ! read the file
      call f%read(filename, header_row=1, status_ok=status_ok)
      ! get the header and type info
      call f%get_header(header,status_ok)
      call f%variable_types(itypes,status_ok)

      ! get some data
      call f%get(1, time, status_ok)
      call f%get(2, outcome, status_ok)
      call f%get(3, group, status_ok)
      ! destroy the file
      call f%destroy()

      

    end subroutine read_event_data
    
end module mod_io