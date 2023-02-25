module mod_io
    use iso_fortran_env, only:stdin => input_unit, stdout => output_unit, &
                              real32, real64, int32
    
    use mod_events, only: event
    use csv_module

    implicit none
    private
    public :: read_event_data, write_event_data, write_sim_data
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
    
    subroutine write_event_data(filename, array)
      character(len=*),intent(in) :: filename
      real,intent(in) :: array(:, :)
      type(csv_file) :: f
      logical :: status_ok
      integer(int32) :: i, n

      ! set optional inputs:
      call f%initialize(verbose = .true.)
      ! open the file
      call f%open(filename, n_cols = size(array, 2), status_ok=status_ok)
      ! add header
      call f%add([filename])
      call f%next_row()
      ! add some data:
      n = size(array, 1)
      do i = 1, n
        call f%add(array(i, :), real_fmt='(F20.15)')
        call f%next_row()
      end do
      ! finished
      call f%close(status_ok)
      
    end subroutine

    subroutine write_sim_data(filename, time, outcome, group)
      character(len=*),intent(in) :: filename
      real(real32),intent(in) :: time(:)
      integer(int32),intent(in) :: outcome(:)
      integer(int32),intent(in) :: group(:)

      type(csv_file) :: f
      logical :: status_ok
      integer(int32) :: i, n
      ! set optional inputs:
      call f%initialize(verbose = .true.)
      ! open the file
      call f%open(filename, n_cols = 3, status_ok=status_ok)
      ! add some data:
      n = size(time, 1)
      do i = 1, n
        call f%add([time(i)], real_fmt='(F20.15)')
        call f%add([outcome(i), group(i)])
        call f%next_row()
      end do
      ! finished
      call f%close(status_ok)
      
    end subroutine write_sim_data
end module mod_io