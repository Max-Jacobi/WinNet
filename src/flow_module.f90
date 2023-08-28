!> @file flow_module.f90
!!
!! The error file code for this file is ***W20***.
!! @brief Module \ref flow_module for calculating reaction flows
!!

!> Provides subroutines to calculate reaction flows
!!
!! @author  Christian Winteler
!! @date    07.10.10
!!
!! \b Edited:
!!           - 03.04.18, M. Jacobi   , Rewrote the module, the flows are now
!!                                     calculated with the help of the Jacobian
!!           - 22.01.21, M. Reichert , added more comments
!!           - 25.08.23, M. Jacobi   , cumulativ flow integration
!! .
#include "macros.h"
module flow_module
  use error_msg_class,  only: raise_exception
  use global_class,     only: net_size, ineu, ipro, ihe4, flow_type
  use pardiso_class,    only: jind, vals
  use global_class,     only: isotope_type, isotope
  use parameter_class,  only: flow_every, cum_flow_every
#ifdef USE_HDF5
  use parameter_class,  only: h_flow_every, h_cum_flow_every
#endif
  implicit none

  real(r_kind), parameter :: flow_limit = 1e-99           !< ignore smaller flows in output
  integer                 :: mode                         !< 1 for flows 2 for cum_flows 3 for both
                                                          !< will be set automatically on init
  integer :: flow_size                                    !< size of flow array
  real(r_kind)  :: cum_dt                                 !< cumulated dt
  type(flow_type), dimension(:), allocatable, Target :: flows     !> flows
  type(flow_type), dimension(:), allocatable, Target :: cum_flows !> cumulative flows


  ! arrays for output
  integer                                 :: output_n_flows !< number of flows in the output
  integer, dimension(:), allocatable      :: output_n_in    !< neutron number of in isotopes
  integer, dimension(:), allocatable      :: output_p_in    !< proton number of in isotopes
  integer, dimension(:), allocatable      :: output_n_out   !< neutron number of out isotopes
  integer, dimension(:), allocatable      :: output_p_out   !< proton number of out isotopes
  real(r_kind), dimension(:), allocatable :: output_y_in    !< abundances of in isotopes
  real(r_kind), dimension(:), allocatable :: output_y_out   !< abundances of out isotopes
  real(r_kind), dimension(:), allocatable :: output_flow    !< flows
  real(r_kind)                            :: output_cum_dt  !< cumulative time step

  !
  ! Public and private fields and methods of the module
  !
  public:: &
      flow_init, flowcalc, flowprint, flowsort

contains



!>
!! Initialise flow subroutine
!!
!! This subroutine counts the number of possible flows
!! and allocates the \ref flows array.
!!
!! \b Edited:
!!          - 11.01.14
!!          - 03.04.18, M. Jacobi
!!          - 25.08.23, M. Jacobi cumulativ flow integration
subroutine flow_init()
   implicit none

   integer :: i, j, ij, ji
   integer :: zi, zj, ni, nj, ai, aj, da

   INFO_ENTRY("flow_init")

   ! set mode for convenience later
   mode = 0
   if ((flow_every .gt. 0) &
#ifdef USE_HDF5
       .or. (h_flow_every .gt. 0) &
#endif
       ) mode = 1
   if ((cum_flow_every .gt. 0) &
#ifdef USE_HDF5
       .or. (h_cum_flow_every .gt. 0) &
#endif
       ) mode = mode + 2

   cum_dt = 0
   ! loop first to find the number of flows
   flow_size = 0
   do i = 1, net_size-1
      ai = isotope(i)%mass
      do j = i+1, net_size
         aj = isotope(j)%mass
         da = aj - ai

         ! exclude flows from projectiles to products but
         ! include fission flows which should have ai and da > 4
         if ((abs(da) > ai) .and. (ai <= 4)) cycle

         ij = jind(i, j)
         ji = jind(j, i)
         if ((ij .eq. 0) .and. (ji .eq. 0)) cycle

         flow_size = flow_size + 1
      end do
   end do

   ! allocate flows
   select case (mode)
   case (1)
      allocate(flows(flow_size))
   case (2)
      allocate(cum_flows(flow_size))
   case (3)
      allocate(flows(flow_size))
      allocate(cum_flows(flow_size))
   case default
      call raise_exception("unknown flow mode (this should never happen)", "flow_init")
   end select


   ! loop again to fill the index arrays
   flow_size = 0
   do i = 1, net_size-1
      ni = isotope(i)%n_nr
      zi = isotope(i)%p_nr
      ai = isotope(i)%mass
      do j = i+1, net_size
         nj = isotope(j)%n_nr
         zj = isotope(j)%p_nr
         aj = isotope(j)%mass

         da = aj - ai

         ! exclude flows from projectiles to products but
         ! include fission flows which should have ai and da > 4
         if ((abs(da) > ai) .and. (ai <= 4)) cycle

         ij = jind(i, j)
         ji = jind(j, i)

         if ((ij .eq. 0) .and. (ji .eq. 0)) cycle

         flow_size = flow_size + 1

         select case (mode)
            case (1)
               flows(flow_size)%i = i
               flows(flow_size)%j = j
               flows(flow_size)%ij = ij
               flows(flow_size)%ji = ji
               flows(flow_size)%zi = zi
               flows(flow_size)%zj = zj
               flows(flow_size)%ni = ni
               flows(flow_size)%nj = nj
               flows(flow_size)%fl = 0
            case (2)
               cum_flows(flow_size)%i = i
               cum_flows(flow_size)%j = j
               cum_flows(flow_size)%ij = ij
               cum_flows(flow_size)%ji = ji
               cum_flows(flow_size)%zi = zi
               cum_flows(flow_size)%zj = zj
               cum_flows(flow_size)%ni = ni
               cum_flows(flow_size)%nj = nj
               cum_flows(flow_size)%fl = 0
            case (3)
               flows(flow_size)%i = i
               flows(flow_size)%j = j
               flows(flow_size)%ij = ij
               flows(flow_size)%ji = ji
               flows(flow_size)%zi = zi
               flows(flow_size)%zj = zj
               flows(flow_size)%ni = ni
               flows(flow_size)%nj = nj
               flows(flow_size)%fl = 0

               cum_flows(flow_size)%i = i
               cum_flows(flow_size)%j = j
               cum_flows(flow_size)%ij = ij
               cum_flows(flow_size)%ji = ji
               cum_flows(flow_size)%zi = zi
               cum_flows(flow_size)%zj = zj
               cum_flows(flow_size)%ni = ni
               cum_flows(flow_size)%nj = nj
               cum_flows(flow_size)%fl = 0
         end select
      end do
   end do

   INFO_EXIT("flow_init")

end subroutine flow_init



!>
!! Flow calculation from jacobian. It is calculated with the help of the Jacobian.
!! \f[
!! F_{ij} = (J_{ji} \times Y_j - J_{ij} \times Y_i) \Delta t
!! \f]
!!
!! At each timestep the flows are cumulatively added to the flows and/or cum_flows arrays.
!! On each output iteration the flows array and the cumulative timestep is reset to zero.
!!
!! @note Using the jacobian directly has the advantage
!!       that the flow will be correct if the calculation
!!       is correct. In previous versions, the flow
!!       was not calculated by using the jacobian.
!!
!! \b Edited:
!!          - 03.04.18, M. Jacobi
!!          - 25.08.23, M. Jacobi cumulativ flow integration
!! .
subroutine flowcalc(Y, dt)
   implicit none

   ! MJ: these could in principle be used from single_zone_vars
   real(r_kind), dimension(:), intent(in)  :: Y  !< abundances
   real(r_kind), intent(in)                :: dt !< time step

   real(r_kind) :: fl
   integer :: i, j, ij, ji, n

   ! update cumulative time step
   cum_dt = cum_dt + dt

   do n = 1, flow_size
      if (mode.eq.1) then
         i = flows(n)%i
         j = flows(n)%j
         ij = flows(n)%ij
         ji = flows(n)%ji
      else
         i = cum_flows(n)%i
         j = cum_flows(n)%j
         ij = cum_flows(n)%ij
         ji = cum_flows(n)%ji
      end if

      ! add up momentary flows in temporary variable
      fl = 0
      if (ij.ne.0) fl = - vals(ij) * Y(i) * dt
      if (ji.ne.0) fl = fl + vals(ji) * Y(j) * dt

      ! add flow only if it is not abnormally large
      ! (can sometimnes happen in the beginning)
      if (abs(fl) < 1e99) then
      select case (mode)
      case (1)
         flows(n)%fl = flows(n)%fl + fl
      case (2)
         cum_flows(n)%fl = cum_flows(n)%fl + fl
      case (3)
         flows(n)%fl = flows(n)%fl + fl
         cum_flows(n)%fl = cum_flows(n)%fl + fl
      end select
      endif
   end do

   INFO_EXIT("flowcalc")

end subroutine flowcalc


!>
!! Sort flows and remove zero flows to prepare them for output.
!! Updates the output_* arrays.
!! If cumulative = .true., the cumultive flows are used,
!! otherwise the momentary flows are used.
!! In the second case, the flows array and cum_dt are reset to zero.
!!
!!
!! \b Edited:
!!         - 25.08.23, M. Jacobi cumulativ flow integration
!! .
subroutine flowsort(cumulative)
   use file_handling_class
   use single_zone_vars, only: Y
   implicit none

   logical, intent(in)   :: cumulative !< if true, cumulative flows are used

   integer               :: i, j, n, m
   real(r_kind)          :: fl
   type(flow_type), dimension(:), pointer :: flow_ptr

   INFO_ENTRY("flowsort")

   if (cumulative) then
      flow_ptr => cum_flows
   else
      flow_ptr => flows
   end if

   ! count flows in output and set output_n_flows
   output_n_flows = 0
   do n = 1, flow_size
      fl = flow_ptr(n)%fl
      ! if flow is used it is normalized to the cumulative dt
      if (.not. cumulative) fl = fl / cum_dt

      ! cycle flows that are to small
      if (abs(fl) < flow_limit) cycle
      output_n_flows = output_n_flows + 1
   end do

   ! allocate output arrays
   if (allocated(output_n_in)) then
      deallocate(output_n_in)
      deallocate(output_p_in)
      deallocate(output_n_out)
      deallocate(output_p_out)
      deallocate(output_y_in)
      deallocate(output_y_out)
      deallocate(output_flow)
   end if

   allocate(output_n_in(output_n_flows))
   allocate(output_p_in(output_n_flows))
   allocate(output_n_out(output_n_flows))
   allocate(output_p_out(output_n_flows))
   allocate(output_y_in(output_n_flows))
   allocate(output_y_out(output_n_flows))
   allocate(output_flow(output_n_flows))

   output_cum_dt = cum_dt

   output_n_flows = 0
   do n = 1, flow_size
      fl = flow_ptr(n)%fl
      if (.not. cumulative) fl = fl / cum_dt

      ! cycle flows that are to small
      if (abs(fl) < flow_limit) cycle
      output_n_flows = output_n_flows + 1

      ! set output arrays
      if (fl > 0) then
         output_n_in(output_n_flows) = flow_ptr(n)%ni
         output_p_in(output_n_flows) = flow_ptr(n)%zi
         output_n_out(output_n_flows) = flow_ptr(n)%nj
         output_p_out(output_n_flows) = flow_ptr(n)%zj
         output_y_in(output_n_flows) = Y(flow_ptr(n)%i)
         output_y_out(output_n_flows) = Y(flow_ptr(n)%j)
         output_flow(output_n_flows) = fl
      else ! revert flow
         output_n_in(output_n_flows) = flow_ptr(n)%nj
         output_p_in(output_n_flows) = flow_ptr(n)%zj
         output_n_out(output_n_flows) = flow_ptr(n)%ni
         output_p_out(output_n_flows) = flow_ptr(n)%zi
         output_y_in(output_n_flows) = Y(flow_ptr(n)%j)
         output_y_out(output_n_flows) = Y(flow_ptr(n)%i)
         output_flow(output_n_flows) = -fl
      endif
   end do

   ! reset cumulative timestep and flows
   if (.not. cumulative) then
      cum_dt = 0
      do n = 1, flow_size
         flows(n)%fl = 0
      end do
   endif

   INFO_EXIT("flowsort")

end subroutine flowsort

!>
!! Output reaction flows to a file
!!
!! An example of this file may look like:
!!\file{
!! time    temp    dens
!! 1.03895957612263E-01   7.19136097013393E+00   1.40977753502083E+06
!!  nin     zin     yin    nout    zout    yout    flow
!! 2   1   4.81807892321990E-08   1   1   2.13994533749120E-06   0.00000000000000E+00
!! 1   2   1.26489216252989E-09   1   1   2.13994533749120E-06   0.00000000000000E+00
!! 1   2   1.26489216252989E-09   2   1   4.81807892321990E-08   1.58426675189734E-10
!! 4   2   9.86495465952053E-13   3   3   2.15833022688002E-11   8.53665754802155E-13
!! ...}
!!
!! \b Edited:
!!         - 11.01.14
!!         - 03.04.18, M. Jacobi
!!         - 25.08.23, M. Jacobi cumulative flow integration
!! .
subroutine flowprint(t, t9, dens, cnt, cumulative)
   use global_class, only: isotope_type, isotope
   use file_handling_class
   implicit none

   real(r_kind), intent(in)  :: t                !< time [s]
   real(r_kind), intent(in)  :: t9               !< temperature [GK]
   real(r_kind), intent(in)  :: dens             !< density [g/cm3]
   integer, intent(in)       :: cnt              !< flow snapshot counter
   logical, intent(in)       :: cumulative       !< if true, cumulative flows are used

   integer      :: flowunit
   integer      :: n
   character*50 :: flowfile

   INFO_ENTRY("flowprint")

   ! prepare flows for output
   call flowsort(cumulative)

   ! determine filename
   if (cumulative) then
      write(flowfile, '(a13, i4.4, a4)')'flow/cum_flow_', cnt, '.dat'
   else
      write(flowfile, '(a13, i4.4, a4)')'flow/flow_', cnt, '.dat'
   endif

   flowunit= open_outfile (adjustl(flowfile))

   ! write header
   write(flowunit, '(4a23)') 'time', 'temp', 'dens', 'cum. timestep'
   write(flowunit, '(4es23.14)') t, t9, dens, output_cum_dt
   write(flowunit, '(2(2a5, a23),2a23)') 'nin', 'zin', 'yin', &
       'nout', 'zout', 'yout', 'flow [abundance/s]'

   ! write flows
   do n = 1, output_n_flows
      write(flowunit,'(2(2i5,es23.14E3),2es23.14E3)') &
          output_n_in(n), output_p_in(n), output_y_in(n), &
          output_n_out(n), output_p_out(n), output_y_out(n), &
          output_flow(n)
   end do

   call close_io_file(flowunit, flowfile)

   INFO_EXIT("flowprint")

end subroutine flowprint

end module flow_module
