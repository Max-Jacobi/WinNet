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
!!           - 25.08.23, M. Jacobi   , cumulative flow integration
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

  real(r_kind), parameter :: flow_limit = 1d-99           !< ignore smaller flows in output
  integer :: flow_size                                    !< size of flow array
  real(r_kind)  :: cum_dt                                 !< dt since last momentary flow output
  real(r_kind)  :: total_dt                               !< total integrated dt, never reset
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
!!          - 25.08.23, M. Jacobi cumulative flow integration
subroutine flow_init()
   implicit none

   integer :: i, j, ij, ji
   integer :: zi, zj, ni, nj, ai, aj, da
   type(flow_type) :: f

   INFO_ENTRY("flow_init")

   cum_dt = 0
   total_dt = 0
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
   if ((flow_every .gt. 0) &
#ifdef USE_HDF5
       .or. (h_flow_every .gt. 0) &
#endif
       ) allocate(flows(flow_size))
   if ((cum_flow_every .gt. 0) &
#ifdef USE_HDF5
       .or. (h_cum_flow_every .gt. 0) &
#endif
       ) allocate(cum_flows(flow_size))

   if ((.not. allocated(flows)) .and. (.not. allocated(cum_flows))) then
      call raise_exception("flow_init called without any flow output enabled", "flow_init")
   end if


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

         f%i = i
         f%j = j
         f%ij = ij
         f%ji = ji
         f%zi = zi
         f%zj = zj
         f%ni = ni
         f%nj = nj
         f%fl = 0

         if (allocated(flows))     flows(flow_size) = f
         if (allocated(cum_flows)) cum_flows(flow_size) = f
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
!!          - 25.08.23, M. Jacobi cumulative flow integration
!! .
subroutine flowcalc(Y, dt)
   implicit none

   ! MJ: these could in principle be used from single_zone_vars
   real(r_kind), dimension(:), intent(in)  :: Y  !< abundances
   real(r_kind), intent(in)                :: dt !< time step

   real(r_kind) :: fl
   integer :: i, j, ij, ji, n
   logical :: momentary, cumulative

   INFO_ENTRY("flowcalc")

   momentary  = allocated(flows)
   cumulative = allocated(cum_flows)

   ! update cumulative time steps
   cum_dt = cum_dt + dt
   total_dt = total_dt + dt

   do n = 1, flow_size
      if (momentary) then
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
      ! (can sometimes happen in the beginning)
      if (abs(fl) < 1d99) then
         if (momentary)  flows(n)%fl = flows(n)%fl + fl
         if (cumulative) cum_flows(n)%fl = cum_flows(n)%fl + fl
      endif
   end do

   INFO_EXIT("flowcalc")

end subroutine flowcalc


!>
!! Sort flows and remove zero flows to prepare them for output.
!! Updates the output_* arrays.
!! If cumulative = .true., the cumulative flows are used,
!! otherwise the momentary flows are used.
!! In the second case, the flows array and cum_dt are reset to zero.
!!
!!
!! \b Edited:
!!         - 25.08.23, M. Jacobi cumulative flow integration
!! .
subroutine flowsort(cumulative)
   use file_handling_class
   use single_zone_vars, only: Y
   implicit none

   logical, intent(in)   :: cumulative !< if true, cumulative flows are used

   integer               :: n
   real(r_kind)          :: fl
   type(flow_type), dimension(:), pointer :: flow_ptr

   INFO_ENTRY("flowsort")

   if (cumulative) then
      flow_ptr => cum_flows
   else
      flow_ptr => flows
   end if

   ! count flows in output and set output_n_flows
   ! ponytail: ASCII and HDF5 momentary output share one accumulator/window;
   ! if both are enabled, whichever outputs first resets the other's average.
   ! Guarding cum_dt > 0 keeps the second call NaN-free (it writes an empty list).
   output_n_flows = 0
   do n = 1, flow_size
      fl = flow_ptr(n)%fl
      ! momentary flows are normalized to the accumulated dt
      if ((.not. cumulative) .and. (cum_dt > 0)) fl = fl / cum_dt

      ! cycle flows that are too small
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

   ! cumulative flows integrate over the whole run, momentary ones
   ! over the window since the last momentary output
   if (cumulative) then
      output_cum_dt = total_dt
   else
      output_cum_dt = cum_dt
   end if

   output_n_flows = 0
   do n = 1, flow_size
      fl = flow_ptr(n)%fl
      if ((.not. cumulative) .and. (cum_dt > 0)) fl = fl / cum_dt

      ! cycle flows that are too small
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
      write(flowfile, '(a, i4.4, a)')'flow/cum_flow_', cnt, '.dat'
   else
      write(flowfile, '(a, i4.4, a)')'flow/flow_', cnt, '.dat'
   endif

   flowunit= open_outfile (adjustl(flowfile))

   ! write header
   write(flowunit, '(4a23)') 'time', 'temp', 'dens', 'cum. timestep'
   write(flowunit, '(4es23.14)') t, t9, dens, output_cum_dt
   if (cumulative) then
      write(flowunit, '(2(2a5, a23),2a23)') 'nin', 'zin', 'yin', &
          'nout', 'zout', 'yout', 'flow [abundance]'
   else
      write(flowunit, '(2(2a5, a23),2a23)') 'nin', 'zin', 'yin', &
          'nout', 'zout', 'yout', 'flow [abundance/s]'
   endif

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
