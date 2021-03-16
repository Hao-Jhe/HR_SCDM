
module ct07_ozone_mod

!---------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use            fms_mod, only: error_mesg, FATAL, file_exist,       &
                              open_namelist_file, set_domain,      &
                              read_data, check_nml_error,          &
                              mpp_pe, mpp_root_pe, close_file,     &
                              write_version_number, stdlog,        &
                              uppercase,&  
                              mpp_clock_id,mpp_clock_begin,mpp_clock_end,CLOCK_COMPONENT,& 
                              WARNING

use   time_manager_mod, only: time_type, get_time, Time_for_hs

use   diag_manager_mod, only: register_diag_field, send_data

use  field_manager_mod, only: MODEL_ATMOS, parse
use tracer_manager_mod, only: query_method, get_number_tracers, &
                              get_tracer_index, NO_TRACER
use   interpolator_mod, only: interpolate_type, interpolator_init, &
                              interpolator, interpolator_end, &
                              CONSTANT, INTERP_WEIGHTED_P, INTERP_LINEAR_P
use   zenith_angle_mod, only: zenith_angle

implicit none
private

!------ Cariolle & Teyssedre (2007) ----------------------------------

  public :: ct07_ozone, ct07_ozone_init, calc_column

   type(interpolate_type),save  ::  a1_interp
   type(interpolate_type),save  ::  a2_interp
   type(interpolate_type),save  ::  a3_interp
   type(interpolate_type),save  ::  a4_interp
   type(interpolate_type),save  ::  a5_interp
   type(interpolate_type),save  ::  a6_interp
   type(interpolate_type),save  ::  a7_interp
   type(interpolate_type),save  ::  a8_interp
   type(interpolate_type),save  ::  a9_interp

   type(interpolate_type),save  ::  ozone_interp

!------ namlist ------------------------------------------------------

!  defult option

   logical :: no_ozone = .false.

   character(len=256) :: ozone_option = 'ct07_ozone' ! 'ct07_ozone' / 'from_file'
   character(len=256) :: ozone_file   = 'O3'

   logical :: anthropogenic_forcing = .true.         ! option only for ozone_option = 'cto7_ozone'
   logical :: cold_tracer_forcing   = .true.         ! option only for anthropogenic_forcing = '.true.'
   logical :: cold_tracer_output    = .true.
   logical :: ozone_tendency_output = .true.
   logical :: total_ozone_output    = .true.

!---------------------------------------------------------------------

   namelist /ct07_ozone_nml/ no_ozone,                              &
                             ozone_option,                          &
                             ozone_file,                            &
                             anthropogenic_forcing,                 &
                             cold_tracer_forcing,                   &
                             cold_tracer_output,                    &
                             ozone_tendency_output,                 &
                             total_ozone_output 

!---------------------------------------------------------------------

   character(len=128) :: version='$Id: ct07_ozone.f90, 2019/03/23 hjh $'
   character(len=128) :: tagname='$Name: riga_O3_SW_control $'

   integer            :: id_o3
   integer            :: id_to3
   integer            :: id_do3dt
   integer            :: id_rtnd_do3, id_rtnd_dtemp, id_rtnd_dto3, id_rtnd_anthropogenic
   integer            :: id_cold_tracer
   real               :: missing_value = -1.e10
   character(len=14)  :: mod_name = 'ct07_ozone'


   logical :: module_is_initialized = .false.


contains

!#####################################################################

 subroutine ct07_ozone_init ( axes, Time, lonb, latb, external_ozone )

   integer, intent(in)                        :: axes(4)
   type(time_type), intent(in)                :: Time
   real, intent(in), optional, dimension(:,:) :: lonb, latb
   logical, intent(out)                       :: external_ozone

!---------------------------------------------------------------------

    integer :: unit, io, ierr

! --- read namelist ---

 #ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=ct07_ozone_nml, iostat=io)
     ierr = check_nml_error(io, 'ct07_ozone_nml')
 #else
     if (file_exist('input.nml')) then
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
           read  (unit, nml=ct07_ozone_nml, iostat=io, end=10)
           ierr = check_nml_error (io, 'ct07_ozone_nml')
        enddo
  10    call close_file (unit)
     endif
 #endif

 if ( trim(ozone_option) == 'ct07_ozone' .or. trim(ozone_option) == 'from_file' ) then

! --- write version info and namelist to log file ---

    call write_version_number (version,tagname)
    if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=ct07_ozone_nml)
    
    if (no_ozone) return

 else
   call error_mesg( 'ct07_ozone_nml', trim(ozone_option)//' is not valid for ozone_option', FATAL)
 endif
  
! --- register diagnostic fields ----

    id_o3 = register_diag_field ( mod_name, 'o3', axes(1:3), Time,   &
                                  'ozone', 'mol/mol',                &
                                  missing_value=missing_value        )

    if ( total_ozone_output ) then
      id_to3 = register_diag_field ( mod_name, 'to3', axes(1:2), Time, &
                                     'total ozone', 'DU',              &
                                     missing_value=missing_value       ) 
    endif
   
    if ( ozone_tendency_output ) then
                   id_do3dt = register_diag_field ( mod_name, 'do3dt', axes(1:3), Time, &
                                                    'ozone tendency', 'mol mol-1 s-1',  & 
                                                    missing_value=missing_value         ) 
                id_rtnd_do3 = register_diag_field ( mod_name, 'rtnd_do3', axes(1:3), Time, &
                                                    'rtnd_do3', 'mol mol-1 s-1',           &
                                                    missing_value=missing_value            )
              id_rtnd_dtemp = register_diag_field ( mod_name, 'rtnd_dtemp', axes(1:3),Time, &
                                                    'rtnd_dtemp', 'mol mol-1 s-1',          &
                                                    missing_value=missing_value             )
               id_rtnd_dto3 = register_diag_field ( mod_name, 'rtnd_dto3', axes(1:3), Time, &
                                                    'rtnd_dto3', 'mol mol-1 s-1',           &
                                                    missing_value=missing_value             )
      id_rtnd_anthropogenic = register_diag_field ( mod_name,'rtnd_anthropogenic', axes(1:3), Time, &
                                                    'rtnd_anthropogenic', 'mol mol-1 s-1',           &
                                                    missing_value=missing_value                     )
    endif
 
    if ( cold_tracer_output ) then 
      id_cold_tracer = register_diag_field ( mod_name, 'cold_tracer', axes(1:3), Time, &
                                             'cold_tracer', 'none',                    &
                                             missing_value=missing_value               )
    endif

!   ---- set Cariolle & Teyssedre (2007) coefficient ----

    if ( trim(ozone_option) == 'ct07_ozone' ) then
      call interpolator_init (a1_interp, 'A1.nc', lonb, latb, data_names=(/'A1'/), data_out_of_bounds=(/CONSTANT/))
      call interpolator_init (a2_interp, 'A2.nc', lonb, latb, data_names=(/'A2'/), data_out_of_bounds=(/CONSTANT/))
      call interpolator_init (a3_interp, 'A3.nc', lonb, latb, data_names=(/'A3'/), data_out_of_bounds=(/CONSTANT/))
      call interpolator_init (a4_interp, 'A4.nc', lonb, latb, data_names=(/'A4'/), data_out_of_bounds=(/CONSTANT/))
      call interpolator_init (a5_interp, 'A5.nc', lonb, latb, data_names=(/'A5'/), data_out_of_bounds=(/CONSTANT/))
      call interpolator_init (a6_interp, 'A6.nc', lonb, latb, data_names=(/'A6'/), data_out_of_bounds=(/CONSTANT/))
      call interpolator_init (a7_interp, 'A7.nc', lonb, latb, data_names=(/'A7'/), data_out_of_bounds=(/CONSTANT/))
      call interpolator_init (a8_interp, 'A8.nc', lonb, latb, data_names=(/'A8'/), data_out_of_bounds=(/CONSTANT/))
      call interpolator_init (a9_interp, 'A9.nc', lonb, latb, data_names=(/'A9'/), data_out_of_bounds=(/CONSTANT/))
      external_ozone = .false.
    endif

!   ---- set ozone from file ----

    if ( trim(ozone_option) == 'from_file' ) then
      call interpolator_init (ozone_interp, trim(ozone_file)//'.nc', lonb, latb, data_names=(/trim(ozone_file)/), data_out_of_bounds=(/CONSTANT/))
      external_ozone = .true.
    endif


    module_is_initialized = .true.
    
 end subroutine ct07_ozone_init

!#######################################################################

 subroutine ct07_ozone ( is, ie, js, je, dt, Time, lon, lat, p_half, p_full,      &
                         tm, rm, rdt, ozone_init, ct07_coeff, count_int )

!-----------------------------------------------------------------------

   integer, intent(in)                          :: is, ie, js, je
   integer, intent(in)                          :: count_int
      real, intent(in)                          :: dt
 type(time_type), intent(in)                    :: Time
      real, intent(in),    dimension(:,:      ) :: lon, lat
      real, intent(in),    dimension(:,:,:    ) :: p_half, p_full
      real, intent(in),    dimension(:,:,:    ) :: tm
      real, intent(inout), dimension(:,:,:,:  ) :: rm
      real, intent(inout), dimension(:,:,:,:  ) :: rdt
      real, intent(inout), dimension(:,:,:,:  ) :: ozone_init
      real, intent(inout), dimension(:,:,:,:,:) :: ct07_coeff

!-----------------------------------------------------------------------

      real, dimension(size(rm,1),size(rm,2))            :: to3
      real, dimension(size(rm,1),size(rm,2),size(rm,3)) :: rst, rtnd, rtnd_anthropogenic
      real, dimension(size(rm,1),size(rm,2),size(rm,3)) :: rst_future
      real, dimension(size(tm,1),size(tm,2),size(tm,3)) :: tst
      real, dimension(size(rm,1),size(rm,2),size(rm,3)) :: r_column
      real, dimension(size(rm,1),size(rm,2),size(rm,3)) :: ozone

      real, dimension(size(rm,1),size(rm,2))            :: cos_theta   ! cos(zenith_angle)
      real, dimension(size(rm,1),size(rm,2),size(rm,3)) :: cold_tracer, cold_tracertnd
      real, dimension(size(rm,1),size(rm,2),size(rm,3)) :: a1, a2, a3, a4, a5, a6, a7, a8, a9

      real, dimension(size(rm,1),size(rm,2),size(rm,3)) :: rtnd_do3, rtnd_dtemp, rtnd_dto3
      real, dimension(size(rm,1),size(rm,2),size(rm,3)) :: do3, dtemp, dto3

      integer :: i, j, k, n
      integer :: num_tracers
      integer :: n_ozone
      integer :: n_cold_tracer
      integer :: days, seconds
      logical :: used

      real    :: du2mol = 2.69e+16  ! [DU] = 2.69e+16 [molecules cm-2]

!---------------------------------------------------------------------

      type(time_type) :: Time2
      integer         :: idays_init, t_days_ind
      integer         :: days_of_year = 365

!---------------------------------------------------------------------

    if (no_ozone) return

    if ( .not. module_is_initialized) call error_mesg ('ct07_ozone','ct07_ozone_init has not been called', FATAL)

!---------------------------------------------------------------------

! === read ozone from file ===

  if ( trim(ozone_option) == 'from_file' ) then

    if ( count_int == 1 ) then
      do idays_init = 1, days_of_year
        call Time_for_hs ( Time2, idays_init )
        call get_ozone ( Time2, p_half, ozone )
        ozone_init(:,:,:,idays_init) = ozone
      enddo
    endif
    call get_time ( Time, seconds, days )
    t_days_ind = mod(days,365)+1
    ozone = ozone_init(:,:,:,t_days_ind)

!   replace ozone in tracer 

    call get_number_tracers( MODEL_ATMOS, num_tracers=num_tracers )
    n_ozone = get_tracer_index( MODEL_ATMOS, 'ozone' )

    if ( num_tracers == size(rdt,4) ) then
      rm(:,:,:,n_ozone) = ozone
    else
      call error_mesg( 'ct07_ozone','size(rdt,4) not equal to num_tracers', FATAL )
    endif

!   discard ozone above 0.1 hPa of the model

    do k = 1, size(rm,3)
      if ( p_full(1,1,k) < 10 ) then
        ozone(:,:,k) = 0.
      endif
    enddo

!   calculate cumulative ozone at each level

    call calc_column ( p_half, p_full, ozone, r_column, to3 )

!   output ozone

    if (id_o3 > 0) used = send_data ( id_o3, ozone, Time, is, js )
  
    if ( total_ozone_output ) then
      if ( id_to3 > 0 ) used = send_data ( id_to3, to3, Time, is, js )
    endif

! === use parametrized ozone ===

  else if ( trim(ozone_option) == 'ct07_ozone' ) then

!   get coefficient for ozone scheme

    if ( count_int == 1 ) then
      do idays_init = 1, days_of_year
        call Time_for_hs( Time2, idays_init )
        call get_ct07_coeff( Time2, p_half, p_full, a1, a2, a3, a4, a5, a6, a7, a8, a9 )
        ct07_coeff(:,:,:,idays_init,1) = a1
        ct07_coeff(:,:,:,idays_init,2) = a2
        ct07_coeff(:,:,:,idays_init,3) = a3
        ct07_coeff(:,:,:,idays_init,4) = a4
        ct07_coeff(:,:,:,idays_init,5) = a5
        ct07_coeff(:,:,:,idays_init,6) = a6
        ct07_coeff(:,:,:,idays_init,7) = a7
        ct07_coeff(:,:,:,idays_init,8) = a8
        ct07_coeff(:,:,:,idays_init,9) = a9
      enddo
    endif
    call get_time ( Time, seconds, days )
    t_days_ind = mod(days,365)+1
    a1 = ct07_coeff(:,:,:,t_days_ind,1)
    a2 = ct07_coeff(:,:,:,t_days_ind,2) 
    a3 = ct07_coeff(:,:,:,t_days_ind,3)
    a4 = ct07_coeff(:,:,:,t_days_ind,4)
    a5 = ct07_coeff(:,:,:,t_days_ind,5)
    a6 = ct07_coeff(:,:,:,t_days_ind,6)
    a7 = ct07_coeff(:,:,:,t_days_ind,7)
    a8 = ct07_coeff(:,:,:,t_days_ind,8)
    a9 = ct07_coeff(:,:,:,t_days_ind,9)

!   read ozone from tracer field
    
    call get_number_tracers( MODEL_ATMOS, num_tracers=num_tracers )
    n_ozone       = get_tracer_index( MODEL_ATMOS, 'ozone' )
    n_cold_tracer = get_tracer_index( MODEL_ATMOS, 'cold_tracer' )

    if ( num_tracers == size(rdt,4) ) then 
      rst         = rm(:,:,:,n_ozone)
      cold_tracer = rm(:,:,:,n_cold_tracer)
      tst         = tm
    else
      call error_mesg( 'ct07_ozone','size(rdt,4) not equal to num_tracers', FATAL )
    endif

    where ( rst < 0. ) rst = 0.  ! ozone should be >= 0.

!   discard ozone above 0.1 hPa of the model

    do k = 1, size(rm,3)
      if ( p_full(1,1,k) < 10 ) then
        rst(:,:,k) = 0.
      endif
    enddo

!   calculate cumulative ozone at each level

    call calc_column ( p_half, p_full, rst, r_column, to3 )
    r_column = r_column*du2mol  ! [DU] to [molecules cm-2]

!   calclate anthropogenic forcing

    rtnd_anthropogenic = 0.

    if ( anthropogenic_forcing ) then
      call zenith_angle ( Time, lon, lat, cos_theta )
!     (a) T < 195K in sunlit area
      if ( .not. cold_tracer_forcing ) then 
        do i = 1, size(rm,1)
          do j = 1, size(rm,2)
            do k = 1, size(rm,3)
              if ( tst(i,j,k) <= 195. .and. cos_theta(i,j) > 0. ) then
                rtnd_anthropogenic(i,j,k) = a8(i,j,k)*rst(i,j,k)
              endif
            enddo
          enddo
        enddo
!     (b) Cold tracer in sunlit area
      else 
        where ( cold_tracer > 1. ) 
          cold_tracer = 1.
        end where
        where ( cold_tracer < 0. )
          cold_tracer = 0.
        end where
        do i = 1, size(rm,1)
          do j = 1, size(rm,2)
            do k = 1, size(rm,3)
              if ( cos_theta(i,j) > 0. ) then
                rtnd_anthropogenic(i,j,k) = (a8(i,j,k)*cold_tracer(i,j,k)*(195./tst(i,j,k))**4.5)*rst(i,j,k)
              endif
            enddo
          enddo
        enddo
      endif
    endif

!   calculate ozone tendency
 
    do3        = rst-a3
    dtemp      = tst-a5
    dto3       = r_column-a7 
    rtnd_do3   = a2*do3
    rtnd_dtemp = a4*dtemp
    rtnd_dto3  = a6*dto3

    rtnd = a1 + rtnd_do3 + rtnd_dtemp + rtnd_dto3 + rtnd_anthropogenic

      !discard ozone tendency above 0.1 hPa of the model
      do k = 1, size(rm,3)
        if ( p_full(1,1,k) < 10 ) then
          rtnd(:,:,k) = 0.
        endif
      enddo

      ! correction for possible negative ozone in the future
      rst_future = rst+dt*rtnd
      where ( rst_future < 0. )
        rtnd = -rst/dt + 1.e-8/dt
      end where

!   update ozone tendency

    rdt(:,:,:,n_ozone) = rdt(:,:,:,n_ozone) + rtnd

!   update cold tracer tendency

    cold_tracertnd = 0.
    if ( cold_tracer_forcing ) then
      call cold_tracer_tendency ( cold_tracer, tst, a9, cold_tracertnd )  
      rdt(:,:,:,n_cold_tracer) = rdt(:,:,:,n_cold_tracer) + cold_tracertnd
    endif

!   ---- set output ----

    if ( id_o3 > 0) used = send_data ( id_o3, rst, Time, is, js )

    if ( total_ozone_output ) then
      if ( id_to3 > 0 ) used = send_data ( id_to3, to3, Time, is, js )
    endif

    if ( ozone_tendency_output ) then
      if ( id_do3dt > 0 )              used = send_data ( id_do3dt, rtnd, Time, is, js )
      if ( id_rtnd_do3 > 0 )           used = send_data ( id_rtnd_do3, rtnd_do3, Time, is, js )
      if ( id_rtnd_dtemp > 0 )         used = send_data ( id_rtnd_dtemp, rtnd_dtemp, Time, is, js )
      if ( id_rtnd_dto3 > 0 )          used = send_data ( id_rtnd_dto3, rtnd_dto3, Time, is, js )
      if ( id_rtnd_anthropogenic > 0 ) used = send_data ( id_rtnd_anthropogenic, rtnd_anthropogenic, Time, is, js )  
    endif

    if ( cold_tracer_output ) then
      if ( id_cold_tracer > 0 ) used = send_data ( id_cold_tracer, cold_tracer, Time, is, js )
    endif

  endif

 end subroutine ct07_ozone

!#######################################################################

 subroutine get_ozone ( Time, p_half, ozone )

    type(time_type), intent(in)            ::  Time
    real, intent(in),    dimension(:,:,:)  ::  p_half
    real, intent(inout), dimension(:,:,:)  ::  ozone

!-----------------------------------------------------------------------

    call interpolator (ozone_interp, Time, p_half, ozone, trim(ozone_file))

 end subroutine get_ozone

!#######################################################################

 subroutine get_ct07_coeff ( Time, p_half, p_full, a1, a2, a3, a4, a5, a6, a7, a8, a9 )

!-----------------------------------------------------------------------

    type(time_type), intent(in)            ::  Time
    real, intent(in),    dimension(:,:,:)  ::  p_half, p_full
    real, intent(inout), dimension(:,:,:)  ::  a1, a2, a3, a4, a5, a6, a7, a8, a9

!-----------------------------------------------------------------------

    real :: du2mol = 2.69e+16  ! [DU] = 2.69e+16 [molecules cm-2]

!-----------------------------------------------------------------------

    call interpolator (a1_interp, Time, p_half, a1, 'A1')
    call interpolator (a2_interp, Time, p_half, a2, 'A2')
    call interpolator (a3_interp, Time, p_half, a3, 'A3')
    call interpolator (a4_interp, Time, p_half, a4, 'A4')
    call interpolator (a5_interp, Time, p_half, a5, 'A5')
    call interpolator (a6_interp, Time, p_half, a6, 'A6')
    !call interpolator (a7_interp, Time, p_half, a7, 'A7')
    call interpolator (a8_interp, Time, p_half, a8, 'A8')
    call interpolator (a9_interp, Time, p_half, a9, 'A9')

!   ---- calculate a7 from a3 ----

    call calc_column ( p_half, p_full, a3, a7 )
    a7 = a7*du2mol  ! [DU] to [molecules cm-2]

 end subroutine get_ct07_coeff

!#######################################################################

 subroutine calc_column ( p_half, p_full, r, r_column, to3 )

!-----------------------------------------------------------------------

!         --------- p_half(1)
!
!  ro3(1) ========= p_full(1)
!                             dp_above = p_half(2) - p_full(1)    
!         --------- p_half(2) 
!                             dp_below = p_full(2) - p_half(2)
!  ro3(2) ========= p_full(2) 
!
!  r_column(2) = r_column(1) + ro3(1)*dp_above/g + ro3(2)*dp_below/g
    
!-----------------------------------------------------------------------

    real, intent(in),    dimension(:,:,:) :: p_half, p_full, r

    real, intent(inout), dimension(:,:,:) :: r_column

    real, intent(inout), dimension(:,:), optional :: to3

!-----------------------------------------------------------------------

    real, dimension(size(r,1),size(r,2)) :: dp_above, dp_below

    real :: rv2rm  = 48./28.97           ! O3 volume mixing ratio to mass mixing ratio
    real :: pstp   = 1.E5                ! [Pa], Note that p_full/p_half are in [Pa]
    real :: tstp   = 273.15              ! [K]  
    real :: Ro3    = 8.31*(1.e+3)/48.    ! specific gas constant of ozone [J K-1 kg-1]
    real :: g      = 9.81                ! [m s-2]
    real :: density_o3

    integer :: k

!-----------------------------------------------------------------------

!   ---- Calculate total ozone [DU] ----
     
    density_o3 = pstp/Ro3/tstp   ! [kg m-3]

    do k = 1, size(r,3)
      if (k == 1) then
        dp_below(:,:)   = p_full(:,:,k) - p_half(:,:,k)
        r_column(:,:,k) = (r(:,:,k)*rv2rm)*dp_below(:,:)/g
      else
        dp_above(:,:)   = p_half(:,:,k) - p_full(:,:,k-1)
        dp_below(:,:)   = p_full(:,:,k) - p_half(:,:,k)
        r_column(:,:,k) = r_column(:,:,k-1) + (r(:,:,k-1)*rv2rm)*dp_above(:,:)/g &
                                            + (r(:,:,k)*rv2rm)*dp_below(:,:)/g
      endif
    enddo

    r_column = (r_column*1.e+5/density_o3)

    if (present(to3)) then 
      to3(:,:) = r_column(:,:,size(r,3))
    endif

 end subroutine calc_column

!#######################################################################

 subroutine cold_tracer_tendency ( cold_tracer, tst, a9, cold_tracertnd )

!-----------------------------------------------------------------------

   real, intent(in),  dimension(:,:,:) :: cold_tracer, tst, a9
   real, intent(out), dimension(:,:,:) :: cold_tracertnd 

!-----------------------------------------------------------------------

   integer :: i, j, k

   real :: tau1 = 1./(6.*60.*60.)  ! 6 hours, a typical time scale for chlorine activation

!-----------------------------------------------------------------------

   do i = 1, size(cold_tracer,1)
     do j = 1, size(cold_tracer,2)
       do k = 1, size(cold_tracer,3)
         if ( tst(i,j,k) > 195. ) then
           cold_tracertnd(i,j,k) = a9(i,j,k)*cold_tracer(i,j,k)
         else
           cold_tracertnd(i,j,k) = tau1*(1.-cold_tracer(i,j,k))
         endif
       enddo
     enddo
   enddo

 end subroutine cold_tracer_tendency

!#######################################################################

end module ct07_ozone_mod
