
module ozone_forcing_mod

!-----------------------------------------------------------------------

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
use tracer_manager_mod, only: query_method, get_number_tracers
use   interpolator_mod, only: interpolate_type, interpolator_init, &
                              interpolator, interpolator_end, &
                              CONSTANT, INTERP_WEIGHTED_P, INTERP_LINEAR_P

use     ct07_ozone_mod, only: calc_column
use   zenith_angle_mod, only: zenith_angle, earth_sun_dis

implicit none
private

!-----------------------------------------------------------------------

  public :: ozone_forcing, ozone_forcing_init

   type(interpolate_type),save  ::  ground_albedo_interp

!------ namlist --------------------------------------------------------

!   --- ozone forcing option ---

   logical :: no_ozone = .false.
   logical :: no_ozone_forcing = .false.

   character(len=256) :: ozone_sw_option = 'ls74_sw'  
   logical            :: ozone_sw_ouput = .true.

   character(len=256) :: albedo_file = 'ALBEDO'

!-----------------------------------------------------------------------

   namelist /ozone_forcing_nml/ no_ozone,                          &
                                no_ozone_forcing,                  &
                                ozone_sw_option,                   &
                                ozone_sw_ouput,                    &
                                albedo_file

!-----------------------------------------------------------------------

   character(len=128) :: version='$Id: ozone_forcing.f90, 2019/03/23 hjh $'
   character(len=128) :: tagname='$Name: riga_O3_SW_control $'

   integer            :: id_tdt_sw
   integer            :: id_s
   integer            :: id_absorb
   integer            :: id_cos_theta 
   real               :: missing_value = -1.e10
   character(len=14)  :: mod_name = 'ozone_forcing'

   logical :: module_is_initialized = .false.


contains

!#######################################################################

 subroutine ozone_forcing_init ( axes, Time, lonb, latb )
    
            integer, intent(in) :: axes(4)
    type(time_type), intent(in) :: Time
    real, intent(in), optional, dimension(:,:) :: lonb, latb

!-----------------------------------------------------------------------

    integer :: unit, io, ierr

!   ---- read namelist ----

    #ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=ozone_forcing_nml, iostat=io)
        ierr = check_nml_error(io, 'ozone_forcing_nml')
    #else
        if (file_exist('input.nml')) then
           unit = open_namelist_file ( )
           ierr=1; do while (ierr /= 0)
              read  (unit, nml=ozone_forcing_nml, iostat=io, end=10)
              ierr = check_nml_error (io, 'ozone_forcing_nml')
           enddo
  10       call close_file (unit)
        endif
    #endif

    if ( trim(ozone_sw_option) == 'ls74_sw' ) then

!   ---- write version info and namelist to log file ----

       call write_version_number (version,tagname)
       if ( mpp_pe() == mpp_root_pe() ) write ( stdlog(), nml=ozone_forcing_nml )

       if ( no_ozone ) return
       if ( no_ozone_forcing ) return

    else
      call error_mesg( 'ozone_forcing_nml', trim(ozone_sw_option)//' is not valid for ozone_sw_option', FATAL)
    endif

!   ---- register diagnostic fields -----

    if ( ozone_sw_ouput ) then
         id_tdt_sw = register_diag_field ( mod_name, 'tdt_sw', axes(1:3), Time,         &
                                           'temperature tendency due to SW', 'deg/sec', &
                                           missing_value=missing_value                  )
         id_absorb = register_diag_field ( mod_name, 'absorb_sw', axes(1:3), Time, &
                                           'absorbed SW due to ozone', '%',        &
                                           missing_value=missing_value             )
      id_cos_theta = register_diag_field ( mod_name, 'costheta', axes(1:2), Time, &
                                           'cos(zenith_angle)', '1',              &
                                           missing_value=missing_value            )
              id_s = register_diag_field ( mod_name, 's', axes(1:2), Time, &
                                           'insolation at TOA', 'W/m^2',   &
                                           missing_value=missing_value     )
    endif

!   --- set ground albedo ---

    if ( trim(ozone_sw_option) == 'ls74_sw' ) then
      call interpolator_init ( ground_albedo_interp, trim(albedo_file)//'.nc', lonb, latb, data_names=(/trim(albedo_file)/), data_out_of_bounds=(/CONSTANT/))
    endif
 

    module_is_initialized = .true.

 end subroutine ozone_forcing_init

!#######################################################################

 subroutine ozone_forcing ( is, ie, js, je, dt, Time, lon, lat, p_half, p_full, &
                            tm, rm, tdt )

!-----------------------------------------------------------------------

 use  tracer_manager_mod, only: get_tracer_index, NO_TRACER

!-----------------------------------------------------------------------
 
   integer, intent(in)                        :: is, ie, js, je
      real, intent(in)                        :: dt
 type(time_type), intent(in)                  :: Time
      real, intent(in),    dimension(:,:    ) :: lon, lat
      real, intent(in),    dimension(:,:,:  ) :: p_half, p_full
      real, intent(in),    dimension(:,:,:  ) :: tm
      real, intent(in),    dimension(:,:,:,:) :: rm
      real, intent(inout), dimension(:,:,:  ) :: tdt

!-----------------------------------------------------------------------

   real, dimension(size(rm,1),size(rm,2),size(rm,3)) :: rst
   real, dimension(size(rm,1),size(rm,2),size(rm,3)) :: ttnd

   integer :: num_tracers
   integer :: n_ozone
   logical :: used

   real, dimension(size(rm,1),size(rm,2))            :: s
   real, dimension(size(rm,1),size(rm,2))            :: cos_theta  ! cos(zenith_angle)
   real, dimension(size(rm,1),size(rm,2),size(rm,3)) :: absorb

   integer :: i, j, k

!-----------------------------------------------------------------------

    if ( no_ozone ) return
    if ( no_ozone_forcing ) return

    if ( .not. module_is_initialized ) call error_mesg ('ozone_forcing','ozone_forcing_init has not been called', FATAL)

!-----------------------------------------------------------------------

!   --- ozone forcing for lacis & hansen (1974) sw calculation ---

    if ( ozone_sw_option == 'ls74_sw' ) then

!     extract ozone from tracer
      call get_number_tracers ( MODEL_ATMOS, num_tracers=num_tracers)
      n_ozone = get_tracer_index( MODEL_ATMOS, 'ozone' )
      if ( num_tracers == size(rm,4) ) then
        rst = rm(:,:,:,n_ozone)
      else
        call error_mesg( 'ozone_forcing','size(rm,4) not equal to num_tracers',FATAL )
      endif

      where ( rst < 0. ) rst = 0. ! ozone should be >= 0.

!     discard ozone above 0.1 hPa of the model
      do k = 1, size(rm,3) 
        if ( p_full(1,1,k) < 10 ) then
          rst(:,:,k) = 0.
        endif
      enddo

!     lacis & hansen (1974) ozone absorption
      call ls74_sw ( Time, lon, lat, p_half, p_full, tm, rst, ttnd, s, absorb, cos_theta )

!     discard temperature tendency above 0.7 hPa of the model
      do k = 1, size(rm,3) 
        if ( p_full(1,1,k) < 70 ) then
          ttnd(:,:,k) = 0.
        endif
      enddo

      tdt = tdt + ttnd

    endif  

! --- set output ---

    if ( ozone_sw_ouput ) then
      if ( id_tdt_sw > 0 )    used = send_data ( id_tdt_sw, ttnd, Time, is, js )
      if ( id_s > 0 )         used = send_data ( id_s, s, Time, is, js )
      if ( id_absorb > 0 )    used = send_data ( id_absorb, absorb, Time, is, js )
      if ( id_cos_theta > 0 ) used = send_data ( id_cos_theta, cos_theta, Time, is, js )
    endif

!-----------------------------------------------------------------------

 end subroutine ozone_forcing

!#######################################################################

 subroutine ls74_sw ( Time, lon, lat, p_half, p_full, tm, r, ttnd, s, absorb, cos_theta )

!-----------------------------------------------------------------------

    type(time_type), intent(in)         :: Time
    real, intent(in),  dimension(:,:  ) :: lon, lat
    real, intent(in),  dimension(:,:,:) :: p_half
    real, intent(in),  dimension(:,:,:) :: p_full
    real, intent(in),  dimension(:,:,:) :: tm
    real, intent(in),  dimension(:,:,:) :: r
    real, intent(out), dimension(:,:,:) :: ttnd

    real, intent(out), dimension(:,:)   :: s
    real, intent(out), dimension(:,:)   :: cos_theta  ! cos(zenith_angle)
    real, intent(out), dimension(:,:,:) :: absorb

!-----------------------------------------------------------------------

    real, dimension(size(r,1),size(r,2)) :: mu0
    real                                 :: es_dis     ! square of ratio of mean earth-Sun distance and earth-Sun distance

    real, dimension(size(r,1),size(r,2)) :: M          ! magnification factor M
    real                                 :: Mb         ! magnification factor Mb
    real, dimension(size(r,1),size(r,2)) :: Rb         ! albedo of the reflecting region
    real, dimension(size(r,1),size(r,2)) :: Rab        ! effective albedo of the lower atmosphere
    real, dimension(size(r,1),size(r,2)) :: Rg         ! ground reflectivity
    real                                 :: Rabs       ! spherical albedo of the reflecting region

    real, dimension(size(r,1),size(r,2),size(r,3))   :: r_column   ! overhead ozone at each level
    real, dimension(size(r,1),size(r,2),size(r,3)-1) :: ul, ull, ut
    real, dimension(size(r,1),size(r,2),size(r,3)-1) :: xl, xll, xls, xlls
    real, dimension(size(r,1),size(r,2))             :: rl, rll, rls, rlls
    real, dimension(size(r,1),size(r,2),size(r,3)-1) :: A      
    real, dimension(size(r,1),size(r,2))             :: dp
    real, dimension(size(r,1),size(r,2))             :: dair
    real, dimension(size(r,1),size(r,2),size(r,3)+1) :: z_half 
    real, dimension(size(r,1),size(r,2))             :: dz

    integer :: i, j, k

    real, parameter :: s0    = 1361. ! mean solar constant
    real, parameter :: g0    = 9.81  ! gravity
    real, parameter :: cp    = 1004. ! specific heat at constant pressure
    real, parameter :: dtime = 1.    ! sec
    real, parameter :: pi    = acos(-1.)    
    real, parameter :: p0    = 101325.        ! surface pressure
    real, parameter :: ha    = 287.*255./9.81 ! scale height


!-----------------------------------------------------------------------
  
!   --- calculate solar insolation ---

    call zenith_angle  ( Time, lon, lat, cos_theta )
    call earth_sun_dis ( Time, es_dis )
  
    s(:,:) = s0*es_dis
    where ( cos_theta < 0. ) s = 0. ! insolation = 0. after dark

!   --- calculate magnification factor for the slant path and refraction M
!       & magnification factor for diffuse upward radiation Mb             ---

    mu0 = cos_theta
    M  = 35./sqrt(1224.*(mu0**2.)+1.)
    Mb = 1.9

!   ---  calculate albedo of the reflecting region ---

!     albedo of the lower atmosphere due to Rayleigh scattering for clear skies
      Rab  = 0.219/(1.+0.816*mu0)
      Rabs = 0.144

!     ground surface reflectivity
      call get_ground_albedo ( Time, Rg )

!     atmospheric albedo for clear skies
      Rb = Rab + (1.-Rab)*(1.-Rabs)*Rg/(1.-Rabs*Rg)

!   --- calculate fraction of the total solar flux absorbed in each layer of the atmosphere ---

!     convert column ozone from [DU] to [cm]      
      call calc_column ( p_half, p_full, r, r_column )
      r_column = r_column*1.e-3 

!     total ozone amount [cm] above the lth & (l+1)th layer
      ul  (:,:,1:size(r,3)-1) = r_column(:,:,1:size(r,3)-1)
      ull (:,:,1:size(r,3)-1) = r_column(:,:,2:size(r,3))

!     total ozone amount above reflecting layer (ground for clear skies)
      do k = 1, size(r,3)-1
        ut (:,:,k) = r_column(:,:,size(r,3)) 
      enddo

!     ozone amount traversed by the direct solar beam in reaching the lth & (l+1)th layer
      do k = 1, size(r,3)-1
        xl  (:,:,k) =  ul(:,:,k)*M(:,:)
        xll (:,:,k) = ull(:,:,k)*M(:,:)
      enddo
!     ozone amount traversed by the diffuse radiation illuminating the lth & (l+1)th layer from below
      do k = 1, size(r,3)-1
        xls (:,:,k) = ut(:,:,k)*M(:,:) + Mb*(ut(:,:,k)- ul(:,:,k))
        xlls(:,:,k) = ut(:,:,k)*M(:,:) + Mb*(ut(:,:,k)-ull(:,:,k))
      enddo

!     calculate frequency-integrated ozone absorption & total absorption in the lth layer
      do k = 1, size(r,3)-1
        call ozone_sw_absorption ( xl(:,:,k)  , rl )
        call ozone_sw_absorption ( xll(:,:,k) , rll )
        call ozone_sw_absorption ( xls(:,:,k) , rls )
        call ozone_sw_absorption ( xlls(:,:,k), rlls )
        A(:,:,k) = mu0(:,:)*(rll(:,:)-rl(:,:) + Rb(:,:)*(rls(:,:)-rlls(:,:)))
        
      enddo

      do i = 1, size(r,1)
        do j = 1, size(r,2)
          if ( cos_theta(i,j) < 0 ) then  ! no absorption after dark
            A(i,j,:) = 0.                 ! also note that A should be > 0. since there is
          endif                           ! no shortwave cooling
        enddo
      enddo
      absorb = 0. 
      do k = 1, size(r,3)-1
        absorb(:,:,k) = A(:,:,k)
      enddo

!   --- calculate heating rate due to ozone absorption ---
    
    ttnd = 0.
    do k = 1, size(r,3)-1 
      dp(:,:) = p_half(:,:,k+1) - p_half(:,:,k) 
      ttnd(:,:,k) = (s(:,:)*g0*A(:,:,k))/(cp*dp(:,:)*dtime)
    enddo

 end subroutine ls74_sw

!#######################################################################

 subroutine ozone_sw_absorption ( x, ratio ) 

    real, intent(in),  dimension(:,:) :: x
    real, intent(out), dimension(:,:) :: ratio  ! fraction of absorption

    integer :: i, j 

!-----------------------------------------------------------------------

    do i = 1, size(x,1)
      do j = 1, size(x,2)
        ratio(i,j) = 0.02118*x(i,j)/(1.+0.042*x(i,j)+0.000323*(x(i,j)**2.)) + &  !Chappuis band (visible)
                     1.082*x(i,j)/((1.+138.6*x(i,j))**0.805) + 0.0658*x(i,j)/(1.+(103.6*x(i,j))**3.)  !ultraviolet
      enddo
    enddo

 end subroutine ozone_sw_absorption

!#######################################################################

 subroutine get_ground_albedo ( Time, Rg )

    type(time_type), intent(in)          ::  Time
    real, intent(inout), dimension(:,:)  ::  Rg

!-----------------------------------------------------------------------
 
    call interpolator ( ground_albedo_interp, Time, Rg, trim(albedo_file) )

 end subroutine get_ground_albedo

!#######################################################################

end module ozone_forcing_mod
