#!/bin/csh
#SBATCH --time=48:00:00 
#SBATCH --nodes=3
#SBATCH --ntasks=32
#SBATCH --account=YOUR ACCOUNT
#SBATCH --partition=YOUR SERVER NAME
#SBATCH --mail-user=YOUR OWN MAIL
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH -o /YOUR OWN PTAH/test.o
#SBATCH -e /YOUR OWN PATH/test.e

set echo

  set name = FILENAME

  set npes       = 32
  set work       = /YOUR OWN PATH FOR WORK/${name}
  set root       = /YOUR OWN PATH FOR SCDM
  set executable = $root/fms_hs/compile/obj/wr_hs.x
  set archive    = /YOUR OWN PATH FOR ARCHIVE
  set outputDir  = $archive/$name
  set scriptName = $root/scripts/${name}
  set initCond   = ""   # cold start if "" or reload_commands if existent, or "cpio" name of a restart file

  # prescribed teq and tau files
  set equilibrium_option = 'from_file'
  set temp_file          = $root/data/teq_ite_drag1.35_cycle29.nc
  set tau_file           = $root/data/tau_ctr.nc

  # prescribed drag file
  set drag_file_name = $root/data/drag_1.35.txt
  
  # ozone setup
  set no_ozone         = .false.
  set ozone_option     = 'ct07_ozone'  # option: 'ct07_ozone' or 'from_file'
  set ozone_file       = ''            # given ozone file name if ozone_option = 'from_file'
  set no_ozone_forcing = .false.
  set ozone_sw_option  = 'ls74_sw'
  set albedo_file      = $root/ct07_coeff/INPUT_MERRA2_ALBEDO_1980_2018.nc

  # list of days for each continous integration
  set dayslist = ( 365 365 365 365 365 365 365 365 365 365 365 365 365 365 365 )

  alias time_stamp $root/time_stamp.csh

########################################################################
#---------- only highly premeditated user changes below here -----------
########################################################################

  set runsPerScript = $#dayslist
  echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
  echo $name
  echo ${npes}pe_${runsPerScript}run
  echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

  set irun = 1
  set ireload = 1
  set reload_file = $outputDir/ascii/reload_commands
  if ( -d $outputDir )  then
    if ( -f $reload_file ) then
      if ( -r $reload_file ) then
        echo "**************************************************"
        echo "RESTART from reload_file" $reload_file
        echo "**************************************************"
        source $reload_file
      else
        unset echo
        echo "ERROR: reload file is not readable: $reload_file"
        exit 1
      endif
    endif
  else
    mkdir -p $outputDir
  endif

  if ( ! -d $work ) mkdir -p $work 
  unset echo
  if ( ! -e $work ) then
    echo "ERROR: working directory could not be created: $work"
    exit 1
  else if ( ! -d $work ) then
    echo "ERROR: $work exists, but is not a directory."
    exit 1
  else if ( ! -r $work ) then
    echo "ERROR: working directory is not readable: $work"
    exit 1
  else if ( ! -w $work ) then
    echo "ERROR: working directory is not writable: $work"
    exit 1
  endif
  set echo

##############################################################################
  cd $work
##############################################################################

  rm -rf * >& /dev/null			# rm *res
  if ( ! -d INPUT ) mkdir INPUT 
  mkdir RESTART

  unset echo
  if ( ! -e INPUT ) then
    echo "ERROR: input directory could not be created: $work/INPUT/"
    exit 1
  else if ( ! -d INPUT ) then
    echo "ERROR: $work/INPUT/ exists, but is not a directory."
    exit 1
  else if ( ! -r INPUT ) then
    echo "ERROR: input directory is not readable: $work/INPUT/"
    exit 1
  else if ( ! -w INPUT ) then
    echo "ERROR: input directory is not writable: $work/INPUT/"
    exit 1
  endif
  if ( ! -e RESTART ) then
    echo "ERROR: restart directory could not be created: $work/RESTART/"
    exit 1
  else if ( ! -d RESTART ) then
    echo "ERROR: $work/RESTART/ exists, but is not a directory."
    exit 1
  else if ( ! -r RESTART ) then
    echo "ERROR: restart directory is not readable: $work/RESTART/"
    exit 1
  else if ( ! -w RESTART ) then
    echo "ERROR: restart directory is not writable: $work/RESTART/"
    exit 1
  endif
  set echo

##############################################################################
  cd INPUT
##############################################################################

  if ( "$initCond" != "" ) then
    cpio -iv  < $initCond
    #if (-e dynamics.res) mv dynamics.res bgrid_prog_var.res	# tjr
  endif

  # ct07_ozone coefficient    # hjh
  ln -s $root/ct07_coeff/A1.nc                     A1.nc
  ln -s $root/ct07_coeff/A2.nc                     A2.nc
  ln -s $root/ct07_coeff/A3_MERRA2_O3_1980_2018.nc A3.nc  # 3-D ozone climatology 
  ln -s $root/ct07_coeff/A4.nc                     A4.nc
  ln -s $root/ct07_coeff/A5_MERRA2_T_1980_2018.nc  A5.nc  # 3-D temperature climatology
  ln -s $root/ct07_coeff/A6.nc                     A6.nc
  ln -s $root/ct07_coeff/A7.nc                     A7.nc
  ln -s $root/ct07_coeff/A8.nc                     A8.nc
  ln -s $root/ct07_coeff/A9.nc                     A9.nc

  # ground albedo   # hjh
  ln -s $albedo_file ALBEDO.nc 

  # input ozone file   # hjh
  ln -s $ozone_file O3.nc

  ln -s $tau_file  tau.nc
  ln -s $temp_file temp.nc

##############################################################################
  cd $work
##############################################################################

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# DIAG TABLE
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#"file_name", output_freq, "output_units", format, "time_units", "long_name",
#
#output_freq:  > 0  output frequency in "output_units"
#              = 0  output frequency every time step
#              =-1  output frequency at end of run
#
#output_units = units used for output frequency
#               (years, months, days, minutes, hours, seconds)
#
#time_units   = units used to label the time axis
#               (days, minutes, hours, seconds)
#
#  FORMAT FOR FIELD ENTRIES (not all input values are used)
#
#"module_name", "field_name", "output_name", "file_name" "time_sampling", time_avg, "other_opts", packing
#
#time_avg = .true. or .false.
#
#packing  = 1  double precision
#         = 2  float
#         = 4  packed 16-bit integers
#         = 8  packed 1-byte (not tested?)

  cat >> diag_table <<EOF
$name
0 0 0 0 0 0
#output files
"atmos_daily",    24, "hours", 1, "days", "time",
"atmos_monthly", 720, "hours", 1, "days", "time",
"atmos_average",  -1, "hours", 1, "days", "time",
#diagnostic field entries.
  "dynamics",      "bk",                  "bk",                  "atmos_daily",    "all", .true., "none", 2,
  "dynamics",      "pk",                  "pk",                  "atmos_daily",    "all", .true., "none", 2,
  "dynamics",      "ps",                  "ps",                  "atmos_daily",    "all", .true., "none", 2,
  "dynamics",      "height",              "height",              "atmos_daily",    "all", .true., "none", 2,
  "dynamics",      "zsurf",               "zsurf",               "atmos_daily",    "all", .true., "none", 2,
  "dynamics",      "ucomp",               "ucomp",               "atmos_daily",    "all", .true., "none", 2,
  "dynamics",      "vcomp",               "vcomp",               "atmos_daily",    "all", .true., "none", 2,
  "dynamics",      "temp",                "temp",                "atmos_daily",    "all", .true., "none", 2,
  "dynamics",      "omega",               "omega",               "atmos_daily",    "all", .true., "none", 2,
  "ct07_ozone",    "o3",                  "o3",                  "atmos_daily",    "all", .true., "none", 2,   
  "ct07_ozone",    "to3",                 "to3",                 "atmos_daily",    "all", .true., "none", 2,   
  "ct07_ozone",    "do3dt",               "do3dt",               "atmos_daily",    "all", .true., "none", 2,   
  "ct07_ozone",    "rtnd_do3",            "rtnd_do3",            "atmos_daily",    "all", .true., "none", 2,   
  "ct07_ozone",    "rtnd_dtemp",          "rtnd_dtemp",          "atmos_daily",    "all", .true., "none", 2,   
  "ct07_ozone",    "rtnd_dto3",           "rtnd_dto3",           "atmos_daily",    "all", .true., "none", 2,   
  "ct07_ozone",    "rtnd_anthropogenic",  "rtnd_anthropogenic",  "atmos_daily",    "all", .true., "none", 2,   
  "ct07_ozone",    "cold_tracer",         "cold_tracer",         "atmos_daily",    "all", .true., "none", 2,
  "ozone_forcing", "tdt_sw",              "tdt_sw",              "atmos_daily",    "all", .true., "none", 2, 
  "ozone_forcing", "absorb_sw",           "absorb_sw",           "atmos_daily",    "all", .true., "none", 2, 
  "ozone_forcing", "costheta",            "costheta",            "atmos_daily",    "all", .true., "none", 2,
  "ozone_forcing", "s",                   "s",                   "atmos_daily",    "all", .true., "none", 2, 
EOF

# no tracers, however it may be that still one tracer is written, tjr 05/05/06
# "dynamics",    "tracer1",        "tracer1",        "atmos_average",  "all", .true.,  "none", 2,
# "dynamics",    "tracer2",        "tracer2",        "atmos_average",  "all", .true.,  "none", 2,

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FIELD TABLE
# if omitted no tracers are written (tjr)
# at least one field table must be supplied - else crash at the end of run - tjr 05/05/06
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  cat >> field_table <<EOF
#"TRACER", "atmos_mod",                "sphum"
#          "holefilling",              "off"
#          "numerical_representation", "grid"
#         "advect_vert",              "finite_volume_parabolic"/
#"TRACER", "atmos_mod",                "age_of_air"
#          "holefilling",              "off"
#          "numerical_representation", "grid"
#          "advect_vert",              "van_leer_linear"/
#"TRACER", "atmos_mod",                "pv"
#          "numerical_representation", "grid"
#          "advect_vert",              "van_leer_linear"/
#"TRACER", "atmos_mod",                "apv"
#          "numerical_representation", "grid"/
#"TRACER", "atmos_mod",                "methane"
#          "holefilling",              "off"
#          "numerical_representation", "grid"
#          "advect_vert",              "van_leer_linear"/
"TRACER", "atmos_mod",                "ozone"
          "holefilling",              "off"
          "numerical_representation", "grid"
          "advect_vert",              "van_leer_linear"/
"TRACER", "atmos_mod",                "cold_tracer"
          "holefilling",              "off"
          "numerical_representation", "grid"
          "advect_vert",              "van_leer_linear"/
#"TRACER", "atmos_mod",                "tracer1"
#          "holefilling",              "off"
#          "numerical_representation", "grid"
#          "advect_vert",              "van_leer_linear"/
EOF

#########################################################################
#########################################################################
  while ($irun <= $runsPerScript)
#########################################################################
#########################################################################

    set days = $dayslist[$irun]

    cat > input.nml <<EOF
 &diag_manager_nml
    mix_snapshot_average_fields=.false.,
    debug_diag_manager=.true.
/

 &fms_io_nml
    threading_write = 'single',
    fileset_write = 'single'
/

 &fms_nml
         print_memory_usage=.true.,
         domains_stack_size = 400000
/

 &hs_forcing_nml
    no_forcing = .false.,
    trsink=0.,
    surface_forcing_input = .false.,
    do_conserve_energy    = .true.,
    equilibrium_tau_option = '$equilibrium_option',
    equilibrium_t_option   = '$equilibrium_option',
    pv_sat_flag   =  .false.,
    sat_only_flag =  .false.,
    pv_gamma = 4.e-3
    pv_phi0 = 50.,
    pv_dphi = -10.,
    sc_flag = .false.,
    sponge_flag =  .true.,
    sponge_pbottom = 5.e1,
    sponge_tau_days = 0.5,
    sigma_strat1 = 0.15,
    sigma_strat2 = 0.095,
    drag_file = '$drag_file_name',
/

 &ct07_ozone_nml
    no_ozone = $no_ozone,
    ozone_option = '$ozone_option',
    ozone_file = 'O3',
    anthropogenic_forcing = .true.,
    cold_tracer_forcing = .true.,
    cold_tracer_output = .true.,
    ozone_tendency_output = .true.,
    total_ozone_output = .true.,
/

 &ozone_forcing_nml
    no_ozone = $no_ozone,
    no_ozone_forcing = $no_ozone_forcing,
    ozone_sw_option = '$ozone_sw_option',
    ozone_sw_ouput = .true.,
    albedo_file = 'ALBEDO',
/

 &main_nml
         days   = $days,
         dt_atmos = 800
/

 &spectral_dynamics_nml
    damping_option          = 'resolution_dependent',
    damping_order           = 4,
    damping_coeff           = 1.15741e-4,
    eddy_sponge_coeff       = 0.00000e+6,
    zmu_sponge_coeff        = 0.00000e+5,
    zmv_sponge_coeff        = 0.00000e+5,
    do_mass_correction      =.true.,
    do_energy_correction    =.true.,
    do_water_correction     =.false.,
    use_virtual_temperature =.false.,
    vert_advect_uv          = 'second_centered',
    vert_advect_t           = 'second_centered',
    longitude_origin        = 0.,
    robert_coeff            = .04,
    alpha_implicit          = .5,
    reference_sea_level_press=1.e5,
    lon_max                 = 128,
    lat_max                 = 64,
    num_levels              = 40,
    num_fourier             = 42,
    num_spherical           = 43,
    fourier_inc             = 1,
    triang_trunc            =.true.,
    vert_coord_option       = 'input'
    scale_heights           = 11.0,
    surf_res                = 0.5,
    exponent                = 3.0
    ocean_topog_smoothing   = 0.0,
/

 &spectral_init_cond_nml
        topography_option = 'interpolated',
        topog_file_name = '/YOUR OWN PATH FOR SCDM/data/navy_topography.data.nc',
        topog_field_name = 'zdat'
/

 &topography_nml
        topog_file = '/YOUR OWN PATH FOR SCDM/data/navy_topography.data.nc',
/

 &gaussian_topog_nml
    height(1) = 4000,
    olon(1) = 0, 
    olat(1) = 45,
    wlon(1) = 50,
    wlat(1) = 20, 
    rlon(1) = 0, 
    rlat(1) = 0,
    height(2) = 4000,
    olon(2) = 180, 
    olat(2) = 45,
    wlon(2) = 50,
    wlat(2) = 20, 
    rlon(2) = 0, 
    rlat(2) = 0,
/

 &vert_coordinate_nml
        pk(1)=0,bk(1)=0,
        pk(2)=1.9201,bk(2)=0,
        pk(3)=4.7775,bk(3)=0,
        pk(4)=10.325,bk(4)=0,
        pk(5)=20.133,bk(5)=0,
        pk(6)=36.284,bk(6)=0,
        pk(7)=61.443,bk(7)=0,
        pk(8)=98.954,bk(8)=0,
        pk(9)=152.9,bk(9)=0,
        pk(10)=228.08,bk(10)=0,
        pk(11)=330.42,bk(11)=0,
        pk(12)=466.6,bk(12)=0,
        pk(13)=644.22,bk(13)=0,
        pk(14)=872.41,bk(14)=0,
        pk(15)=1161.2,bk(15)=0,
        pk(16)=1520.9,bk(16)=0,
        pk(17)=1965.7,bk(17)=0,
        pk(18)=2508.8,bk(18)=0,
        pk(19)=3166.4,bk(19)=0,
        pk(20)=3954.7,bk(20)=0,
        pk(21)=4892,bk(21)=0,
        pk(22)=5625.4,bk(22)=0.0037013,
        pk(23)=6387,bk(23)=0.009005,
        pk(24)=7163.2,bk(24)=0.016314,
        pk(25)=7933.7,bk(25)=0.0261,
        pk(26)=8665.8,bk(26)=0.038875,
        pk(27)=9328.2,bk(27)=0.055237,
        pk(28)=9894.4,bk(28)=0.07595,
        pk(29)=10310,bk(29)=0.10175,
        pk(30)=10519,bk(30)=0.13348,
        pk(31)=10468,bk(31)=0.17219,
        pk(32)=10085,bk(32)=0.21897,
        pk(33)=9286.4,bk(33)=0.27495,
        pk(34)=7988.8,bk(34)=0.34166,
        pk(35)=6085.8,bk(35)=0.42044,
        pk(36)=3464.7,bk(36)=0.51291,
        pk(37)=0,bk(37)=0.6209,
        pk(38)=0,bk(38)=0.7025,
        pk(39)=0,bk(39)=0.7925,
        pk(40)=0,bk(40)=0.8914,
        pk(41)=0,bk(41)=1
/

EOF

######################################################
#   --- run the model ---
######################################################

    echo loop_$irun/$runsPerScript
    echo run

    module load intel impi netcdf-c netcdf-f
    mpirun -np $npes $executable > fms.out

    if ( $status != 0 ) then
      wait
      set MPI_FAIL
      echo "ERROR: in mpirun, loop $irun" 
    endif

########################################################
  cd $work
########################################################

#   --- generate date for file names ---
    set begindate = `time_stamp -bf digital`				# 00010101
    if ( $begindate == "" ) set begindate = tmp`date '+%j%H%M%S'`
    set enddate = `time_stamp -ef digital`				# 00010201
    if ( $enddate == "" ) set enddate = tmp`date '+%j%H%M%S'`

    if ( -f time_stamp.out ) rm -f time_stamp.out

#   --- save ascii output files to local disk ---
    if ( ! -d $outputDir/ascii ) mkdir -p $outputDir/ascii
    foreach out (`ls *.out`)
      mv $out $outputDir/ascii/$begindate.$out
    end

########################################################
    cd RESTART
########################################################

    set resfiles = `ls *.res *.nc`
    if ( $#resfiles > 0 ) then
      set restart_file = $outputDir/restart/$enddate.cpio
      if ( ! -d $restart_file:h ) mkdir -p $restart_file:h
      unset echo
      if ( ! -e $restart_file:h ) then
	echo "ERROR: restart directory could not be created: $restart_file:h"
        exit 1
      else if ( ! -d $restart_file:h ) then
	echo "ERROR: $restart_file:h exists, but is not a directory."
        exit 1
      else if ( ! -w $restart_file:h ) then
	echo "ERROR: restart directory is not writable: $restart_file:h"
        exit 1
      endif
      set echo

      cp ../input.nml .
      cp ../*_table .
      set files = ( $resfiles input.nml *_table )
      ls $files | cpio -o > $restart_file

      if ( $irun < $runsPerScript ) then
         rm -f ../INPUT/*.res 
         mv -f *.res  ../INPUT
         mv -f *.nc  ../INPUT
      endif
    endif

########################################################
    cd $work
########################################################

#   --- save netcdf and data files ---
    if ( ! -d $outputDir/history ) mkdir -p $outputDir/history
    unset echo
    if ( ! -e $outputDir/history ) then
      echo "ERROR: history directory could not be"
      echo "       created: $outputDir/history"
      exit 1
    else if ( ! -d $outputDir/history ) then
      echo "ERROR: $outputDir/history exists, but is not a directory."
      exit 1
    else if ( ! -w $outputDir/history ) then
      echo "ERROR: history directory is not"
      echo "       writable: $outputDir/history"
      exit 1
    endif
    set echo

    foreach suffix ( nc data )
      foreach file ( `ls *.$suffix *.$suffix.????` ) 
         mv $file $begindate.$file
      end
      set files = `ls $begindate.*.$suffix $begindate.*.$suffix.????`
      if ( $#files > 0 ) then
         set arfile = $outputDir/history/$begindate.$suffix.0000.cpio
         if ( ! -e $arfile:h ) mkdir -p $arfile:h
         ls $files | cpio -o > $arfile
         rm -f $files
      endif
    end

#   --- terminate script if mpirun crashed ---
    if ( $?MPI_FAIL ) then
      unset echo
      exit 1
    endif

#   --- terminate script if there are no restart files ---
    if ( $#resfiles == 0 ) then
      unset echo
      echo "ERROR: no restart files exist, loop $irun" 
      exit 1
    endif

#   --- write new reload information ---
    @ irun++
    echo "set initCond     =  $restart_file"  > $reload_file

#########################################################################
#########################################################################
end # runsPerScript
#########################################################################
#########################################################################

echo end_of_run
echo "NOTE: Natural end-of-script for $scriptName."
exit 0
  
