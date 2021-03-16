# HR_SCDM

Simplified Chemistry-Dynamical Model

We extend the GFDL's idealized general circulation model (GCM) by incporating a linear ozone scheme and a parameterization method for the shortwave ozone heating in the model. We call the new model as simplified chemistry-dynamical model (SCDM), which allows more efficient simulations of ozone-circulation interaction with a much lower computational cost. Our model is based on the improving version of the idealized GCM proposed by Wu and Reichler (2018), who optimized the Newtonian relaxation term by introducing a zonally asymmetric equilibrium temperature (Teq) profile. We use their D7 model setup with surface drag being 1.35 1/day and with the Newtonian relaxation time scale (Tau) as in Jucker et al. (2014). Since the shortwave ozone effet is internally described by the SCDM, we update the Teq profile used in Wu and Reichler (2018) to avoid accounting additional shortwave ozone heating from the Newtonian term. Therefore, the overall diabatic heating in the SCMD becomes the sum of the Newtonian relaxztion term and the shortwave ozone heating. A detailed description for the Teq correction procedure and an initial evaluation of the SCDM ozone simulation are described in Hong and Reichler (2021). 

This code is based on WR_simpleGCM (https://github.com/ZhengWinnieWu/WR_simpleGCM.git). The following program and modules were modified:
• atmos_model
• hs_forcing_mod
• atmosphere_mod
• time_manager_mod

The new code reads in twelve monthly Teq and Tau profiles during initialization and then uses linear interpolation to calculate daily fields. The interpolation is done only once during initialization to reduce the actual run time. Additional modification concerns the inclusion of a more flexible surface drag. The Rayleigh damping (surface drag) in the original code is defined in module “hs_forcing_mod”, which is a fixed number. In our modification, we introduce the Rayleigh damping as a two-dimensional matrix (longitude-latitude), and allow the code to read in the Rayleigh damping data files. Certain new variables are read in through the name-list:
• equilibrium_option = 'from_file',
• drag_file_name = 'name of drag data file.txt',
• temp_file = 'name of Teq data file.nc',
• tau_file = 'name of Tau data file.nc'.

--- Main reference ---

Cariolle, D. and H. Teyssèdre (2007): A revised linear ozone photochemistry parameterization for use in transport and general circulation models: multi-annual simulations, Atmos. Chem. Phys., 7, 2183–2196, doi:10.5194/acp-7-2183-2007.

Held, I. M. and M. J. Suarez (1994): A proposal for the intercomparison of the dynamical cores of atmospheric general circulation models, Bull. Amer. Meteor. Soc., 75, 1825–1830.

Hong, H.-J. and T. Reichler (2021): A simplified chemistry-dynamical model, J. Adv. Model. Earth Sys. (in preparation)

Jucker, M., S. Fueglistaler, and G. K. Vallis (2014): Stratospheric sudden warmings in an idealized GCM, J. Geophys. Res. Atmos., 119, 11054–11064, doi:10.1002/2014JD022170.

Lacis, A. A. and J. Hansen (1974): A parameterization for the absorption of solar radiation in the earth’s atmosphere, J. Atmos. Sci., 31, 118–133, doi:10.1175/1520-0469(1974)031%3C0118:APFTAO%3E2.0.CO;2.

Wu, Z. and T. Reichler (2018): Towards a More Earth-like Circulation in Idealized Models, J. Adv. Model. Earth Sys., 30(24), 10101-10116, doi:10.1029/2018MS001356.
