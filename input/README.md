# Ressources

### Input data from DIVA and CIAM
The exact usage of DIVA and CIAM data in FRISIA is reported in Ramme et al. (in preparation).
Here we upload pre-processed data that are based on CIAM and DIVA data, which were taken from the [CIAM GitHub repository](https://github.com/raddleverse/MimiCIAM.jl) (Wong et al., 2022). Original DIVA data are coming from Vafeides et al. (2008) and Hinkel and Klein (2009).

The fingerprint file for calculating local SLR is taken from the CIAM GitHub repository based on Slangen et al. (2014).
The coastal GDP and population file is a pre-processed file based on IIASA SSP data (Riahi et al., 2017; Cuaresma, 2017; Samir and Lutz, 2017).

### Temperature and ocean heat content data
This repository contains a subset of the FaIR v2.1.1 data (calibration 1.1.0) (Smith et al., 2018; Leach et al., 2021) used to force FRISIAv1.0 in the analysis of Ramme et al. (in preparation).


### CO2 emissions data
These are global CO2 emissions data, aggregated from files used to force Earth System Models in CMIP6. These data include emissions from fossil burning and cement production, but not those from land-use changes.
CO2 emissions are used for calculating the expected SLR in FRISIA.

Historical gridded data are from Hoesly et al. (2018). The description of CO2 emissions in future scenarios for CMIP6 is described in Rogelj et al. (2018) and Gidden et al. (2019).
The extensions of emissions data beyond 2100 are based on Meinshausen et al. (2020).


#### References

Cuaresma, J. C.: Income projections for climate change research: A framework based on human capital dynamics, Global Environmental Change, 42, 226–236, https://doi.org/10.1016/j.gloenvcha.2015.02.012, 2017

Gidden, M. J., Riahi, K., Smith, S. J., Fujimori, S., Luderer, G., Kriegler, E., van Vuuren, D. P., van den Berg, M., Feng, L., Klein, D., Calvin, K., Doelman, J. C., Frank, S., Fricko, O., Harmsen, M., Hasegawa, T., Havlik, P., Hilaire, J., Hoesly, R., Horing, J., Popp, A., Stehfest, E., and Takahashi, K.: Global emissions pathways under different socioeconomic scenarios for use in CMIP6: a dataset of harmonized emissions trajectories through the end of the century, Geosci. Model Dev. Dev., 12, 1443-1475, https://doi.org/10.5194/gmd-12-1443-2019, 2019

Hinkel, J. and Klein, R. J.: Integrating knowledge to assess coastal vulnerability to sea-level rise: The development of the DIVA tool, Global Environmental Change, 19, 384–395, https://doi.org/10.1016/j.gloenvcha.2009.03.002, 2009

Hoesly, R. M., Smith, S. J., Feng, L., Klimont, Z., Janssens-Maenhout, G., Pitkanen, T., Seibert, J. J., Vu, L., Andres, R. J., Bolt, R. M., Bond, T. C., Dawidowski, L., Kholod, N., Kurokawa, J.-I., Li, M., Liu, L., Lu, Z., Moura, M. C. P., O'Rourke, P. R., and Zhang, Q.: Historical (1750–2014) anthropogenic emissions of reactive gases and aerosols from the Community Emissions Data System (CEDS), Geosci. Model Dev., 11, 369–408, https://doi.org/10.5194/gmd-11-369-2018, 2018

Leach, N. J., Jenkins, S., Nicholls, Z., Smith, C. J., Lynch, J., Cain, M., Walsh, T., Wu, B., Tsutsui, J., and Allen, M. R.: FaIRv2.0.0: A generalized impulse response model for climate uncertainty and future scenario exploration, Geoscientific Model Development, 14, 3007–3036, https://doi.org/10.5194/gmd-14-3007-2021, 2021

Meinshausen, M., Nicholls, Z. R. J., Lewis, J., Gidden, M. J., Vogel, E., Freund, M., Beyerle, U., Gessner, C., Nauels, A., Bauer, N., Canadell, J. G., Daniel, J. S., John, A., Krummel, P. B., Luderer, G., Meinshausen, N., Montzka, S. A., Rayner, P. J., Reimann, S., Smith, S. J., van den Berg, M., Velders, G. J. M., Vollmer, M. K., and Wang, R. H. J.: The shared socio-economic pathway (SSP) greenhouse gas concentrations and their extensions to 2500, Geosci. Model Dev., 13, 3571–3605, https://doi.org/10.5194/gmd-13-3571-2020, 2020

Riahi, K., van Vuuren, D. P., Kriegler, E., Edmonds, J., O’Neill, B. C., Fujimori, S., Bauer, N., Calvin, K., Dellink, R., Fricko, O., Lutz, W., Popp, A., Cuaresma, J. C., KC, S., Leimbach, M., Jiang, L., Kram, T., Rao, S., Emmerling, J., Ebi, K., Hasegawa, T., Havlik, P., Humpenöder, F., Da Silva, L. A., Smith, S., Stehfest, E., Bosetti, V., Eom, J., Gernaat, D., Masui, T., Rogelj, J., Strefler, J., Drouet, L., Krey, V., Luderer, G., Harmsen, M., Takahashi, K., Baumstark, L., Doelman, J. C., Kainuma, M., Klimont, Z., Marangoni, G., Lotze-Campen, H., Obersteiner, M., Tabeau, A., and Tavoni, M.: The Shared Socioeconomic Pathways and their energy, land use, and greenhouse gas emissions implications: An overview, Global Environmental Change, 42, 153–168, https://doi.org/10.1016/j.gloenvcha.2016.05.009, 2017

Rogelj, J., Popp, A., Calvin, K.V., Luderer, G., Emmerling, J., Gernaat, D., Fujimori, S., Strefler, J., Hasegawa, T., Marangoni, G., Krey, V., Kriegler, E., Riahi, K., van Vuuren, D.P., Doelman, J., Drouet, L., Edmonds, J., Fricko, O., Harmsen, M., Havlik, P., Humpenöder, F., Stehfest, E., Tavoni, M., Scenarios towards limiting global mean temperature increase below 1.5 °C. Nature Climate Change 8, 2018, 325-332. https://doi.org/10.1038/s41558-018-0091-3

Samir, K. and Lutz, W.: The human core of the shared socioeconomic pathways: Population scenarios by age, sex and level of education for all countries to 2100, Global Environmental Change, 42, 181–192, https://doi.org/10.1016/j.gloenvcha.2014.06.004, 2017

Slangen, A. B. A., Carson, M., Katsman, C. A., van de Wal, R. S. W., Köhl, A., Vermeersen, L. L. A., & Stammer, D. (2014). Projecting twenty-first century regional sea-level changes. Climatic Change, 124(1–2), 317–332. https://doi.org/10.1007/s10584-014-1080-9

Smith, C. J., Forster, P. M., Allen, M., Leach, N., Millar, R. J., Passerello, G. A., and Regayre, L. A.: FAIR v1.3: a simple emissions-based impulse response and carbon cycle model, Geoscientific Model Development, 11, 2273–2297, https://doi.org/10.5194/gmd-11-2273-2018, 2018

Vafeidis, A. T., Nicholls, R. J., McFadden, L., Tol, R. S. J., Hinkel, J., Spencer, T., Grashoff, P. S., Boot, G., and Klein, R. J. T.: A New Global Coastal Database for Impact and Vulnerability Analysis to Sea-Level Rise, Journal of Coastal Research, 24, 917–924, https://doi.org/10.2112/06-0725.1, 2008

Wong, T. E., Bakker, A. M. R., Ruckert, K., Applegate, P., Slangen, A. B. A., and Keller, K.: BRICK v0.2, a simple, accessible, and transparent model framework for climate and regional sea-level projections, Geoscientific Model Development, 10, 2741–2760, https://doi.org/10.5194/gmd-10-2741-2017, 2017
