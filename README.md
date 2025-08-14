# FRISIA

FRISIA is the *Feedback-based knowledge Repository for Intergrated assessments of Sea level rise Impacts and Adaptation*. Its first published version 1.0 is fully described in [Ramme et al. (preprint)](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-1875/). FRISIAv1.0 is developed as part of the new integrated assessment model [FRIDA](https://github.com/metno/WorldTransFRIDA) within the [WorldTrans project](https://worldtrans-horizon.eu/) of the Horizon Europe research and innovation programs under grant agreement no. 101081661.

## Quick start
Setting up the enviroment using conda (e.g. [miniconda](https://docs.anaconda.com/miniconda/) or [miniforge](https://github.com/conda-forge/miniforge)). This requires a python version >= 3.8.3 (versions beyond 3.10.10 have not been tested).


```
conda env create -f environment.yml
conda activate frisia
```

The example script runs the global aggregation of FRISIA for the scenarios SSP1-1.9, SSP2-4.5 and SSP5-8.5, using 3 ensemble members from [FaIR](https://github.com/OMS-NetZero/FAIR) (17th, 50th and 83rd percentile) for each scenario.

```python
python EXAMPLE_runFRISIA.py
```

## User manual

### Running the model
Here you find instructions for a set of experiments one might want to conduct with FRISIA. You can also find simple examples for running FRISIA in this repository:

```
EXAMPLE_calculateSLR.py
EXAMPLE_runFRISIA.py
```

A more detailed analysis can be found in [this data repository](https://doi.org/10.5281/zenodo.15249065) on Zenodo.

#### Creating time series of sea level rise
The first step to running FRISIA is to generate time series of sea level rise. Here we do this using the example input from FaIR that come with the FRISIA reository. For more info on loading the input data have a look at EXAMPLE scripts.

```python
SLRModel = globalSLRModel(T, OHC_change, population=population, sy=1850, ey=2200)   # create the instance of the SLR model
SLRModel.integrate()                                                                # Integrate
SLRModel.align(2010)                                                                # Shift the SLR components to be 0 in 2010.
```

The inputs to the model should be one-dimensional arrays covering time from the start to the end year. Now, the total SLR and the individual components are available via

```python
SLR_total  = SLRModel.getSLRTotal()
SLR_thermo = SLRModel.getSLRThermo()
SLR_MG     = SLRModel.getSLRMG()
SLR_LWS    = SLRModel.getSLRLWS()
SLR_GIS    = SLRModel.getSLRGIS()
SLR_AIS    = SLRModel.getSLRAIS()
```

If we wanted to add the possibility of higher rates of sea level at certain levels of warming ("tipping behaviour"), we can create the instance of the SLRModel with the ```include_tipping=True``` option in the call to create the instance.
 
#### Running the impacts and adaptation model
The typical way to run the FRISIA impact and adaptation model is the following call:

```python
Model = SLRImpactModel(SLR, T, CO2emis, coastal_pop, coastal_gdp,
                       ey=2200, ndim=1, version='DIVA_global',
		       include_SLR_components=True, SLR_components=SLR_components)
Model.integrate()
```

The first line creates the instance of the impact model and the second line integrates the model from the start year 2010 until the user defined end year. 

Please take into account:

- the input timeseries for SLR, T and CO2emis have to be one-dimensional arrays of the global mean values that cover exactly the range over which the model integrates (e.g. here 191 timesteps from 2010 to 2200).
- SLR should be given as an anomaly since the start year, whereas T is the anomaly with respect to the pre-industrial mean (1850-1900). CO2emis are fossil emissions of CO2 in units of Gt C per year.
- the input timeseries for coastal population and coastal GDP have to be two-dimensional arrays, where the first dimension is the number of coastal regions (n=1 in this example) and the second dimension is time. These inputs are available from the input folder for the different default aggregation types coming with FRISIA. Have a look at the EXAMPLE script on how to load these.
- Currently allowed aggregation types are `DIVA_global` (ndim=1), `DIVA_bipolar` (ndim=2) and `DIVA_regional` (ndim=7) (see Ramme et al., 2025), with the respective inputs available here. If new aggregation types are added, these also require the pre-processed inputs to be generated (see corresponding section below).
- We have activated the `include_SLR_components` option, which means that instead of using the global mean sea level rise (SLR) to calculate coastal SLR (scaled linearly with region-specific pre-processed factors), we use the individual components of SLR and regionally averaged wheights for pre-processed input files to create local SLR for each coastal region, which is slightly more precise and allows for a larger range of uncertainty if the components vary between runs.

Output can be used from the model by either referring to the respective variable of the model directly

```python
Nflooded_byRegion = np.copy(SLRImpactModel.annual_people_flooded)
```

or by calling one of the designated output functions (e.g. for getting the sum over all coastal regions)

```python
Nflooded_total = SLRImpactModel.getTotalPeopleFlooded()
```

Currently available output functions can be found at the end of the model source code.

#### Varying the uncertainty parameters

So far we have used the default values for each parameter within FRISIA. However, most parameters have a range of possible values. To maximize on the capabilities of FRISIA to provide ranges of possible outputs, one would normally vary all the uncertainty parameters independently and create an ensemble of runs. The ensemble should then ideally also cover uncertainty over the inputs, most importantly for the temperature time series. In the example below T and OHC_change are now 2-dimensions arrays with the first dimension being time and the second dimension being the ensemble member:

```python
# creating the SLRModel instance once with activated randomize option is enough 
SLRModel = globalSLRModel(T[:,0], OHC_change[0,:], population=population, sy=1850, ey=2200, randomize=True)

for i in range(N):
    SLRModel.T_anomaly = T[:,i]
    SLRModel.OHC_change = OHC_change[:,i]
    SLRModel.reset_SLR()
    SLRModel.integrate()
    SLRModel.align(2010, silent=True)

    SLR = SLRModel.getSLRTotal()[160:] # use only values from 2010 onwards here
    SLR_components = [SLRModel.getSLRThermo()[160:], SLRModel.getSLRLWS()[160:], SLRModel.getSLRLMG()[160:],
                      SLRModel.getSLRGIS()[160:], SLRModel.getSLRAIS()[160:]]
	
    Model = SLRImpactModel(SLR, T[:,i], CO2emis, coastal_pop, coastal_gdp, ey=2200, ndim=1, version='DIVA_global',
		       include_SLR_components=True, SLR_components=SLR_components, randomize=True)
    Model.integrate()

    ...
```

For certain types of analysis it might make sense to use the same set of parametes for two model ensembles that will be compared afterwards. This is possible with specific functions to get and set the uncertainty parameters from the model. In the following example we once run a model ensemble without adaptation and then use the same parameters in a scenario with protection as adaptation strategy.

```python
for i in range(N):
    ...
    Model = SLRImpactModel(SLR, T[:,i], CO2emis, coastal_pop, coastal_gdp, ey=2200, ndim=1, version='DIVA_global',
		       include_SLR_components=True, SLR_components=SLR_components, randomize=True)
    Model.integrate()
    uncertainty_parameters = Model.getUncertaintyParameters()

    ...

    Model = SLRImpactModel(SLR, T[:,i], CO2emis, coastal_pop, coastal_gdp, ey=2200, ndim=1, version='DIVA_global',
		       include_SLR_components=True, SLR_components=SLR_components, randomize=True)

    Model.willingness_to_invest_in_fp = 1.0
    Model.setUncertaintyParameters(uncertainty_parameters)
    Model.integrate()
    
    ...
```

#### Running adaptation scenarios
As can be seen from the example above activating protection as adaptation strategy is possible by simply setting

```python
Model.willingness_to_invest_in_fp = 1.0
```

after creating the model instance, but before calling ```Model.integrate()```. Accordingly, retreat can be activated by setting

```python
Model.willingness_to_retreat = 1.0
```

In FRISIA version 1.0 we advise the user to not set both options at the same time, because this will lead to unrealistic costs as coastal regions would retreat and build protection at the same time. Future versions of the model are planned, which will allow to include time-dependent changes in adaptation and hybrid adaptation strategies.

#### Activate additional feedback

The default settings of FRISIA v1.0 are to turn off the additional feedback implemented in FRISIA compared to CIAM, to make model results comparable. Additional feedback can be activated by setting the following model switches after creating the model instance, but before integrating the model.

- ```Model.asset_feedback_switch = 1.0```: This will allow storm surge damages to reduce coastal asset values.
- ```Model.people_feedback_switch = 1.0```: This will allow storm surge fatalities to reduce coastal population.
- ```Model.include_reduced_growth = True```: This will include the feedback that coastal assets grow slower as there is an increase in expected exposure to storm surges under SLR.
- ```Model.include_fp_investment_cap = True```: This will limit the amount of money annually available for flood protection expanses (including maintenance costs) to a fraction of coastal GDP.
- ```Model.include_GDP_effect = True```: This will allow reductions in coastal assets and populations from the feedbacks above to impact also coastal GDP and GDP per capita, which has secondary effect like storm surge resilience and money available for flood protection.
- ``` Model.include_retreat_exposure_reduction = True```: This will let the retreat of assets or people reduce the annual exposure to storm surges.
- ``` Model.include_initial_fp = False```: This is for counterfactual scenarios in which the initial flood protection height that is taken from the DIVA dataset is set to 0 (default is True, i.e. the inital value from DIVA is used)


## Creating new aggregation level
If you want to create a new type of aggregation for FRISIA, e.g. separating coastal segments by income, the `processing/createNewAggregation.ipynb` notebook is the place to start. It is set up in a way to produce all input data for the already existing "DIVA_regional" aggregation. See the text in the beginning of the notebook on how to update the file.

Note that running this notebook exactly as it is will overwrite the input files for the "DIVA_regional" aggregation. The resulting files will not be bit-identical to what is commited in the repository, as the fit functions might find different fit parameters. The overall behaviour of the new fit parameters should however be the same. 


## Coupling of FRISIA to FRIDA
FRISIA is coupled to the [FRIDA model](https://github.com/metno/WorldTransFRIDA) in its "DIVA_bipolar" aggregation type. The coupled version has all the additional feedbacks activated. The coupling was done as a complete reimplementation of FRISIA in Stella Architect, the modelling language of FRIDA. The version used in FRIDA v2.1 is FRISIA v1.0. FRIDA will not be continuously updated with newer versions of FRISIA. Only larger updates that affect the feedback to other components of FRIDA will be made.


## Model description
For a full documentation of the model implementation, we refer the user to [Ramme et al., 2025](https://doi.org/10.5194/egusphere-2025-1875).

### FRISIA sea-level rise module
The sea-level rise (SLR) component of FRISIA is largely based on MAGICC (Nauels et al., 2017) and BRICK (Wong et al., 2017). It calculates SLR from five different components (thermosteric, land water storage, mountain glacier, Greenland ice sheet, Antarctic ice sheet). 

#### Thermosteric sea-level rise 
Thermosteric SLR follows ocean heat content changes approximately linearly with the rate factor r (m/YJ) (r is "EEH" in Marti et al., 2022). We here use a default value of r=0.11 (uncertainty range: 0.10 - 0.12).
Reference values I found are:
- 0.145 (Marti et al.2022)
- 0.12 (Levitus et al., 2012)
- 0.15 (Church et al., 2011)
- 0.105 (BRICK v0.3)
- 0.1199 from IPCC estimates
- 0.115 from MPIESM fit
- 0.0975 from fitting OBS data of OHC to thermoSLR

#### Land water storage
Very little consistent data exist on how this will evolve under climate change. At the same time it is largely a function of socio-economic parameters. Here we provide two different options:

1. A linear increase of SLR_LWS as in BRICKv0.3 (Wong, 2017).
2. A linear relationship between (dSLR_LWS)/dt and population.

Both versions are calibrated to historical data and IPCC estimates. The FRISIA version implemented in FRIDA uses groundwater anomaly and installed hydropower capacity from FRIDA to calculate this component.

#### Mountain glaciers
We calculate this component as in Perette et al. (2013), with calibrated parameters and uncertainty ranges.

#### Greenland ice sheet
Calculation as in MAGICC6 (Nauels, 2017), with recalibrated parameters and uncertainty ranges.

#### Antarctic ice sheet
Caclulation as in BRICK (Wong et al., 2017), which employs the DAIS model (Shaffer, 2014), and with recalibrated parameters and uncertainty ranges.


### FRISIA impacts and adaptation module
The impacts and adaptation module is a substantially modified, aggregated version of the [CIAM model](https://github.com/raddleverse/MimiCIAM.jl) (Diaz, 2016, Wong et al., 2022). It incorporates information from coastal segments by aggregating them over specific regions or segment types. The segment information it uses are from CIAM and originally stemming from the DIVA database (Vafeides et al., 2008). For any further information we refer the user to the FRISIA documentation paper.

## How to cite:
Please cite the following paper for any reference to FRISIA.

Ramme, L., Blanz, B., Wells, C., Wong, T. E., Schoenberg, W., Smith, C., and Li, C.: Feedback-based sea level rise impact modelling for integrated assessment models with FRISIAv1.0, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-1875, 2025


## References:
- Church, J. A., N. J. White, L. F. Konikow, C. M. Domingues, J. G. Cogley, E. Rignot, J. M. Gregory, M. R. van den Broeke, A. J. Monaghan, and I. Velicogna (2011), Revisiting the Earth's sea-level and energy budgets from 1961 to 2008, Geophys. Res. Lett., 38, L18601, doi:10.1029/2011GL048794
- Diaz, D. B.: Estimating global damages from sea level rise with the Coastal Impact and Adaptation Model (CIAM), Climatic Change, 137, 143–156, https://doi.org/10.1007/s10584-016-1675-4, 2016
- Levitus, S., et al. (2012), World ocean heat content and thermosteric sea level change (0–2000 m), 1955–2010, Geophys. Res. Lett., 39, L10603, doi:10.1029/2012GL051106
- Marti, F., Blazquez, A., Meyssignac, B., Ablain, M., Barnoud, A., Fraudeau, R., Jugier, R., Chenal, J., Larnicol, G., Pfeffer, J., Restano, M., and Benveniste, J.: Monitoring the ocean heat content change and the Earth energy imbalance from space altimetry and space gravimetry, Earth Syst. Sci. Data, 14, 229–249, https://doi.org/10.5194/essd-14-229-2022, 2022
- Nauels, A., Meinshausen, M., Mengel, M., Lorbacher, K., and Wigley, T. M. L.: Synthesizing long-term sea level rise projections – the MAGICC sea level model v2.0, Geosci. Model Dev., 10, 2495–2524, https://doi.org/10.5194/gmd-10-2495-2017, 2017.
- Perrette, M., Landerer, F., Riva, R., Frieler, K., and Meinshausen, M.: A scaling approach to project regional sea level rise and its uncertainties, Earth Syst. Dynam., 4, 11–29, https://doi.org/10.5194/esd-4-11-2013, 2013
- Shaffer, G.: Formulation, calibration and validation of the DAIS model (version 1), a simple Antarctic ice sheet model sensitive to variations of sea level and ocean subsurface temperature, Geosci. Model Dev., 7, 1803–1818, https://doi.org/10.5194/gmd-7-1803-2014, 2014
- Vafeidis, A. T., Nicholls, R. J., McFadden, L., Tol, R. S. J., Hinkel, J., Spencer, T., Grashoff, P. S., Boot, G., and Klein, R. J. T.: A New Global Coastal Database for Impact and Vulnerability Analysis to Sea-Level Rise, Journal of Coastal Research, 24, 917–924, https://doi.org/10.2112/06-0725.1, 2008
- Wong, T. E., Bakker, A. M. R., Ruckert, K., Applegate, P., Slangen, A. B. A., and Keller, K.: BRICK v0.2, a simple, accessible, and transparent model framework for climate and regional sea-level projections, Geosci. Model Dev., 10, 2741–2760, https://doi.org/10.5194/gmd-10-2741-2017, 2017
- Wong, T. E., Ledna, C., Rennels, L., Sheets, H., Errickson, F. C., Diaz, D., and Anthoff, D.: Sea Level and Socioeconomic Uncertainty Drives High-End Coastal Adaptation Costs, Earth’s Future, 10, e2022EF003 061, https://doi.org/10.1029/2022EF003061, 2022
