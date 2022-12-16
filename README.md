# Modelica district heating optimization MPC
This is a Modelica package for district heating and building heating systems modelling, optimization and control. The main contribution of this package are:
1) A fast and effective approach to model, optimize and control district heating and building heating systems;
2) Year-level measured data for a campus district heating system and a university building in Norway is opened for testing and research.

The package includes the following parts:
1) Package for MPC of building heating systems, [MPC building heating system](https://github.com/haoranli1988/Modelica-district-heating-optimization-MPC/tree/main/MPC%20building%20heating%20system). It applies MPC to buildings with space heating and ventilation systems. The package includes an RC model, radiator and ventilation system model components and an MPC framework. The package also contains the measured data on a university building in Norway for testing.   
2) Package for optimizing district heating systems, [Optimization DH system](https://github.com/haoranli1988/Modelica-district-heating-optimization-MPC/tree/main/Optimization%20DH%20system). It applies up to year-level optimal operation for district heating systems. The package includes model components of the main substation, data centre waste heat recovery, network and buildings, and an optimization framework. It also contains measured data on a university campus district heating system in Norway for testing.
3) Package for MPC district heating systems, [MPC DH system](https://github.com/haoranli1988/Modelica-district-heating-optimization-MPC/tree/main/MPC%20DH%20system). It applies MPC to district heating systems. The package uses system models and data from [Optimization DH system](https://github.com/haoranli1988/Modelica-district-heating-optimization-MPC/tree/main/Optimization%20DH%20system). It contains an MPC framework for district heating systems.

The package is developed by [Haoran Li](https://www.linkedin.com/in/haoran-li-4397311ba/) and [Juan Hou](https://www.linkedin.com/in/juan-hou-4a54a622a/) during their study at the Norwegian University of Science and Technology (NTNU). If you have any questions regarding the package, feel free to contact us at lihaoran198811@gmail.com and juan.hou2022@gmail.com. For questions regarding the case study and measured data, please contact natasa.nord@ntnu.no.   

## How to cite the package
For the package for MPC of building heating systems, please cite [Hou J, Li H, Nord N, Huang G. Model predictive control under weather forecast uncertainty for HVAC systems in university buildings. Energy and Buildings. 2022 Feb 15;257:111793](https://www.sciencedirect.com/science/article/pii/S037877882101077X).

For the package for optimizing district heating systems, please cite [Li H, Hou J, Tian Z, Hong T, Nord N, Rohde D. Optimize heat prosumers' economic performance under current heating price models by using water tank thermal energy storage. Energy. 2022 Jan 15;239:122103](https://www.sciencedirect.com/science/article/pii/S0360544221023513).  

For the package for MPC of DH system, please cite [Hou J. Investigation of Model Predictive Control Application in District Heating Systems with Distributed Sources](https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/3016473).  

## Publications using the package
[Hou J, Li H, Nord N. Nonlinear model predictive control for the space heating system of a university building in Norway. Energy. 2022 Aug 15;253:124157](https://www.sciencedirect.com/science/article/pii/S036054422201060X).

[Li H, Hou J, Hong T, Nord N. Distinguish between the economic optimal and lowest distribution temperatures for heat-prosumer-based district heating systems with short-term thermal energy storage. Energy. 2022 Jun 1;248:123601](https://www.sciencedirect.com/science/article/pii/S0360544222005047).

[Li H. Economic optimization for heat prosumer-based district heating systems in unidirectional heating markets](https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/2837506).

[Hou J. Investigation of Model Predictive Control Application in District Heating Systems with Distributed Sources](https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/3016473).

## Acknowledgements
The package is the output of the Norwegian research council project under FRIPRO/FRINATEK program Understanding Behaviour of District Heating Systems Integrating Distributed Sources (project number 262707). In addition, the maintenance of the package receives funding from the EU research project Climate Positive Circular Communities (ARV) (grant agreement ID: 101036723).
