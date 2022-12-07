within ;
package Plants
  package Components
    model Building
      // Note 1: building envelope= wall + roof
      // Note 2: load shifting only conducted by radiator system
      // Initial states

      // Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws "Supply water temperature (degC)";
      input Modelica.SIunits.MassFlowRate mw "Mass flow rate of water (kg/s)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Dif "Temperature difference (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr "Temperature of return water (degC)";
      // Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Real Ten_init=13;
      parameter Real Tia_init=13;
      parameter Real Tma_init=13;
      parameter Real Cen=4.5e10 "Thermal capacity of building envelop (J/K)";
      parameter Real Cia=1.3e9 "Thermal capacity of indoor air (J/K)";
      parameter Real Cma=2.9e9 "Thermal capacity of inetrior thermal medium, i.e. interior walls and furniture (J/K)";
      parameter Real Roe=1.03 "The thermal resistance between building envelope and outdoor air (m2K/W)";
      parameter Real Rie=1.18 "The thermal resistance between buidling envelope and indoor air (m2K/W)";
      parameter Real Rmi=0.35 "The thermal resistance between indoor air and interior walls and furniture (m2K/W)";
      parameter Real Rwin=0.48 "The thermal resistance between indoor and outdoor air through windows (m2K/W)";
      parameter Real Aen=1.8e5 "The area of building envelope (m2)";
      parameter Real Ama=4.2e5 "The area of interior walls and furniture (m2)";
      parameter Real Awin=4.1e4 "The area of building windowns (m2)";
      // Manipulated variable: heat flow, all in postive values
      Modelica.Blocks.Interfaces.RealInput qra "Heat flow rate from Radiator to Indoor air (W)";
      // Input variables: weather
      Modelica.Blocks.Interfaces.RealInput Toa;
      // Input variables: heat flow, all in postive values
      Modelica.Blocks.Interfaces.RealInput qin "Heat flow rate of internal heat gain (W)";
      Modelica.Blocks.Interfaces.RealInput qinf "Heat flow rate of infiltration (W)";
      Modelica.Blocks.Interfaces.RealInput qve "Heat flow rate of ventilation system (W)";
      Modelica.Blocks.Interfaces.RealInput qdhw "Heat flow rate of DHW system (W)";
      // Variables
      Real qload "Heating load of building(W)";
      Real qsh "Heat flow rate of SH system (W)";
      Real qoe "Heat flow rate from outdoor air to building envelope (W)";
      Real qie "Heat flow rate from indoor air to building envelope (W)";
      Real qwin "Heat flow rate from outdoor air to indoor air through windows (W)";
      Real qmi "Heat flow rate from interior walls and furnitures to indoor air (W)";
      // State variables
      Real Ten( nominal=10) "Temperature of buidling envelop (Centigrade)";
      Real Tia( nominal=20) "Temperature of indoor air (Centigrade)";
      Real Tma( nominal=20) "Temperature of interior walls and furniture (Centigrade)";
    initial equation
      Ten=Ten_init;
      Tia=Tia_init;
      Tma=Tma_init;
    equation
      // Building
      qoe=Aen*(Toa-Ten)/Roe;
      qie=Aen*(Tia-Ten)/Rie;
      qwin=Awin*(Toa-Tia)/Rwin;
      qmi=Ama*(Tma-Tia)/Rmi;
      Cen*der(Ten)=qoe+qie;
      Cia*der(Tia)=-qie+qmi+qwin-qinf+qra+qin "vantilation system has no impact";
      Cma*der(Tma)=-qmi;
      // Heating system
      qsh=qra+qve;
      qload=qsh+qdhw;
      // Connection
      Dif=Tws-Twr;
      qload=cw*mw*Dif;
    end Building;

    model Infiltration
      // Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity ca=1002 "Specific heat capacity at constant volume of air (J/(kg.K))";
      // Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tai "Temperature of inlet air (degC)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tao "Temperature of outlet air (degC)";
      input Modelica.SIunits.MassFlowRate ma "Mass flow rate of air (kg/s)";
      output Modelica.SIunits.HeatFlowRate qinf "Total Heat flow rate from Infiltration to indoor air (W)";
    equation
      qinf=ca*ma*(Tao-Tai);
    end Infiltration;

    model Radiator
      // Initial states
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tra_init[N]=fill(40,N);
      // Parameters
      parameter Integer N=10 "Number of sections for radiator";
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Real Scra=1.2*7380 "Number of radiators";
      parameter Real a=7.6187 "characteristic coefficients of the radiator";
      parameter Real b=1.2831 "characteristic coefficients of the radiator";
      parameter Modelica.SIunits.HeatCapacity Cwra=16000 "Heat capacity of water in radiator (= cp*m)(J/K)";
      // Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws "Water supply tmeperature of radiator (degC)";
      input Modelica.SIunits.MassFlowRate mw "Mass flow rate of water (kg/s)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tia "Temperature of indoor air (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr "Water return tmeperature of radiator (degC)";
      output Modelica.SIunits.HeatFlowRate qra_sum "Total Heat flow rate from radiator to indoor air (W)";
      // Variables
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tra[N] "Temperature of water flow and radiator of each section (degC)";
      Modelica.SIunits.HeatFlowRate qra[N] "Heat flow rate from radiator (each section) to indoor air (W)";
    initial equation
      Tra=Tra_init;
    equation
      qra_sum=Scra*sum(qra);
      qra[1]=a/N*max(0,Tra[1]-Tia)^b;
      Cwra/N*der(Tra[1])=cw*(mw/Scra)*(Tws-Tra[1])-qra[1];
      for i in 2:(N-1) loop
          qra[i]=a/N*max(0,Tra[i]-Tia);
          Cwra/N*der(Tra[i])=cw*(mw/Scra)*(Tra[i-1]-Tra[i])-qra[i];
      end for;
      qra[end]=a/N*max(0,Tra[end]-Tia)^b;
      Cwra/N*der(Tra[end])=cw*(mw/Scra)*(Tra[end-1]-Tra[end])-qra[end];
      Twr=Tra[end];
    end Radiator;

    model AHU
      // Initial states
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Ta_init[N]=fill(20,N);
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tm_init[N]=fill(40,N);
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tw_init[N]=fill(40,N);
      // Parameters
      parameter Integer N=10 "Number of sections for AHU";
      parameter Integer nt=10 "Number of tubes for AHU";
      parameter Modelica.SIunits.Diameter dt=0.0117 "Diameter of tube for AHU";
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Modelica.SIunits.SpecificHeatCapacity ca=1002 "Specific heat capacity at constant volume of air (J/(kg.K))";
      parameter Real Scra=470 "Number of AHU";
      parameter Modelica.SIunits.Area Aw=0.404 "The area of water side (m2)";
      parameter Modelica.SIunits.Area Aa=5.181 "The area of air side (m2)";
      parameter Modelica.SIunits.Efficiency Eps=0.6 "Eps of heat recovery";
      parameter Modelica.SIunits.HeatCapacity Caahu=483 "Heat capacity of air in AHU (= cp*m)(J/K)";
      parameter Modelica.SIunits.HeatCapacity Cmahu=4047 "Heat capacity of metal wall in AHU (= cp*m)(J/K)";
      parameter Modelica.SIunits.HeatCapacity Cwahu=4836 "Heat capacity of water in AHU (= cp*m)(J/K)";
      // Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tia "Temperature of indoor air (degC)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Toa "Temperature of outdoor air (degC)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws "Temperature of supply water (degC)";
      input Modelica.SIunits.MassFlowRate mw "Mass flow rate of water (kg/s)";
      input Modelica.SIunits.MassFlowRate ma "Mass flow rate of air (kg/s)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tao "Temperature of outlet air (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr "Temperature of return water (degC)";
      output Modelica.SIunits.HeatFlowRate qahu_sum "Total Heat flow rate from AHU to indoor air (W)";
      // Variables
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tai "Temperature of inlet air (degC)";
      Modelica.SIunits.HeatFlowRate qw2m[N] "Heat flow rate from water to metal wall of each section (W)";
      Modelica.SIunits.HeatFlowRate qm2a[N] "Heat flow rate from meatal wall to air of each section (W)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tm[N] "Temperature of metal wall of each section (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tw[N] "Temperature of water of each section (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Ta[N] "Temperature of outlet air of each section (degC)";
      Modelica.SIunits.CoefficientOfHeatTransfer Uw2m[N] "U value between water and metal wall (W/(m2.K))";
      Modelica.SIunits.CoefficientOfHeatTransfer Um2a[N] "U value between metal wall and air (W/(m2.K))";
    initial equation
      Ta=Ta_init;
      Tm=Tm_init;
      Tw=Tw_init;
    equation
      Eps=(Tai-Toa)/(Tia-Toa+1E-6);
      qahu_sum=Scra*sum(qm2a);
      Uw2m[1]=(5.823+1.153*10^(-1)*Tw[1]-1.48*10^(-4)*Tw[1]^2)*((mw/Scra)/(nt*dt^2.25))^0.8;
      Um2a[1]=60.37+140.3*exp((28.62-Tm[1])/9.796);
      qw2m[1]=Uw2m[1]*Aw/N*(Tw[1]-Tm[1]);
      qm2a[1]=Um2a[1]*Aa/N*(Tm[1]-Ta[1]);
      Caahu/N*der(Ta[1])=qm2a[1]-ca*(ma/Scra/N)*(Ta[1]-Tai);
      Cmahu/N*der(Tm[1])=qw2m[1]-qm2a[1];
      Cwahu/N*der(Tw[1])=cw*(mw/Scra)*(Tws-Tw[1])-qw2m[1];
      for i in 2:N loop
          Uw2m[i]=(5.823+1.153*10^(-1)*Tw[i]-1.48*10^(-4)*Tw[i]^2)*((mw/Scra)/(nt*dt^2.25))^0.8;
          Um2a[i]=60.37+140.3*exp((28.62-Tm[i])/9.796);
          qw2m[i]=Uw2m[i]*Aw/N*(Tw[i]-Tm[i]);
          qm2a[i]=Um2a[i]*Aa/N*(Tm[i]-Ta[i]);
          Caahu/N*der(Ta[i])=qm2a[i]-ca*(ma/Scra/N)*(Ta[i]-Tai);
          Cmahu/N*der(Tm[i])=qw2m[i]-qm2a[i];
          Cwahu/N*der(Tw[i])=cw*(mw/Scra)*(Tw[i-1]-Tw[i])-qw2m[i];
      end for;
      Twr=Tw[end];
      Tao=sum(Ta)/N;
    end AHU;

    model HExchanger_SH
      // Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Real Scra=23 "Number of heat exchangers";
      parameter Modelica.SIunits.Area A=13.2 "The area of heat exchanger wall (m2)";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer U=6000 "U value of heat exchanger wall (W/(m2.K))";
      // Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws1 "Temperature of supply water in primary side (degC)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr2 "Temperature of return water in secondary side (degC)";
      input Modelica.SIunits.MassFlowRate mw1 "Mass flow rate of water in primary side (kg/s)";
      input Modelica.SIunits.MassFlowRate mw2 "Mass flow rate of water in secondary side (kg/s)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr1 "Temperature of return water in primary side (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws2 "Temperature of supply water in secondary side (degC)";
      // Variables
      Real NTU "NTU of heat exchanger";
      Modelica.SIunits.Efficiency Eps "Eps of heat exchanger";
      Modelica.SIunits.EntropyFlowRate Cw1 "Heat capacity rate of water flow in primary side";
      Modelica.SIunits.EntropyFlowRate Cw2 "Heat capacity rate of water flow in secondary side";
      Modelica.SIunits.EntropyFlowRate Cw_min "Heat capacity rate of water flow in Min side";
      Modelica.SIunits.EntropyFlowRate Cw_max "Heat capacity rate of water flow in Max side";
      Real Cwr "Cw_min/Cw_max";
      Modelica.SIunits.HeatFlowRate qhx "Heat flow rate of heat exchanger (W)";
      Modelica.SIunits.HeatFlowRate qhx_sum "Total Heat flow rate of heat exchanger (W)";
    equation
      qhx_sum=Scra*qhx;
      Cw1=cw*(mw1/Scra)+1E-6;
      Cw2=cw*(mw2/Scra)+1E-6;
      Cw_min=min(Cw1,Cw2);
      Cw_max=max(Cw1,Cw2);
      Cwr=Cw_min/Cw_max;
      NTU=U*A/Cw_min;
      qhx=Cw1*(Tws1-Twr1);
      qhx=Cw2*(Tws2-Twr2);
      Eps=qhx/Cw_min/(Tws1-Twr2);
      Eps=(1-exp(-NTU*(1-Cwr)))/(1-Cwr*exp(-NTU*(1-Cwr)));
    end HExchanger_SH;

    model HExchanger_DHW
      // Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Real Scra=23 "Number of heat exchangers";
      parameter Modelica.SIunits.Area A=1.2 "The area of heat exchanger wall (m2)";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer U=6000 "U value of heat exchanger wall (W/(m2.K))";
      // Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws1 "Temperature of supply water in primary side (degC)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr2 "Temperature of return water in secondary side (degC)";
      input Modelica.SIunits.MassFlowRate mw1 "Mass flow rate of water in primary side (kg/s)";
      input Modelica.SIunits.MassFlowRate mw2 "Mass flow rate of water in secondary side (kg/s)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr1 "Temperature of return water in primary side (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws2 "Temperature of supply water in secondary side (degC)";
      // Variables
      Real NTU "NTU of heat exchanger";
      Modelica.SIunits.Efficiency Eps "Eps of heat exchanger";
      Modelica.SIunits.EntropyFlowRate Cw1 "Heat capacity rate of water flow in primary side";
      Modelica.SIunits.EntropyFlowRate Cw2 "Heat capacity rate of water flow in secondary side";
      Modelica.SIunits.EntropyFlowRate Cw_min "Heat capacity rate of water flow in Min side";
      Modelica.SIunits.EntropyFlowRate Cw_max "Heat capacity rate of water flow in Max side";
      Real Cwr "Cw_min/Cw_max";
      Modelica.SIunits.HeatFlowRate qhx "Heat flow rate of heat exchanger (W)";
      Modelica.SIunits.HeatFlowRate qhx_sum "Total Heat flow rate of heat exchanger (W)";
    equation
      qhx_sum=Scra*qhx;
      Cw1=cw*(mw1/Scra)+1E-6;
      Cw2=cw*(mw2/Scra)+1E-6;
      Cw_min=min(Cw1,Cw2);
      Cw_max=max(Cw1,Cw2);
      Cwr=Cw_min/Cw_max;
      NTU=U*A/Cw_min;
      qhx=Cw1*(Tws1-Twr1);
      qhx=Cw2*(Tws2-Twr2);
      Eps=qhx/Cw_min/(Tws1-Twr2);
      Eps=(1-exp(-NTU*(1-Cwr)))/(1-Cwr*exp(-NTU*(1-Cwr)));
    end HExchanger_DHW;

    model BuildingAndSystem
      Plants.Components.Building Building(Ten_init=Ten_init,Tia_init=Tia_init,Tma_init=Tma_init);
      Plants.Components.Radiator Radiator(N=N,Tra_init=Tra_init);
      Plants.Components.Infiltration Infiltration;
      Plants.Components.AHU AHU(N=N,Ta_init=Ta_init,Tm_init=Tm_init,Tw_init=Tw_init);
      Plants.Components.HExchanger_SH HExchanger_SH;
      Plants.Components.HExchanger_DHW HExchanger_DHW;
      // Global initial states
      parameter Integer N=10 "Number of sections for AHU, Radiator, and water tank";
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Ten_init=10;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tia_init=20;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tma_init=20;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tra_init[N]=fill(40,N);
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Ta_init[N]=fill(20,N);
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tm_init[N]=fill(40,N);
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tw_init[N]=fill(40,N);
      // Global parameters
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tgrou=5 "Temperature of ground (degC)";
      // Global input variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Toa "Temperature of outdoor air (degC)";
      input Modelica.SIunits.HeatFlowRate qin "Internal heat gain internal heat gains due to lighting, equipment, and occupants (W)";
      input Modelica.SIunits.MassFlowRate mwdhw "Water Mass flow rate of tap water (kg/s)";
      input Modelica.SIunits.MassFlowRate maahu "Air Mass flow rate of AHU (kg/s)";
      input Modelica.SIunits.MassFlowRate mainf "Air Mass flow rate of Infiltration (kg/s)";
      // Manipulated variable
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws1 "Temperature of supply water in primary side (degC)";
      input Modelica.SIunits.MassFlowRate mw1sh "Mass flow rate of water in primary side of space heating system (kg/s)";
      input Modelica.SIunits.MassFlowRate mw1dhw "Mass flow rate of water in primary side of domestic hot water system (kg/s)";
      input Modelica.SIunits.MassFlowRate mwra "Water Mass flow rate of radiator (kg/s)";
      input Modelica.SIunits.MassFlowRate mwahu "water Mass flow rate of AHU (kg/s)";
      // Global output variables
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr1sh "Temperature of return water in primary side of space heating system (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr1dhw "Temperature of return water in primary side of domestic hot water system (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr1 "Temperature of return water in primary side (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws2sh "Temperature of supply water in secondary side of space heating system (degC)";
      Modelica.SIunits.MassFlowRate mw1 "Mass flow rate of water in primary side (kg/s)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tia "Temperature of indoor air (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Taoahu "Temperature of outlet air of AHU (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twodhw "Temperature of outlet water of domestic hot water (degC)";
      Modelica.SIunits.HeatFlowRate qload "Heating load of space heating and domestic hot water system (W)";
    equation
      // Level HeatExcahnger_SH
      HExchanger_SH.Tws1=Tws1;
      HExchanger_SH.mw1=mw1sh;
      HExchanger_SH.Twr2*HExchanger_SH.mw2=Radiator.Twr*Radiator.mw+AHU.Twr*AHU.mw;
      HExchanger_SH.mw2=Radiator.mw+AHU.mw;
      // Level HeatExcahnger_DHW
      HExchanger_DHW.Tws1=Tws1;
      HExchanger_DHW.mw1=mw1dhw;
      HExchanger_DHW.Twr2=Tgrou;
      HExchanger_DHW.mw2=mwdhw;
      // Level Building
      Building.Toa=Toa;
      Building.qin=qin;
      Building.qra=Radiator.qra_sum;
      Building.qve=AHU.qahu_sum;
      Building.qinf=Infiltration.qinf;
      Building.ma=maahu;
      // Level Infiltration
      Infiltration.Tai=Toa;
      Infiltration.Tao=Building.Tia;
      Infiltration.ma=mainf;
      // Level Radiator
      Radiator.Tws=HExchanger_SH.Tws2;
      Radiator.mw=mwra;
      Radiator.Tia=Building.Tia;
      // Level AHU
      AHU.Tia=Tia;
      AHU.Toa=Toa;
      AHU.Tws=HExchanger_SH.Tws2;
      AHU.mw=mwahu;
      AHU.ma=maahu;
      // Addition outputs
      Tia=Building.Tia;
      Taoahu=AHU.Tao;
      Twr1sh=HExchanger_SH.Twr1;
      Tws2sh=HExchanger_SH.Tws2;
      HExchanger_DHW.Twr1=Twr1dhw;
      HExchanger_DHW.Tws2=Twodhw;
      mw1*Twr1=mw1sh*Twr1sh+mw1dhw*Twr1dhw;
      mw1=mw1sh+mw1dhw;
      qload=Building.qload+HExchanger_DHW.qhx_sum;
    end BuildingAndSystem;

    model ThermalTank
      // Initial states
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt_ini=25;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt_init[N]=fill(Twt_ini,N);
      // Parameters
      parameter Integer N=10 "Number of sections for radiator";
      parameter Modelica.SIunits.Volume Vwt "Volume of water tank thermal storage (m3)";
      parameter Modelica.SIunits.Height Hwt=15 "Height of water tank thermal storage (m)";
      parameter Modelica.SIunits.Diameter Dwt=2*(Vwt/Hwt/3.14)^0.5 "Diameter of water tank thermal storage (m)";
      parameter Modelica.SIunits.Area Awtl=3.14*Dwt*Hwt "Lateral Area of water tank thermal storage (m2)";
      parameter Modelica.SIunits.Area Awtb=3.14*(Dwt/2)^2 "Up or Down Base Area of water tank thermal storage (m2)";
      parameter Modelica.SIunits.HeatCapacity Cwt=cw*Vwt*rho "Heat capacity of water in radiator (= cp*m)(J/K)";
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tgrou=5 "Temperature of ground (degC)";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer Uwt=1.2 "U value between water tank thermal storage and ground (W/(m2.K))";
      parameter Modelica.SIunits.ThermalConductivity keff=4.1 "Effective ThermalConductivity of water tank thermal storage between each section (W/(m.K))";
      parameter Modelica.SIunits.Density rho=995.6  "Density of water (kg/m3)";
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      // Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws1 "Water supply tmeperature of heat source side (degC)";
      input Modelica.SIunits.MassFlowRate mw1 "Mass flow rate of water of heat source side (kg/s)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr2 "Water return tmeperature of heat user side (degC)";
      input Modelica.SIunits.MassFlowRate mw2 "Mass flow rate of water of heat user side (kg/s)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws2 "Water supply tmeperature of heat user side (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr1 "Water return tmeperature of heat source side (degC)";
      // Variables
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt[N] "Temperature of water in water tank thermal storage of each section (degC)";
      Modelica.SIunits.HeatFlowRate qloss[N] "Heat loss flow rate from water tank thermal storage (each section) to ground (W)";
      Modelica.SIunits.HeatFlowRate qlosswt "Heat loss flow rate from water tank thermal storage to ground (W)";
      Modelica.SIunits.ThermalConductivity kwt "ThermalConductivity of water tank thermal storage between each section considers inversion effect(W/(m.K))";
      Modelica.SIunits.Heat Qstore "Stored heat (J)";
    initial equation
      Twt=Twt_init;
    equation
      Qstore=Cwt*(sum(Twt)/N-Twt_ini);
      qlosswt=sum(qloss);
      qloss[1]=Uwt*(Awtl/N+Awtb)*(Twt[1]-Tgrou);
      kwt=keff*(1+exp(5*(Twt[1]-Twt[end])));
      Cwt/N*der(Twt[1])=cw*mw1*(Twt[2]-Twt[1])+cw*mw2*(Twr2-Twt[1])-qloss[1]+kwt*Awtb/(Hwt/N)*(Twt[2]-Twt[1]);
      for i in 2:(N-1) loop
          qloss[i]=Uwt*Awtl/N*(Twt[i]-Tgrou);
          Cwt/N*der(Twt[i])=cw*mw1*(Twt[i+1]-Twt[i])+cw*mw2*(Twt[i-1]-Twt[i])-qloss[i]+kwt*Awtb/(Hwt/N)*(Twt[i+1]-Twt[i])+kwt*Awtb/(Hwt/N)*(Twt[i-1]-Twt[i]);
      end for;
      qloss[end]=Uwt*Awtl/N*(Twt[end]-Tgrou);
      Cwt/N*der(Twt[end])=cw*mw1*(Tws1-Twt[end])+cw*mw2*(Twt[end-1]-Twt[end])-qloss[end]+kwt*Awtb/(Hwt/N)*(Twt[end-1]-Twt[end]);
      Tws2=Twt[end];
      Twr1=Twt[1];
    end ThermalTank;

    model Substation
      // Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      // Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws2 "Temperature of supply water in secondary side (degC)";
      input Modelica.SIunits.MassFlowRate mw2 "Mass flow rate of water in secondary side (kg/s)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr2 "Temperature of return water in secondary side (degC)";
      output Modelica.SIunits.HeatFlowRate qhx "Heat flow rate of heat exchanger (W)";
    equation
      qhx=cw*mw2*(Tws2-Twr2);
    end Substation;

    model DataCenter
      // Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Modelica.SIunits.MassFlowRate mwdc=18.7 "Mass flow rate of water into data center (kg/s)";
      // Connection variables
      input Modelica.SIunits.HeatFlowRate qdc "Heat feedin flow rate of data center (W)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twidc "Inlet Temperature of data center (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twodc "Outlet Temperature of data center (degC)";
    equation
      qdc=cw*mwdc*(Twodc-Twidc);
    end DataCenter;

    model BuildingLoadData
      // Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      // Connection variables
      input Modelica.SIunits.HeatFlowRate qload "Heat supply flow rate of DH system (W)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws "Supply water temperature (degC)";
      input Modelica.SIunits.MassFlowRate mw "Mass flow rate of water (kg/s)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Dif "Temperature difference (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr "Temperature of return water (degC)";
    equation
      Dif=Tws-Twr;
      qload=cw*mw*Dif;
    end BuildingLoadData;

    model HPData
      // Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Modelica.SIunits.MassFlowRate mwevap=36.5 "Mass flow rate of water of Evaporator (kg/s)";
      parameter Modelica.SIunits.MassFlowRate mwcond=18.7 "Mass flow rate of water of Condenser (kg/s)";
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tevapin=11.2 "Inlet water temperature of Evaporator (degC)";
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tevapout=7.0 "Outlet water temperature of Evaporator (degC)";
      parameter Real a=-1.5021;
      parameter Real b=10.5052;
      parameter Real c=-16.9627;
      parameter Real d=24.3869;
      parameter Real e=-5.7348;
      parameter Real f=-9.3339;
      // Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tctin "Inlet water temperature of Condenser (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tcondout "Outlet water temperature of Condenser (degC)";
      output Modelica.SIunits.HeatFlowRate qloss "Heat flow rate loss of Cooling tower (W)";
      output Modelica.SIunits.Power P "Power of HP (W)";
      output Modelica.SIunits.HeatFlowRate qevap "Heat flow rate of Evaporator (W)";
      output Modelica.SIunits.HeatFlowRate qcond "Heat flow rate of Condenser (W)";
      output Real COP "Heating COP of HP (1)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tcondin "Inlet water temperature of Condenser (degC)";
    equation
      qloss=max(0,cw*mwcond*(Tctin-60));
      Tcondin=min(60,Tctin);
      P=1000*(a*Tevapin+b*Tevapout+c*Tcondin+d*Tcondout+e*mwevap+f*mwcond);
      qcond=qevap+P;
      qevap=cw*mwevap*(Tevapin-Tevapout);
      qcond=cw*mwcond*(Tcondout-Tcondin);
      COP=qcond/P;
    end HPData;

    model PipeTDrop
      // Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Modelica.SIunits.Density rho=995.6  "Density of water (kg/m3)";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer K=2.2 "Coefficient Of HeatTransfer for pipe (W/(m2.K))";
      parameter Modelica.SIunits.Length L=1500 "Length Of pipe (m)";
      parameter Modelica.SIunits.Diameter Do=0.2731 "Outer Diameter of pipe (m)";
      parameter Modelica.SIunits.Diameter Di=0.2223 "Inner Diameter of pipe (m)";
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tgrou=5 "Temperature of ground (degC)";
      // Variables
      Modelica.SIunits.Velocity v "Water Velocity inside pipe (m/s)";
      Modelica.SIunits.HeatFlowRate qloss "Heat loss flow rate from pipe to ground (W)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tdrop "Water temperature drop of pipe (degC)";
      // Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tin "Inlet water temperature of pipe (degC)";
      input Modelica.SIunits.MassFlowRate mw "Mass flow rate of water of heat source side (kg/s)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tout "Outlet water temperature of pipe (degC)";
    equation
      mw=3.14*(Di/2)^2*v*rho;
      Tdrop=4*K*L*Do*(Tin-Tgrou)/(v*Di^2*rho*cw);
      Tout=Tin-Tdrop;
      qloss=cw*mw*Tdrop;
    end PipeTDrop;

    model PipeHeatLoss
      // Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Modelica.SIunits.ThermalConductivity Kgrou=1.5 "ThermalConductivity of ground (W/(m.K))";
      parameter Modelica.SIunits.ThermalConductivity Kinsu=0.03 "ThermalConductivity of insulation (W/(m.K))";
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tgrou=5 "Temperature of ground (degC)";
      parameter Modelica.SIunits.Length L=1500 "Length Of pipe (m)";
      parameter Modelica.SIunits.Diameter d=0.273 "Outer Diameter of pipe (m)";
      parameter Modelica.SIunits.Diameter D=d+0.05*2 "Outer Diameter of insulation (m)";
      parameter Modelica.SIunits.Length s=1.2 "Distance between two pipes (m)";
      parameter Modelica.SIunits.Length h=1.2 "Distance between cnetral pipe and ground surface (m)";
      parameter Modelica.SIunits.ThermalInsulance Ri=(d/2/Kinsu)*log(D/d) "The thermal resistance of insulation (m2.K/W)";
      parameter Modelica.SIunits.ThermalInsulance Rg=(d/2/Kgrou)*log(4*h/d) "The thermal resistance of ground (m2.K/W)";
      parameter Modelica.SIunits.ThermalInsulance Rc=(d/2/Kgrou)*log(((2*h/s)^2+1)^0.5) "The thermal resistance of coinciding (m2.K/W)";
      // Variables
      Modelica.SIunits.HeatFlowRate qlosss "Heat loss flow rate from supply pipe (W)";
      Modelica.SIunits.HeatFlowRate qlossr "Heat loss flow rate from return pipe (W)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tdifs "Temperature difference between supply pipe and ground (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tdifr "Temperature difference between return pipe and ground (degC)";
      // Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tins "Inlet water temperature of supply pipe (degC)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tinr "Inlet water temperature of return pipe (degC)";
      input Modelica.SIunits.MassFlowRate mw "Mass flow rate of water in the pipe supply = return (kg/s)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Touts "Outlet water temperature of supply pipe (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Toutr "Outlet water temperature of return pipe (degC)";
    equation
      Tdifs=Tins-Tgrou;
      Tdifr=Tinr-Tgrou;
      qlosss=L*3.14*d*((Rg+Ri)*Tdifs-Rc*Tdifr)/((Rg+Ri)^2-Rc^2);
      qlossr=L*3.14*d*((Rg+Ri)*Tdifr-Rc*Tdifs)/((Rg+Ri)^2-Rc^2);
      qlosss=cw*mw*(Tins-Touts);
      qlossr=cw*mw*(Tinr-Toutr);
    end PipeHeatLoss;

    package Validation
      model Building_V
        parameter Modelica.SIunits.Time T=3600 "Time Constant";
        parameter Real k(unit="1")=1 "Gain";
        Real u "Connector of Real input signal";
        Real y "Connector of Real output signal";
        Plants.Components.Building Building;
        Modelica.SIunits.HeatFlowRate Error_Tia;
        Modelica.SIunits.HeatFlowRate qload_meas;
        Modelica.SIunits.HeatFlowRate Error_qload;
        Modelica.Blocks.Sources.CombiTimeTable Inputs(
          tableOnFile=true,
          table=fill(
              0.0,
              0,
              16),
          tableName="Data",
          fileName="C:/Work/Models/InputData/BuildingAndSystem_V.mat") "Time, 1 Toa, 2 Tiaref, 3 Tws1, 4 Twr1, 5 qMS, 6 qDC, 7 qin, 8 qen, 9 qinf, 10 qve, 11 qDHW, 12 mw1, 13 mainf, 14 maahu, 15 mwdhw
"         annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      equation
        Building.Toa=Inputs.y[1];
        Building.qin=Inputs.y[7];
        Building.qve=Inputs.y[10];
        Building.qinf=Inputs.y[9];
        Building.ma=Inputs.y[14];
        Error_Tia=Inputs.y[2]-Building.Tia;
        u=min(8.5E6,20E6*max(0,min(1,Error_Tia)));
        der(y) = (k*u - y)/T;
        Building.qra=y;
        qload_meas=max(0,Inputs.y[8]+Inputs.y[9]+Inputs.y[10]-Inputs.y[7]);
        Error_qload=qload_meas-Building.qload;
        annotation (experiment(StopTime=31536000, Interval=3600));
      end Building_V;

      model Radiator_V
        Plants.Components.Radiator Radiator;
      equation
        Radiator.Tws=75;
        Radiator.mw=200;
        Radiator.Tia=20;
        annotation (experiment(StopTime=86400, Interval=3600));
      end Radiator_V;

      model AHU_V
        Plants.Components.AHU AHU;
        Modelica.SIunits.HeatFlowRate Error_qve;
        Modelica.Blocks.Sources.CombiTimeTable Inputs(
          tableOnFile=true,
          table=fill(
                    0.0,
                    0,
                    5),
          tableName="Data",
          fileName="C:/Work/Models/InputData/AHU_V.mat")
          "Time, 1 Toa, 2 Tiaref, 3 qve, 4 maahu"
          annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
      equation
        AHU.Tia=Inputs.y[2];
        AHU.Toa=Inputs.y[1];
        AHU.Tws=70;
        AHU.mw=min(1,max(0,AHU.Tia-AHU.Tao))*300;
        AHU.ma=Inputs.y[4];
        Error_qve=Inputs.y[3]-AHU.qahu_sum;
        annotation (experiment(StopTime=31536000, Interval=3600));
      end AHU_V;

      model Infiltration_V
        Plants.Components.Infiltration Infiltration;
      equation
        Infiltration.ma=70;
        Infiltration.Tai=-15;
        Infiltration.Tao=20;
        annotation (experiment(StopTime=86400, Interval=3600));
      end Infiltration_V;

      model HExchanger_SH_V
        Plants.Components.HExchanger_SH HExchanger_SH;
      equation
        HExchanger_SH.Tws1=80;
        HExchanger_SH.Twr2=40;
        HExchanger_SH.mw1=80;
        HExchanger_SH.mw2=160;
        annotation (experiment(StopTime=86400, Interval=3600));
      end HExchanger_SH_V;

      model HExchanger_DHW_V
        parameter Modelica.SIunits.Time T=1 "Time Constant";
        parameter Real k(unit="1")=1 "Gain";
        Real u "Connector of Real input signal";
        Real y( start=0.6) "Connector of Real output signal";
        Plants.Components.HExchanger_DHW HExchanger_DHW;
        parameter Real Kp=6 "Kp of PID";
        Real Err( start=50) "Error of PID";
        Modelica.Blocks.Sources.CombiTimeTable Inputs(
          tableOnFile=true,
          table=fill(
              0.0,
              0,
              16),
          tableName="Data",
          fileName="C:/Work/Models/InputData/BuildingAndSystem_V.mat") "Time, 1 Toa, 2 Tiaref, 3 Tws1, 4 Twr1, 5 qMS, 6 qDC, 7 qin, 8 qen, 9 qinf, 10 qve, 11 qDHW, 12 mw1, 13 mainf, 14 maahu, 15 mwdhw
"         annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      equation
        Err=55-HExchanger_DHW.Tws2;
        HExchanger_DHW.Tws1=65;
        HExchanger_DHW.Twr2=5;
        u=Kp*min(1,max(0,Err));
        der(y) = (k*u - y)/T;
        HExchanger_DHW.mw1=y;
        HExchanger_DHW.mw2=Inputs.y[15];
        annotation (experiment(StopTime=31536000, Interval=3600));
      end HExchanger_DHW_V;

      model ThermalTank_V
        Plants.Components.ThermalTank ThermalTank(Vwt=1200);
      equation
        ThermalTank.mw1=150;
        ThermalTank.mw2=150;
        ThermalTank.Tws1=if time<8640000/2 then 40 else 80;
        ThermalTank.Twr2=20;
        annotation (experiment(StopTime=8640000, Interval=3600));
      end ThermalTank_V;

      model Subsation_V
        Plants.Components.Substation Substation;
      equation
        Substation.Tws2=90;
        Substation.Twr2=50;
        Substation.mw2=80;
        annotation (experiment(StopTime=86400, Interval=3600));
      end Subsation_V;

      model DataCenter_V
        Plants.Components.DataCenter DataCenter;
      equation
        DataCenter.qdc=1E6;
        DataCenter.Twidc=if time<86400/2 then 40 else 70;
        annotation (experiment(StopTime=86400, Interval=3600));
      end DataCenter_V;

      model BuildingAndSystem_V
        Plants.Components.BuildingAndSystem BuildingAndSystem;
        parameter Real Kpdhw=20 "Kp of PID of DHW system";
        Real Errdhw( start=0) "Error of PID of DHW system";
        parameter Real Kpra=2000 "Kp of PID of Radiator system";
        Real Erra( start=0) "Error of PID of Radiator system";
        parameter Real Kpahu=1000 "Kp of PID of AHU system";
        Real Errahu( start=0) "Error of PID of AHU system";
        Real qloadmeas "Measured heating load";
        Real Errqload( start=0) "Error of simulated heating load";
        Real Errqdhw( start=0) "Error of simulated heating load";
        Real Errqra( start=0) "Error of simulated heating load";
        Real Errqve( start=0) "Error of simulated heating load";
        Modelica.Blocks.Sources.CombiTimeTable Inputs(
          tableOnFile=true,
          table=fill(
              0.0,
              0,
              16),
          tableName="Data",
          fileName="C:/Work/Models/InputData/BuildingAndSystem_V.mat") "Time, 1 Toa, 2 Tiaref, 3 Tws1, 4 Twr1, 5 q_MS, 6 qDC, 7 qin, 8 qen, 9 qinf, 10 q_ve, 11 q_DHW, 12 mw1, 13 mainf, 14 maahu, 15 mwdhw
"         annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      equation
        BuildingAndSystem.Toa=Inputs.y[1];
        BuildingAndSystem.qin=Inputs.y[7];
        BuildingAndSystem.mwdhw=Inputs.y[15];
        BuildingAndSystem.maahu=Inputs.y[14];
        BuildingAndSystem.mainf=Inputs.y[13];
        BuildingAndSystem.Tws1=min(100,max(80,80-2*Inputs.y[1]));
        BuildingAndSystem.mw1sh=0.8*(BuildingAndSystem.mwra+BuildingAndSystem.mwahu);
        Errdhw=55-BuildingAndSystem.Twodhw;
        Erra=Inputs.y[2]-BuildingAndSystem.Tia;
        Errahu=Inputs.y[2]-BuildingAndSystem.Taoahu;
        BuildingAndSystem.mw1dhw=min(10,Kpdhw*max(1e-6,min(1,Errdhw)));
        BuildingAndSystem.mwra=min(250,Kpra*max(1e-6,min(1,Erra)));
        BuildingAndSystem.mwahu=min(150,Kpahu*max(1e-6,min(1,Errahu)));
        qloadmeas=Inputs.y[8]+Inputs.y[9]+Inputs.y[10]+Inputs.y[11]-Inputs.y[7];
        Errqload=qloadmeas-BuildingAndSystem.qload;
        Errqdhw=Inputs.y[11]-BuildingAndSystem.HExchanger_DHW.qhx_sum;
        Errqra=Inputs.y[8]+Inputs.y[9]-Inputs.y[7]-BuildingAndSystem.Building.qra;
        Errqve=Inputs.y[10]-BuildingAndSystem.Building.qve;
        annotation (experiment(StopTime=31536000, Interval=3600));
      end BuildingAndSystem_V;

      model BuildingLoadData_V
        Plants.Components.BuildingLoadData BuildingData;
        Real Error_Twr;
        Real Error_qdh;
        Modelica.Blocks.Sources.CombiTimeTable Inputs(
          tableOnFile=true,
          table=fill(
              0.0,
              0,
              6),
          tableName="Data",
          fileName="C:/Work/Models/InputData/BuildingData_V.mat")      "Time, 1 Ts, 2 qdh, 3 mw, 4 Toa, 5 Twr"   annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      equation
        BuildingData.Tws=Inputs.y[1];
        BuildingData.mw=2000*min(1,max(1E-5,Error_qdh/1e5));
        BuildingData.Toa=Inputs.y[4];
        Error_Twr=Inputs.y[5]-BuildingData.Twr;
        Error_qdh=Inputs.y[2]-BuildingData.qload;
        annotation (experiment(StopTime=31536000, Interval=3600));
      end BuildingLoadData_V;

      model HPData_V
        Plants.Components.HPData HPData;
        Modelica.Blocks.Sources.CombiTimeTable Inputs(
          tableOnFile=true,
          table=fill(
              0.0,
              0,
              2),
          tableName="Data",
          fileName="C:/Work/Models/InputData/HPData_V.mat")      "Time, 1 Tcoolin"   annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      equation
        HPData.Tcoolin=Inputs.y[1];
        annotation (experiment(StopTime=31536000, Interval=3600));
      end HPData_V;

      model PipeTDrop_V
        Plants.Components.PipeTDrop PipeTDrop;
      equation
        PipeTDrop.Tin=90;
        PipeTDrop.mw=20;
        annotation (experiment(StopTime=86400, Interval=3600));
      end PipeTDrop_V;

      model PipeHeatLoss_V
        Plants.Components.PipeHeatLoss Pipe;
      equation
        Pipe.Tins=80;
        Pipe.Tinr=60;
        Pipe.mw=10;
        annotation (experiment(StopTime=86400, Interval=3600));
      end PipeHeatLoss_V;
    end Validation;

  end Components;

  package Models
    partial model CampusDH
      // Note 1: BuildingAndSystem.mw1>DataCenter.mwdc
      Plants.Components.Substation MS;
      Plants.Components.ThermalTank WT(N=N,Twt_init=Twt_init);
      Plants.Components.BuildingAndSystem Bu(N=N,Ten_init=Ten_init,Tia_init=Tia_init,Tma_init=Tma_init,Tra_init=Tra_init,Ta_init=Ta_init,Tm_init=Tm_init,Tw_init=Tw_init);
      Plants.Components.DataCenter DC;
      // Global initial states
      parameter Integer N=10 "Number of sections for AHU, Radiator, and water tank";
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Ten_init=10;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tia_init=20;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tma_init=20;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tra_init[N]=fill(40,N);
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Ta_init[N]=fill(20,N);
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tm_init[N]=fill(40,N);
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tw_init[N]=fill(40,N);
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt_init[N]=fill(60,N);
      // Global parameters
      parameter Real Kp_DHW=20 "Kp of PID of DHW system";
      parameter Real Kp_Ra=2000 "Kp of PID of Radiator system";
      parameter Real Kp_AHU=1000 "Kp of PID of AHU system";
      parameter Modelica.SIunits.MassFlowRate m_DC=18.7 "Water Mass flow rate of data center (kg/s)";
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_DHW=55 "Seeting Temperature of DHW (degC)";
      // Global input variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Toa "Temperature of outdoor air (degC)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tiaref "Reference Temperature of indoor air (degC)";
      input Modelica.SIunits.HeatFlowRate qin "Internal heat gain internal heat gains due to lighting, equipment, and occupants (W)";
      input Modelica.SIunits.MassFlowRate mwdhw "Water Mass flow rate of tap water (kg/s)";
      input Modelica.SIunits.MassFlowRate maahu "Air Mass flow rate of AHU (kg/s)";
      input Modelica.SIunits.MassFlowRate mainf "Air Mass flow rate of Infiltration (kg/s)";
      input Modelica.SIunits.HeatFlowRate qdc "Heat feedin flow rate of data center (W)";
      // Manipulated variable
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MSs "Temperature of supply water in secondary side of substation (degC)";
      input Modelica.SIunits.MassFlowRate m_MS "Mass flow rate of water in secondary side of substation (kg/s)";
      input Modelica.SIunits.MassFlowRate m_SH "Mass flow rate of water in primary side of space heating system (kg/s)";
      // Global variables
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Err_DHW( start=0) "Error of PID of DHW system";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Err_Ra( start=0) "Error of PID of Radiator system";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Err_AHU( start=0) "Error of PID of AHU system";
      // Global output variables
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tia "Temperature of indoor air (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twodhw "Temperature of outlet water of domestic hot water (degC)";
      output Modelica.SIunits.HeatFlowRate qsub "Heat flow rate of substation (W)";
      output Modelica.SIunits.HeatFlowRate qloss "Heat loss flow rate from cooling tower (W)";
    equation
      // Level Substation
      MS.Tws2=T_MSs;
      MS.mw2=m_MS;
      MS.Twr2=WT.Twr1;
      // Level ThermalTank
      WT.Tws1=MS.Tws2;
      WT.mw1=MS.mw2;
      WT.Twr2*WT.mw2=Bu.Twr1*(Bu.mw1-DC.mwdc)+DC.Twodc*DC.mwdc;
      WT.mw2=Bu.mw1;
      // Level BuildingAndSystem
      Bu.Toa=Toa;
      Bu.qin=qin;
      Bu.mwdhw=mwdhw;
      Bu.maahu=maahu;
      Bu.mainf=mainf;
      Bu.Tws1=WT.Tws2;
      Bu.mw1sh=m_SH;
      Err_DHW=T_DHW-Bu.Twodhw;
      Err_Ra=Tiaref-Bu.Tia;
      Err_AHU=Tiaref-Bu.Taoahu;
      Bu.mw1dhw=min(10,Kp_DHW*max(1e-6,min(1,Err_DHW)));
      Bu.mwra=min(250,Kp_Ra*max(1e-6,min(1,Err_Ra)));
      Bu.mwahu=min(150,Kp_AHU*max(1e-6,min(1,Err_AHU)));
      // Level DataCenter
      DC.qdc=qdc;
      DC.Twidc=Bu.Twr1;
      // Level outputs
      Tia=Bu.Tia;
      Twodhw=Bu.Twodhw;
      qsub=MS.qhx;
      qloss=DC.qloss;
    end CampusDH;

    model CampusDHData
      Plants.Components.Substation MS;
      Plants.Components.ThermalTank WT(Vwt=Vwt);
      Plants.Components.BuildingLoadData Bu;
      Plants.Components.HPData DC;
      Plants.Components.PipeHeatLoss Pipe;
      // Global initial states
      parameter Modelica.SIunits.Volume Vwt "Volume of water tank thermal storage (m3)";
      // Global input variables
      input Modelica.SIunits.HeatFlowRate qload "Building heating load (W)";
      // Manipulated variable
      input Modelica.SIunits.MassFlowRate m_MS "Mass flow rate of water in secondary side of substation (kg/s)";
      input Modelica.SIunits.MassFlowRate m_Bu "Mass flow rate of water building (kg/s)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_s "Temperature of supply water in secondary side of substation (degC)";
      // Global variables
      Modelica.SIunits.HeatFlowRate qMS "Heat flow rate of substation (W)";
      Modelica.SIunits.Power Pdc "Power of data center (W)";
      Modelica.SIunits.HeatFlowRate qlosswt "Heat loss flow rate from water tank (W)";
      Modelica.SIunits.HeatFlowRate qlosspipe "Heat loss flow rate from pipes (W)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_r "Temperature of return water in secondary side of substation (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_Bu_s "Temperature of supply water of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_Bu_r "Temperature of return water of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_WT_s "Temperature of supply water of water tank at the user side (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_WT_r "Temperature of return water of water tank at the user side (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Dif "Temperature difference of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt[WT.N] "Temperature of water in water tank thermal storage of each section (degC)";
    equation
      // Level Substation
      MS.Tws2=T_MS_s;
      MS.mw2=m_MS;
      MS.Twr2=WT.Twr1;
      // Level ThermalTank
      WT.Tws1=MS.Tws2;
      WT.mw1=MS.mw2;
      WT.mw2=Pipe.mw;
      if Bu.mw>=DC.mwcond then
         Pipe.Tinr*Pipe.mw=Bu.Twr*(Bu.mw-DC.mwcond)+DC.Tcondout*DC.mwcond;
         DC.Tctin=Bu.Twr;
      else
         Pipe.Tinr=DC.Tcondout;
         DC.Tcondout*(DC.mwcond-Bu.mw)+Bu.Twr*Bu.mw=DC.Tctin*DC.mwcond;
      end if;
      WT.Twr2=Pipe.Toutr;
      // Level BuildingAndSystem
      Pipe.Tins=WT.Tws2;
      Bu.Tws=Pipe.Touts;
      m_Bu=Bu.mw;
      Pipe.mw=Bu.mw;
      // Level outputs
      qMS=MS.qhx;
      Pdc=DC.P;
      qload=Bu.qload;
      qlosswt=WT.qlosswt;
      qlosspipe=Pipe.qlosss+Pipe.qlossr;
      T_MS_r=MS.Twr2;
      T_Bu_s=Bu.Tws;
      T_Bu_r=Bu.Twr;
      Twt=WT.Twt;
      Dif=Bu.Dif;
      T_WT_r=WT.Twr2;
      T_WT_s=WT.Tws2;
    end CampusDHData;

    model CampusDHDataIF
      Plants.Components.Substation MS;
      Plants.Components.ThermalTank WT(Vwt=Vwt);
      Plants.Components.BuildingLoadData Bu;
      Plants.Components.HPData DC;
      Plants.Components.PipeTDrop PipeS;
      Plants.Components.PipeTDrop PipeR;
      // Global initial states
      parameter Modelica.SIunits.Volume Vwt "Volume of water tank thermal storage (m3)";
      // Global input variables
      input Modelica.SIunits.HeatFlowRate qload "Building heating load (W)";
      // Manipulated variable
      input Modelica.SIunits.MassFlowRate m_MS "Mass flow rate of water in secondary side of substation (kg/s)";
      input Modelica.SIunits.MassFlowRate m_Bu "Mass flow rate of water building (kg/s)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_s "Temperature of supply water in secondary side of substation (degC)";
      // Global variables
      Modelica.SIunits.HeatFlowRate qMS "Heat flow rate of substation (W)";
      Modelica.SIunits.Power Pdc "Power of data center (W)";
      Modelica.SIunits.HeatFlowRate qlosswt "Heat loss flow rate from water tank (W)";
      Modelica.SIunits.HeatFlowRate qlosspipe "Heat loss flow rate from pipes (W)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_r "Temperature of return water in secondary side of substation (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_Bu_s "Temperature of supply water of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_Bu_r "Temperature of return water of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_WT_s "Temperature of supply water of water tank at the user side (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_WT_r "Temperature of return water of water tank at the user side (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Dif "Temperature difference of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt[WT.N] "Temperature of water in water tank thermal storage of each section (degC)";
    equation
      // Level Substation
      MS.Tws2=T_MS_s;
      MS.mw2=m_MS;
      MS.Twr2=WT.Twr1;
      // Level ThermalTank
      WT.Tws1=MS.Tws2;
      WT.mw1=MS.mw2;
      WT.mw2=PipeR.mw;
      if Bu.mw>=DC.mwcond then
         PipeR.Tin*PipeR.mw=Bu.Twr*(Bu.mw-DC.mwcond)+DC.Tcondout*DC.mwcond;
         DC.Tcondin=Bu.Twr;
      else
         PipeR.Tin=DC.Tcondout;
         DC.Tcondout*(DC.mwcond-Bu.mw)+Bu.Twr*Bu.mw=DC.Tcondin*DC.mwcond;
      end if;
      PipeR.mw=Bu.mw;
      WT.Twr2=PipeR.Tout;
      // Level BuildingAndSystem
      PipeS.Tin=WT.Tws2;
      Bu.Tws=PipeS.Tout;
      m_Bu=Bu.mw;
      PipeS.mw=Bu.mw;
      // Level outputs
      qMS=MS.qhx;
      Pdc=DC.P;
      qload=Bu.qload;
      qlosswt=WT.qlosswt;
      qlosspipe=PipeS.qloss+PipeR.qloss;
      T_MS_r=MS.Twr2;
      T_Bu_s=Bu.Tws;
      T_Bu_r=Bu.Twr;
      Twt=WT.Twt;
      Dif=Bu.Dif;
      T_WT_r=WT.Twr2;
      T_WT_s=WT.Tws2;
    end CampusDHDataIF;

    model CampusDHDataSimple
      Plants.Components.Substation MS1;
      Plants.Components.Substation MS2;
      Plants.Components.ThermalTank WT(Vwt=Vwt);
      Plants.Components.Building Bu;
      Plants.Components.DataCenter DC;
      Plants.Components.PipeHeatLoss Pipe;
      // Global initial states
      parameter Modelica.SIunits.Volume Vwt "Volume of water tank thermal storage (m3)";
      // Global input variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Toa;
      input Modelica.SIunits.HeatFlowRate qin "Heat flow rate of internal heat gain (W)";
      input Modelica.SIunits.HeatFlowRate qinf "Heat flow rate of infiltration (W)";
      input Modelica.SIunits.HeatFlowRate qve "Heat flow rate of ventilation system (W)";
      input Modelica.SIunits.HeatFlowRate qdhw "Heat flow rate of DHW system (W)";
      // Manipulated variable
      input Modelica.SIunits.HeatFlowRate qdc "Heat feedin flow rate of data center (W)";
      input Modelica.SIunits.HeatFlowRate qra "Heat flow rate from Radiator to Indoor air (W)";
      input Modelica.SIunits.MassFlowRate m_MS1 "Mass flow rate of water in secondary side of the substation for charging (kg/s)";
      input Modelica.SIunits.MassFlowRate m_MS2 "Mass flow rate of water in secondary side of the substation for campus (kg/s)";
      input Modelica.SIunits.MassFlowRate m_Bu "Mass flow rate of water building (kg/s)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_s1 "Temperature of supply water in secondary side of the substation for charging (degC)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_s2 "Temperature of supply water in secondary side of the substation for campus (degC)";
      // Global variables
      Modelica.SIunits.HeatFlowRate qload "Heat flow rate of the building (W)";
      Modelica.SIunits.HeatFlowRate qMS "Heat flow rate of the substation for charging (W)";
      Modelica.SIunits.HeatFlowRate qMS1 "Heat flow rate of the substation for charging (W)";
      Modelica.SIunits.HeatFlowRate qMS2 "Heat flow rate of the substation for campus (W)";
      Modelica.SIunits.HeatFlowRate qlosswt "Heat loss flow rate from water tank (W)";
      Modelica.SIunits.HeatFlowRate qlosspipe "Heat loss flow rate from pipes (W)";
      Modelica.SIunits.Heat Qstore "Stored heat (J)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tia "Temperature of indoor air (Centigrade)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_r1 "Temperature of return water in secondary side of the substation for charging (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_r2 "Temperature of return water in secondary side of the substation for campus (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_Bu_s "Temperature of supply water of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_Bu_r "Temperature of return water of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_WT_s "Temperature of supply water of water tank at the user side (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_WT_r "Temperature of return water of water tank at the user side (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Dif "Temperature difference of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt[WT.N] "Temperature of water in water tank thermal storage of each section (degC)";
    equation
      // Level Substation 1
      MS1.Tws2=T_MS_s1;
      MS1.mw2=m_MS1;
      MS1.Twr2=WT.Twr1;
      // Level Substation 2
      MS2.Tws2=T_MS_s2;
      MS2.mw2=m_MS2;
      MS2.Twr2=Pipe.Toutr;
      // Level ThermalTank
      WT.Tws1=MS1.Tws2;
      WT.mw1=MS1.mw2;
      WT.mw2=Pipe.mw-MS2.mw2;
      Pipe.Tinr*Pipe.mw=Bu.Twr*(Bu.mw-DC.mwdc)+DC.Twodc*DC.mwdc;
      DC.Twidc=Bu.Twr;
      WT.Twr2=Pipe.Toutr;
      // Level BuildingAndSystem
      Pipe.Tins*Pipe.mw=WT.Tws2*WT.mw2+MS2.Tws2*MS2.mw2;
      Bu.Tws=Pipe.Touts;
      m_Bu=Bu.mw;
      Pipe.mw=Bu.mw;
      // Level outputs
      qMS=qMS1+qMS2;
      qMS1=MS1.qhx;
      qMS2=MS2.qhx;
      Toa=Bu.Toa;
      qin=Bu.qin;
      qinf=Bu.qinf;
      qve=Bu.qve;
      qdhw=Bu.qdhw;
      qload=Bu.qload;
      Tia=Bu.Tia;
      qra=Bu.qra;
      qdc=DC.qdc;
      qlosswt=WT.qlosswt;
      qlosspipe=Pipe.qlosss+Pipe.qlossr;
      Qstore=WT.Qstore;
      T_MS_r1=MS1.Twr2;
      T_MS_r2=MS2.Twr2;
      T_Bu_s=Bu.Tws;
      T_Bu_r=Bu.Twr;
      Twt=WT.Twt;
      Dif=Bu.Dif;
      T_WT_r=WT.Twr2;
      T_WT_s=WT.Tws2;
    end CampusDHDataSimple;

    model CampusDHDataSimpleIF
      Plants.Components.Substation MS;
      Plants.Components.ThermalTank WT(Vwt=Vwt);
      Plants.Components.BuildingLoadData Bu;
      Plants.Components.DataCenter DC;
      Plants.Components.PipeHeatLoss Pipe;
      // Global initial states
      parameter Modelica.SIunits.Volume Vwt "Volume of water tank thermal storage (m3)";
      // Global input variables
      input Modelica.SIunits.HeatFlowRate qload "Building heating load (W)";
      input Modelica.SIunits.HeatFlowRate qdc "Heat feedin flow rate of data center (W)";
      // Manipulated variable
      input Modelica.SIunits.MassFlowRate m_MS "Mass flow rate of water in secondary side of substation (kg/s)";
      input Modelica.SIunits.MassFlowRate m_Bu "Mass flow rate of water building (kg/s)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_s "Temperature of supply water in secondary side of substation (degC)";
      // Global variables
      Modelica.SIunits.HeatFlowRate qMS "Heat flow rate of substation (W)";
      Modelica.SIunits.HeatFlowRate qlosswt "Heat loss flow rate from water tank (W)";
      Modelica.SIunits.HeatFlowRate qlosspipe "Heat loss flow rate from pipes (W)";
      Modelica.SIunits.HeatFlowRate qlossct "Heat loss flow rate from cooling tower (W)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_r "Temperature of return water in secondary side of substation (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_Bu_s "Temperature of supply water of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_Bu_r "Temperature of return water of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_WT_s "Temperature of supply water of water tank at the user side (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_WT_r "Temperature of return water of water tank at the user side (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Dif "Temperature difference of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt[WT.N] "Temperature of water in water tank thermal storage of each section (degC)";
    equation
      // Level Substation
      MS.Tws2=T_MS_s;
      MS.mw2=m_MS;
      MS.Twr2=WT.Twr1;
      // Level ThermalTank
      WT.Tws1=MS.Tws2;
      WT.mw1=MS.mw2;
      WT.mw2=Pipe.mw;
      if Bu.mw>=DC.mwdc then
         Pipe.Tinr*Pipe.mw=Bu.Twr*(Bu.mw-DC.mwdc)+DC.Twodc*DC.mwdc;
         DC.Twidc=Bu.Twr;
      else
         Pipe.Tinr=DC.Twodc;
         DC.Twodc*(DC.mwdc-Bu.mw)+Bu.Twr*Bu.mw=DC.Twidc*DC.mwdc;
      end if;
      WT.Twr2=Pipe.Toutr;
      // Level BuildingAndSystem
      Pipe.Tins=WT.Tws2;
      Bu.Tws=Pipe.Touts;
      m_Bu=Bu.mw;
      Pipe.mw=Bu.mw;
      // Level outputs
      qMS=MS.qhx;
      qload=Bu.qload;
      qdc=DC.qdc;
      qlosswt=WT.qlosswt;
      qlossct=DC.qloss;
      qlosspipe=Pipe.qlosss+Pipe.qlossr;
      T_MS_r=MS.Twr2;
      T_Bu_s=Bu.Tws;
      T_Bu_r=Bu.Twr;
      Twt=WT.Twt;
      Dif=Bu.Dif;
      T_WT_r=WT.Twr2;
      T_WT_s=WT.Tws2;
    end CampusDHDataSimpleIF;

    model Simu
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tiaref;
      extends Plants.Models.CampusDHDataSimple;
    end Simu;

    package Validation
      model CampusDH_V
        Plants.Models.CampusDH CampusDH;
        Modelica.Blocks.Sources.CombiTimeTable Inputs(
          tableOnFile=true,
          table=fill(
              0.0,
              0,
              16),
          tableName="Data",
          fileName="C:/Work/Models/InputData/BuildingAndSystem_V.mat") "Time, 1 Toa, 2 Tiaref, 3 Tws1, 4 Twr1, 5 q_MS, 6 qDC, 7 qin, 8 qen, 9 qinf, 10 q_ve, 11 q_DHW, 12 mw1, 13 mainf, 14 maahu, 15 mwdhw
"         annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      equation
        CampusDH.Toa=Inputs.y[1];
        CampusDH.Tiaref=Inputs.y[2];
        CampusDH.qin=Inputs.y[7];
        CampusDH.mwdhw=Inputs.y[15];
        CampusDH.maahu=Inputs.y[14];
        CampusDH.mainf=Inputs.y[13];
        CampusDH.qdc=Inputs.y[6];
        CampusDH.m_MS=200;
        CampusDH.T_MSs=90;
        CampusDH.m_SH=200;
        annotation (experiment(StopTime=31536000, Interval=3600));
      end CampusDH_V;

      model CampusDHData_V
        Plants.Models.CampusDHData CampusDHData(Vwt=12000);
      equation
        CampusDHData.qload=4e6;
        CampusDHData.T_MS_s=90;
        CampusDHData.m_MS=50;
        CampusDHData.m_Bu=100;
        annotation (experiment(StopTime=31536000, Interval=3600));
      end CampusDHData_V;

      model Simu_V
        Plants.Models.Simu Simu(Vwt=12000);
      equation
        Simu.Toa=0;
        Simu.qload=12E6;
        Simu.qdc=1E6;
        Simu.m_MS1=80;
        Simu.m_MS2=20;
        Simu.m_Bu=100;
        Simu.T_MS_s1=80;
        Simu.T_MS_s2=80;
        annotation (experiment(StopTime=31536000, Interval=3600));
      end Simu_V;
    end Validation;

    package Simulations
      model Simulations
        Plants.Models.Simulations.Plant Plant;
        parameter Real Kp=2000 "Kp of PID";
        parameter Modelica.SIunits.Time T=60 "Time Constant";
        parameter Real k(unit="1")=1 "Gain";
        //Real Err(start=1) "Error of PID";
        //Real u "Connector of Real input signal";
        //Real y(start=50) "Connector of Real output signal";
        Modelica.Blocks.Sources.CombiTimeTable Inputs(
          tableOnFile=true,
          table=fill(
              0.0,
              0,
              10),
          tableName="Data",
          fileName="C:/Work/Models/Outputs/InputsForDymola.mat")      "Time, 1 Toa, 2 qload_ref, 3 m_MS[0], 4 m_MS[1], 5 m_MS[2], 6 m_MS[3], 7 m_MS[4], 8 m_MS[5], 9 m_MS[6]"   annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      equation
        Plant.Toa=Inputs.y[1];
        Plant.qload_ref=Inputs.y[2];
        Plant.m_MS=Inputs.y[9];
        annotation (experiment(StopTime=31536000, Interval=3600));
      end Simulations;

      model Plant
        Plants.Models.Simu Plant;
        Real dif_load;
        parameter Real Nom_dif_load=1E4;
        Modelica.Blocks.Sources.CombiTimeTable Inputs(
          tableOnFile=true,
          table=fill(
              0.0,
              0,
              10),
          tableName="Data",
          fileName="C:/Work/Models/Outputs/InputsForDymola.mat")      "Time, 1 Toa, 2 qload_ref, 3 m_MS[0], 4 m_MS[1], 5 m_MS[2], 6 m_MS[3], 7 m_MS[4], 8 m_MS[5], 9 m_MS[6]"   annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      equation
        Plant.Toa=Inputs.y[1];
        Plant.qload_ref=Inputs.y[2];
        Plant.m_MS=Inputs.y[9];
        dif_load=Plant.qload_ref-Plant.qload;
        Plant.m_Bu=2000*min(1,max(1E-2,dif_load/Nom_dif_load));
        annotation (experiment(StopTime=31536000, Interval=3600));
      end Plant;
    end Simulations;
  end Models;

  annotation (uses(Modelica(version="3.2.3")));
end Plants;
