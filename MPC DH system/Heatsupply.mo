within ;
package Heatsupply
  package Components
    model DataCentreHP
      // Parameters
      parameter Real a=-1.5417;
      parameter Real b=10.5450;
      parameter Real c=-16.9801;
      parameter Real d=24.4038;
      parameter Real e=-5.7424;
      parameter Real f=-9.3203;
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Modelica.SIunits.MassFlowRate mwevap=36.5 "Mass flow rate of water of Evaporator (kg/s)";
      parameter Modelica.SIunits.MassFlowRate mwcond=18.7 "Mass flow rate of water of Condenser (kg/s)";
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tevapin=11 "Inlet water temperature of Evaporator (degC)";

      // Variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tre "Return temperature from campus district heating system (degC)";
      input Modelica.SIunits.Power P "Power of HP (W)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tevapout "Outlet water temperature of Evaporator (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tcondout "Outlet water temperature of Condenser (degC)";
      output Modelica.SIunits.HeatFlowRate qloss "Heat flow rate loss of Cooling tower (W)";
      output Modelica.SIunits.HeatFlowRate qevap "Heat flow rate of Evaporator (W)";
      output Modelica.SIunits.HeatFlowRate qcond "Heat flow rate of Condenser (W)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tcondin "Inlet water temperature of Condenser (degC)";

    equation
      qloss=max(0,cw*mwcond*(Tre-60));
      Tcondin=min(60,Tre);
      P=1000*(a*Tevapin+b*Tevapout+c*Tcondin+d*Tcondout+e*mwevap+f*mwcond);
      P=qcond-qevap;
      qevap=cw*mwevap*(Tevapin-Tevapout);
      qcond=cw*mwcond*(Tcondout-Tcondin);

    end DataCentreHP;

    model MainSubstation
      // Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      // Variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws2 "Temperature of supply water in secondary side (degC)";
      input Modelica.SIunits.MassFlowRate mw2 "Mass flow rate of water in secondary side (kg/s)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr2 "Temperature of return water in secondary side (degC)";
      output Modelica.SIunits.HeatFlowRate qhx "Heat flow rate of heat exchanger (W)";
    equation
      qhx=cw*mw2*(Tws2-Twr2);
    end MainSubstation;

    model WaterTank
      //Initial states
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt1_init=52;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt2_init=54;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt3_init=56;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt4_init=58;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt5_init=60;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt6_init=62;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt7_init=64;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt8_init=66;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt9_init=68;
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt10_init=70;
      //Parameters
      parameter Integer N=10 "Number of nodes used in the spatial approximation";
      parameter Modelica.SIunits.Volume Vwt=900 "Volume of water tank thermal storage (m3)";
      parameter Modelica.SIunits.Height Hwt=16 "Height of water tank theraml storage (m)";
      parameter Modelica.SIunits.Diameter Dwt=2*(Vwt/3.14/Hwt)^0.5 "Diameter of water tank thermal storage (m)";
      parameter Modelica.SIunits.Area Awt1=3.14*Dwt*Hwt "Surface area except top and down of water tank (m2)";
      parameter Modelica.SIunits.Area Awt2=3.14*(Dwt/2)^2 "Top or down surface area of water tank (m2)";
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Modelica.SIunits.Density rho=995.6  "Density of water (kg/m3)";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer Uwt=1.2 "U value between water tank thermal storage and ground (W/(m2.K))";
      parameter Modelica.SIunits.ThermalConductivity keff=4.1 "Effective Thermal Conductivity of water tank thermal storage between each section (W/(m.K))";
      parameter Modelica.SIunits.HeatCapacity Cwt=cw*Vwt*rho "Heat capacity of water in water tank (= cp*m)(J/K)";
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tgrou=5 "Temperature of ground (degC)";
      //Connection variables
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws1 "Water supply tmeperature of heat source side (degC)";
      input Modelica.SIunits.MassFlowRate mw1 "Mass flow rate of water of heat source side (kg/s)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr2 "Water return tmeperature of heat user side (degC)";
      input Modelica.SIunits.MassFlowRate mw2 "Mass flow rate of water of heat user side (kg/s)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws2 "Water supply tmeperature of heat user side (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr1 "Water return tmeperature of heat source side (degC)";
      //State Variables
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt[N]( nominal=fill(Twt1_init,N)) "Temperature of water in water tank thermal storage of each section (degC)";
      //Variables
      Modelica.SIunits.HeatFlowRate qloss[N] "Heat loss flow rate from water tank thermal storage (each section) to ground (W)";
      Modelica.SIunits.HeatFlowRate qlosswt "Heat loss flow rate from water tank thermal storage to ground (W)";
      Modelica.SIunits.HeatFlowRate qwt "Heat flow rate from water tank (W)";
      Modelica.SIunits.ThermalConductivity kwt "ThermalConductivity of water tank thermal storage between each section considers inversion effect(W/(m.K))";
      Modelica.SIunits.Heat Qstore "Stored heat (J)";
    initial equation
      Twt[1]=Twt1_init;
      Twt[2]=Twt2_init;
      Twt[3]=Twt3_init;
      Twt[4]=Twt4_init;
      Twt[5]=Twt5_init;
      Twt[6]=Twt6_init;
      Twt[7]=Twt7_init;
      Twt[8]=Twt8_init;
      Twt[9]=Twt9_init;
      Twt[10]=Twt10_init;
    equation
      Qstore=Cwt*(sum(Twt)/N-(Twt1_init+Twt2_init+Twt3_init+Twt4_init+Twt5_init+Twt6_init+Twt7_init+Twt8_init+Twt9_init+Twt10_init)/N);
      qlosswt=sum(qloss);
      qloss[1]=Uwt*(Awt1/N+Awt2)*(Twt[1]-Tgrou);
      qwt=cw*mw2*(Tws2-Twr2);
      kwt=keff*(1+exp(5*(Twt[1]-Twt[end])));
      Cwt/N*der(Twt[1])=cw*mw1*(Twt[2]-Twt[1])+cw*mw2*(Twr2-Twt[1])-qloss[1]+kwt*Awt2*(Twt[2]-Twt[1])/(Hwt/N);
      for i in 2:(N-1) loop
        qloss[i]=Uwt*Awt1/N*(Twt[i]-Tgrou);
        Cwt/N*der(Twt[i])=cw*mw1*(Twt[i+1]-Twt[i])+cw*mw2*(Twt[i-1]-Twt[i])-qloss[i]+kwt*Awt2/(Hwt/N)*(Twt[i+1]-Twt[i])+kwt*Awt2/(Hwt/N)*(Twt[i-1]-Twt[i]);
      end for;
      qloss[end]=Uwt*Awt1/N*(Twt[end]-Tgrou);
      Cwt/N*der(Twt[end])=cw*mw1*(Tws1-Twt[end])+cw*mw2*(Twt[end-1]-Twt[end])-qloss[end]+kwt*Awt2*(Twt[end-1]-Twt[end])/(Hwt/N);
      Tws2=Twt[end];
      Twr1=Twt[1];

    end WaterTank;

    model Building
      //Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      //Connection variables
      input Modelica.SIunits.HeatFlowRate qload "Heat supply flow rate of DH system (W)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tws "Supply water temperature (degC)";
      input Modelica.SIunits.MassFlowRate mw "Mass flow rate of water (kg/s)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Dif "Temperature difference (degC)";
      output Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twr "Temperature of return water (degC)";
    equation
      Dif=Tws-Twr;
      qload=cw*mw*Dif;
    end Building;

    model PipeHeatLoss
      //Parameters
      parameter Modelica.SIunits.SpecificHeatCapacity cw=4184 "Specific heat capacity at constant volume of water (J/(kg.K))";
      parameter Modelica.SIunits.ThermalConductivity Kgrou=1.5 "ThermalConductivity of ground (W/(m.K))";
      parameter Modelica.SIunits.ThermalConductivity Kinsu=0.03 "ThermalConductivity of insulation (W/(m.K))";
      parameter Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tgrou=5 "Temperature of ground (degC)";
      parameter Modelica.SIunits.Length L=1500*10 "Length Of pipe (m)";
      parameter Modelica.SIunits.Diameter d=0.273 "Outer Diameter of pipe (m)";
      parameter Modelica.SIunits.Diameter D=d+0.05*2 "Outer Diameter of insulation (m)";
      parameter Modelica.SIunits.Length s=1.2 "Distance between two pipes (m)";
      parameter Modelica.SIunits.Length h=1.2 "Distance between cnetral pipe and ground surface (m)";
      parameter Modelica.SIunits.ThermalInsulance Ri=(d/2/Kinsu)*log(D/d) "The thermal resistance of insulation (m2.K/W)";
      parameter Modelica.SIunits.ThermalInsulance Rg=(d/2/Kgrou)*log(4*h/d) "The thermal resistance of ground (m2.K/W)";
      parameter Modelica.SIunits.ThermalInsulance Rc=(d/2/Kgrou)*log(((2*h/s)^2+1)^0.5) "The thermal resistance of coinciding (m2.K/W)";
      //Variables
      Modelica.SIunits.HeatFlowRate qlosss "Heat loss flow rate from supply pipe (W)";
      Modelica.SIunits.HeatFlowRate qlossr "Heat loss flow rate from return pipe (W)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tdifs "Temperature difference between supply pipe and ground (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Tdifr "Temperature difference between return pipe and ground (degC)";
      //Connection variables
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

    model Circulationpump
      //Parameters
      parameter Modelica.SIunits.Density rho=995.6  "Density of water (kg/m3)";
      parameter Real S=48437500 "resistance friction coefficient of pipeline network (Pa/(m3/s)^2)";
      parameter Real effi=0.7 "efficiency of circulation pump";
      //Variables
      Modelica.SIunits.VolumeFlowRate G "Volume flow rate of circulation pump (m3/s)";
      Modelica.SIunits.Pressure H "hydraulic resistance of the pipeline network (Pa)";
      //Connection Variables
      input Modelica.SIunits.MassFlowRate mw "Mass flow rate of circulation pump (kg/s)";
      //Output variables
      output Modelica.SIunits.Power P "Power of circulation pump (W)";
    equation
      G=mw/rho;
      H=S*G^2;
      P=H*G/effi;

    end Circulationpump;
  end Components;

  package Model

    model CampusDH
      Heatsupply.Components.MainSubstation MS1;
      Heatsupply.Components.MainSubstation MS2;
      Heatsupply.Components.WaterTank WT;
      Heatsupply.Components.Building Bu;
      Heatsupply.Components.DataCentreHP DC;
      Heatsupply.Components.PipeHeatLoss Pipe;
      Heatsupply.Components.Circulationpump Pump;
      //Global input variables
      input Modelica.SIunits.HeatFlowRate qload "Building heating load (W)";
      //Manipulated variable
      input Modelica.SIunits.MassFlowRate m_MS1 "Mass flow rate of water in secondary side of the substation for charging (kg/s)";
      input Modelica.SIunits.MassFlowRate m_MS2 "Mass flow rate of water in secondary side of the substation for campus (kg/s)";
      input Modelica.SIunits.MassFlowRate m_WT "Mass flow rate of water tank in user side (kg/s)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_s1 "Temperature of supply water in secondary side of the substation for charging (degC)";
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_s2 "Temperature of supply water in secondary side of the substation for campus (degC)";
      input Modelica.SIunits.Power P_DC "Power of heat pump in data centre";
      //Global variables
      Modelica.SIunits.Power P_pump "Power of circulation pump";
      Modelica.SIunits.HeatFlowRate qMS "Heat flow rate of the main substation  (W)";
      Modelica.SIunits.HeatFlowRate qMS1 "Heat flow rate of the substation for charging (W)";
      Modelica.SIunits.HeatFlowRate qMS2 "Heat flow rate of the substation for campus (W)";
      Modelica.SIunits.HeatFlowRate qwt "Heat flow rate of the water tank for campus (W)";
      Modelica.SIunits.HeatFlowRate qlosswt "Heat loss flow rate from water tank (W)";
      Modelica.SIunits.HeatFlowRate qlosspipe "Heat loss flow rate from pipes (W)";
      Modelica.SIunits.Heat Qstore "Stored heat (J)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_r1 "Temperature of return water in secondary side of the substation for charging (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_MS_r2 "Temperature of return water in secondary side of the substation for campus (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_Bu_s "Temperature of supply water of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_Bu_r "Temperature of return water of building (degC)";
      Modelica.SIunits.MassFlowRate m_Bu "Mass flow rate of building (kg/s)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_WT_s "Temperature of supply water of water tank at the user side (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_WT_r "Temperature of return water of water tank at the user side (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_evapout "Outlet temperature of evaprator in data centre (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Dif "Temperature difference of building (degC)";
      Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Twt[WT.N] "Temperature of water in water tank thermal storage of each section (degC)";
    equation
      //Level Substation 1
      MS1.Tws2=T_MS_s1;
      MS1.mw2=m_MS1;
      MS1.Twr2=WT.Twr1;
      //Level Substation 2
      MS2.Tws2=T_MS_s2;
      MS2.mw2=m_MS2;
      MS2.Twr2=Pipe.Toutr;
      Pump.mw=m_Bu;
      //Level Water Tank
      WT.Tws1=MS1.Tws2;
      WT.mw1=MS1.mw2;
      WT.mw2=m_WT;
      WT.qwt=qwt;
      WT.mw2=Pipe.mw-MS2.mw2;
      Pipe.Tinr*Pipe.mw=Bu.Twr*(Bu.mw-DC.mwcond)+DC.Tcondout*DC.mwcond;
      DC.Tre=Bu.Twr;
      WT.Twr2=Pipe.Toutr;
      //Level Building
      m_Bu=Bu.mw;
      Pipe.Tins*Pipe.mw=WT.Tws2*WT.mw2+MS2.Tws2*MS2.mw2;
      Bu.Tws=Pipe.Touts;
      Pipe.mw=Bu.mw;
      //Level outputs
      qMS=qMS1+qMS2;
      qMS1=MS1.qhx;
      qMS2=MS2.qhx;
      qload=Bu.qload;
      P_DC=DC.P;
      T_evapout=DC.Tevapout;
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
      P_pump=Pump.P;

    end CampusDH;

    model Simu
      input Modelica.SIunits.Conversions.NonSIunits.Temperature_degC Toa;
      extends Heatsupply.Model.CampusDH;
    end Simu;
  end Model;

  package Simulation

    model Simulation
      Heatsupply.Model.Simu Simu;
      //PI control
      parameter Real Tevapout_ref=6.5 "Reference value for the outlet water temperature of evaprator";
      parameter Real k1=1E6 "Gain of controller";
      parameter Real Ti1=1E3 "Time constant of integrator block (s)";
      parameter Real k2=1E3 "Gain of controller";
      parameter Real Ti2=1E3 "Time constant of integrator block (s)";
      Real dif_Tevapout "Temperature difference between the outlet water temperature of evaprator and its reference value";
      Real x "Derivative";
      Real T_Bu_r_ref "Reference signal of return temperature of builiding";
      Real Tdif_Bu_r "Temperature difference between return temperature of building and its reference value";
      Real y "Derivative";
    equation
      Simu.Toa=1.2;
      Simu.qload=5759215;
      Simu.m_MS1=0;
      Simu.m_WT=40;
      Simu.T_MS_s1=80;
      Simu.T_MS_s2=68;
      dif_Tevapout=Simu.T_evapout-Tevapout_ref;
      der(x)=dif_Tevapout;
      T_Bu_r_ref=Simu.T_Bu_s-(0.7017*Simu.T_Bu_s-27.1661);
      Tdif_Bu_r=T_Bu_r_ref-Simu.T_Bu_r;
      der(y)=Tdif_Bu_r;
      Simu.P_DC=min(max(k1*dif_Tevapout+Ti1*x,100000),600000);
      Simu.m_MS2=min(60,max(20,k2*Tdif_Bu_r+Ti2*y));
    end Simulation;
  end Simulation;
  annotation (uses(Modelica(version="3.2.2")));
end Heatsupply;
