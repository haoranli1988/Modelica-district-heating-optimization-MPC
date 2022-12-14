##########################################################################################  
### 1 Import Python and JModelica.org
# 1.1 Import from Python
import os.path
import csv
import numpy as N
import matplotlib.pyplot as plt
import pandas as pd
from collections import OrderedDict
# 1.2 Import from JModelica.org
from pymodelica import compile_fmu
from pyfmi import load_fmu
from pyjmi import transfer_optimization_problem, get_files_path
from pyjmi.optimization.casadi_collocation import BlockingFactors
from pyjmi.symbolic_elimination import BLTOptimizationProblem, EliminationOptions
from pyfmi.common.plotting import plot_gui
from pyjmi.optimization.casadi_collocation import ExternalData

##########################################################################################  
### 2 Global Setting
# 2.1  MPC settting
Sim_period=600                              # time step for optimization window shifting (s)
Opt_period=3600                             # Resolution for control action (s)
horizon=12                                  # Optimization horizon (hour) 
days=1                                      # MPC control period (1-7 dags)
predict_time=Opt_period*horizon             # Optimization horizon (s)    
final_time=days*24*3600                     # MPC control period (s) 
n_e_per_sample=1                            # Collocation elements / sample
n_e=n_e_per_sample*horizon                  # Total collocation elements
number_samp_tot=int(final_time/Sim_period)  # Total number of samples to do
# 2.2 Input varibels from 14/12/2016  00:00:00 Index 384  to 21/12/2016  00:00:00
Index_start=384
csv_data_Toa = pd.read_csv("C:/GitHub/MPC/Measurement data.csv")
csv_data_Toa_fore = pd.read_csv("C:/GitHub/MPC/Forecast data.csv")
csv_data_Qin = pd.read_csv("C:/GitHub/MPC/Internal heat gain.csv")
csv_data_Hv = pd.read_csv("C:/GitHub/MPC/Mechanical ventilation heat transfer coefficient.csv")
csv_data_Tiaref = pd.read_csv("C:/GitHub/MPC/Indoor air temperature reference.csv")
csv_data_P_EDC = pd.read_csv("C:/GitHub/MPC/Heat price energy demand.csv")
Timesamples = csv_data_Toa.Time[Index_start:-1]-csv_data_Toa.Time[Index_start]
Toa=csv_data_Toa.T_meas[Index_start:-1]
Toa_fore=csv_data_Toa_fore.T_fore[Index_start:-1]
Qin=csv_data_Qin.Q_in[Index_start:-1]
Hv=csv_data_Hv.Hv[Index_start:-1]
Tiaref=csv_data_Tiaref.Tiaref[Index_start:-1]
P_EDC=csv_data_P_EDC.P_EDC[Index_start:-1]
Timesamples.reset_index(drop=True, inplace=True)
Toa.reset_index(drop=True, inplace=True)
Toa_fore.reset_index(drop=True, inplace=True)
Qin.reset_index(drop=True, inplace=True)
Hv.reset_index(drop=True, inplace=True)
Tiaref.reset_index(drop=True, inplace=True)
P_EDC.reset_index(drop=True, inplace=True)

##########################################################################################  
### 2  Define MPC problem
# 2.1 Initial guess trajectories
# 2.1.1 Get simulation model
init_fmu = compile_fmu("Building.Plant","C:/GitHub/MPC/OptModel.mop", compiler_options={'generate_html_diagnostics': True})
init_model = load_fmu(init_fmu)
# 2.1.2 Set manipulated variables
init_model.set('Ts',50)
init_model.set('mr',0.5)
# 2.1.3 Set input variables
init_model.set('Toa',-2.4)
init_model.set('Qin',91000)
init_model.set('Hv',1000)
# 2.1.4 Set initial states
Ten_ini=10
Tia_ini=20
Tma_ini=20
Tra1_ini=50
Tra2_ini=45
Tra3_ini=40
Tra4_ini=35
Tra5_ini=30
init_model.set('Ten_init', Ten_ini)
init_model.set('Tia_init', Tia_ini)
init_model.set('Tma_init', Tma_ini)
init_model.set('Tra1_init', Tra1_ini)
init_model.set('Tra2_init', Tra2_ini)
init_model.set('Tra3_init', Tra3_ini)
init_model.set('Tra4_init', Tra4_ini)
init_model.set('Tra5_init', Tra5_ini)
# 2.1.5 Compute initial guess trajectories
init_res = init_model.simulate(start_time=0., final_time=predict_time)

# 2.2  Define simulation problem
sim_fmu = compile_fmu("Building.Plant","C:/GitHub/MPC/OptModel.mop", compiler_options={'generate_html_diagnostics': True})
sim_model = load_fmu(sim_fmu)

# 2.3 Define the optimal control problem
# 2.3.1 Compile and load optimization problem
op = transfer_optimization_problem("Building.BuildingMPC","C:/GitHub/MPC/OptModel.mop", compiler_options={'generate_html_diagnostics': True, 'equation_sorting': True})
# 2.3.2 Create blocking factors with quadratic penalty and bound on 'T_ras' and 'm_ra'
bf_list = [1]*horizon
factors = {'Ts': bf_list, 'mr': bf_list}
du_quad_pen = {'Ts': 100, 'mr': 100}
du_bounds = {'Ts': 5, 'mr': 0.5}
bf = BlockingFactors(factors,du_bounds,du_quad_pen)
# 2.3.3 Set collocation options
opt_opts = op.optimize_options()
opt_opts['n_e'] = n_e
opt_opts['n_cp'] = 1
opt_opts['init_traj'] = init_res
opt_opts['nominal_traj'] = init_res
#opt_opts['blocking_factors'] = bf

# 2.4 Lists for logging 
# 2.4.1 Global variable
t_MPC = []
# 2.4.2 Manipulated variable
Ts_MPC = []
mr_MPC = []
# 2.4.3 Input variables
Toa_MPC = []
Qin_MPC = []
Hv_MPC = []
# 2.4.4 Variable
Trs_MPC = []
Trr_MPC = []
ma_MPC = []
Qload_MPC = []
Qsh_MPC = []
Qra_MPC = []
Qve_MPC = []
Qinf_MPC = []
Qen_MPC = []
Qmaen_MPC = []
Qma_MPC = []
# 2.4.5 State variables
Ten_MPC = []
Tia_MPC = []
Tma_MPC = []
Tra1_MPC = []
Tra2_MPC = []
Tra3_MPC = []
Tra4_MPC = []
Tra5_MPC = []
# 2.4.6 Refrence for MPC
Tiaref_MPC = []
P_EDC_MPC = []

# 2.5 MPC loops
for k in range(number_samp_tot):
    # 2.5.1 Prepare input for optimization
    t0 = k*Sim_period
    tf = k*Sim_period+predict_time  
    data_input_opt_select = OrderedDict()
    data_input_opt_select['Toa'] = N.vstack([Timesamples[0:horizon+1],Toa[int(t0/Opt_period):int(tf/Opt_period)+1]])
    data_input_opt_select['Qin'] = N.vstack([Timesamples[0:horizon+1],Qin[int(t0/Opt_period):int(tf/Opt_period)+1]])
    data_input_opt_select['Hv'] = N.vstack([Timesamples[0:horizon+1],Hv[int(t0/Opt_period):int(tf/Opt_period)+1]])
    data_input_opt_select['Tiaref'] = N.vstack([Timesamples[0:horizon+1],Tiaref[int(t0/Opt_period):int(tf/Opt_period)+1]])
    data_input_opt_select['P_EDC'] = N.vstack([Timesamples[0:horizon+1],P_EDC[int(t0/Opt_period):int(tf/Opt_period)+1]])
    opt_opts['external_data'] = ExternalData(eliminated=data_input_opt_select)
    # 2.5.2 Solve optimization problem
    res = op.optimize(options=opt_opts)
    # 2.5.3 Get optimal Manipulated variable
    Ts = res['Ts'][0]
    mr = res['mr'][0]
    # 2.5.4 Reset the model and set the new initial states before simulating
    sim_model.set('Ts',float(Ts))
    sim_model.set('mr',float(mr))  
    data_sim = N.vstack(( Timesamples[int(t0/Opt_period)], Toa[int(t0/Opt_period)], Qin[int(t0/Opt_period)], Hv[int(t0/Opt_period)]))
    data_input_sim = (['Toa', 'Qin', 'Hv'],N.transpose(data_sim))
    # 2.5.5 Simulate using optimal inputs
    sim_res = sim_model.simulate(start_time= k*Sim_period, final_time=(k+1)*Sim_period, input=data_input_sim)
    # 2.5.6 Logging
    # 2.5.6.1 Global variable
    t_MPC.append(sim_res['time'][-1])
    # 2.5.6.2 Manipulated variable
    mr_MPC.append(sim_res['mr'][-1])
    Ts_MPC.append(sim_res['Ts'][-1])
    # 2.5.6.3 Input variables
    Toa_MPC.append(sim_res['Toa'][-1])
    Qin_MPC.append(sim_res['Qin'][-1])
    Hv_MPC.append(sim_res['Hv'][-1])
    # 2.5.6.4 Variable
    Trs_MPC.append(sim_res['Trs'][-1])
    Trr_MPC.append(sim_res['Trr'][-1])
    ma_MPC.append(sim_res['ma'][-1])
    Qload_MPC.append(sim_res['Qload'][-1])
    Qsh_MPC.append(sim_res['Qsh'][-1])
    Qra_MPC.append(sim_res['Qra'][-1])
    Qve_MPC.append(sim_res['Qve'][-1])
    Qinf_MPC.append(sim_res['Qinf'][-1])
    Qen_MPC.append(sim_res['Qen'][-1])
    Qmaen_MPC.append(sim_res['Qmaen'][-1])
    Qma_MPC.append(sim_res['Qma'][-1])
    # 2.5.6.5 State variables
    Ten_MPC.append(sim_res['Ten'][-1])
    Tia_MPC.append(sim_res['Tia'][-1])
    Tma_MPC.append(sim_res['Tma'][-1])
    Tra1_MPC.append(sim_res['Tra[1]'][-1])
    Tra2_MPC.append(sim_res['Tra[2]'][-1])
    Tra3_MPC.append(sim_res['Tra[3]'][-1])
    Tra4_MPC.append(sim_res['Tra[4]'][-1])
    Tra5_MPC.append(sim_res['Tra[5]'][-1])
    # 2.5.6.6 Refrence for MPC
    Tiaref_MPC.append(res['Tiaref'][0])
    P_EDC_MPC.append(res['P_EDC'][0])
    # 2.5.6.7 Generate random noise
    noise_MPC = 1
    # 2.5.6.8 Extract staets at the end of the simulation
    Ten = sim_res['Ten'][-1]*noise_MPC
    Tia = sim_res['Tia'][-1]*noise_MPC
    Tma = sim_res['Tma'][-1]*noise_MPC
    Tra1 = sim_res['Tra[1]'][-1]*noise_MPC
    Tra2 = sim_res['Tra[2]'][-1]*noise_MPC
    Tra3 = sim_res['Tra[3]'][-1]*noise_MPC
    Tra4 = sim_res['Tra[4]'][-1]*noise_MPC
    Tra5 = sim_res['Tra[5]'][-1]*noise_MPC
    # 2.5.6.9 Set initial conditions for the next step
    op.set('Ten_init', float(Ten))
    op.set('Tia_init', float(Tia))
    op.set('Tma_init', float(Tma))
    op.set('Tra1_init', float(Tra1))
    op.set('Tra2_init', float(Tra2))
    op.set('Tra3_init', float(Tra3))
    op.set('Tra4_init', float(Tra4))
    op.set('Tra5_init', float(Tra5))
    sim_model.reset ()
    sim_model.set('Ten_init', float(Ten))
    sim_model.set('Tia_init', float(Tia))
    sim_model.set('Tma_init', float(Tma))
    sim_model.set('Tra1_init', float(Tra1))
    sim_model.set('Tra2_init', float(Tra2))
    sim_model.set('Tra3_init', float(Tra3))
    sim_model.set('Tra4_init', float(Tra4))
    sim_model.set('Tra5_init', float(Tra5))

##########################################################################################  
### 3 Plot
# Air Temperature
plt.close('Air Temperature')
plt.figure('Air Temperature')
plt.plot(t_MPC,Tia_MPC)
plt.plot(t_MPC,Toa_MPC)
plt.legend(('Tia','Toa'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')
plt.title('Air Temperature')
plt.show()

# Water temperature
plt.close('Water temperature')
plt.figure('Water temperature')
plt.plot(t_MPC,Ts_MPC)
plt.plot(t_MPC,Trr_MPC)
plt.legend(('Ts','Trr'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')
plt.title('Water Temperature')
plt.show()

# State temperature
plt.close('State temperature')
plt.figure('State temperature')
plt.plot(t_MPC,Tia_MPC)
plt.plot(t_MPC,Ten_MPC)
plt.plot(t_MPC,Tma_MPC)
plt.legend(('Tia','Ten','Tma'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')
plt.title('State Temperature')
plt.show()

# Mass flow rate
plt.close('Mass flow')
plt.figure('Mass flow')
plt.plot(t_MPC,mr_MPC)
plt.legend(('mr'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Mass flow (kg/s)')
plt.title('Mass flow')
plt.show()

# Heat flow rate
plt.close('Heat flow')
plt.figure('Heat flow')
plt.plot(t_MPC,Qin_MPC)
plt.plot(t_MPC,[i * -1 for i in Qma_MPC])
plt.plot(t_MPC,Qsh_MPC)
plt.plot(t_MPC,Qra_MPC)
plt.plot(t_MPC,[i * -1 for i in Qmaen_MPC])
plt.legend(('Qin','Qma','Qsh','Qra','Qmaen'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Q (W)')
plt.title('Heat flow rate')
plt.show()