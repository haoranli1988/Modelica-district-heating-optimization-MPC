##########################################################################################  
### 1 Import Python and JModelica.org
# 1.1 Import from Python
import csv
import numpy as np
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
# 2.1 Paths
input_paths = "C:/GitHub/MPC DH system/Inputs.csv"
file_paths = ("C:/GitHub/MPC DH system/MPC.mop", "C:/GitHub/MPC DH system/Heatsupply.mo")
# 2.2 Setting and inputs
N=10                                          # Sections for water tank
Days=1                                        # MPC control period (1-7 dags)
Sim_period=3600                               # time step for optimization window shifting (s)
Opt_period=3600                               # Resolution for control action (s)
horizon=12                                    # Optimization horizon (hour)
predict_time=Opt_period*horizon               # Optimization horizon (s)
final_time=Days*24*3600                       # MPC control period (s)
n_e=horizon                                   # Number of finite elements
n_cp=1                                        # Number of collocation points in each element
number_samp_tot=int(final_time/Sim_period)    # Total number of samples to do 
# 2.3 Read input file
Index_start=3432
csv_data = pd.read_csv(input_paths)
Timesamples = csv_data.Time[Index_start:-1]-csv_data.Time[Index_start]
Toa=csv_data.Toa[Index_start:-1]
qload=csv_data.qload[Index_start:-1]
P_ele=csv_data.P_ele[Index_start:-1]
P_edc=csv_data.P_edc[Index_start:-1]
Timesamples.reset_index(drop=True, inplace=True)
Toa.reset_index(drop=True, inplace=True)
qload.reset_index(drop=True, inplace=True)
P_ele.reset_index(drop=True, inplace=True)
P_edc.reset_index(drop=True, inplace=True)

##########################################################################################  
### 3  Define MPC problem
# 3.1 Initial guess trajectories
# 3.1.1 Get simulation model
init_fmu = compile_fmu("Package.Simu",file_paths,compiler_options={'generate_html_diagnostics': True})
init_model = load_fmu(init_fmu)
# 3.1.2 Set input variables
init_model.set('qload',3E6)
init_model.set('Toa',4)
# 3.1.3 Set manipulated variables
init_model.set('T_MS_s1',60)
init_model.set('m_MS1',10)
init_model.set('T_MS_s2',80)
init_model.set('m_MS2',40)
init_model.set('m_WT',0)
init_model.set('P_DC',3E5)
# 3.1.4 Set initial states
Twt1_ini=60;
Twt2_ini=60;
Twt3_ini=60;
Twt4_ini=60;
Twt5_ini=60; 
Twt6_ini=60;
Twt7_ini=60;
Twt8_ini=60;
Twt9_ini=60;
Twt10_ini=60;
init_model.set('WT.Twt1_init',Twt1_ini)
init_model.set('WT.Twt2_init',Twt2_ini)
init_model.set('WT.Twt3_init',Twt3_ini)
init_model.set('WT.Twt4_init',Twt4_ini)
init_model.set('WT.Twt5_init',Twt5_ini)
init_model.set('WT.Twt6_init',Twt6_ini)
init_model.set('WT.Twt7_init',Twt7_ini)
init_model.set('WT.Twt8_init',Twt8_ini)
init_model.set('WT.Twt9_init',Twt9_ini)
init_model.set('WT.Twt10_init',Twt10_ini)
# 3.1.5 Compute initial guess trajectories
init_res = init_model.simulate(start_time=0., final_time=final_time)

# 3.2 Define simulation problem
sim_fmu = compile_fmu("Package.Simu",file_paths,compiler_options={'generate_html_diagnostics': True})
sim_model = load_fmu(sim_fmu)

# 3.3 Define the optimal control problem
# 3.3.1 Compile and load optimization problem
op = transfer_optimization_problem("Package.Opt", file_paths, compiler_options={'generate_html_diagnostics': True, 'equation_sorting': True})
# 3.3.2 Create blocking factors with quadratic penalty and bound on 'm_MS1', 'm_MS2', 'm_WT', 'T_MS_s1', 'T_MS_s2', and 'P_DC'
bf_list = [1]*horizon
factors = {'m_MS1': bf_list, 'm_MS2': bf_list, 'm_WT': bf_list, 'T_MS_s1': bf_list, 'T_MS_s2': bf_list, 'P_DC': bf_list}
du_quad_pen = {'m_MS1': 0.5, 'm_MS2': 0.001, 'm_WT': 1, 'T_MS_s1': 0.001, 'T_MS_s2': 0.001, 'P_DC': 0.000000001}
du_bounds = {'m_MS1': 5, 'm_MS2': 10, 'm_WT': 5, 'T_MS_s1': 20, 'T_MS_s2': 20, 'P_DC': 50000}
bf = BlockingFactors(factors,du_bounds,du_quad_pen)
# 3.3.3 Set collocation options
opt_opts = op.optimize_options()
opt_opts['n_e'] = n_e
opt_opts['n_cp'] = n_cp
opt_opts['init_traj'] = init_res
opt_opts['nominal_traj'] = init_res
opt_opts['blocking_factors'] = bf
opt_opts['IPOPT_options']['tol'] = 1e-4
opt_opts['IPOPT_options']['acceptable_tol'] = 1e-3
opt_opts['IPOPT_options']['acceptable_iter'] = 10
opt_opts['IPOPT_options']['acceptable_obj_change_tol'] = 1
opt_opts["IPOPT_options"]["max_iter"] = 300

# 3.4 Lists for logging 
# 3.4.1 Global variable
t_MPC = []
# 3.4.2 Manipulated variable
m_MS1_MPC = []
m_MS2_MPC = []
m_WT_MPC = []
T_MS_s1_MPC = []
T_MS_s2_MPC = []
P_DC_MPC = []
# 3.4.3 Input variables
Toa_MPC = []
qload_MPC = []
# 3.4.4 Variable
qMS_MPC = []
qMS1_MPC = []
qMS2_MPC = []
qwt_MPC = []
qlosswt_MPC = []
qlosspipe_MPC = []
Qstore_MPC = []
T_MS_r1_MPC = []
T_MS_r2_MPC = []
T_Bu_s_MPC = []
T_Bu_r_MPC = []
Dif_MPC = []
T_WT_s_MPC = []
T_WT_r_MPC = []
T_evapout_MPC = []
Tcondin_MPC = []
Tcondout_MPC = []
qcond_MPC = []
# 3.4.5 State variables
Twt1_MPC = []
Twt2_MPC = []
Twt3_MPC = []
Twt4_MPC = []
Twt5_MPC = []
Twt6_MPC = []
Twt7_MPC = []
Twt8_MPC = []
Twt9_MPC = []
Twt10_MPC = []
# 3.4.6 Energy price
P_edc_MPC = []
P_ele_MPC = []

# 3.5 MPC loops
for k in range(number_samp_tot):
    if k == 0:
       op.set('WT.Twt1_init',Twt1_ini)
       op.set('WT.Twt2_init',Twt2_ini)
       op.set('WT.Twt3_init',Twt3_ini)
       op.set('WT.Twt4_init',Twt4_ini)
       op.set('WT.Twt5_init',Twt5_ini)
       op.set('WT.Twt6_init',Twt6_ini)
       op.set('WT.Twt7_init',Twt7_ini)
       op.set('WT.Twt8_init',Twt8_ini)
       op.set('WT.Twt9_init',Twt9_ini)
       op.set('WT.Twt10_init',Twt10_ini)
       sim_model.set('WT.Twt1_init',Twt1_ini)
       sim_model.set('WT.Twt2_init',Twt2_ini) 
       sim_model.set('WT.Twt3_init',Twt3_ini)
       sim_model.set('WT.Twt4_init',Twt4_ini) 
       sim_model.set('WT.Twt5_init',Twt5_ini)
       sim_model.set('WT.Twt6_init',Twt6_ini)
       sim_model.set('WT.Twt7_init',Twt7_ini) 
       sim_model.set('WT.Twt8_init',Twt8_ini)
       sim_model.set('WT.Twt9_init',Twt9_ini) 
       sim_model.set('WT.Twt10_init',Twt10_ini)
    # 3.5.1 Prepare input for optimization
    t0 = k*Sim_period
    tf = k*Sim_period+predict_time  
    data_input_opt_select = OrderedDict()
    data_input_opt_select['Toa'] = np.vstack([Timesamples[0:horizon],Toa[int(t0/Opt_period):int(tf/Opt_period)]])
    data_input_opt_select['qload'] = np.vstack([Timesamples[0:horizon],qload[int(t0/Opt_period):int(tf/Opt_period)]])
    data_input_opt_select['P_edc'] = np.vstack([Timesamples[0:horizon],P_edc[int(t0/Opt_period):int(tf/Opt_period)]])
    data_input_opt_select['P_ele'] = np.vstack([Timesamples[0:horizon],P_ele[int(t0/Opt_period):int(tf/Opt_period)]])
    opt_opts['external_data'] = ExternalData(eliminated=data_input_opt_select)
    # 3.5.2 Solve optimization problem
    res = op.optimize(options=opt_opts)
    # 3.5.3 Get optimal Manipulated variable
    T_MS_s1 = res['T_MS_s1'][0]
    m_MS1 = res['m_MS1'][0]
    T_MS_s2 = res['T_MS_s2'][0]
    m_MS2 = res['m_MS2'][0]
    m_WT = res['m_WT'][0]
    P_DC = res['P_DC'][0]
    # 3.5.4 Reset the model and set the new initial states before simulating
    sim_model.set('T_MS_s1',float(T_MS_s1))
    sim_model.set('m_MS1',float(m_MS1))
    sim_model.set('T_MS_s2',float(T_MS_s2))
    sim_model.set('m_MS2',float(m_MS2))
    sim_model.set('m_WT',float(m_WT))
    sim_model.set('P_DC',float(P_DC))
    data_sim = np.vstack(( Timesamples[int(t0/Opt_period)], Toa[int(t0/Opt_period)], qload[int(t0/Opt_period)]))
    data_input_sim = (['Toa', 'qload'],np.transpose(data_sim))
    # 3.5.5 Simulate using optimal inputs
    sim_res = sim_model.simulate(start_time= k*Sim_period, final_time=(k+1)*Sim_period, input=data_input_sim)
    # 3.5.6 Logging
    # 3.5.6.1 Global variable
    t_MPC.append(sim_res['time'][-1])
    # 3.5.6.2 Manipulated variable
    m_MS1_MPC.append(sim_res['m_MS1'][-1])
    m_MS2_MPC.append(sim_res['m_MS2'][-1])
    m_WT_MPC.append(sim_res['m_WT'][-1])
    T_MS_s1_MPC.append(sim_res['T_MS_s1'][-1])
    T_MS_s2_MPC.append(sim_res['T_MS_s2'][-1])
    P_DC_MPC.append(sim_res['P_DC'][-1])
    # 3.5.6.3 Input variables
    Toa_MPC.append(sim_res['Toa'][-1])
    qload_MPC.append(sim_res['qload'][-1])
    # 3.5.6.4 Variable
    qMS_MPC.append(sim_res['qMS'][-1])
    qMS1_MPC.append(sim_res['qMS1'][-1])
    qMS2_MPC.append(sim_res['qMS2'][-1])
    qwt_MPC.append(sim_res['qwt'][-1])
    qlosswt_MPC.append(sim_res['qlosswt'][-1])
    qlosspipe_MPC.append(sim_res['qlosspipe'][-1])
    Qstore_MPC.append(sim_res['Qstore'][-1])
    T_MS_r1_MPC.append(sim_res['T_MS_r1'][-1])
    T_MS_r2_MPC.append(sim_res['T_MS_r2'][-1])
    T_Bu_s_MPC.append(sim_res['T_Bu_s'][-1])
    T_Bu_r_MPC.append(sim_res['T_Bu_r'][-1])
    Dif_MPC.append(sim_res['Dif'][-1])
    T_WT_s_MPC.append(sim_res['T_WT_s'][-1])
    T_WT_r_MPC.append(sim_res['T_WT_r'][-1])
    T_evapout_MPC.append(sim_res['T_evapout'][-1])
    Tcondin_MPC.append(sim_res['DC.Tcondin'][-1])
    Tcondout_MPC.append(sim_res['DC.Tcondout'][-1])
    qcond_MPC.append(sim_res['DC.qcond'][-1])
    # 3.5.6.5 State variables
    Twt1_MPC.append(sim_res['WT.Twt[1]'][-1])
    Twt2_MPC.append(sim_res['WT.Twt[2]'][-1])
    Twt3_MPC.append(sim_res['WT.Twt[3]'][-1])
    Twt4_MPC.append(sim_res['WT.Twt[4]'][-1])
    Twt5_MPC.append(sim_res['WT.Twt[5]'][-1])
    Twt6_MPC.append(sim_res['WT.Twt[6]'][-1])
    Twt7_MPC.append(sim_res['WT.Twt[7]'][-1])
    Twt8_MPC.append(sim_res['WT.Twt[8]'][-1])
    Twt9_MPC.append(sim_res['WT.Twt[9]'][-1])
    Twt10_MPC.append(sim_res['WT.Twt[10]'][-1])
    # 3.5.6.6 Energy price
    P_edc_MPC.append(res['P_edc'][0])
    P_ele_MPC.append(res['P_ele'][0])
    # 3.5.6.7 Generate random noise
    noise_MPC = 1
    # 3.5.6.8 Extract staets at the end of the simulation
    Twt1 = sim_res['WT.Twt[1]'][-1]*noise_MPC
    Twt2 = sim_res['WT.Twt[2]'][-1]*noise_MPC
    Twt3 = sim_res['WT.Twt[3]'][-1]*noise_MPC
    Twt4 = sim_res['WT.Twt[4]'][-1]*noise_MPC
    Twt5 = sim_res['WT.Twt[5]'][-1]*noise_MPC
    Twt6 = sim_res['WT.Twt[6]'][-1]*noise_MPC
    Twt7 = sim_res['WT.Twt[7]'][-1]*noise_MPC
    Twt8 = sim_res['WT.Twt[8]'][-1]*noise_MPC
    Twt9 = sim_res['WT.Twt[9]'][-1]*noise_MPC
    Twt10 = sim_res['WT.Twt[10]'][-1]*noise_MPC
    # 3.5.6.9 Set initial conditions for the next step
    op.set('WT.Twt1_init',float(Twt1))
    op.set('WT.Twt2_init',float(Twt2))
    op.set('WT.Twt3_init',float(Twt3))
    op.set('WT.Twt4_init',float(Twt4))
    op.set('WT.Twt5_init',float(Twt5))
    op.set('WT.Twt6_init',float(Twt6))
    op.set('WT.Twt7_init',float(Twt7))
    op.set('WT.Twt8_init',float(Twt8))
    op.set('WT.Twt9_init',float(Twt9))
    op.set('WT.Twt10_init',float(Twt10))
    sim_model.reset ()
    sim_model.set('WT.Twt1_init',float(Twt1))
    sim_model.set('WT.Twt2_init',float(Twt2))
    sim_model.set('WT.Twt3_init',float(Twt3))
    sim_model.set('WT.Twt4_init',float(Twt4))
    sim_model.set('WT.Twt5_init',float(Twt5))
    sim_model.set('WT.Twt6_init',float(Twt6))
    sim_model.set('WT.Twt7_init',float(Twt7))
    sim_model.set('WT.Twt8_init',float(Twt8))
    sim_model.set('WT.Twt9_init',float(Twt9))
    sim_model.set('WT.Twt10_init',float(Twt10))

##########################################################################################  
### 4 Plot
# Outdoor temperature
plt.close('Outdoor temperature')
plt.figure('Outdoor temperature')
plt.plot(t_MPC,Toa_MPC)
plt.legend(('Toa'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Outdoor temperature (Centigrade)')
plt.title('Outdoor temperature')
plt.show()

# MS1 Water temperature
plt.close('MS1 water temperature')
plt.figure('MS1 water temperaturee')
plt.plot(t_MPC,T_MS_s1_MPC)
plt.plot(t_MPC,T_MS_r1_MPC)
plt.legend(('T_MS_s1','T_MS_r1'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')
plt.title('MS1 water temperature')
plt.show()

# MS2 Water temperature
plt.close('MS2 water temperature')
plt.figure('MS2 water temperaturee')
plt.plot(t_MPC,T_MS_s2_MPC)
plt.plot(t_MPC,T_MS_r2_MPC)
plt.legend(('T_MS_s2','T_MS_r2'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')
plt.title('MS2 water temperature')
plt.show()

# Building water temperature
plt.close('Building water temperature')
plt.figure('Building water temperature')
plt.plot(t_MPC,T_Bu_s_MPC)
plt.plot(t_MPC,T_Bu_r_MPC)
plt.plot(t_MPC,Dif_MPC)
plt.legend(('T_Bu_s','T_Bu_r','Dif'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')
plt.title('Building water temperature')
plt.show()

# Water tank temperature
plt.close('Water tank temperature')
plt.figure('Water tank temperature')
plt.plot(t_MPC,T_WT_s_MPC)
plt.plot(t_MPC,T_WT_r_MPC)
plt.legend(('T_WT_s','T_WT_r'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')
plt.title('Water tank temperature')
plt.show()

# Condensor temperature
plt.close('Condensor temperature')
plt.figure('Condensor temperature')
plt.plot(t_MPC,Tcondin_MPC)
plt.plot(t_MPC,Tcondout_MPC)
plt.legend(('Tcondin','Tcondout'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')
plt.title('Condensor temperature')
plt.show()

# Outlet temperature of evaprator
plt.close('Outlet temperature of evaprator')
plt.figure('Outlet temperature of evaprator')
plt.plot(t_MPC,T_evapout_MPC)
plt.legend(('T_evapout'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')
plt.title('Outlet temperature of evaprator')
plt.show()

# State temperature
plt.close('State temperature')
plt.figure('State temperature')
plt.plot(t_MPC,Twt1_MPC)
plt.plot(t_MPC,Twt2_MPC)
plt.plot(t_MPC,Twt3_MPC)
plt.plot(t_MPC,Twt4_MPC)
plt.plot(t_MPC,Twt5_MPC)
plt.plot(t_MPC,Twt6_MPC)
plt.plot(t_MPC,Twt7_MPC)
plt.plot(t_MPC,Twt8_MPC)
plt.plot(t_MPC,Twt9_MPC)
plt.plot(t_MPC,Twt10_MPC)
plt.legend(('Twt1','Twt2','Twt3','Twt4','Twt5','Twt6','Twt7','Twt8','Twt9','Twt10'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')
plt.title('State Temperature')
plt.show()

# Mass flow rate
plt.close('Mass flow')
plt.figure('Mass flow')
plt.plot(t_MPC,m_MS1_MPC)
plt.plot(t_MPC,m_MS2_MPC)
plt.plot(t_MPC,m_WT_MPC)
plt.legend(('m_MS1','m_MS2','m_WT'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Mass flow (kg/s)')
plt.title('Mass flow')
plt.show()

# Heat flow rate
plt.close('Heat flow')
plt.figure('Heat flow')
plt.plot(t_MPC,qload_MPC)
plt.plot(t_MPC,qMS_MPC)
plt.plot(t_MPC,qcond_MPC)
plt.plot(t_MPC,qwt_MPC)
plt.legend(('qload','qMS','qcond','qwt'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Q (W)')
plt.title('Heat flow rate')
plt.show()

# Electricity of HP in data centre
plt.close('Electricity of HP')
plt.figure('Electricity of HP')
plt.plot(t_MPC,P_DC_MPC)
plt.legend(('P_DC'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Power(W)')
plt.title('Electricity of HP')
plt.show()

# Stored heat rate
plt.close('Stored heat rate')
plt.figure('Stored heat rate')
plt.plot(t_MPC,Qstore_MPC)
plt.legend(('Qstore'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Q (W)')
plt.title('Stored heat rate')
plt.show()