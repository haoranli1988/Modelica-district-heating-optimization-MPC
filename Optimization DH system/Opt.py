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
# 2.1 Setting for Paths
input_paths = "C:/GitHub/Optimization/Inputs.csv"
file_paths = ("C:/GitHub/Optimization/OpFile.mop",
              "C:/GitHub/Optimization/Plants.mo")
# 2.2 Setting for JModelica Optimization  
Days=30                              # Horizon (1-365 days), remember to set OpFile "Days" also
Act=24                               # Resolution (control actions per day)
n_cp=1                               # Resolution (grid for optimization per action)
final_time=Days*24*3600              # Setting for JModelica
n_e=Days*Act                         # Setting for JModelica
# 2.3 Setting for Modelica Model
# 2.3.1 Sections for water tank
N=10                                 
# 2.3.2 Get Modelica model
init_fmu = compile_fmu("Package.Simu",file_paths, compiler_options={'generate_html_diagnostics': True})
init_model = load_fmu(init_fmu)
# 2.3.3 Set Input variables
init_model.set('Toa',5)
init_model.set('Tiaref',20)
init_model.set('qin',2E5)
init_model.set('qinf',9E5)
init_model.set('qve',2E6)
init_model.set('qdhw',5E5)
init_model.set('qdc',1E6)
# 2.3.4 Set Manipulated variables
init_model.set('T_MS_s1',80)
init_model.set('m_MS1',20)
init_model.set('T_MS_s2',80)
init_model.set('m_MS2',80)
init_model.set('m_Bu',100)
# 2.3.5 Set parameters
Vwt=1700
init_model.set('Vwt',Vwt)
# 2.3.6 Get inputs
csv_data = pd.read_csv(input_paths)
Timesamples = csv_data.Time
Toa=csv_data.Toa
Tiaref=csv_data.Tiaref
qin=csv_data.qin
qinf=csv_data.qinf
qve=csv_data.qve
qdhw=csv_data.qdhw
# 2.3.7 Compute initial guess trajectories
init_traj = init_model.simulate(start_time=0., final_time=final_time)  

##########################################################################################  
### 3 Mian Code - Optimization
# 3.1 Define the optimization problem
op = transfer_optimization_problem("Package.Opt", file_paths, compiler_options={'generate_html_diagnostics': True, 'equation_sorting': True})
op.set('Vwt', Vwt)
# 3.2 Set optimization and simulation options
opt_opts = op.optimize_options()
opt_opts['result_mode'] = 'mesh_points'
opt_opts['n_e'] = n_e
opt_opts['n_cp'] = n_cp
# 3.3 Prepare input for optimization
data_input_opt_select = OrderedDict()
data_input_opt_select['Toa'] = np.vstack([Timesamples,Toa])
data_input_opt_select['Tiaref'] = np.vstack([Timesamples,Tiaref])
data_input_opt_select['qin'] = np.vstack([Timesamples,qin])
data_input_opt_select['qinf'] = np.vstack([Timesamples,qinf])
data_input_opt_select['qve'] = np.vstack([Timesamples,qve])
data_input_opt_select['qdhw'] = np.vstack([Timesamples,qdhw])
opt_opts['external_data'] = ExternalData(eliminated=data_input_opt_select)
opt_opts['IPOPT_options']['tol'] = 1e-3
opt_opts["IPOPT_options"]["max_iter"] = 1000
op = BLTOptimizationProblem(op)  
# 3.4 Warm start
opt_opts['init_traj'] = init_traj
opt_opts['nominal_traj'] = init_traj   
# 3.5 Solve the optimization problem
opt_res = op.optimize(options=opt_opts)

##########################################################################################  
### 4 Log data and Plot
# 4.1 Log data
opt_time=opt_res['time'].flatten()
opt_qload=opt_res['qload'].flatten()
opt_qdc=opt_res['qdc'].flatten()
opt_qMS1=opt_res['qMS1'].flatten()
opt_qMS2=opt_res['qMS2'].flatten()
opt_T_MS_s1=opt_res['T_MS_s1'].flatten()
opt_T_MS_r1=opt_res['T_MS_r1'].flatten()
opt_T_MS_s2=opt_res['T_MS_s2'].flatten()
opt_T_MS_r2=opt_res['T_MS_r2'].flatten()
opt_T_Bu_s=opt_res['T_Bu_s'].flatten()
opt_T_Bu_r=opt_res['T_Bu_r'].flatten()
opt_T_WT_s=opt_res['T_WT_s'].flatten()
opt_T_WT_r=opt_res['T_WT_r'].flatten()
opt_m_MS1=opt_res['m_MS1'].flatten()
opt_m_MS2=opt_res['m_MS2'].flatten()
opt_m_Bu=opt_res['m_Bu'].flatten()
# 4.2 Plot
# 4.2.1 Heat demand and supply
plt.figure('Heat demand and supply')
plt.plot(opt_time,opt_qload)
plt.plot(opt_time,opt_qdc)
plt.plot(opt_time,opt_qMS1)
plt.plot(opt_time,opt_qMS2)
plt.legend(('heat demand- building','heat supply- DC','heat supply- MS1','heat supply- MS2'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Heat (W)')
plt.title('Heat demand vs Supply')
plt.show()
# 4.2.2 Operating temperature
plt.figure('Operating temperature')
plt.plot(opt_time,opt_T_MS_s1)
plt.plot(opt_time,opt_T_MS_r1)
plt.plot(opt_time,opt_T_MS_s2)
plt.plot(opt_time,opt_T_MS_r2)
plt.plot(opt_time,opt_T_WT_s)
plt.plot(opt_time,opt_T_WT_r)
plt.plot(opt_time,opt_T_Bu_s)
plt.plot(opt_time,opt_T_Bu_r)
plt.legend(('supply temperature- MS1','return temperature- MS1','supply temperature- MS2','return temperature- MS2','supply temperature- WT','return temperature- WT','supply temperature- building','return temperature- building'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (C)')
plt.title('Water temperature')
plt.show()
# 4.2.3 Operating flow rate
plt.figure('Operating flow rate')
plt.plot(opt_time,opt_m_MS1)
plt.plot(opt_time,opt_m_MS2)
plt.plot(opt_time,opt_m_Bu)
plt.legend(('flow- MS1','flow- MS2','flow- building'))
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Flow rate (kg/s)')
plt.title('Water flow rate')
plt.show()