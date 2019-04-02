"""2D Hartmann flow example

Uses vector potential form of MHD; incompressible Navier-Stokes, no further approximations.


"""
import time
import os
import dedalus.public as de
from dedalus.extras import flow_tools
from dedalus.tools import post
import numpy as np

import logging
logger = logging.getLogger(__name__)

Lx = 8.
Ly = 2.
nx = 64
ny = 64

Ha = 1.
Re = 1.
Rm = 1.
Pi = 1.
tau = 0.1

stop_time = 10
data_dir = "scratch"

x = de.Fourier('x', nx, interval=[0,Lx], dealias=3/2)
y = de.Chebyshev('y', ny, interval=[-Ly/2, Ly/2], dealias=3/2)

domain = de.Domain([x,y],grid_dtype='float')

hartmann = de.IVP(domain, variables=['vx', 'vy', 'Az', 'p', 'vx_y', 'vy_y', 'Az_y'])
hartmann.parameters['Ha'] = Ha
hartmann.parameters['Re'] = Re
hartmann.parameters['Rm'] = Rm
hartmann.parameters['Pi'] = Pi # pressure gradient driving flow
hartmann.parameters['Lx'] = Lx
hartmann.parameters['Ly'] = Ly
hartmann.parameters['tau'] = tau
hartmann.substitutions['Bx'] = "Az_y"
hartmann.substitutions['By'] = "-dx(Az) + 1."
hartmann.substitutions['Lap(A, Ay)'] = "dx(dx(A)) + dy(Ay)"
hartmann.substitutions['Jz'] = "-Lap(Az, Az_y)"
hartmann.substitutions['Avg_x(A)'] = "integ(A,'x')/Lx"

# Navier Stokes 
hartmann.add_equation("dt(vx) + dx(p) - Lap(vx, vx_y)/Re = -vx*dx(vx) - vy*dy(vx) - Ha**2/(Re*Rm) * Jz*By - Pi*(exp(-t/tau) - 1)")
hartmann.add_equation("dt(vy) + dy(p) - Lap(vy, vy_y)/Re = -vx*dx(vy) - vy*dy(vy) + Ha**2/(Re*Rm) * Jz*Bx")

# div(v) = 0
hartmann.add_equation("dx(vx) + vy_y = 0")

# Az
hartmann.add_equation("dt(Az) - Lap(Az, Az_y)/Rm = vx*By - vy*Bx")

# first order form
hartmann.add_equation("dy(Az) - Az_y = 0")
hartmann.add_equation("dy(vx) - vx_y = 0")
hartmann.add_equation("dy(vy) - vy_y = 0")

# boundary conditions
hartmann.add_bc("left(vx) = 0")
hartmann.add_bc("right(vx) = 0")
hartmann.add_bc("left(vy) = 0")
hartmann.add_bc("right(vy) = 0", condition="(nx != 0)")
hartmann.add_bc("right(p) = 0", condition="(nx == 0)")
hartmann.add_bc("left(Az) = 0")
hartmann.add_bc("right(Az) = 0")

# build solver
solver = hartmann.build_solver(de.timesteppers.MCNAB2)
logger.info("Solver built")

# Integration parameters
solver.stop_sim_time = stop_time
solver.stop_wall_time = 5*24*60.*60
solver.stop_iteration = np.inf
dt = 1e-3

# Initial conditions are zero by default in all fields

# Analysis
analysis_tasks = []
check = solver.evaluator.add_file_handler(os.path.join(data_dir,'checkpoints'), wall_dt=3540, max_writes=50)
check.add_system(solver.state)
analysis_tasks.append(check)

snap = solver.evaluator.add_file_handler(os.path.join(data_dir,'snapshots'), sim_dt=1e-1, max_writes=200)
snap.add_task("Bx", scales=1)
snap.add_task("By", scales=1)
snap.add_task("Az", scales=1)
snap.add_task("vx", scales=1)
snap.add_task("vy", scales=1)

analysis_tasks.append(snap)

integ = solver.evaluator.add_file_handler(os.path.join(data_dir,'integrals'), sim_dt=1e-3, max_writes=200)
integ.add_task("Avg_x(vx)", name='<vx>_x', scales=1)
integ.add_task("Avg_x(vy)", name='<vy>_x', scales=1)
integ.add_task("Avg_x(Bx)", name='<Bx>_x', scales=1)
integ.add_task("Avg_x(By)", name='<By>_x', scales=1)
analysis_tasks.append(integ)

timeseries = solver.evaluator.add_file_handler(os.path.join(data_dir,'timeseries'), sim_dt=1e-2)
timeseries.add_task("0.5*integ(vx**2 + vy**2)",name='Ekin')
timeseries.add_task("0.5*integ(Bx**2 + By**2)",name='Emag')
analysis_tasks.append(timeseries)

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=1)
flow.add_property("0.5*(vx**2 + vy**2)", name='Ekin')

try:
    logger.info('Starting loop')
    start_run_time = time.time()

    while solver.ok:
        solver.step(dt)
        if (solver.iteration-1) % 1 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('Max E_kin = %17.12e' %flow.max('Ekin'))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %f' %(end_run_time-start_run_time))


logger.info('beginning join operation')
for task in analysis_tasks:
    logger.info(task.base_path)
    post.merge_analysis(task.base_path)

