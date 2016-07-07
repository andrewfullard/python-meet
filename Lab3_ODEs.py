# -*- coding: utf-8 -*-
"""
Week 3: Ordinary Differential Equations
Mark Hannam
"""

from pylab import *

# Function to define Euler step. 
# Assume no explicit time dependence in the solution
def EulerStep(yd, dt, rhs):
    dy = rhs(yd)
    return (yd + dy*dt)

# Define ODE solver
# This assumes a single ODE (not a system), and takes as input
#    the initial time t0
#    the final time Tmax, 
#    the initial value y0, 
#    the number of integration steps Npts
#    the integration method to use (Euler is currently the only choice)
#    the right-hand-side function "rhs".
# The function therefore solves 
#    y' = rhs(y), with y(t0) = y0, up to Tmax with Npts steps
# The funtion returns the times and the values as separate arrays
def solveODE(t0, y0, Tmax, Npts, rhs, method = 'Euler'):
    dt = Tmax/Npts
    # Make sure to ask for Npts+1 arrays, so that we go all the way to Tmax.
    tdat = arange(t0,Tmax+dt,dt)
    yvec = zeros(Npts+1)
    yvec[0] = y0
    for i in range(0,Npts):
        if (method == 'Euler'):
            yvec[i+1]= EulerStep(yvec[i],dt,rhs)
    return tdat, yvec
  
# RHS for the ODE
#    y' = y  
def f1(yd):
    return yd
    
# RHS for the ODE
#    y' = -y
def f2(yd):
    return -yd


# 2. Solve y' = y with y(0) = 1, up to t=1.
# The analytic solution is y(t) = exp(t).
#
# Show both the numercal and analytic solutions for some coarse resolution,
# as a demonstration that we get qualitatively the right solution.
# Use two time resolutions to show that the numerical solution gets closer
# to the analytic solution as the time step is reduced. 

times1, sol1 = solveODE(0., 1., 1., 10, f1, 'Euler')
times2, sol2 = solveODE(0., 1., 1., 40, f1, 'Euler')
figure(1)
clf()
plot(times1, sol1, 'rx',label='Numerical solution, dt = 1/10')
plot(times2, sol2, 'bx',label='Numerical solution, dt = 1/40')
# Construct the analytic solution
analytic = exp(times2)
plot(times2, analytic,'b-', label='Analytic solution, exp(t)')
title('Solution of y\' = y')
xlabel('t')
ylabel('y')
legend(loc='upper left')

#
# Demonstrate first-order convergence with respect to the analytic solution 
# by comparing solutions with 100 and 200 points. 

N1 = 100
times1, sol1 = solveODE(0., 1., 1., N1, f1, 'Euler')
times2, sol2 = solveODE(0., 1., 1., 2*N1, f1, 'Euler')

analytic1 = exp(times1)
analytic2 = exp(times2)

# Differences between the numerical and analytic solutions.
# The factor of two reflects our expectation of first-order convergence. 
diff1 = sol1 - analytic1
diff2 = 2.*(sol2 - analytic2)

figure(2)
clf()
plot(times1,diff1, 'r-', label='dt=1/100 error')
plot(times2,diff2, 'b-', label='dt=1/200 error, \n scaled for 1st-order convergence')
title('Convergence test against analytic solution')
xlabel('t')
ylabel('Error in y(t)')
legend(loc='lower left')
# The figure clearly shows that the difference between the numerical and
# analytic solutions falls linarly with the timestep, indicating first-order
# convrgence.
 
 
# 3. Demonstrate self convergence

# Function to construct differences for convergence test. 
def ConvergenceTest(Tmax, N, method, rhs):
    # Produce solution, with time step successively reduced by factor of 2
    times1, sol1 = solveODE(0., 1., 1., N, rhs, method)
    times2, sol2 = solveODE(0., 1., 1., 2*N, rhs, method)
    times3, sol3 = solveODE(0., 1., 1., 4*N, rhs, method)
    # Make sure to choose points so that the times line up between solutions
    # produced with different time steps
    diff1 = sol1 - sol2[::2]
    diff2 = sol2[::2] - sol3[::4]
    return times1, diff1, diff2
    
# Use base time steps 1/40, 1/80 and 1/160 for convergence tests. 
t1, d1, d2 = ConvergenceTest(1., 40, 'Euler', f1)
t2, d3, d4 = ConvergenceTest(1., 80, 'Euler', f1)
t3, d5, d6 = ConvergenceTest(1., 160, 'Euler', f1)

# Self-convergence test with dt = 1/40, 1/80 and 1/160.
figure(3)
clf()
plot(t1,d1,'r-', label = '(dt=1/40) - (dt=1/80)')
plot(t1,2.*d2,'b-', label = '(dt=1/80) - (dt=1/160) \n scaled for 1st-order convergence')
title('Self-convergence test for y\' = y')
xlabel('t')
ylabel('Differences in y(t)')
legend(loc='lower left')
# The curves lie almost on top of each other, verifying almost perfect first-order
# convergence. 

# Look at the ratios of the differences; if they are close to 2, then we
# have first-order convergence.
figure(4)
clf()
# Don't include the first point, because the differences are zero there,
# and we don't want to divide 0/0 !
plot(t1[1:],d1[1:]/d2[1:], 'r-', label='dt = 1/40, 1/80, 1/160')
plot(t2[1:],d3[1:]/d4[1:], 'b-', label='dt = 1/80, 1/160, 1/320')
plot(t3[1:],d5[1:]/d6[1:], 'g-', label='dt = 1/160, 1/320, 1/640')
title('Self-convergence test (as ratios) for y\' = y')
xlabel('t')
ylabel('Ratios of differences in y(t)')
legend(loc='lower left')
# The convergence rate gets closer to exactly first-order as the time
# resolution is improved. If the code is correct, and the resolution good
# enough, the method should asymptote to perfect first-order convergence 
# as dt -> 0, and this plot suggests that this does indeed happen.  
    
    
# 4. Solve y' = -y
# The right-hand-side (RHS) for this problem is given by the function f2
times1, sol1 = solveODE(0., 1., 16., 2, f2, 'Euler')
times2, sol2 = solveODE(0., 1., 16., 4, f2, 'Euler')
times3, sol3 = solveODE(0., 1., 16., 8, f2, 'Euler')
times4, sol4 = solveODE(0., 1., 16., 16, f2, 'Euler')
times5, sol5 = solveODE(0., 1., 16., 32, f2, 'Euler')

# Look at stability by choosing a range of large time steps
figure(5)
clf()
plot(times1,sol1, 'r-', label='dt = 8')
plot(times2,sol2, 'g-', label='dt = 4')
plot(times3,sol3, 'b-', label='dt = 2')
plot(times4,sol4, 'c-', label='dt = 1')
plot(times5,sol5, 'm-', label='dt = 0.5')
title('Stability study of y\' = -y')
xlabel('t')
ylabel('y(t)')
legend(loc='upper left')
# The analytic solution is exp(-t), i.e., exponential decay
# The numerical solution does not show the correct qualitative behaviuor 
# (the solution does not monotonically decay) for the choices dt = (8, 4, 2).
# This is consistent with our expectation that the system is unstable if 
# |1 - dt| > 1. 

# The largest time step for which |1 - dt| < 1 is dt = 1. 
# Now check that the dt =1 and dt = 0.5 solution show the correct qualitative 
# behaviour (i.e., they *are* stable).
# Do this on a separate plot, because it is not clear in Fig. 5.
figure(6)
clf()
plot(times4,sol4, 'c-', label='dt = 1')
plot(times5,sol5, 'm-', label='dt = 0.5')
title('Stability study of y\' = -y \n (stable timestep choices)')
xlabel('t')
ylabel('y(t)')
legend(loc='upper left')
# We see that dt = 1 is stable, but extemely inaccurate: The initial condition
# is y(0) = 1. The next step will be y(1) = 1 - dt*1 = 0. The next step will 
# be y(2) = 0 - dt*0 = 0, and the solution will remain at zero for all 
# successive steps. 
#
# The choice dt = 0.5 looks more like a standard exponential decay. 


