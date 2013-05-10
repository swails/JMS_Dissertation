"""
This program uses matplotlib to make the various plots included in my
dissertation.
"""
from __future__ import division

from argparse import ArgumentParser
from matplotlib import pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D
import numpy as np

def _evaluate(x, func):
   result = np.ndarray(x.shape)
   for i, val in enumerate(x):
      result[i] = func(val)
   return result

def timestep():
   """ Creates the time step plot """
   fig = plt.figure(1, figsize=(8,5))
   dumfig = plt.figure(2, figsize=(1,1))

   TIMESTEP = 1.0

   # Sinusoidal function
   def get_func(k):
      """ Returns the force function and the position function """
      return lambda x: -k * x, lambda x: np.cos(np.sqrt(k) * x)

   # Get arrow locations
   def arrow_locations(du, t0, x0, v0, dt, stop):
      head, tail = [], []
      t = t0
      x = x0
      v = v0
      while t < stop:
         head.append((t, x))
         v += du(x) * dt
         dx = v * dt
         tail.append((dt, dx))
         x += dx
         t += dt
      return head, tail

   # Generate the subplots and capture the axis instances
   ax1 = fig.add_subplot(122, xlim=(0,10), ylim=(-1.5,2))
   ax2 = fig.add_subplot(121, xlim=(0,10), ylim=(-1.5,2))

   # Add a grid
   ax1.grid(lw=1)
   ax2.grid(lw=1)

   # Generate the plot data
   du1, f1 = get_func(0.5)
   du2, f2 = get_func(1.5)

   x1 = np.arange(0, 11, 0.1)
   y1 = f1(x1)
   x2 = np.arange(0, 11, 0.1)
   y2 = f2(x2)
   
   # Apply the axis labels
   ax1.set_xlabel('Time', family='sans-serif', fontdict={'fontsize' : 20})
   ax2.set_ylabel('Position', family='sans-serif', fontdict={'fontsize' : 20})
   ax2.set_xlabel('Time', family='sans-serif', fontdict={'fontsize' : 20})
   ax1.set_yticklabels([])
#  ax2.set_ylabel('Position', family='sans-serif', fontdict={'fontsize' : 20})

   # Plot the data
   pl1, = ax1.plot(x1, y1, 'r-', lw=2, label='Exact Position')
   pl2, = ax2.plot(x2, y2, 'r-', lw=2, label='Exact Position')

   # Draw the arrows
   heads, tails = arrow_locations(du1, 0.0, 1.0, 0.0, TIMESTEP, 10.0)
   ax1.arrow(heads[0][0], heads[0][1], tails[0][0], tails[0][1], color='k',
             shape='full', head_width=0.10, length_includes_head=True)
   for i in range(1, len(heads)):
      head, tail = heads[i], tails[i]
      ax1.arrow(head[0], head[1], tail[0], tail[1], color='k', shape='full',
                head_width=0.10, length_includes_head=True)
   heads, tails = arrow_locations(du2, 0.0, 1.0, 0.0, TIMESTEP, 10.0)
   ax2.arrow(heads[0][0], heads[0][1], tails[0][0], tails[0][1], color='k',
             shape='full', head_width=0.10, length_includes_head=True)
   for i in range(1, len(heads)):
      head, tail = heads[i], tails[i]
      ax2.arrow(head[0], head[1], tail[0], tail[1], color='k', shape='full',
                head_width=0.10, length_includes_head=True)

   pr, = dumfig.add_subplot(111).plot(x1, y1, 'k-')
   ax1.legend((pr, pl1), ('Integrated Position', 'Exact Position'), loc='best')
   ax2.legend((pr, pl2), ('Integrated Position', 'Exact Position'), loc='best')

   # Annotate the graph

   fig.savefig('TimeStepDemo.ps')

def hydrogen():
   """ Makes a plot of parametrizing the H-H molecule bond """
   HAR_TO_KCAL = 627.509469 # conversion factor
   # From Kolos and Wolniewicz
   xdata = np.asarray([0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,
               0.85, 0.90, 1.00, 1.10, 1.2, 1.3, 1.35, 1.4, 1.401, 1.4011, 1.41,
               1.45, 1.50, 1.60, 1.70, 1.80, 1.90, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5,
               2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.6, 3.7, 7])
   ydata = np.asarray([-0.1202028, -0.3509282, -0.5266270, -0.6627583,
              -0.7696253, -0.8543531, -0.9220185, -0.9763287, -1.0200487,
              -1.0836362, -1.1245331, -1.1500512, -1.1649294, -1.1723414,
              -1.1739581, -1.1744472, -1.1744699, -1.1744701, -1.1744701,
              -1.1744556, -1.1740513, -1.1728492, -1.1685773, -1.1624521,
              -1.1550616, -1.1468425, -1.1381236, -1.1291528, -1.1201190,
              -1.1111659, -1.1024035, -1.0939149, -1.0857627, -1.0779927,
              -1.0706404, -1.0637259, -1.0572607, -1.0512547, -1.0457057,
              -1.0406020, -1.0359419, -1.0278471, -1.0243742, -1.0])
   xdata2 = np.arange(0, 10, 0.1)
   # Convert to kcal/mol
   ydata = ydata * HAR_TO_KCAL

   alpha = 1.011
   xeq = 1.401
   de = 0.1744701
   def convert(func):
      return lambda x: HAR_TO_KCAL * func(x)
   # Set up the fitting functions
   @convert
   def quadratic(x):
      return alpha ** 2 * 2 * de * (x - xeq) ** 2
   @convert
   def quartic(x):
      return de * (alpha ** 2 - alpha ** 3 * (x - xeq) + 7 / 12 * alpha ** 4 *
                  (x - xeq) ** 2) * (x - xeq) ** 2
   @convert
   def morse(x):
      return de * (1 - np.exp(alpha * (xeq - x))) ** 2


   fig = plt.figure(3, figsize=(8,5))

   ax = fig.add_subplot(111)
   # Set up options
   ax.grid(lw=1)
   ax.set_xlim((0,6))
   ax.set_ylim((-750,-500))
   
   # Plot the actual data
   realdat, = ax.plot(xdata, ydata, linestyle='-', marker='o', markersize=5,
                      linewidth=1, color='k')
   quaddat, = ax.plot(xdata2, quadratic(xdata2)+min(ydata), 
                      linestyle='--', color='b', lw=2)
   quardat, = ax.plot(xdata2, quartic(xdata2)+min(ydata), 
                      linestyle='-.', color='r', lw=2)
   mordat, = ax.plot(xdata2, morse(xdata2)+min(ydata), 
                     linestyle=':', color='g', lw=3)
   
   # Labels and legends
   ax.set_xlabel('Nuclear Separation ($\AA$)', family='sans-serif',
                 fontdict={'fontsize' : 16})
   ax.set_ylabel('Energy (kcal/mol)', family='sans-serif', 
                 fontdict={'fontsize' : 16})

   ax.legend((realdat, quaddat, quardat, mordat),
             ('Exact Potential', 'Quadratic Potential', 'Quartic Potential',
              'Morse Potential'), loc=1)

   fig.savefig('HydrogenMoleculeBond.ps')

def LennardJones():
   """ Makes the lennard-jones plot """
   fig = plt.figure(4, figsize=(8,5))

   ax = fig.add_subplot(111)

   # Set up axes
   ax.set_xlim((3,6))
   ax.set_ylim((-0.5,1.0))
   ax.set_xlabel('Nuclear Separation ($\AA$)', family='sans-serif',
                 fontdict={'fontsize' : 16})
   ax.set_ylabel(r'Energy ($kcal$  $mol^{-1}$)', family='sans-serif',
                 fontdict={'fontsize' : 16})
   RMIN = 3.816
   EPS = 0.1094
   SIG = RMIN / 2 ** (1/6)
   ACOEF = EPS * RMIN ** 12
   BCOEF = 2.0 * EPS * RMIN ** 6

   lj_attr = lambda x: - BCOEF / x ** 6
   lj_repu = lambda x: ACOEF / x ** 12
   lj = lambda x: lj_attr(x) + lj_repu(x)

   xdata = np.arange(0.2,10,0.05)
   ax.grid(lw=1)

   real, = ax.plot(xdata, lj(xdata), color='k', lw=2)
   attr, = ax.plot(xdata, lj_attr(xdata), linestyle='--', lw=2)
   repu, = ax.plot(xdata, lj_repu(xdata), linestyle='-.', lw=2)
   axis, = ax.plot([0,10], [0,0], color='k', lw=1)

   ax.legend((real, attr, repu),
             ('LJ Potential', 'Attractive Part', 'Repulsive Part'), loc=1)

   # Add attributes
   r1 = ax.arrow(RMIN, -EPS, 0, EPS - 0.5, linestyle='dashdot', shape='full',
            length_includes_head=True, color='k', head_width=0.03)
   r2 = ax.arrow(RMIN, -EPS, -RMIN + 3, 0, linestyle='dashdot', shape='full',
            length_includes_head=True, color='k', head_width=0.03)
   r3 = ax.arrow(SIG, 0, 0, -0.5, linestyle='dashdot', shape='full',
            length_includes_head=True, color='k', head_width=0.03)

   ax.annotate('$R_{min,i,j}$', (RMIN,-0.5), xytext=(4.58, -0.32), size=20,
               arrowprops={'arrowstyle':'fancy', 'fc':'0.6', 'ec':'none',
                           'patchB':r1})
   ax.annotate(r'$-\varepsilon_{i,j}$', (3.0, -EPS), xytext=(3.10,-0.3),
               size=20, arrowprops={'arrowstyle':'fancy', 'fc':'0.6',
                                    'ec':'none', 'patchB':r2})
   ax.annotate(r'$\sigma_{i,j}$', (SIG, -0.5), xytext=(4.2, -0.3), size=20,
               arrowprops={'arrowstyle':'fancy', 'fc':'0.6', 'ec':'none',
                           'patchB':r3})
   fig.savefig('LennardJones.ps')

def cmap(fname):
   """ Generates the cmap plot """
   KB = 0.00199 # kcal/mol / K
   TEMP = 300   # K
   ydata, xdata = np.loadtxt(fname).transpose()
   hist, hedge, vedge = np.histogram2d(xdata, ydata, bins=(100,100),
                              range=((-180,180), (-180,180)))
   hist /= np.max(hist)
   # Now turn it into a free energy
   hist = - KB * TEMP * np.log(hist)
   # Now time to generate the heat map
   fig = plt.figure(1, figsize=(8,5))

   ax = fig.add_subplot(111)
   
   ax.set_xlim((-180,180))
   ax.set_ylim((-180,180))
   ax.set_xlabel(r'$\phi_1$', family='sans-serif', fontdict={'fontsize' : 20})
   ax.set_ylabel(r'$\phi_2$', family='sans-serif', fontdict={'fontsize' : 20})

   img = ax.imshow(hist, extent=(-180,180,-180,180))
   fig.colorbar(img).set_label(r'$\Delta G(\phi_1, \phi_2)$')

   fig.savefig('cmap.png')

def distance_dielectric():
   """ Creates the plot for a distance-dependent dielectric """
   def eps_eff(eps, S):
      """ Returns a distance-dependent dielectric function of r """
      return lambda r: eps - (eps-1) / 2 * (r*r*S*S + 2*r*S + 2) * np.exp(-r*S)
   
   xdata = np.arange(0, 35, 0.1)
   s0_1 = eps_eff(78.5, 0.2)
   s0_5 = eps_eff(78.5, 0.5)
   s1 = eps_eff(78.5, 1.0)

   fig = plt.figure(5, figsize=(8,5))

   ax = fig.add_subplot(111)

   ax.grid(linewidth=1)
   ax.set_xlim((0,30))
   ax.set_ylim((0,80))
   ax.set_xlabel(r"Distance ($\AA$)", family='sans-serif',
                 fontdict={'fontsize' : 20})
   ax.set_ylabel(r"Effective Dielectric ($\varepsilon_{eff}$)",
                 family='sans-serif', fontdict={'fontsize' : 20})

   pl1, = ax.plot(xdata, s0_1(xdata), linewidth=3, linestyle='-', color='b')
   pl2, = ax.plot(xdata, s0_5(xdata), linewidth=3, linestyle='--', color='r')
   pl3, = ax.plot(xdata, s1(xdata), linewidth=3, linestyle='-.', color='g')

   ax.legend((pl1, pl2, pl3), ('S = 0.2', 'S = 0.5', 'S = 1.0'), loc=4)

   fig.savefig('DistanceDielectric.ps')

def cutoff_effects():
   
   CUTOFF = 16.0
   SWITCH_START = 8.0

   fig = plt.figure(6, figsize=(16,10))

   # Create all of the axes
   trunce = fig.add_subplot(231, xlim=(1,20), ylim=(-100,5))
   truncf = fig.add_subplot(234, xlim=(1,20), ylim=(-10,100))
   switche = fig.add_subplot(232, xlim=(1,20), ylim=(-100,5))
   switchf = fig.add_subplot(235, xlim=(1,20), ylim=(-10,100))
   shifte = fig.add_subplot(233, xlim=(1,20), ylim=(-100,5))
   shiftf = fig.add_subplot(236, xlim=(1,20), ylim=(-10,100))

   # Make a list for handling them easily
   grids = [trunce, truncf, switche, switchf, shifte, shiftf]

   # Show the grid
   for grid in grids: grid.grid(lw=1)

   xd = np.arange(1, 20, 0.01)

   # The various potentials
   def no_cutoff(x): return -332.0522017 / x
   def hard_cutoff(x):
      if x <= CUTOFF: return no_cutoff(x)
      return 0.0
   def _switch(x):
      return ((CUTOFF**2 - x*x)**2 * (CUTOFF**2+2*x*x-3*SWITCH_START**2)) / (
               CUTOFF**2-SWITCH_START**2)**3
   def switch(x):
      if x < SWITCH_START: return no_cutoff(x)
      if x < CUTOFF:
         return no_cutoff(x) * _switch(x)
      return 0.0
   def shift(x):
      if x < CUTOFF:
         return no_cutoff(x) * (1 - (x/CUTOFF)**2)**2
      return 0.0

   # The derivatives of the various potentials
   def dno_cutoff(x): return 332.0522017 / (x*x)
   def dhard_cutoff(x):
      if x <= CUTOFF: return dno_cutoff(x)
      return 0.0
   def dswitch(x):
      # Do the numerical derivative here, since it's simpler (and I know
      # switch(x) is continuous)
      return (switch(x+1e-8) - switch(x-1e-8)) / 2e-8
   def dshift(x):
      if x < CUTOFF:
         return dno_cutoff(x)*(1-(x/CUTOFF)**2)**2 + no_cutoff(x)*(
                2*(1-(x/CUTOFF)**2)*(-2*(x/CUTOFF**2)))
      return 0.0

   # Define the font dictionaries
   titlefont = {'fontsize' : 25, 'family' : 'sans-serif'}
   axisfont = {'fontsize' : 20, 'family' : 'sans-serif'}

   # Plot everything
   trunce.plot(xd, _evaluate(xd,no_cutoff), 'k-', lw=3, label='No Cutoff')
   trunce.plot(xd, _evaluate(xd,hard_cutoff), 'r--', lw=2, label='Hard Cutoff')
   truncf.plot(xd, _evaluate(xd,dno_cutoff), 'k-', lw=3, label='No Cutoff')
   truncf.plot(xd, _evaluate(xd,dhard_cutoff), 'r--', lw=2, label='Hard Cutoff')

   switche.plot(xd, _evaluate(xd,no_cutoff), 'k-', lw=3, label='No Cutoff')
   switche.plot(xd, _evaluate(xd,switch), 'r--', lw=2, label='Switching Function')
   switchf.plot(xd, _evaluate(xd,dno_cutoff), 'k-', lw=3, label='No Cutoff')
   switchf.plot(xd, _evaluate(xd,dswitch), 'r--', lw=2, label='Switching Function')

   shifte.plot(xd, _evaluate(xd,no_cutoff), 'k-', lw=3, label='No Cutoff')
   shifte.plot(xd, _evaluate(xd,shift), 'r--', lw=2, label='Shifting Function')
   shiftf.plot(xd, _evaluate(xd,dno_cutoff), 'k-', lw=3, label='No Cutoff')
   shiftf.plot(xd, _evaluate(xd,dshift), 'r--', lw=2, label='Shifting Function')

   # Create the force subplots to zoom in on important parts
   subtrunc = plt.axes([1/6.3, 0.25, 0.5/3, 0.2])
   subtrunc.set_xlim((15,20))
   subtrunc.set_ylim((-1,2))
   subtrunc.grid(lw=1)
   subtrunc.plot(xd, _evaluate(xd, dno_cutoff), 'k-', lw=3)
   subtrunc.plot(xd, _evaluate(xd,dhard_cutoff), 'r--', lw=2)
   subtrunc.set_xticklabels(())

   subswitch = plt.axes([1/6.7+1/3, 0.25, 0.5/3, 0.2])
   subswitch.set_xlim((7, 20))
   subswitch.set_ylim((-1,10))
   subswitch.grid(lw=1)
   subswitch.plot(xd, _evaluate(xd, dno_cutoff), 'k-', lw=3)
   subswitch.plot(xd, _evaluate(xd, dswitch), 'r--', lw=3)
   subswitch.set_xticklabels(())

   subshift = plt.axes([1/7.2+2/3, 0.25, 0.5/3, 0.2])
   subshift.set_xlim((5, 20))
   subshift.set_ylim((-1,15))
   subshift.grid(lw=1)
   subshift.plot(xd, _evaluate(xd, dno_cutoff), 'k-', lw=3)
   subshift.plot(xd, _evaluate(xd, dshift), 'r--', lw=3)
   subshift.set_xticklabels(())

   # Now draw the rectangles
   truncrect = patches.Rectangle((15,-1), 5, 3, fill=False, lw=1, ls='solid',
                                 ec='k')
   switchrect = patches.Rectangle((7,-1), 13, 11, fill=False, lw=1, ls='solid', 
                                 ec='k')
   shiftrect = patches.Rectangle((5,-1), 15, 16, fill=False, lw=1, ls='solid', 
                                 ec='k')

   # Now add the necessary lines
   ar1 = Line2D([15,7.206], [2,35.433], c='k', ls='--')
   ar2 = Line2D([20,18.842], [2,35.433], c='k', ls='--')
   ar3 = Line2D([7,7.51277], [10, 35.433], c='k', ls='--')
   ar4 = Line2D([20,19.094], [10,35.433], c='k', ls='--')
   ar5 = Line2D([5,7.65601], [15,35.433], c='k', ls='--')
   ar6 = Line2D([20,19.2368], [15,35.433], c='k', ls='--')

   truncf.add_patch(truncrect); truncf.add_line(ar1); truncf.add_line(ar2)
   switchf.add_patch(switchrect); switchf.add_line(ar3); switchf.add_line(ar4)
   shiftf.add_patch(shiftrect); shiftf.add_line(ar5); shiftf.add_line(ar6)
   

   # Show the legend
   for grid in [trunce, shifte, switche]:
      grid.legend(loc='lower right')

   # Set titles on the energy plots
   trunce.set_title('A) Hard Cutoff', fontdict=titlefont)
   switche.set_title('B) Switching Function', fontdict=titlefont)
   shifte.set_title('C) Shifting Function', fontdict=titlefont)

   # Set axis labels on all force plots (x-axis) and the 2 trunc plots (y-axis)
   trunce.set_ylabel(r'Energy (kcal mol$^{-1}$)',
                     fontdict=axisfont)
   truncf.set_ylabel(r'Force (kcal mol$^{-1}\,\AA^{-1}$)',
                     fontdict=axisfont)
   truncf.set_xlabel(r'Distance ($\AA$)', fontdict=axisfont)
   switchf.set_xlabel(r'Distance ($\AA$)', fontdict=axisfont)
   shiftf.set_xlabel(r'Distance ($\AA$)', fontdict=axisfont)

   # Tight layout
   fig.tight_layout()

   fig.savefig('Cutoff.png')

def ewald():
   """ Plots each charge with its neutralizing charge distribution """
   ALPHA = 5.0
   def neutralizing(q, xeq, x):
      x -= xeq
#     return -q*ALPHA**3/(np.pi**1.5)*np.exp(-ALPHA*ALPHA*x*x)
      return -q*np.sqrt(2)*np.pi**-0.5*np.exp(-ALPHA*ALPHA*x*x)

   axisfont = {'fontsize' : 25, 'family' : 'sans-serif'}

   fig = plt.figure(7, figsize=(8,5))
   ax = fig.add_subplot(111)
   ax.set_xlim((0,10))
   ax.set_ylim((-1.1,1.1))
   ax.set_ylabel('Charge ($e^{-}$)', fontdict=axisfont)
   ax.set_xlabel('Distance ($\\AA$)', fontdict=axisfont)
   ax.grid(lw=1)

   # Plot the first charge with neutralizing distribution
   xdata = np.arange(1,2.01,0.05)
   ydata = _evaluate(xdata, lambda x: neutralizing(1, 1.5, x))

   chdi, = ax.plot(xdata, ydata, 'b-', lw=2)
   ch, = ax.plot([1.5,1.5], [0,1], 'r-', lw=2)

   # Plot the second charge with neutralizing distribution (at 2.1)
   xdata = np.arange(1.6, 2.61, 0.05)
   ydata = _evaluate(xdata, lambda x: neutralizing(-1, 2.1, x))

   ax.plot(xdata, ydata, 'b-', lw=2)
   ax.plot([2.1,2.1], [0,-1], 'r-', lw=2)

   # Plot the third charge with neutralizing distribution (at 3.2)
   xdata = np.arange(2.7, 3.71, 0.05)
   ydata = _evaluate(xdata, lambda x: neutralizing(-0.6, 3.2, x))

   ax.plot(xdata, ydata, 'b-', lw=2)
   ax.plot([3.2,3.2], [0,-0.6], 'r-', lw=2)

   # Plot the fourth charge with neutralizing distribution (at 5.0)
   xdata = np.arange(4.5, 5.51, 0.05)
   ydata = _evaluate(xdata, lambda x: neutralizing(0.8, 5.0, x))

   ax.plot(xdata, ydata, 'b-', lw=2)
   ax.plot([5.0,5.0], [0,0.8], 'r-', lw=2)

   # Plot the fifth charge with neutralizing distribution (at 6.3)
   xdata = np.arange(5.8, 6.81, 0.05)
   ydata = _evaluate(xdata, lambda x: neutralizing(0.7, 6.3, x))

   ax.plot(xdata, ydata, 'b-', lw=2)
   ax.plot([6.3,6.3], [0,0.7], 'r-', lw=2)

   # Plot the sixth charge with neutralizing distribution (at 8.6)
   xdata = np.arange(8.1, 9.11, 0.05)
   ydata = _evaluate(xdata, lambda x: neutralizing(-0.9, 8.6, x))

   ax.plot(xdata, ydata, 'b-', lw=2)
   ax.plot([8.6,8.6], [0,-0.9], 'r-', lw=2)

   ax.plot([0,10], [0,0], 'k-', lw=2)

   ax.legend((ch, chdi), ('Point Charge', 'Neutralizing Distribution'),
             loc='best')
   fig.tight_layout()
   fig.savefig('Ewald.ps')

def softcore():
   """ Create the hard core LJ figure """
   fig = plt.figure(8, figsize=(16,10))

   ax1 = fig.add_subplot(221)
   ax2 = fig.add_subplot(222)
   ax3 = fig.add_subplot(223)
   ax4 = fig.add_subplot(224)
   # Set up axes
   for ax in (ax1, ax2, ax3, ax4):
      ax.set_xlim((0,6))
      ax.set_ylim((-0.5,2.0))
      ax.set_xlabel('Nuclear Separation ($\AA$)', family='sans-serif',
                    fontdict={'fontsize' : 18})
      ax.set_ylabel(r'Energy ($kcal$  $mol^{-1}$)', family='sans-serif',
                    fontdict={'fontsize' : 18})
      ax.grid(lw=1)
   RMIN = 3.816
   EPS = 0.1094
   SIG = RMIN / 2 ** (1/6)
#  ACOEF = EPS * RMIN ** 12
#  BCOEF = 2.0 * EPS * RMIN ** 6

   def vdw(x, alpha, lam):
      return 4*EPS*(1-lam)*(1/(alpha*lam + (x/SIG)**6)**2 - 
               1/(alpha*lam + (x/SIG)**6))

   xdata = np.arange(-0.2,10,0.05)

   for ax, lam in zip((ax1,ax2,ax3,ax4), (0.0, 0.20, 0.80, 1.00)):
      ax.set_title(r'$\lambda = %.2f$' % lam,
                    fontdict={'fontsize' : 18})
      al1, = ax.plot(xdata,_evaluate(xdata,lambda x:vdw(x,0.1,lam)),color='k',lw=2)
      al2, = ax.plot(xdata,_evaluate(xdata,lambda x:vdw(x,0.5,lam)),color='b',lw=2)
      al3, = ax.plot(xdata,_evaluate(xdata,lambda x:vdw(x,1.0,lam)),color='r',lw=2)
      al4, = ax.plot(xdata,_evaluate(xdata,lambda x:vdw(x,2.0,lam)),color='g',lw=2)
      axis, = ax.plot([0,10], [0,0], color='k', lw=1)

      ax.legend((al1, al2, al3, al4, ),#al5),
                (r'$\alpha = 0.1$',
                 r'$\alpha = 0.5$',
                 r'$\alpha = 1.0$',
                 r'$\alpha = 2.0$',
                ), loc=1)

   fig.tight_layout()
   fig.savefig('SoftCore.png')

def hardcore():
   """ Create figure showing hard cores of disappearing atoms """
   fig = plt.figure(9, figsize=(8,5))

   ax = fig.add_subplot(111)
   # Set up axes
   ax.set_xlim((0,6))
   ax.set_ylim((-0.5,2.0))
   ax.set_xlabel('Nuclear Separation ($\AA$)', family='sans-serif',
                 fontdict={'fontsize' : 16})
   ax.set_ylabel(r'Energy ($kcal$  $mol^{-1}$)', family='sans-serif',
                 fontdict={'fontsize' : 16})
   ax.grid(lw=1)

   RMIN = 3.816
   EPS = 0.1094
   SIG = RMIN / 2 ** (1/6)
#  ACOEF = EPS * RMIN ** 12
#  BCOEF = 2.0 * EPS * RMIN ** 6

   def vdw(x, lam):
      xos6 = (x/SIG)**6
      return 4*EPS*(1-lam)*(1/(xos6*xos6) - 1/(xos6))

   xdata = np.arange(0.1,10,0.05)
   x2 = np.arange(-1,10,0.1)

#  ax.set_title(r'$\lambda = %.2f$' % lam)
   al1, = ax.plot(xdata,_evaluate(xdata,lambda x:vdw(x,0)),color='k',lw=2)
   al2, = ax.plot(xdata,_evaluate(xdata,lambda x:vdw(x,0.5)),color='b',lw=2)
   al3, = ax.plot(xdata,_evaluate(xdata,lambda x:vdw(x,0.9)),color='r',lw=2)
   al4, = ax.plot(xdata,_evaluate(xdata,lambda x:vdw(x,.99)),color='g',lw=2)
   al5, = ax.plot(x2,_evaluate(x2,lambda x:vdw(x,1.0)),color='m',lw=2)
   axis, = ax.plot([0,10], [0,0], color='k', lw=1)

   ax.legend((al1, al2, al3, al4, al5),
             (r'$\lambda = 0.0$',
              r'$\lambda = 0.5$',
              r'$\lambda = 0.9$',
              r'$\lambda = 0.99$',
              r'$\lambda = 1.0$',
             ), loc=1)

   fig.tight_layout()
   fig.savefig('HardCore.ps')
#  plt.show()

def umbrella():
   """ Creates the umbrella sampling PMF example """
   fig = plt.figure(10, figsize=(8,5))

   ax = fig.add_subplot(111)

   # Set up axes
   ax.set_xlim((0,8))
   ax.set_ylim((0,10))

   ax.set_xlabel('Reaction Coordinate', family='sans-serif',
                 fontdict={'fontsize' : 16})
   ax.set_ylabel(r'Energy ($k_B T$)', family='sans-serif',
                 fontdict={'fontsize' : 16})
   ax.grid(lw=1)

   def shift(func):
      SHIFT = -2
      def new_func(x):
         return func(x-SHIFT)
      return new_func

   @shift
   def pmf(x): return 3*(np.cos(1.5*x)+1)

   def bias(k, x, center):
      return 0.5*k*(x-center)*(x-center) + pmf(center)

   def biased(k, x, center):
      return pmf(x) + 0.5*k*(x-center)*(x-center)

   xdata = np.arange(0,8,0.01)

   ftp = 4 / 3 * np.pi - 2# four-thirds pi - 2
   K = 5
   pl1, = ax.plot(xdata, _evaluate(xdata, pmf), 'k-', lw=3)
   pl2, = ax.plot(xdata, _evaluate(xdata, lambda x: bias(K, x, 6.3)), 'r-.',
                   lw=1)
   ydata = _evaluate(xdata, lambda x: biased(K, x, 6.3))
   min1 = ydata.min()
   ydata -= min1
   pl3, = ax.plot(xdata, ydata, 'r--', lw=2)
   pl4, = ax.plot(xdata, _evaluate(xdata, lambda x: bias(K, x, ftp)), 'b-.',
                   lw=1)
   ydata = _evaluate(xdata, lambda x: biased(K, x, ftp))
   min2 = ydata.min()
   ydata -= min2
   pl5, = ax.plot(xdata, ydata, 'b--', lw=2)

   # Create the annotations
   ax.annotate('Unbiased\nPMF', xy=(1.2, pmf(1.2)), xytext=(0.2,6.1),
           xycoords='data', textcoords='data', arrowprops=dict(arrowstyle='->'))
   ax.annotate('Umbrella 1', xy=(2.6, bias(K, 2.6, ftp)), xytext=(2.01,8.2),
           xycoords='data', textcoords='data', arrowprops=dict(arrowstyle='->'))
   ax.annotate('Biased\nPMF 1', xy=(2.2, biased(K, 2.2, ftp)-min2),
           xytext=(2.2,2.2), xycoords='data', textcoords='data',
           arrowprops=dict(arrowstyle='->'))
   ax.annotate('Umbrella 2', xy=(6.9, bias(K, 6.9, 6.3)), xytext=(6.01,8.2),
           xycoords='data', textcoords='data', arrowprops=dict(arrowstyle='->'))
   ax.annotate('Biased\nPMF 2', xy=(6.5, biased(K, 6.5, 6.3)-min1),
           xytext=(6.2,2.2), xycoords='data', textcoords='data',
           arrowprops=dict(arrowstyle='->'))

#  plt.show()
   fig.tight_layout()
   fig.savefig('FreeEnergyProfile.ps')

if __name__ == '__main__':
   """ Determine which plots to make """
   parser = ArgumentParser()
   group = parser.add_argument_group('Chapter 1 Figures')
   group.add_argument('--time-step', dest='timestep', default=False,
                      action='store_true', help='Create the time step plot')
   group.add_argument('--hydrogen-atom', dest='hydrogen', default=False,
                      action='store_true', help='''Create the plot for the
                      hydrogen atom.''')
   group.add_argument('--lennard-jones', dest='lj', default=False,
                      action='store_true', help='''Create the plot for the
                      Lennard-Jones potential.''')
   group.add_argument('--cmap', dest='cmap', default=None, metavar='FILE',
                      help='Generate a CMAP plot with the given input file.')
   group = parser.add_argument_group('Chapter 2 Figures')
   group.add_argument('--distance-dielectric', dest='distdiel', default=False,
                      action='store_true', help='''Create the plot for the
                      distance-dependent dielectric figure.''')
   group.add_argument('--cutoff-effects', dest='cutoff', default=False,
                      action='store_true', help='''Create the plot showing the
                      effect of hard cutoffs and shifting/switching
                      functions''')
   group.add_argument('--ewald', dest='ewald', default=False,
                      action='store_true', help='''Create the figure
                      demonstrating how the Ewald sum works''')
   group.add_argument('--hard-core', dest='hardcore', default=False,
                      action='store_true', help='''Create the figure showing how
                      a vanishing atom has a hard vdW core.''')
   group.add_argument('--soft-core', dest='softcore', default=False,
                      action='store_true', help='''Create the figure showing
                      what soft cores are.''')
   group.add_argument('--umbrella', dest='umbrella', default=False,
                      action='store_true', help='''Create the figure showing the
                      effect of umbrellas for umbrella sampling''')
   
   opt = parser.parse_args()

   if opt.timestep:
      timestep()
   if opt.hydrogen:
      hydrogen()
   if opt.lj:
      LennardJones()
   if opt.cmap is not None:
      cmap(opt.cmap)
   if opt.distdiel:
      distance_dielectric()
   if opt.cutoff:
      cutoff_effects()
   if opt.ewald:
      ewald()
   if opt.hardcore:
      hardcore()
   if opt.softcore:
      softcore()
   if opt.umbrella:
      umbrella()
