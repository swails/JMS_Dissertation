"""
This program uses matplotlib to make the various plots included in my
dissertation.
"""
from __future__ import division

from argparse import ArgumentParser
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats.kde import gaussian_kde

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
   x1 = xdata[7:20]
   y1 = ydata[7:20] - min(ydata)
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
