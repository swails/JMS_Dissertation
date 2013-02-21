"""
This program uses matplotlib to make the various plots included in my
dissertation.
"""
from __future__ import division

from argparse import ArgumentParser
from matplotlib import pyplot as plt
import numpy as np


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

if __name__ == '__main__':
   """ Determine which plots to make """
   parser = ArgumentParser()
   parser.add_argument('--time-step', dest='timestep', default=False,
                       action='store_true', help='Create the time step plot')

   
   opt = parser.parse_args()

   if opt.timestep:
      timestep()

