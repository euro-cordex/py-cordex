"""contour plot example script
"""

import cordex.plot.animation as cxani
# pyplot needed to plot the dataset, but animation only needed much further down.
from matplotlib import pyplot as plt, animation
# This is needed to display graphics calculated outside of jupyter notebook
from IPython.display import HTML, display


# We need to create a function that updates the values for the colormesh, as well as the title.
figu = plt.figure(figsize=(12,6))
infile ='data/video.nc'
caxt,var = cxani.contour_video(infile,'t2m')
def animate(frame):
    caxt.set_array(var[frame,:,:].values.flatten())
    plt.title("Time = " + str(var.coords['time'].values[frame])[:13])
 

# Finally, we use the animation module to create the animation.
ani = animation.FuncAnimation(
    figu,             # figure
    animate,         # name of the function above
    frames=6,       # Could also be iterable or list
    interval=800     # ms between frames
)

#HTML(ani.to_jshtml()) 

#save it as .mp4
ani.save('results/Python_Animation_04.mp4')
#display(HTML("<video controls><source src='Python_Animation_04.mp4' type='video/mp4'></#video>"))