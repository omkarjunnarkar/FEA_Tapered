import pandas as pd
import matplotlib.pyplot as plt

strain=pd.read_csv("strain.csv")
stress=pd.read_csv("stress.csv")
u=pd.read_csv("displacement.csv")

fig, axs = plt.subplots(1,2,figsize=(25,25))
axs[0].set(xlabel = 'Strain')
axs[0].set(ylabel = 'Stress')
axs[1].set(xlabel = 'Time Step')
axs[1].set(ylabel = 'Displacement of Node')
axs[0].plot(strain,stress)
axs[1].plot(u)
plt.show()