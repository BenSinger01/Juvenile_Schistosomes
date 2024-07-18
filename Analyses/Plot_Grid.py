## Plot grids
## Benjamin John Singer, July 2024.
## Code to create grid figures, e.g. figures 2 and 3 in Singer et al. 2024

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import sys

plt.rcParams['font.size'] = '10'
plt.rcParams["font.family"] = "Times New Roman"

epgpwp = int(sys.argv[1])
juvenile_duration = int(sys.argv[2])
seasonality = int(sys.argv[3])
t_start = int(sys.argv[4])

static_force = int(sys.argv[5])
dyn_type = sys.argv[6]
coverage = int(sys.argv[7])
nonadherence = int(sys.argv[8])

n_years = int(sys.argv[9]) # How many years of annual MDA? (0 to 5)
delay = int(sys.argv[10]) # How many days after last MDA round are outcomes measured? (0 or 5)

eff_as = [100,95,90]
eff_js = [0,50,90,95,100]

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(6.5,4), sharex=True, sharey=True)
ims = []
grid = np.zeros((3,5))

for settn,sett in enumerate(['Low','High']):
    for valn,val in enumerate(['Prevalence 3dx2s','Nonlight Prevalence 3dx2s']):
        for i,eff_a in enumerate(eff_as):
            for j,eff_j in enumerate(eff_js):
                if seasonality==0:
                    if (epgpwp==4)&(juvenile_duration==6):
                        data = pd.read_csv('Data/MDA_sims/Grid_sims/'+sett+'/'+str(eff_a)+'a'+str(eff_j)+'j_'+str(static_force)+dyn_type+'_'+str(coverage)+'cov'+str(nonadherence)+'.csv')
                    elif (epgpwp!=4)&(juvenile_duration==6):
                        data = pd.read_csv('Data/MDA_sims/Grid_sims/'+sett+'_'+str(epgpwp)+'epgpwp/'+str(eff_a)+'a'+str(eff_j)+'j_'+str(static_force)+dyn_type+'_'+str(coverage)+'cov'+str(nonadherence)+'.csv')
                    elif (epgpwp==4)&(juvenile_duration!=6):
                        data = pd.read_csv('Data/MDA_sims/Grid_sims/'+sett+'_'+str(juvenile_duration)+'juve/'+str(eff_a)+'a'+str(eff_j)+'j_'+str(static_force)+dyn_type+'_'+str(coverage)+'cov'+str(nonadherence)+'.csv')
                else:
                    data = pd.read_csv('Data/MDA_sims/Grid_sims/'+sett+'_'+str(seasonality)+'seasonal'+str(t_start)+'/'+str(eff_a)+'a'+str(eff_j)+'j_'+str(static_force)+dyn_type+'_'+str(coverage)+'cov'+str(nonadherence)+'.csv')
                if val=='Nonlight Prevalence 3dx2s':
                    data['Nonlight Prevalence 3dx2s'] = data['Prevalence 3dx2s']*(100-data['Light 3dx2s'])/100
                value_mean = data.loc[(data['t']==n_years*52+delay),['t',val]].pivot_table(values=val,index='t',aggfunc=np.mean)
                grid[i,j] = float(value_mean.iloc[0])
        im = axes[valn,settn].imshow(grid, cmap='viridis', interpolation='nearest')
        ims.append(im)
        axes[valn,settn].set_xticks(np.arange(len(eff_js)))
        axes[valn,settn].set_yticks(np.arange(len(eff_as)))
        axes[valn,settn].set_xticklabels(eff_js)
        axes[valn,settn].set_yticklabels(eff_as)

# separate colorbars for each panel
fig.colorbar(ims[0], ax=axes[0,0])
fig.colorbar(ims[1], ax=axes[1,0])
fig.colorbar(ims[2], ax=axes[0,1],label='Prevalence (%)')
fig.colorbar(ims[3], ax=axes[1,1],label='Moderate and heavy\nprevalence (%)')

# label axes
axes[0,0].set_ylabel('Adult efficacy (%)')
axes[1,0].set_ylabel('Adult efficacy (%)')
axes[1,0].set_xlabel('Juvenile efficacy (%)')
axes[1,1].set_xlabel('Juvenile efficacy (%)')

# add labels for settings
for ax, col in zip(axes[0], ['Low','High']):
    ax.set_title(col+" endemicity")

# save figure
if delay==0:
    dl = 'pre'
else:
    dl = 'post'
plt.savefig('Figures/grid_plot_y'+str(n_years+1)+'_'+dl+str(static_force)+dyn_type+'_'+str(coverage)+'cov'+str(nonadherence)+'_static_'+str(seasonality)+'seasonal'+str(t_start)+'_'+str(epgpwp)+'epgpwp_'+str(juvenile_duration)+'juve.svg',transparent=True)
