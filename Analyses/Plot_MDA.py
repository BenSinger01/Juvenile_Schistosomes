## Plot MDA
## Benjamin John Singer, July 2024.
## Code to plot prevalence over time, e.g. figure 1 in Singer et al. 2024

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import sys

pd.set_option('display.max_rows', None)

# font formatting
plt.rcParams['font.size'] = '10'
plt.rcParams["font.family"] = "Times New Roman"

epgpwp = int(sys.argv[1])
juvenile_duration = int(sys.argv[2])
seasonality = int(sys.argv[3])
t_start = int(sys.argv[4])
dyn_type = sys.argv[5]

coverage = int(sys.argv[6])
nonadherence = int(sys.argv[7])
eff_scale = int(sys.argv[8])
wash = int(sys.argv[9])

nyears = int(sys.argv[10])

# generate plots
fig,axes=plt.subplots(3,3,sharey=True,sharex=True,figsize=(6.5,6.5))
for statn in range(3):
	stat = [100,50,0][statn]
	statname = ['Static','Semi-dynamic','Dynamic'][statn]
	for settn in range(3):
		sett = ['Low','Moderate','High'][settn]
		for stratn in range(5):
			strat = ['Single','NovelA','NovelB','NovelC','Double'][stratn]
			followup_week = [5,5,5,5,11][stratn]
			stratname = ['1Rx PZQ','Novel A','Novel B','Novel C','2Rx PZQ'][stratn]
			color = ['#648FFF','#FFB000','#DC267F','#FF832B','#785EF0'][stratn]

			# load data
			if seasonality==0:
				if (epgpwp==4)&(juvenile_duration==6):
					extension = ''
				elif (epgpwp!=4)&(juvenile_duration==6):
					extension = '_'+epgpwp+'epgpwp'
				elif (epgpwp==4)&(juvenile_duration!=6):
					extension = '_'+juvenile_duration+'juve'
			else:
				extension = seasonality+'seasonal'+t_start

			data = pd.read_csv('Data/MDA_sims/'+sett+extension+'/'+strat+'_'+str(stat)+dyn_type+'_'+str(coverage)+'cov'+str(nonadherence)+'_'+str(eff_scale)+'scale_'+str(wash)+'wash_Mean.csv',index_col=False)
			if (strat=='Single')|(strat=='Double'):
				data_up = pd.read_csv('Data/MDA_sims/'+sett+extension+'/'+strat+'_'+str(stat)+dyn_type+'_'+str(coverage)+'cov'+str(nonadherence)+'_'+str(eff_scale)+'scale_'+str(wash)+'wash_Upper.csv',index_col=False)
				data_lo = pd.read_csv('Data/MDA_sims/'+sett+extension+'/'+strat+'_'+str(stat)+dyn_type+'_'+str(coverage)+'cov'+str(nonadherence)+'_'+str(eff_scale)+'scale_'+str(wash)+'wash_Lower.csv',index_col=False)
			else:
				data_up = data_lo = data

			median_prev = data.loc[(data['t']<=nyears*52+followup_week),['t','Strategy','Prevalence 3dx2s']].pivot_table(values='Prevalence 3dx2s',index='t',columns='Strategy',aggfunc=np.median)
			bot_prev = data_up.loc[(data_up['t']<=nyears*52+followup_week),['t','Strategy','Prevalence 3dx2s']].pivot_table(values='Prevalence 3dx2s',index='t',columns='Strategy',aggfunc=lambda x: x.quantile(.25))
			top_prev = data_lo.loc[(data_lo['t']<=nyears*52+followup_week),['t','Strategy','Prevalence 3dx2s']].pivot_table(values='Prevalence 3dx2s',index='t',columns='Strategy',aggfunc=lambda x: x.quantile(.75))

			axes[statn,settn].plot(median_prev,label=stratname,c=color,linewidth=1)

			median_juves = data.loc[(data['t']<=nyears*52),['t','Strategy','Juvenile Prevalence']].pivot_table(values='Juvenile Prevalence',index='t',columns='Strategy',aggfunc=np.median)
			axes[statn,settn].plot(median_juves,alpha=0.5,linestyle=':',c=color,linewidth=1)

			print(sett,stratname,statname)
			print("{:.1f}".format(median_prev[strat].iloc[-1])+'% ('+"{:.1f}".format(bot_prev[strat].iloc[-1])+'–'+"{:.1f}".format(top_prev[strat].iloc[-1])+')')
		axes[0,settn].set_title(sett+' endemicity')
	axes[statn,2].set_ylabel(statname,fontsize=13)
	axes[statn,2].yaxis.set_label_position('right')
handles, labels = plt.gca().get_legend_handles_labels()
order = [0,4,1,2,3]
axes[0,0].legend([handles[idx] for idx in order],[labels[idx] for idx in order])
axes[1,0].set_ylabel(r'Prevalence ($\%$)')
axes[0,0].set_xticks(range(0,(nyears+1)*52,52),labels=range(1,nyears+2))
axes[2,1].set_xlabel('Year')
plt.tight_layout()
# save figure
plt.savefig('Figures/MDA_'+str(wash)+'wash_'+str(coverage)+'cov'+str(nonadherence)+'_'+dyn_type+'_'+str(seasonality)+'seasonal'+str(t_start)+'_'+str(epgpwp)+'epgpwp_'+str(juvenile_duration)+'juve_'+str(eff_scale)+'scale.svg',transparent=True)
