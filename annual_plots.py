from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

pd.set_option('display.max_rows', None)

# font formatting
plt.rcParams['font.size'] = '10'
plt.rcParams["font.family"] = "Times New Roman"


dyn_type = 'no_saturation'
epgpwp = 4
# choose coverage
cov = 75
nonadherence = 10
# length of time
nyears = 4
# followup?
followup = True
# juves=8


# generate plots
# for cov in [90]:
fig,axes=plt.subplots(3,3,sharey=True,sharex=True,figsize=(6.5,6.5))
for statn in range(3):
	stat = [100,50,0][statn]
	statname = ['Static','Semi-dynamic','Dynamic'][statn]
	for settn in range(3):
		sett = ['Low'
		,'Moderate'
		,'High'][settn]
		for stratn in range(5):
			strat = ['single','novel1','novelC','novelB','double'][stratn]
			followup_week = [5,5,5,5,11][stratn]
			stratname = ['1Rx PZQ','Novel A','Novel B','Novel C','2Rx PZQ'][stratn]
			color = ['#648FFF','#FFB000','#DC267F','#FF832B','#785EF0'][stratn]
			# load data

			data = pd.read_csv('Tables_smallsample/'+strat+'_'+sett+'_'+str(stat)+dyn_type+'_'+str(cov)+'cov'+str(nonadherence)+'_'+str(epgpwp)+'epgpwp_Mean.csv',index_col=False)
			data_up = pd.read_csv('Tables_smallsample/'+strat+'_'+sett+'_'+str(stat)+dyn_type+'_'+str(cov)+'cov'+str(nonadherence)+'_'+str(epgpwp)+'epgpwp_Upper.csv',index_col=False)
			data_lo = pd.read_csv('Tables_smallsample/'+strat+'_'+sett+'_'+str(stat)+dyn_type+'_'+str(cov)+'cov'+str(nonadherence)+'_'+str(epgpwp)+'epgpwp_Lower.csv',index_col=False)

			print(sett,stratname, statname)
			# print(data.loc[(data_lo['t']==nyears*52+followup*followup_week)&(data_lo['Setting']==sett)&(data_lo['Coverage']==cov)&(data_lo['Static Force']==stat)].shape)

			median_prev = data.loc[(data['t']<=nyears*52+followup*followup_week)&(data['Setting']==sett)&(data['Coverage']==cov)&(data['Static Force']==stat),['t','Strategy','Prevalence 3dx2s']].pivot_table(values='Prevalence 3dx2s',index='t',columns='Strategy',aggfunc=np.median)
			bot_prev = data_up.loc[(data_up['t']<=nyears*52+followup*followup_week)&(data_up['Setting']==sett)&(data_up['Coverage']==cov)&(data_up['Static Force']==stat),['t','Strategy','Prevalence 3dx2s']].pivot_table(values='Prevalence 3dx2s',index='t',columns='Strategy',aggfunc=lambda x: x.quantile(.25))
			top_prev = data_lo.loc[(data_lo['t']<=nyears*52+followup*followup_week)&(data_lo['Setting']==sett)&(data_lo['Coverage']==cov)&(data_lo['Static Force']==stat),['t','Strategy','Prevalence 3dx2s']].pivot_table(values='Prevalence 3dx2s',index='t',columns='Strategy',aggfunc=lambda x: x.quantile(.75))
			mean_prev = data.loc[(data['t']<=nyears*52+followup*followup_week)&(data['Setting']==sett)&(data['Coverage']==cov)&(data['Static Force']==stat),['t','Strategy','Prevalence 3dx2s']].pivot_table(values='Prevalence 3dx2s',index='t',columns='Strategy',aggfunc=np.mean)

			print("{:.1f}".format(median_prev[strat].iloc[-1])+'% ('+"{:.1f}".format(bot_prev[strat].iloc[-1])+'–'+"{:.1f}".format(top_prev[strat].iloc[-1])+')')
			axes[statn,settn].plot(median_prev[strat],label=stratname,c=color,linewidth=1)

			# median_light = data.loc[(data['t']<=nyears*52)&(data['Setting']==sett)&(data['Coverage']==cov)&(data['Static Force']==stat),['t','Strategy','Light 3dx2s']].pivot_table(values='Light 3dx2s',index='t',columns='Strategy',aggfunc=np.median)
			# nonlight_prev = median_prev[strat] - (median_prev[strat]*median_light[strat])/100
			# axes[statn,settn].plot(nonlight_prev,label=stratname,c=color)

			# median_AMI = data.loc[(data['t']<=nyears*52)&(data['Setting']==sett)&(data['Coverage']==cov)&(data['Static Force']==stat),['t','Strategy','AMI 3dx2s']].pivot_table(values='AMI 3dx2s',index='t',columns='Strategy',aggfunc=np.median)
			# overall_AMI = (median_prev[strat]*median_AMI[strat])/100
			# axes[statn,settn].plot(overall_AMI,label=stratname,c=color)

			# median_GMI = data.loc[(data['t']<=nyears*52)&(data['Setting']==sett)&(data['Coverage']==cov)&(data['Static Force']==stat),['t','Strategy','GMI 3dx2s']].pivot_table(values='GMI 3dx2s',index='t',columns='Strategy',aggfunc=np.median)
			# axes[statn,settn].plot(median_GMI,label=stratname,c=color)

			median_juves = data.loc[(data['t']<=nyears*52)&(data['Setting']==sett)&(data['Coverage']==cov)&(data['Static Force']==stat),['t','Strategy','Juvenile Prevalence']].pivot_table(values='Juvenile Prevalence',index='t',columns='Strategy',aggfunc=np.median)
			axes[statn,settn].plot(median_juves[strat],alpha=0.5,linestyle=':',c=color,linewidth=1)

		axes[0,settn].set_title(sett+' endemicity')
	axes[statn,2].set_ylabel(statname,fontsize=13)
	axes[statn,2].yaxis.set_label_position('right')
handles, labels = plt.gca().get_legend_handles_labels()
order = [0,4,1,2,3]
axes[0,0].legend([handles[idx] for idx in order],[labels[idx] for idx in order])
# axes[0,0].set_ylim([-2.5,54.5])
# axes[0,0].set_yticks([0,15,30,45])
axes[1,0].set_ylabel(r'Prevalence ($\%$)')
# axes[1,0].set_ylabel(r'Moderate and heavy prevalence ($\%$)')
# axes[1,0].set_ylabel(r'Mean intensity (EPG)')
# axes[1,0].set_ylabel(r'Geometric mean intensity (EPG)')
axes[0,0].set_xticks(range(0,(nyears+1)*52,52),labels=range(1,nyears+2))
axes[2,1].set_xlabel('Year')
plt.tight_layout()
# save figure
case = "Mean"
plt.savefig('annual_'+str(cov)+'cov'+str(nonadherence)+'_'+case+'_'+dyn_type+'_'+str(epgpwp)+'epgpwp.svg'
	# ,dpi=300
	,transparent=True
	)
