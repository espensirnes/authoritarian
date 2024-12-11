import pandas as pd


import numpy as np
from matplotlib import pyplot as plt
import paneltime as pt
import os
os.chdir(os.path.dirname(__file__))
import networkdrawing as nw
import markov as mkv
import pickle


from scipy import stats


def execute(desc, fname):
	df_orig = pd.read_csv(fname)
	df, df_lags = create_lags(df_orig)
	#df_sim, df_sim_lags  = generate_paths(190,40, 0.3, 10)
	boot_strap_mean_thresholds = 0.3, 0.7
	if True:
		print('*************** ANALYZING EMPIRICAL DATA *****************************')
		emp = analyze(df, df_lags, 'empirical', boot_strap_mean_thresholds)
	print('*************** ANALYZING SIMULATED DATA *****************************')

	dmpfile = f'results/sim_out_{desc}.dmp'
	arr = []
	if os.path.exists(dmpfile):
		with open(dmpfile, 'rb') as f:
			arr = pickle.load(f)

	analyze_sim(arr,emp, desc)

	for i in range(1000): 
		print(f'Running sim {i}')
		#sim boostrapping:
		df_bs, df_bs_lags = bootstrap(df_orig, df)
		#analysis of sim data:
		d = analyze(df_bs, df_bs_lags, 'simulation', boot_strap_mean_thresholds)
		arr.append(d)
		with open(dmpfile, 'wb') as f:
			pickle.dump(arr, f)
		if len(arr)>=1000:
			break
	analyze_sim(arr,emp, desc)
		
def analyze_sim(bsdata,emp, desc):
	if len(bsdata)==0:
		return
	coefs_cw = []
	coefs_after_cw = []
	markovs = []
	markovs_thr = []
	condv = []
	rel_mid = []
	crackdowns = []
	crackdowns_rat = []
	coefs_chris = []
	coefs_chris = []

	for d in bsdata:
		#coefs_cw.append(d['reg coef']['Cold war'])
		#coefs_after_cw.append(d['reg coef']['Post cold war'])
		if 'Chris' in d:
			coefs_chris.append(d['Chris'])
		markovs.append(d['markov_matrix'])
		markovs_thr.append(d['thresholds'])
		condv.append(d['cond_vol'])
		rel_mid.append(d['cond_vol_rel_mid'])
		crackdowns.append(d['mean crackdown'])
		crackdowns_rat.append(d['crackdown ratio'])

	#handle_regression(coefs_cw, emp['reg coef']['Cold war'], 'Cold war', desc)
	#handle_regression(coefs_cw, emp['reg coef']['Post cold war'], 'Post cold war', desc)
	handle_markov(markovs, emp, markovs_thr, desc)
	handle_crackdown(crackdowns, crackdowns_rat, emp, desc)
	handle_cond_vol(condv, rel_mid, emp, desc)

	coefs_chris = np.array(coefs_chris)

	
	if False:
		print('Regressjon Chris')
		print(np.vstack((est['Chris'], np.mean(coefs_chris,0), np.std(coefs_chris,0))))
		fig, ax = define_plot_reg(desc)
		plot_reg(lambda x: sum(est['Chris'][:4]*np.array([1, x, x**2, x**3])), 'Empirical', ax, fig, '#8B0000' )
		ax.legend()
		fig.show()
		a=0

def cond_vol(df_lags):
	n_groups = 3
	x = np.zeros(n_groups)
	y = np.zeros(n_groups)
	gap = 10/n_groups
	for i in range(n_groups):
		left = [0,3,7][i]
		right = [3,7,10.001][i]
		s = df_lags['oppression_D_L0'][(df_lags['oppression']>=left)*(df_lags['oppression']<right)]
		x[i] = np.std(s)
		y[i] = (left + right)*0.5

	rel_mid = 2* x[1]/(x[0]+x[-1])
	return x ,y, rel_mid
	
def handle_cond_vol(bsarr, rel_mid, est, desc):
	empirical = est['cond_vol']
	y = est['cond_vol_ifphol']
	rm = est['cond_vol_rel_mid']
	rm_s = np.std(rel_mid,0, ddof = 1)/(len(rel_mid)-1)**0.5

	bootstrap = np.mean(bsarr,0)
	s = np.std(bsarr,0, ddof = 1)
	
	fig, ax = plt.subplots()
	ax.bar(y, bootstrap, label='Bootstrap', yerr = 1.96*s)
	ax.bar(y+0.1, empirical, label='Empirical')

	# Adding labels and title
	ax.set_ylabel('Volatility')
	ax.set_xlabel(desc)
	ax.legend()
	ax.set_title(f"Bootstrap and empirical volatility for \nthree intervals of {desc} (0-3, 3-7 and 7-10)]")
	fig.savefig(f"results/volatility.png")
	fig.show()

	
	fig, ax = plt.subplots()

	ax.bar(['Empirical', 'Bootstrap'], np.array([rm, np.mean(rel_mid)]), yerr = 1.96*rm_s)

	# Adding labels and title
	ax.set_ylabel('Volatility')
	ax.set_xlabel(desc)

	ax.legend()
	ax.set_title(f"Relative difference in {desc} at 0-3 and 7-10, relative to volatility at 3-7")


	fig.savefig(f"results/volatility_rel.png")
	fig.show()

def handle_markov(markov, d, markovs_thr, desc):

	e = d[f'markov_matrix {desc}']
	m = np.mean(markov,0)
	s = np.std(markov,0)

	mkv.print_table(e, s, f'Markov transition probabilities empirical {desc}')
	mkv.print_table(m, s, f'Markov transition probabilities bootstrap {desc}')

	nw.draw_network(e,f'Empirical {desc}', s)
	nw.draw_network(m,f'Bootstrap means {desc}', s)

	m = np.mean(markovs_thr,0)
	s = np.std(markovs_thr,0)

	print('Estimated Markov transition thresholds:')
	print(d['thresholds'])
	print('Bootstrapped Markov transition thresholds:')
	print(m)
	print('Bootstrap standard errors Markov transition thresholds:')
	print(s)



def handle_crackdown(crackdowns, crackdowns_rat, d, desc):

	e = d['mean crackdown'].iloc[::-1]
	m = pd.concat(crackdowns, axis=1).mean(axis=1).iloc[::-1]
	s = pd.concat(crackdowns, axis=1).std(axis=1).iloc[::-1]
	ind = np.arange(len(m))
	width = 0.35

	fig, ax = plt.subplots()


	ax.bar(ind - width/2, e.values, width = width, label='Empirical')
	ax.bar(ind + width/2, m.values, width = width, yerr = 1.96*s.values, label='Bootstrap')

	# Adding labels and title
	ax.set_ylabel(desc)
	ax.set_xlabel('Years since crackdown')
	

	# Setting the x-ticks to be in the center of the two sets of bars
	ax.set_xticks(ind)
	ax.set_xticklabels(m.index.str.replace('oppression_L',''), rotation=90)
	fig.subplots_adjust(bottom=0.15)
	fig.show()
	ax.legend(frameon = False)

	fig.savefig(f"results/Crackdown_{desc}.png")

	e = d['crackdown ratio']
	m = np.mean(crackdowns_rat,0)
	s = np.std(crackdowns_rat,0)
	print('Crackdown:')
	print('Per cent of the sample satisfying the criteria')
	print(f'Empirical:{100*e}')
	print(f'Bootstrap:{100*m}')
	print(f'Bootstrap standard error:{100*s}')
	

def handle_regression(a, est,prefix):
	z = 1.96

	a = np.array(a)
	a=a[(a[:,-1]<50)*(a[:,-1]>1e-50)]
	m = np.mean(a, 0)
	s = np.std(a, 0)

	print('Estimated mean, simulated mean and bootstrap standard errors')
	print(np.vstack((est, m,s))[:,:3])
	fig, ax = define_plot_reg()
	plot_reg(lambda x: predfunc(x, est, s, 0), 'Empirical', ax, fig, '#8B0000' )
	plot_reg(lambda x: predfunc(x, m, s, 0), 'Bootstrap', ax, fig, '#6B8E23' )
	plot_reg(lambda x: predfunc(x, est, s, -z), 'Lower bound', ax, fig, 'lightgrey' )
	plot_reg(lambda x: predfunc(x, est, s, z), 'Upper bound', ax, fig, 'lightgrey' )
	ax.legend()
	fig.show()
	fig.savefig(f"results/predicted_change_{prefix}.png")


def predfunc(x, m, s, i):
	return (
		(m[0]+i*s[0])
			+(m[1]+i*s[1])*x
			+(m[2]-i*s[2])*x**2
			)

def define_plot_reg(desc):
	plt.ioff()
	fig, ax = plt.subplots()
	ax.set_xlabel(desc)
	ax.set_ylabel(f'{desc} change')
	return fig, ax

def plot_reg(f, label, ax, fig, color):
	x = np.linspace(0,10)
	y= [f(xi) for xi in x]
	ax.plot(x,y, label=label, color = color)
	
	




def bootstrap(df, df_change):
	"""Creates a new df where changes in oppression are randomly and unsystematically selected from 
	all changes in the dataset"""
	df = df.copy()
	
	#randomizing start value:
	df1 = df.sort_values(by=['country', 'year'])
	df1 = df1.groupby('country').first()
	for i,r in df1.iterrows():
		df.loc[(df['country']==i)*(df['year']==r['year']), 'oppression']=np.random.randint(60)/6

	#bootstrapping:
	changes = df_change['oppression_D']# np.array(df_change[(df_change['oppression']>2)*(df_change['oppression']<8)]['oppression_D'])
	for year in np.sort(df['year'].unique()):
		sel = df['year'] == year
		df.loc[sel,'oppression'] += np.random.choice(changes, size=sum(sel))
		df.loc[sel,'oppression'] = np.minimum(np.maximum(df.loc[sel,'oppression'],0), 10)
	df, df_lags = create_lags(df)
	return df, df_lags
	



	
def analyze(df, df_lags, prefix, thresholds=None):
	d = {}
	sim = prefix =='simulation'
	d['mean crackdown'], d['crackdown ratio'] = analyze_structure(df_lags, prefix, sim)
	#d['reg coef'] = regressions(df, prefix, sim)
	if thresholds is None:
		t0, t1  = mesh(df, 'oppression', sim)
	else:
		t0, t1 = thresholds
	d['thresholds'] = t0, t1 
	markov_matrix, tr, t = mkv.markov(df, [t0,t1], 'oppression')
	d['markov_matrix']  = markov_matrix
	x, y, rel_mid = cond_vol(df_lags)
	d['cond_vol']  = x
	d['cond_vol_ifphol']  = y
	d['cond_vol_rel_mid']  = rel_mid


	#replication of chris:
	#res = pt.execute('totdurny2~ oppression + oppression**2 + oppression**3', df, 
	#		ID = 'country', T = 'year'
	#		)

	#d['Chris'] = res.coef_params
	return d

	

def analyze_structure(df, prefix, sim):
	
	df_crackdown = df[(df['oppression']>=7.0)*(df['oppression_D']>1.0)][[f'oppression_L{i}' for i in range(16)]]
	mean_crackdown = df_crackdown.mean()
	crackdown_ratio = len(df_crackdown)/len(df)
	return mean_crackdown, crackdown_ratio


def generate_paths(N, T, s, maxval):
	# Initialize the paths array with zeros. The first dimension is time, the second is the path index.
	# Start all paths at M/2 to ensure they're within the interval [0, M] at the start.

	paths = np.zeros((T, N))
	paths[0] = np.random.rand(N)*maxval

	
	for t in range(1,T):

		innovations = np.random.normal(0, s, N)
		paths[t] = np.clip(paths[t-1] + innovations, 0, maxval)
		
	country =  np.array([[f'Country {i}']*T for i in range(N)]).T
	year =  np.array([[1960+i]*T for i in range(N)]).T
	data = {'country': country.flatten(), 'oppression':paths.flatten(), 'year':year.flatten(), 'totdurny2':gen_totdurney(N,T, paths)}
	df = data.copy()
	df = df.sort_values(['country', 'year'])
	df, df_lags  = shifts(df)
	return df, df_lags

def shifts(df):
	grouped = df.groupby('country')
	df['oppression_D'] = grouped['oppression'].diff()
	
	df_diff = df.copy().dropna()
	df_diff[f'oppression_L1'] = grouped['oppression'].shift(1)

	for i in range(16):
		df[f'oppression_L{i}'] = grouped['oppression'].shift(i)
	df = df.copy().dropna()
	return df_diff, df

def gen_totdurney(N,T, paths):
	P_CHANGE = 0.28
	a=0
	fpaths = P_CHANGE*(40*paths - 4*paths**2)
	fpaths = fpaths*np.mean(paths)/np.mean(fpaths)

	while np.sum(a)==0:
		a = np.cumsum(np.random.rand(N,T)<fpaths.T*P_CHANGE/T,1)
	b = np.zeros((N,T))
	k=0
	while np.sum(b==0)>0:
		for i in range(len(a)):
			b[i][a[i]==k]=np.sum(a[i]==k)
		k+=1
	return b.T.flatten()


def regressions(df, prefix, sim):
	measure = 'oppression'
	reslts = {}
	for test, heading in [[df['year']>=1991, 'Post cold war'], 
						  [df['year']<1991, 'Cold war']]:
		

		df_sub = df[test].copy()
		res = pt.execute('np.abs(oppression_D)~ oppression + oppression**2 ', df_sub, 
				ID = 'country', T = 'year'
				)
		reslts[heading] = res.coef_params


	return reslts


def create_lags(df):
	# Calculate the one-year change in the 'variable' column
	df = df.copy()
	grouped = df.groupby('country')
	#df['totdurny2_L1'] = df['totdurny2'].shift(1)
	df['oppression_D'] = grouped['oppression'].diff()

	df_diff = df.copy().dropna(subset=['oppression'])
	df_diff[f'oppression_L1'] = grouped['oppression'].shift(1)


	df['oppression_DD'] = grouped['oppression_D'].diff()
	df['oppression_DD_L1'] = grouped['oppression_DD'].shift(1)
	df['oppression_DD_L2'] = grouped['oppression_DD'].shift(2)
	for i in range(10):
		df[f'oppression_D_L{i}'] = grouped['oppression_D'].shift(i)
	for i in range(16):
		df[f'oppression_L{i}'] = grouped['oppression'].shift(i)

	# Drop the first row for each group as it will have a NaN value for the changes
	df = df.dropna(subset=['oppression'])
	
	# Print the resulting DataFrame
	return df_diff, df

def mesh(df, name, sim):
	dumpfile = 'results/thresholds.dmp'
	d = {}
	x = 0.01
	if not sim:
		try:
			with open(dumpfile, 'rb') as f:
				t1, t2  = pickle.load(f)
			return t1, t2
		except:
			pass

	mint = min(df[name][df[name]>0])
	maxt1 = 0.4
	mint2 = 0.6
	maxval = max(df[name])
	maxt = max(df[name][df[name]<maxval])/maxval
	for t1 in np.arange(mint,maxt1,x):
		for t2 in np.arange(mint2,maxt,x):
			markov_matrix, tr, t = mkv.markov(df, [t1,t2], name)
			obj = np.std(np.diag(markov_matrix), ddof=1)
			if not np.isnan(obj): 
				#print(f"{obj}: {np.diag(markov_matrix)}; {(t1,t2)}")
				d[obj] =[t1, t2]

	t1, t2 = d[min(d)]
	markov_matrix, tr, t = mkv.markov(df, [t1,t2], name)
	with open(dumpfile, 'wb') as f:
		pickle.dump((t1, t2), f)
	return t1, t2 

main()