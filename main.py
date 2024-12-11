import pandas as pd

from scipy.optimize import minimize
import numpy as np
from matplotlib import pyplot as plt
import paneltime as pt
import statsmodels.api as sm
import os
os.chdir(os.path.dirname(__file__))
import networkdrawing as nw
import markov
import montecarlo as mt


from scipy import stats




def main():
    df = handle_data()
    sim()
    analyze_structure(df)
    #res = OLS(df, X)

    pres = panel_reg(df)
    
    
    
    t0, t1  = mesh(df, 'ifhpol')
    markov_matrix, tr, t = markov.markov(df, [t0,t1], 'ifhpol')
    mean, std = mt.monte_carlo(df, 1000, [1/3,2/3], 'ifhpol')
    print(mean)
    print(std)
    print(pd.DataFrame(np.round(100*markov_matrix,1)))
    nw.draw_network(markov_matrix)

    markov.save_transitions(df)
    #res = OLS(df, X)

def analyze_structure(df):
    df[(df['ifhpol']<4)*(df['ifhpol_D']<-1)][[f'ifhpol_L{i}' for i in range(16)]].mean().plot(kind='bar')
    t_stat, p_value = stats.ttest_rel(df['ifhpol'], df['ifhpol_L9'])
    print("t-statistic:", t_stat)
    print("p-value:", p_value)

    a=0

def plotting(df, res, measure, coefs, deduct, interval):
    df['change_pred']  = res.predict()['in-sample predicted Y']

    fig, ax = plt.subplots()

    for test, heading in [[df['year']<1991, 'Pre cold war'], 
                 [df['year']>=1991, 'Post cold war']]:
        df_sub = df[test]
        df_sub = df_sub[[measure, 'change_pred']]
        df_sub = df_sub.sort_values(measure)
        df_sub.plot(measure, 'change_pred', ax = ax, label=heading)
    plt.show()
    x = np.linspace(*interval,100)
    plt.plot(x,coefs[0] + coefs[1]*x+ coefs[2]*np.abs(x-deduct)**2)
    plt.title('prediction of normalized_rank and normalized_rank_sq only')
    plt.show()
    a=0
    
def handle_data():
    df = pd.read_csv("datasimplified.csv")
    df = create_annual_rank(df)
    df = create_lags(df)

    df = pd.DataFrame(df)
    df['normalized_rank_sq'] = np.abs(df['normalized_rank']-0.5)**2
    df['demdummy'] = (df['normalized_rank']==1.0)*1.0
    df['semidemdummy'] = (df['normalized_rank']>0.5)*1.0
    return df
    
def panel_reg(df=None):
    pt.options.pqdkm.set([1,1,0,1,1])
    pt.options.fixed_random_time_eff.set(2)
    pt.options.fixed_random_group_eff.set(2)
    pt.options.tolerance.set(0.000001)
    
    if False:#abandoned using rank
        res = pt.execute('abs_one_year_change~normalized_rank+normalized_rank_sq+semidemdummy+demdummy', df, 
                ID = 'country', T = 'year'
                )
        print(res)
        a = res.ll.args.args_names
        plotting(df, res, 'normalized_rank', 
                [a['Intercept'] + a['demdummy'] +a['semidemdummy'], 
                a['normalized_rank'], a['normalized_rank_sq']], 
                0.5, [0,1])


    res = pt.execute('np.abs(ifhpol_D)~ (year<1991) +ifhpol + ifhpol**2 + (ifhpol>5) +(ifhpol==10)', df, 
               ID = 'country', T = 'year'
               )
    print(res)
    a = res.ll.args.args_names
    plotting(df, res,'ifhpol',
              [a['(ifhpol==10)'] +a['(ifhpol>5)'] + a['Intercept'] , 
                a['ifhpol'], a['ifhpol**2']], 
             0, [0,10])

    #replication of chris:
    res = pt.execute('totdurny2~ ifhpol + ifhpol**2 + ifhpol**3', df, 
               ID = 'country', T = 'year'
               )
    print(res)
    res_dd = pt.execute('np.abs(ifhpol_D)~(year<1991) +ifhpol + ifhpol**2 + ifhpol_DD_L1+ifhpol_DD_L2  + (ifhpol>5) +(ifhpol==10)', df, 
               ID = 'country', T = 'year'
               )
    print(res_dd)
    pt.options.pqdkm.set([4,0,0,0,1])
    res_arma = pt.execute('ifhpol_D', df, 
               ID = 'country', T = 'year'
               )
    print(res_arma)
    
    return res
    

def generate_paths(N, T, s, M):
    # Initialize the paths array with zeros. The first dimension is time, the second is the path index.
    # Start all paths at M/2 to ensure they're within the interval [0, M] at the start.

    paths = np.zeros((T, N))
    paths[0] = np.random.rand(N)*10

    
    for t in range(1,T):

        innovations = np.random.normal(0, s, N)
        paths[t] = np.clip(paths[t-1] + innovations, 0, M)
        
    country =  np.array([[f'Country {i}']*T for i in range(N)]).T
    year =  np.array([[1960+i]*T for i in range(N)]).T
    data = {'country': country.flatten(), 'ifhpol':paths.flatten(), 'year':year.flatten(), 'totdurny2':gen_totdurney(N,T, paths)}
    df = pd.DataFrame(data)
    df = df.sort_values(['country', 'year'])
    return df

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


def sim():
    df = generate_paths(190,40, 1, 10)


    grouped = df.groupby('country')
    df['ifhpol_D'] = grouped['ifhpol'].diff()
    for i in range(16):
        df[f'ifhpol_L{i}'] = grouped['ifhpol'].shift(i)
    df = df.dropna()


    analyze_structure(df)



    res = pt.execute('np.abs(ifhpol_D)~ ifhpol + ifhpol**2 ', df, 
               ID = 'country', T = 'year'
               )
    print(res)
    a = res.ll.args.args_names
    plotting(df, res,'ifhpol',
            [a['Intercept'] , 
            a['ifhpol'], a['ifhpol**2']], 
            0, [0,10])
        #replication of chris:
    res = pt.execute('totdurny2~ ifhpol + ifhpol**2 + ifhpol**3', df, 
               ID = 'country', T = 'year'
               )
    print(res)

    t0, t1  = mesh(df, 'ifhpol')
    markov_matrix, tr, t = markov.markov(df, [t0,t1], 'ifhpol')
    a=0
    
def OLS(df, X):
    X = sm.add_constant(df[['normalized_rank', 'normalized_rank_sq', 'semidemdummy', 'demdummy']])
    model = sm.OLS(df['abs_one_year_change'], X)
    results = model.fit()
    df['change_pred'] = results.predict()
    print(results.summary()) 
    return results

def create_annual_rank(df):
    df = pd.DataFrame(df)
    df = df.dropna()
    df = df.sort_values(['year', 'ifhpol'])
    # Group the data by 'year'
    grouped = df.groupby('year')
    # Add a new column 'rank' with the ranking for each identifier within each year
    df['rank'] = grouped['ifhpol'].rank(ascending=True)

    # Calculate the maximum rank for each year
    max_rank = grouped['rank'].transform('max')

    # Calculate the normalized rank for each identifier within each year
    df['normalized_rank'] = df['rank'] / max_rank
    df['max_rank'] = max_rank
    

    df = df.sort_values(['country', 'year'])

    return df

def create_lags(df):
    # Calculate the one-year change in the 'variable' column
    grouped = df.groupby('country')
    df['rank_L1'] = df['rank'].shift(1)
    df['ifhpol_L1'] = df['ifhpol'].shift(1)
    df['totdurny2_L1'] = df['totdurny2'].shift(1)
    df['abs_one_year_change'] = np.abs(grouped['normalized_rank'].diff())
    df['ifhpol_D'] = grouped['ifhpol'].diff()
    df['ifhpol_DD'] = grouped['ifhpol_D'].diff()
    df['ifhpol_DD_L1'] = grouped['ifhpol_DD'].shift(1)
    df['ifhpol_DD_L2'] = grouped['ifhpol_DD'].shift(2)
    for i in range(10):
        df[f'ifhpol_D_L{i}'] = grouped['ifhpol_D'].shift(i)
    for i in range(16):
        df[f'ifhpol_L{i}'] = grouped['ifhpol'].shift(i)
    df['normalized_rank_L1'] = df['rank_L1'] / df['max_rank']

    # Drop the first row for each group as it will have a NaN value for the changes
    df = df.dropna()
    # Print the resulting DataFrame
    return df

def mesh(df, name):
    d = {}
    x = 0.01
    mint = np.min(df[name])
    maxt1 =0.3
    mint2 = 0.9
    for t1 in np.arange(mint,maxt1,x):
        for t2 in np.arange(mint2,1.0,x):
            markov_matrix, tr, t = markov.markov(df, [t1,t2], name)
            obj = np.std(np.diag(markov_matrix), ddof=1)
            if not np.isnan(obj): 
                print(f"{obj}: {np.diag(markov_matrix)}; {(t1,t2)}")
                d[obj] =[t1, t2]

    t1, t2 = d[min(d)]
    markov_matrix, tr, t = markov.markov(df, [t1,t2], name)
    return t1, t2 

main()