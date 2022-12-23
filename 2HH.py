from neuron import h
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from HH import neu
import efel
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, Lasso, Ridge
import copy

from scipy.stats import gaussian_kde
"""
def main():
	

	#samme tid neste uke
	# får ut Aksjonspotentiale DONE
	# lage prior gkabar, gnabar, lage fordelinger av gnabar og gkabar, UNIFORM
	# sammenligne trace, beskrive spikes, FWHM, bunnen, toppen, BRUK EFEL AP1 AMP, AHP1 Depth from peak, AP1 begin voltage, 
	# REJ ABC, trekke fra prior, SE NIC

	# fordel sterrat = neuron, units
	#nrn.rtd.com

	# init
	rng = np.random.default_rng(1) # rng seeded with seed 8

	def gen_uniform(element: str, offset: int = 10, size: int = 1000) -> np.ndarray:
		if element in "gkbar, gbark":
			element = 36
		elif element in "gnabar, gbarna":
			element = 120
		else:
			raise ValueError("'element' must be 'k' or 'Na'.")
		
		return rng.uniform(element - offset*rng.random(), element + offset*rng.random(), size = size)


	def distance_function(x,y):
	    return np.linalg.norm(x-y, ord=2, axis = 1)



	g_k_bar_true = 36 # mS/cm2
	g_na_bar_true = 120 # ms/cm2


	g_k = gen_uniform("k", size=1000)
	g_na = gen_uniform("na",size=1000)

	#g_k = np.linspace(31,41,1000)
	#g_na = np.linspace(116,136, 1000)

	true_t, true_V = neu(g_na_bar_true, g_k_bar_true)
	plt.plot(true_t, true_V)
	plt.show()

	traces = []
	t = np.zeros((true_t.shape[0], true_V.shape[0]))
	V = np.copy(t)

	print(true_t.shape, true_V.shape)

	true_trace = {}
	true_trace['T'] = true_t
	true_trace['V'] = true_V
	true_trace['stim_start'] = [true_t[0]]
	true_trace['stim_end'] = [true_t[-1]]
	traces.append(true_trace)


	for i in range(1000):
		t, V = neu(g_na[i], g_k_bar_true)
		trace = {}
		trace['T'] = t
		trace['V'] = V
		trace['stim_start'] = [t[0]]
		trace['stim_end'] = [t[-1]]
		traces.append(trace)	

	# very slow
	traces_results = efel.getFeatureValues(traces, ['AP1_amp', 'AHP1_depth_from_peak', 'AP1_begin_voltage'], return_list=True)
	# list of dicts to dict of lists
	traces_results = {k: [dic[k] for dic in traces_results] for k in traces_results[0]}

	efelFeatures = list(traces_results.keys())
	

	#Distance function - Absolute difference 


	def rej():
		samples = []
		gnas = []
	
		q = 0.01
	
		distance = np.zeros((len(efelFeatures), g_na.shape[0]))
	
		for idx, name in enumerate(efelFeatures):
	
			distance[idx] = distance_function(traces_results[name][0],
		                         traces_results[name][1:])
		
			epsilon = np.quantile(distance[idx], q)
			
			gnas.append(g_na[distance[idx] <= epsilon])
		
		return gnas


	gnas = rej()
	print(gnas)
	gnas = np.array(gnas)
	gnas = np.mean(gnas, axis=0)
	print(gnas)
	sns.kdeplot(gnas)
	plt.show()
	sns.histplot(g_na)
	plt.show()
#main()
"""

# squeeze transforms shape (1,x,y) -> (x,y)
ths = np.load("[N4000_6]ths_file.npy").squeeze()
s_obs, s_sim = np.load("[N4000_6]stats_file.npy", allow_pickle=True).squeeze()
 

def _linreg_adj(theta, s_obs, s_sim):
	X = copy.deepcopy(s_sim)
	X_obs = copy.deepcopy(s_obs).reshape(1, -1)
	y = np.log(copy.deepcopy(theta))

	if X.ndim == 1:
		X = X.reshape(-1, 1)

	scaler = StandardScaler()
	X = scaler.fit_transform(X)
	X_obs = scaler.transform(X_obs)

	reg_model = LinearRegression(fit_intercept=True)
	reg_model.fit(X, y)

	theta_adj = np.exp(reg_model.predict(X_obs) + y - reg_model.predict(X))
	#data = dict(zip(param_names, np.stack(theta_adj, axis=-1)))
	# return pd.DataFrame.from_dict(data)
	return theta_adj

gna, gk = _linreg_adj(ths, s_obs, s_sim).T


def rmse(y_true, y_pred):
    '''
    Compute Root Mean Square Percentage Error between two arrays.
    '''
    loss = np.sqrt(np.mean(np.square(((y_true - y_pred) / y_true)), axis=0))

    return loss

sns.kdeplot(gna)
plt.title(f"RMSE:{rmse(120, gna)*100:.3f}%, mean={gna.mean():.3f}, median={np.median(gna):.3f}, N*q={len(ths)}")
plt.show()
sns.kdeplot(gk)
plt.title(f"RMSE:{rmse(36, gk)*100:.3f}%, mean={gk.mean():.3f}, median={np.median(gk):.3f}, N*q={len(ths)}")
plt.show()





def hpd_grid(sample, alpha=0.05, roundto=2):
    """Calculate highest posterior density (HPD) of array for given alpha. 
    The HPD is the minimum width Bayesian credible interval (BCI). 
    The function works for multimodal distributions, returning more than one mode
    Parameters
    ----------
    
    sample : Numpy array or python list
        An array containing MCMC samples
    alpha : float
        Desired probability of type I error (defaults to 0.05)
    roundto: integer
        Number of digits after the decimal point for the results
    Returns
    ----------
    hpd: array with the lower 
          
    """
    sample = np.asarray(sample)
    sample = sample[~np.isnan(sample)]
    # get upper and lower bounds
    l = np.min(sample)
    u = np.max(sample)
    density = gaussian_kde(sample)
    x = np.linspace(l, u, 2000)
    y = density.evaluate(x)
    #y = density.evaluate(x, l, u) waitting for PR to be accepted
    xy_zipped = zip(x, y/np.sum(y))
    xy = sorted(xy_zipped, key=lambda x: x[1], reverse=True)
    xy_cum_sum = 0
    hdv = []
    for val in xy:
        xy_cum_sum += val[1]
        hdv.append(val[0])
        if xy_cum_sum >= (1-alpha):
            break
    hdv.sort()
    diff = (u-l)/20  # differences of 5%
    hpd = []
    hpd.append(round(min(hdv), roundto))
    for i in range(1, len(hdv)):
        if hdv[i]-hdv[i-1] >= diff:
            hpd.append(round(hdv[i-1], roundto))
            hpd.append(round(hdv[i], roundto))
    hpd.append(round(max(hdv), roundto))
    ite = iter(hpd)
    hpd = list(zip(ite, ite))
    modes = []
    for value in hpd:
         x_hpd = x[(x > value[0]) & (x < value[1])]
         y_hpd = y[(x > value[0]) & (x < value[1])]
         modes.append(round(x_hpd[np.argmax(y_hpd)], roundto))
    return hpd, x, y, modes

print(hpd_grid(gna))