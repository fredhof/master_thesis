import os
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import time

from neuron import h

import efel


rng = np.random.default_rng(1)

class ABC:
	def __init__(self, N: int = 1000, q: float = 0.05, offset: int = 0,
			np_uniform: bool = True, efel_stats: list = None, use_efel: bool = True):

		self.N = N
		self.q = q
		self.offset = offset
		self.gnabar = 120 #mS
		self.gkbar = 36 #mS

		if np_uniform:
			self.gnas = rng.uniform(self.offset + self.gnabar - self.gnabar/10,
				self.offset + self.gnabar + self.gnabar/10, size=self.N)
			self.gks = rng.uniform(self.offset + self.gkbar - self.gkbar/10,
				self.offset + self.gkbar + self.gkbar/10, size=self.N)
		else:
			self.gnas = np.linspace(self.offset + self.gnabar - self.gnabar/10,
				self.offset + self.gnabar + self.gnabar/10, self.N)
			self.gks = np.linspace(self.offset + self.gkbar - self.gkbar/10,
				self.offset + self.gkbar + self.gkbar/10, self.N)

		self.ths = np.zeros((self.N,2))

		if use_efel is not True:
			raise ValueError(f'Only supports summary statistics from efel.')
		self.efel_stats = efel_stats

		self.traces = []

	def run(self) -> np.ndarray:
		
		# true trace
		self._tracer(gna=self.gnabar, gk=self.gkbar)
		
		#self.nums = int(np.max([len(self.gnas), len(self.gks)]))
		for i in range(self.N):
			pick = rng.choice(np.arange(self.N))
			self.ths[i] = [self.gnas[pick], self.gks[pick]]
			self._tracer(gna=self.gnas[pick], gk=self.gks[pick])
		
		
		stats = self._stats(traces=self.traces)
		
		
		posterior = self._rejecter(stats=stats)

		return posterior
	
	@staticmethod
	def _nrn_HH(gnabar: int = 120, gkbar: int = 36):
		h.load_file('stdrun.hoc') 

		soma = h.Section('soma')
		# soma is a rod
		soma.pt3dadd(60,0,0,30)
		soma.pt3dadd(90,0,0,30)
		soma.insert('hh')
		soma().hh.gnabar = gnabar*1e-3 # mS -> S
		soma().hh.gkbar = gkbar*1e-3


		l = 0.5 # center of soma
		stim = h.IClamp(soma(l))
		#stim.delay = 50
		stim.amp = 0.3 # mA/cm2
		stim.dur = 40 # ms

		vvec = h.Vector().record(soma(l)._ref_v)
		tvec = h.Vector().record(h._ref_t)

		h.finitialize(-65) # mV
		h.continuerun(40) # ms

		return tvec, vvec

	def _tracer(self, gna, gk):
		t, V = self._nrn_HH(gna, gk)
		
		trace = {}
		trace['T'] = t
		trace['V'] = V
		trace['stim_start'] = [t[0]]
		trace['stim_end'] = [t[-1]]

		self.traces.append(trace)

	def _stats(self, traces: list, transform: bool = True, return_list: bool = True) -> dict:
		traces_results = efel.getFeatureValues(traces, self.efel_stats, return_list=True)
		

		# delete traces where we have an empty list (i.e. bad trace)
		for num, i in enumerate(traces_results):
			for j in range(len(self.efel_stats)):
				if not np.array(list(i.values())[j]).size >0:
					del traces_results[num]
					self.ths = np.delete(self.ths, num, axis=0)
		
		# list of dicts to dict of lists
		traces_results = {k: [dic[k] for dic in traces_results] for k in traces_results[0]}
		print(np.array(list(traces_results.values())).shape)
		l = np.array(list(traces_results.values()))
		print(l[0])
		return np.hstack(list(traces_results.values()))
	
	@staticmethod
	def _calc_2D_dist(x, y) -> np.ndarray:
		return np.linalg.norm(x-y, ord=2, axis=1)

	def _rejecter(self, stats):
		th = []
		stat = [stats[0]]

		
		#ths = rng.choice(self.ths.shape[0]-1)

		distance = self._calc_2D_dist(stats[0], stats[1:])
		
		epsilon = np.quantile(distance, self.q)
		th.append(self.ths[distance <= epsilon])
		stat.append(stats[1:][distance <= epsilon])	
		#th.append(self.ths[distance <= epsilon])
		"""for idx, name in enumerate(self.efelFeatures):
	
			distance[idx] = 
		
			epsilon = np.quantile(distance[idx], self.q)
			
			vals.append(val[distance[idx] <= epsilon])
		"""
		return th, stat

	


if __name__ == '__main__':
	start = time.time()
	efel_stats = ['AP1_amp', 'AHP1_depth_from_peak', 'AP1_begin_voltage']
	es2 = ['AP2_amp', 'AHP2_depth_from_peak', 'AP2_begin_voltage']

	o = ABC(N = 4000, q = 0.01, np_uniform=True, efel_stats=efel_stats+es2)
	gnas, stats = o.run()

	# saves output and adds number if already exists
	th_name = 'ths_file'
	stats_name = 'stats_file'
	if os.path.exists(th_name+".npy"):
		n = 0
		pattern = r'[0-9]'
		while os.path.exists(f'{th_name}.npy'):
			# removes number from filename so it doesnt become abc0123456... .npy
			th_name = re.sub(pattern,"",th_name)
			stats_name = re.sub(pattern,"",stats_name)

			th_name = th_name + str(n)
			stats_name = stats_name + str(n)
			n += 1

		np.save(th_name, gnas)
		np.save(stats_name, stats)
	else:
		np.save(th_name, gnas)
		np.save(stats_name, stats)

	print(time.time()-start)
	print('Done.')

"""

def _linreg_adj(self, theta, s_obs, s_sim):
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

"""

# sammenligne med: https://github.com/mackelab/sbi
# https://github.com/eth-cscs/abcpy
# https://github.com/elfi-dev/elfi UIO


# des 9.
# fikse flere parametre, mulig master hvor mange parametre kan vi se på?
# linreg, se kode over
# se bok, Handbook of ABC, se mail


# des 16.
# fikse sånn at man trekker tilfeldig (samtidig) fra gna og gk -> self.ths
	# gjort
# confidence interval
# legg til flere statistikker og test med ridge