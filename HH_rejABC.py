import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from neuron import h

import efel

class ABC:
	def __init__(self, N: int = 1000, q: float = 0.05, offset: int = 0,
			np_uniform: bool = False, use_efel: bool = True):

		self.N = N
		self.q = q
		self.offset = offset
		self.gnabar = 120 #mS
		self.gkbar = 36 #mS

		if np_uniform:
			self.gnas = ...
			self.gks = ...
		else:
			self.gnas = np.linspace(self.offset + self.gnabar - self.gnabar/10,
				self.offset + self.gnabar + self.gnabar/10, self.N)
			self.gks = np.linspace(self.offset + self.gkbar - self.gkbar/10,
				self.offset + self.gkbar + self.gkbar/10, self.N)

		if use_efel is not True:
			raise ValueError(f'Only supports summary statistics from efel.')


		self.traces = []

	def run(self, kind: str) -> np.ndarray:
		if "k" in kind:
			self.gnas = np.ones(len(self.gnas)) * self.gnabar
			val = self.gks
		elif "na" in kind:
			self.gks = np.ones(len(self.gks)) * self.gkbar
			val = self.gnas
		else:
			raise ValueError(f"Supports 'k' and 'na'.")

		# true trace
		self._tracer(gna=self.gnabar, gk=self.gkbar)
		self.nums = int(np.max([len(self.gnas), len(self.gks)]))
		for i in range(self.nums):
			self._tracer(gna=self.gnas[i], gk=self.gks[i])

		stats = self._stats(traces=self.traces)
		posterior = self._rejecter(stats=stats, val=val)

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
		stim.dur = 22 # ms

		vvec = h.Vector().record(soma(l)._ref_v)
		tvec = h.Vector().record(h._ref_t)

		h.finitialize(-65) # mV
		h.continuerun(22) # ms

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
		traces_results = efel.getFeatureValues(traces, ['AP1_amp', 'AHP1_depth_from_peak', 'AP1_begin_voltage'], return_list=True)
		# list of dicts to dict of lists
		traces_results = {k: [dic[k] for dic in traces_results] for k in traces_results[0]}
		self.efelFeatures = list(traces_results.keys())

		return traces_results
	
	@staticmethod
	def _calc_2D_dist(x, y) -> np.ndarray:
		return np.linalg.norm(x-y, ord=2, axis = 1)

	def _rejecter(self, stats, val):
		vals = []
		distance = np.zeros((len(self.efelFeatures), self.nums))
	
		for idx, name in enumerate(self.efelFeatures):
	
			distance[idx] = self._calc_2D_dist(stats[name][0],
		                         stats[name][1:])
		
			epsilon = np.quantile(distance[idx], self.q)
			
			vals.append(val[distance[idx] <= epsilon])
		
		return vals

	


if __name__ == '__main__':
	o = ABC(N = 10000, q = 0.01)
	gnas = o.run('na')
	sns.kdeplot(gnas)
	plt.show()