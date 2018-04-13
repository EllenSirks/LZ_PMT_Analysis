import matplotlib.pyplot as plt
import numpy as np
import csv

class temp_plot:

	def __init__(self, file):

		self.date = []
		self.time = []
		self.temp_A = []
		self.temp_B = []

		i = 0
		with open(file) as f:
			data = csv.reader(f, delimiter=',')
			for row in data:
				if i != 0 and float(row[2]) < 300.:
					self.date.append(row[0])
					self.time.append(row[1])
					self.temp_A.append(float(row[2]))
					self.temp_B.append(float(row[5]))
				i += 1

		self.dates = list(sorted(set(self.date)))


	def hms_to_seconds(self, t):
	    h, m, s = [int(i) for i in t.split(':')]
	    return h + m/60. + s/3600.

	def convert(self):
		start = self.hms_to_seconds(self.time[0])
		self.time = [self.hms_to_seconds(self.time[i]) - start for i in range(len(self.time))]

		for i in range(len(self.time)):
			for j in range(len(self.dates)):
				if self.date[i] == self.dates[j]:
					self.time[i] = self.time[i] + 24*j

		index = 0
		for k in range(len(self.time)):
			if self.time[k] <= 40:
				index = k

		self.time = self.time[:index+2]
		self.temp_A = self.temp_A[:index+2]
		self.temp_B = self.temp_B[:index+2]

	def temp_plot(self):
		self.convert()
		fig = plt.figure()
		plt.plot(self.time, self.temp_A, '-r', label='Temp A')
		plt.plot(self.time, self.temp_B, '-b', label='Temp B')
		plt.axhline(min(self.temp_B), linestyle='--', color='k')
		plt.axhline(max(self.temp_B), linestyle='--', color='k')
		plt.legend(loc='best')
		plt.xlabel('time (hr after start)', fontsize=15)
		plt.ylabel('T (K)', fontsize=15)
		# plt.title('23-11-2017')
		plt.xlim([self.time[0],self.time[-1]])
		plt.show()
		plt.close(fig)
		fig.savefig('/home/ellen/Desktop/MPhys_Project/temp_plots/06122017Cooldown+Warmup_subplot.png')

		diff = [self.temp_B[i] - self.temp_A[i] for i in range(len(self.temp_A))]
		fig = plt.figure()
		plt.plot(self.time, diff, '-', label= r'$T_{B}$ - $T_{A}$')
		plt.axhline(-5., linestyle='--', color='r', label=r'$\pm5K$')
		plt.axhline(0., linestyle='--', color='k')
		plt.axhline(5., linestyle='--', color='r')
		plt.legend(loc='best')
		plt.xlabel('time (hr after start)', fontsize=15)
		plt.ylabel(r'$\Delta$T (K)', fontsize=15)
		# plt.title('23-11-2017')
		plt.xlim([self.time[0],self.time[-1]])
		plt.show()
		plt.close(fig)
		fig.savefig('/home/ellen/Desktop/MPhys_Project/temp_plots/06122017Cooldown+Warmup_diff_subplot.png')

test = temp_plot('/home/ellen/Desktop/MPhys_Project/temp_files/06122017Cooldown+Warmup.csv')
diff = test.temp_plot()

