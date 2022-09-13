import numpy as np

def simulate_measurement(A, t_out, t_in, E_low, sigma_low , E_high, sigma_high,	 n_traces=1, t_tot=2000, acq_rate=1):
	''' Creates measurement traces using the input parameters and returns the maximum of each trace
	  Parameters
	  ----------
	   
	  Dot parameters: 
		  A: float from 0 to 1
			  Proportion of spin down traces
			  
		  t_out: float
			Time constant for a spin up electron to leave the dot (in  mu s)

		  t_in: float
			Time constant for a spin down electron to return to the dot (in	 mu s)
	  
	  Sensor parameters: 
		  E_low: float
			  Epected value of the background
			  
		  sigma_low: float
			  Background noise
			  
		  E_high: float
			  Exoected value at the blip
			  
		  sigma_high: float
			  noise during the blip		  
	  
	  User parameters
		  n_traces: integer
			  Number of simulated traces	 

		  t_tot: float
			Time during which data is acquire (in  mu s)
			
		  acq_rate:
			Acquisition rate in MHz

  '''
  # Testing inputs
  
	if A < 0 or A> 1: raise ValueError(
			"A must be between 0 and 1")
	if type(n_traces) != int:
		raise TypeError("The number of traces should be an *integer*.")
	# and the right values (positive or null)
	if t_out<0 or t_in<0:
		raise ValueError("t_in and t_out should be positive numbers.")
	if t_tot*acq_rate< 3:
		raise ValueError("The readout trace is too short. Try to increase the acquisition rate or the readout time (t_tot).")
	if t_tot< 0 or acq_rate< 0: 
		raise ValueError("The readout time and acquisition rate must be positive numbers.")


	spin_traces = []
	spin_number=np.empty(0)
	n_points_time=int(t_tot*acq_rate)
	
	for i in range(n_traces):
		choose_spin = np.random.binomial(1,A,1) #create 1 with proba A and 0 with proba (1-A)
		spin_number = np.append(spin_number, choose_spin)
		
		if choose_spin == 1 : # spin down because A = <P_down>
		# create a spin down trace i.e. constant E_low with noise sigma_low
			line_spin_down = E_low + np.random.normal(0,sigma_low ,n_points_time)
			spin_traces .append(line_spin_down )
			
			
		if choose_spin == 0 : #spin up
		# create a spin up trace i.e. trace with a blip of E_high 
			line_spin_up =spin_up_trace_sim(t_out, t_in, E_low, sigma_low, E_high, sigma_high, t_tot, acq_rate)		  
			spin_traces .append( line_spin_up )
			
	spin_traces =np.array(spin_traces )
	time=np.arange(0,t_tot, 1/acq_rate)		# sweep time over trace
	return(time,spin_traces,spin_number)	 # spin_traces : array of 1D traces, spin_number: array with their corresponding spin

def spin_up_trace_sim(t_out, t_in, E_low, sigma_low, E_high, sigma_high, t_tot, acq_rate):
	''' Creates a spin up measurement trace
	  Parameters
	  ----------
	   
	  Dot parameters: 
		  t_out: float
			Time constant for a spin up electron to leave the dot (in mu s)

		  t_out: float
			Time constant for a spin down electron to return to the dot (in mu s)
	  
	  Sensor parameters: 
		  E_low: float
			  Epected value of the background
			  
		  sigma_low: float
			  Background noise
			  
		  E_high: float
			  Exoected value at the blip
			  
		  sigma_high: float
			  noise during the blip		  
	  
	  User parameters
		   t_tot: float
			Time during which data is acquire (in mu s)
			
		  acq_rate:
			Acquisition rate in MHz	  
			
	'''	  
	time=np.arange(0,t_tot, 1/acq_rate) # sweep time over trace
	spin_trace_up=np.zeros(len(time)) # spin-up trace
	
	t_leave=np.random.exponential(t_out) #random t_out time
	t_enter=np.random.exponential(t_in) #random t_in time
	
	i_leave = np.where(time<=t_leave)[0][-1] 
	i_enter = np.where(time<=t_enter+t_leave)[0][-1]
	
	
	spin_trace_up[0:i_leave] = E_low + np.random.normal(0,sigma_low,len(time[0:i_leave])) # first part of trace (background)

	if i_enter == i_leave : #means blip time = 0
		spin_trace_up[i_leave:] =  E_low + np.random.normal(0,sigma_low,len(time[i_leave:])) 

	else : 
		spin_trace_up[i_leave:i_enter] = E_high + np.random.normal(0,sigma_low,len(time[i_leave:i_enter])) # blip 
		spin_trace_up[i_enter:] = E_low + np.random.normal(0,sigma_low,len(time[i_enter:])) # last part of the trace (background)

	return spin_trace_up

def moving_avg_filter(spin_traces,time,	 factor, acq_rate):
	''' Applies moving average filter to the generated traces
	  Parameters
	  ----------
		  spin_traces: 2D array
			All data traces	 
			  
		  time: 1D array
			Time sweep during a measurement
			
		  factor: integer
			length of the moving average. When n=20 a moving average is perform on array a every 20 points
		  acq_rate:
			 Acquisition rate in MHz   
	   
	   Returns
	   -------
	   2d array.
		moving average of spin traces
	   
  '''
  
	
	if type(factor) != int:
		raise TypeError("Factor should be an *integer*.")
	# and the right values (positive or null)
	if factor<0:
		raise ValueError("Factor needs to be a positive number.")

	spin_traces_filtered=[]
	time_filtered=[]
	for a in range(len(spin_traces)):
		ret = np.cumsum(spin_traces[a], dtype=float)
		ret[factor:] = ret[factor:] - ret[:-factor] 
		spin_traces_filtered.append(ret[factor - 1:] / factor)
	spin_traces_filtered=np.array(spin_traces_filtered)
	
	
	
	ret = np.cumsum(time, dtype=float)
	ret[factor:] = ret[factor:] - ret[:-factor] 
	time_filtered=(ret[factor - 1:] / factor)
	time_filtered=time_filtered+factor/acq_rate/2 # Initial time is larger due to the filter
	
	return time_filtered, spin_traces_filtered
 
  
 





def get_inter_bins(bins):
    ''' Creates an array with the position of each bin
	  Parameters
	  ----------
		  bins: 1D array
			bins given by the function np.histogram
	'''
	res = []
	for i in range (len(bins)-1):
		res.append( 1/2*(bins[i]+bins[i+1]))
	return(res)

def max_all_traces(time, spin_traces, t_min, t_max):
	''' Calculates the maximum of each trace from t_min to t_max
	  Parameters
	  ----------
		  time: 1D array
			Time sweep during a measurement
		  spin_traces: 2D array
			All data traces		  
		  t_min: float
			start of the measurement trace
		  t_max: float
			end of the measurement trace

	'''
	minim= np.where(time>= t_min)[0][0]
	maxim= np.where(time>= t_max)[0][0]
	max_spin_traces=np.amax(spin_traces[:,minim:maxim],axis=1)
	return max_spin_traces
	

def max_simulated_histogram(time, spin_traces, spin_number, t_min=0, t_max=None, num_bins=100, plot_histo=True):
	from scipy.interpolate import interp1d
	import matplotlib.pyplot as plt
	''' Creates a histogram with the maximum of the spin traces
	  Parameters
	  ----------
		  time: 1D array
			Time sweep during a measurement
		  spin_traces: 2D array
			All data traces		  
		  spin_number: 1D array
			State of each data trace. 0 is spin up and 1 is spin down 
		  t_min: float
			start of the measurement trace
		  t_max: float
			end of the measurement trace
		  bins: integer
			 number of bins in the histogram
		  Plot: bool
			 True if the histogram needs to be visualised
  '''
	if t_max==None:
		t_max=time.max()
		
	n_traces=len(spin_traces) 
	max_spin_traces=max_all_traces(time, spin_traces, t_min, t_max) # maximum of each trace
	range_min=max_spin_traces.min() 
	range_max=max_spin_traces.max()
	print(range_min)  

	# total histogram
	hist_max_ss, bin_edges_max_ss = np.histogram(max_spin_traces,bins=num_bins,range = (range_min,range_max)) 
	bin_axis_ss= get_inter_bins(bin_edges_max_ss)
	normalisation=( bin_edges_max_ss[1]- bin_edges_max_ss[0])*n_traces # for histogram area =1	

      
	if plot_histo:
        
		bin_width=bin_edges_max_ss[1]-bin_edges_max_ss[0]
       
        # Spin up histogram
		spin_up_maximums=max_spin_traces[np.where(spin_number==0)[0]] # repeated
		hist_max_up, bin_edges_max_up = np.histogram(spin_up_maximums,bins=num_bins,range = (range_min,range_max)) 


        # Spin down histogram
		spin_down_maximums=max_spin_traces[np.where(spin_number==1)[0]] # repeated
		hist_max_down, bin_edges_max_down = np.histogram(spin_down_maximums,bins=num_bins,range = (range_min,range_max)) 

		# Plot the histogram 
		plt.figure()
		plt.bar(bin_axis_ss, hist_max_up/normalisation, width=bin_width,label="spin_up",alpha=0.5)
		plt.bar(bin_axis_ss, hist_max_down/normalisation, width=bin_width,label="spin_down",alpha=0.5)
		plt.xlabel("Signal maxima")
		plt.ylabel("Probability")
		plt.title('Histogram of the readout maxima')
		#print('Parameters for the simulation :')
		#print('T_int = %f us, T_in = %f us, T_out = %f us' %(t_int,t_in,t_out))
		#print('A = %f, std = %f'%(A,std_down))
		#plt.title('Simulated histogram as lines')
		f_hist_inter = interp1d(bin_axis_ss, hist_max_ss/normalisation, kind='cubic') # smoothing the histogram
		plt.plot(bin_axis_ss,f_hist_inter(get_inter_bins(bin_edges_max_ss))) 
		plt.legend()
		plt.show()
 
		
	return(bin_axis_ss,hist_max_ss/normalisation)
		







# we can also have it as a class and it gives as the plot as one of its functions
def fidelity(time, spin_traces, spin_number, t_min=0, t_max=None, epsilon=0.005, plot_visibility=True):
	''' Creates a histogram with the maximum of the spin traces
	  Parameters
	  ----------
		  time: 1D array
			Time sweep during a measurement
		  spin_traces: 2D array
			All data traces		  
		  spin_number: 1D array
			State of each data trace. 0 is spin up and 1 is spin down 
		  t_min: float
			start of the measurement trace
		  t_max: float
			end of the measurement trace
		  epsilon: float
			delimitant of the threshol. Possible thesholds go from range_min+epsilon to range_max-epsilon 
		  plot_visibility

  '''
	import matplotlib.pyplot as plt
  # repeated from here

	if t_max==None:
		t_max=time.max()
		
	n_traces=len(spin_traces) 
	max_spin_traces=max_all_traces(time, spin_traces, t_min, t_max) # maximum of each trace
	range_min=max_spin_traces.min() 
	range_max=max_spin_traces.max()
	print(range_min)  
 
 # to here
	signal_threshold = np.linspace(range_min + epsilon, range_max - epsilon,50) # array with all possible thresholds

	fidelity_up = []
	fidelity_down = []

	spin_up_maximums=max_spin_traces[np.where(spin_number==0)[0]]
	spin_down_maximums=max_spin_traces[np.where(spin_number==1)[0]]
	
	for signal in signal_threshold: # calculate the fidelity for each threshold
				integ_up=len(np.where(spin_up_maximums<=signal)[0]) # spin up traces that don't surpass thershold
				integ_down=len(np.where(spin_down_maximums>=signal)[0]) # spin down traces that surpass the thershold

				fidelity_up.append(1-integ_up/n_traces)
				fidelity_down.append(1-integ_down/n_traces)
                
	# Obtain optimal fidelity and visibility
	visibility=np.array(fidelity_up)+np.array(fidelity_down)-1
	integer=np.where(visibility==visibility.max())[0][0]
	fidelity_up_def=np.array(fidelity_up)[integer]
	fidelity_down_def=np.array(fidelity_down)[integer]
	fidelity_tot_def=(fidelity_up_def+fidelity_down_def)/2                
                
	if plot_visibility:
		plt.plot(signal_threshold,fidelity_down, color=(120/255,120/255,120/255), label=r' Spin down fidelity')#Fidelity$_{\downarrow}$')
		plt.plot(signal_threshold,fidelity_up, color= (255/255,66/255,66/255), label=r'Spin up fidelity')#Fidelity$_{\uparrow}$' )
		plt.plot(signal_threshold,np.array(fidelity_up)+np.array(fidelity_down)-1,'--',label='Visbility',c='k')

		plt.legend(loc=1, fontsize=8)
		plt.xlabel("Threshold voltage (V) ", {'color': 'k', 'fontsize': 14})
		plt.ylabel("Fidelity", {'color': 'k', 'fontsize': 14})
		plt.xticks( color='k', size=12)
		plt.yticks( color='k', size=12)
		#plt.savefig("fidelity.png",dpi=600, bbox_inches = "tight", rasterized='True')
		print("max visibility =", visibility.max(),'at', signal_threshold[np.argmax(np.array(fidelity_up)+np.array(fidelity_down)-1)])
		print("max fidelity =", (fidelity_up_def+fidelity_down_def)/2)
	
	
	return fidelity_up_def, fidelity_down_def, fidelity_tot_def





# max_histogram and plot
# fidelity and plot
# __init__.py
# readme.md file
# setup.py file
# license file
