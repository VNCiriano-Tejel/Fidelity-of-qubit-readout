import numpy as np

def simulate_measurement(A, t_out, t_in, E_low, sigma_low , E_high, sigma_high,	 n_traces=1, t_max=2000, acq_rate=1):
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

		  t_max: float
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
	if t_max*acq_rate< 3:
		raise ValueError("The readout trace is too short. Try to increase the acquisition rate or the readout time (t_max).")
	if t_max< 0 or acq_rate< 0: 
		raise ValueError("The readout time and acquisition rate must be positive numbers.")


	spin_traces = []
	spin_number=np.empty(0)
	n_points_time=int(t_max*acq_rate)
	
	for i in range(n_traces):
		choose_spin = np.random.binomial(1,A,1) #create 1 with proba A and 0 with proba (1-A)
		spin_number = np.append(spin_number, choose_spin)
		
		if choose_spin == 1 : # spin down because A = <P_down>
		# create a spin down trace i.e. constant E_low with noise sigma_low
			line_spin_down = E_low + np.random.normal(0,sigma_low ,n_points_time)
			spin_traces .append(line_spin_down )
			
			
		if choose_spin == 0 : #spin up
		# create a spin up trace i.e. trace with a blip of E_high 
			line_spin_up =spin_up_trace_sim(t_out, t_in, E_low, sigma_low, E_high, sigma_high, t_max, acq_rate)		  
			spin_traces .append( line_spin_up )
			
	spin_traces =np.array(spin_traces )
	time=np.arange(0,t_max, 1/acq_rate)		# sweep time over trace
	return(time,spin_traces,spin_number)	 # spin_traces : array of 1D traces, spin_number: array with their corresponding spin

def spin_up_trace_sim(t_out, t_in, E_low, sigma_low, E_high, sigma_high, t_max, acq_rate):
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
		   t_max: float
			Time during which data is acquire (in mu s)
			
		  acq_rate:
			Acquisition rate in MHz	  
			
	'''	  
	time=np.arange(0,t_max, 1/acq_rate) # sweep time over trace
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
		and reduces the measurement bandwidth
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
	#I NEED TO CHANGE HOW TO CALCULATE THE MEAS_BW
	meas_bw=1
	
	ret = np.cumsum(time, dtype=float)
	ret[factor:] = ret[factor:] - ret[:-factor] 
	time_filtered=(ret[factor - 1:] / factor)
	time_filtered=time_filtered+factor/acq_rate/2 # Initial time is larger due to the filter
	
	return time_filtered, spin_traces_filtered, meas_bw
 
  
 




# We can put it as a class and there is a function called plot
def max_histogram(spin_traces, bins=20, Plot=1):
	return 


def simulate_histogram_out(A,range_min,range_max,std_down,std_up, signal_low,signal_high,t_int,t_out,t_in,t_max,n_points_time,n_traces,num_bins = 60,plot_histo = 1): 
    max_ss_up = []
    max_ss_down = []
    max_ss_bg = []
    
    minim= np.where(time>= t_min)[0][0]
    maxim= np.where(time>= t.max())[0][0]
    # Do many traces and get max
    for i in range(n_traces):
        
        # Plot the traces for the first 5 traces 
        #if i<5 :
        #    ss_readout = simulate_ss_out(A,std_down,signal_low,signal_high,t_int,t_out,t_in,t_max,n_points_time, 1)

        #else :
        ss_readout = simulate_ss_out(A,std_down,std_up,signal_low,signal_high,t_int,t_out,t_in,t_max,n_points_time, 0)
            
        if ss_readout[2] == 1 : #choose_spin = 1 so spin is down
            max_ss_down.append(ss_readout[0])
            max_ss_bg.append(ss_readout[1])

        if ss_readout[2] == 0 : #choose_spin = 0 so spin is up
            max_ss_up.append(ss_readout[0])
            max_ss_bg.append(ss_readout[1])
        
        #if i==np.int(n_traces/2):
        #    print('Progress>50%')
    
    ma = np.append(np.array(max_ss_up),np.array(max_ss_down))
    max_ss_bg = np.array(max_ss_bg)


    hist_max_ss, bin_edges_max_ss = np.histogram(ma,bins=num_bins,range = (range_min,range_max))    
    hist_max_ss_bg, bin_edges_max_ss_bg = np.histogram(max_ss_bg, bins=num_bins,range=(range_min,range_max) )
    
    bin_axis_ss= bin_edges_max_ss[0:-1]+( bin_edges_max_ss[1]- bin_edges_max_ss[0])/2
    normalisation=( bin_edges_max_ss[1]- bin_edges_max_ss[0])*n_traces
    
    
    
    
    for a in range(len(sim_distance_rem)):    
        


        maximum_sim_rem[a]=sim_distance_40_rem[a][minim:maxim].max()
    range_min=maximum_sim_rem.min()
    range_max=maximum_sim_rem.max()
    
    
    
    
    #print('number of traces for bg =',len(max_ss_bg),'number of up traces =',len(max_ss_up))
    if plot_histo == 1:
        # Plot the histogram 
        plt.figure()
        plt.hist(ma,bins=num_bins,label='total', alpha = 0.6,range = (range_min,range_max))
        plt.hist(max_ss_bg,bins=num_bins ,label='background', alpha = 0.5,range = (range_min,range_max))
        plt.xlabel("Signal")
        plt.ylabel("Probability")
        plt.title('Simulated histogram')
        print('Parameters for the simulation :')
        print('T_int = %f us, T_in = %f us, T_out = %f us' %(t_int,t_in,t_out))
        print('A = %f, std = %f'%(A,std_down))
        plt.legend()
        plt.show()
        #plt.savefig('1histogram.png')

        #Plot the histogram as lines 
        plt.title('Simulated histogram as lines')
        plt.plot(get_inter_bins(bin_edges_max_ss),hist_max_ss, label='signal')
        plt.plot(get_inter_bins(bin_edges_max_ss_bg),hist_max_ss_bg,label ='background')
        plt.legend()
        plt.xlabel('Distance(V)')
        plt.grid()
        #plt.title('Tint = 2.5us')

    return(bin_axis_ss,hist_max_ss/normalisation
    
    
    
    
        # HISTOGRAM WITH SIMULATION
    maximum_sim_rem=np.zeros(len(sim_distance_40_rem))




# we can also have it as a class and it gives as the plot as one of its functions
def fidelity(histogram, epsilon=0.005):
	return fidelity_up, fidelity_down, fidelity_tot




# max_histogram and plot
# fidelity and plot
# __init__.py
# readme.md file
# setup.py file
# license file
