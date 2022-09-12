import numpy as np

def simulate_measurement(A, t_out, t_in, E_low, sigma_low , E_high, sigma_high,  n_traces=1, t_max=2000, acq_rate=1):
    """ Creates measurement traces using the input parameters and returns the maximum of each trace
      Parameters
      ----------
       
      Dot parameters: 
          A: float from 0 to 1
              Proportion of spin down traces
              
          t_out: float
            Time constant for a spin up electron to leave the dot (in \mu s)

          t_in: float
            Time constant for a spin down electron to return to the dot (in \mu s)
      
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
            Time during which data is acquire (in \mu s)
            
          acq_rate:
            Acquisition rate in MHz

  """
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


    ss = []
    spin_number=np.empty(0)
    n_points_time=int(t_max*acq_rate)
    
    for i in range(n_traces):
        choose_spin = np.random.binomial(1,A,1) #create 1 with proba A and 0 with proba (1-A)
        spin_number = np.append(spin_number, choose_spin)
        
        if choose_spin == 1 : # spin down because A = <P_down>
        # create a spin down trace i.e. constant E_low with noise sigma_low
            line_spin_down = E_low + np.random.normal(0,sigma_low ,n_points_time)
            ss.append(line_spin_down )
            
            
        if choose_spin == 0 : #spin up
        # create a spin up trace i.e. trace with a blip of E_high 
            line_spin_up =spin_up_trace_sim(t_out, t_in, E_low, sigma_low, E_high, sigma_high, t_max, acq_rate)       
            ss.append( line_spin_up )
            
    ss=np.array(ss)
    return(ss,spin_number)  # ss: array of 1D traces, spin_number: array with their corresponding spin


def spin_up_trace_sim(t_out, t_in, E_low, sigma_low, E_high, sigma_high, t_max, acq_rate):
    """ Creates a spin up measurement trace
      Parameters
      ----------
       
      Dot parameters: 
          t_out: float
            Time constant for a spin up electron to leave the dot (in \mu s)

          t_out: float
            Time constant for a spin down electron to return to the dot (in \mu s)
      
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
            Time during which data is acquire (in \mu s)
            
          acq_rate:
            Acquisition rate in MHz   
    """   
    sweep=np.arange(0,t_max, 1/acq_rate) # time sweep over trace
    ss_x=np.zeros(len(sweep)) # spin-up trace
    
    t_leave=np.random.exponential(t_out) #random t_out time
    t_enter=np.random.exponential(t_in) #random t_in time
    
    i_leave = np.where(sweep<=t_leave)[0][-1] 
    i_enter = np.where(sweep<=t_enter+t_leave)[0][-1]
    
    
    ss_x[0:i_leave] = E_low + np.random.normal(0,sigma_low,len(sweep[0:i_leave])) # first part of trace (background)

    if i_enter == i_leave : #means blip time = 0
        ss_x[i_leave:] =  E_low + np.random.normal(0,sigma_low,len(sweep[i_leave:])) 

    else : 
        ss_x[i_leave:i_enter] = E_high + np.random.normal(0,sigma_low,len(sweep[i_leave:i_enter])) # blip 
        ss_x[i_enter:] = E_low + np.random.normal(0,sigma_low,len(sweep[i_enter:])) # last part of the trace (background)

    return ss_x
