from pytest import raises
from .Simulate_measurement import simulate_measurement

from unittest.mock import MagicMock
from pytest import raises, approx

def test_input_sanity():
    """ Check incorrect input do fail """
    
    E_low=MagicMock()
    sigma_low = MagicMock()
    E_high=MagicMock()
    sigma_high=MagicMock()

    with raises(ValueError) as exception: 
        simulate_measurement(-1, 1, 1, E_low, sigma_low , E_high, sigma_high)
    with raises(TypeError) as exception: 
        simulate_measurement(0.5, 1, 1, E_low, sigma_low , E_high, sigma_high, n_traces=4.5)       
    with raises(ValueError) as exception: 
        simulate_measurement(0.5, -1, -1, E_low, sigma_low , E_high, sigma_high)
    with raises(ValueError) as exception: 
        simulate_measurement(0.5, 1, 1, E_low, sigma_low , E_high, sigma_high, t_max=1, acq_rate=2)
    with raises(ValueError) as exception: 
        simulate_measurement(0.5, 1, 1, E_low, sigma_low , E_high, sigma_high, t_max=-1, acq_rate=-2)