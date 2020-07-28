import pytest
from m3h3 import *

def test_add_element(params):
    assert "test" not in params.keys()
    params.add("test", 2)
    assert "test" in params.keys()

def test_set_electro_parameters(params):
    params.set_electro_parameters()

def test_update_electro_parameters(params):
    pass

def test_add_electro_parameters(params):
    pass 

def test_problem_specifications(params):
    problem_specifications =(params["problem_specifications"])
    assert problem_specifications["stimulus"] == None
    assert problem_specifications["applied_current"] == None
    assert problem_specifications["initial_conditions"] == None

def test_problem_specifications_add_stimulus(params):
    problem_specifications =(params["problem_specifications"])
    problem_specifications["stimulus"] = Expression("10*x[0]*t", t = Constant(0.0), degree=1)
    assert problem_specifications["stimulus"][0] == Expression("10*x[0]*t", t = Constant(0.0), degree=1)[0]
    problem_specifications["stimulus"] = Expression("10", t = Constant(0.0), degree = 1)
    assert problem_specifications["stimulus"].str() != Expression("10*x[0]*t", t = Constant(0.0), degree = 1).str()
    assert problem_specifications["stimulus"].str() == Expression("10", t = Constant(0.0), degree = 1).str()

@pytest.fixture 
def params():
    return m3h3.Parameters("M3H3")