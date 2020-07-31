from cbcbeat.splittingsolver import SplittingSolver, BasicSplittingSolver

__all__ = ['ElectroSolver']

class ElectroSolver(SplittingSolver):
  """ Solver for the electrical activity in the heart

  This class inherits all of the functionality from the SplittingSolver 
  from CBCBeat. Decided to do it this way in case it was needed to modify it
  at a later stage. It solves the bidomain and monodomain equations presented
  in Sundnes (2006). For more information, see the documentation for the 
  SplittingSolver in CBCBeat. 

  It solves the bidomain (or monodomain) equations given by: 
  

  """ 
  pass 

