"""
This code is created by Prof. Roberto Roson. You can found more about this and his others projects here  here https://roson.great-site.net/software.html
"""

from sage.all_cmdline import * 
_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_1p0 = RealNumber('1.0'); _sage_const_10 = Integer(10)
from scipy import optimize as op
import functools
import numpy as np 

class Variable: # define indexed variables

  """
  * descr * 
  is a descriptive string (a string, or “a string”) used for illustrating what the variable is about;
  *dim * 

  is an integer number equal to 0, 1, 2, or 3. Zero means a scalar, one a vector of variables, two a
  matrix, etc.;

  * indexlist * 
  is list of sets of strings, consistent with dim. The sets are standard Python sets, to be defined
  before the variable;

  * name * 
  is a string defining the internal name of the variable. What PAMS does after a variable definition
  is creating several SAGE (SimPy) variables, depending on the dimension. The actual name of
  these variables is given by name, followed by the underscore sign (_), followed by one or more
  elements of indexlist (also separated by underscore). For this reason, it is advisable to have a
  name limited to !!!one or a few letters!!!;

  * bounds *
  is a list of two numerical values, expressing lower and upper bounds. They are only considered
  in constrained optimization;

  * inval *
  it is the initial value assigned to the variable(s) in numerical algorithms (optimization and root
  finding). It can be a single numerical value (in this case applied to all sub-variables) or a vector /
  matrix, with the same dimension as in dim.
  
  """

  
  def __init__(self,descr=None,dim=_sage_const_0 ,indexlist=None,name=None,bounds=(None,None),inval=_sage_const_1p0 ):
     "Indexed variable es. v=Variable(descr,dim,indexlist,name,bounds=(None,None),inval=1.0)\nname is a string, possibly one or two characters"
     self.description=descr
     self.dimension=dim
     self.v={}  # gives a Python dictionary, where the key are elements in indexlist and the corresponding values are the internal SAGE (SimPy) variable names;
     self.b={}  # gives a Python dictionary, where the key are elements in indexlist and the corresponding values are couples of lower and upper bounds;
     self.i={}  # gives a Python dictionary, where the key are elements in indexlist and the corresponding values are the initial points used in numerical algorithms
     if dim==_sage_const_0 : # scalar 
         self.v[_sage_const_0 ]=var(name, domain=RR)
         self.b[_sage_const_0 ]=bounds
         self.i[_sage_const_0 ]=inval
     if dim==_sage_const_1 :  # vector
         for i in range(len(indexlist)):
             self.v[indexlist[i]]=var(name+'_'+indexlist[i], domain=RR)
             self.b[indexlist[i]]=bounds
             if type(inval)==type([_sage_const_1p0 ]):
              self.i[indexlist[i]]=inval[i]
             else:
              self.i[indexlist[i]]=inval
     elif dim==_sage_const_2 :  # matriz 
         for i in range(len(indexlist[_sage_const_0 ])):
             for j in range(len(indexlist[_sage_const_1 ])):
                 self.v[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j]]=var(name+'_'+indexlist[_sage_const_0 ][i]+'_'+indexlist[_sage_const_1 ][j], domain=RR)
                 self.b[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j]]=bounds
                 if type(inval)==type(np.array([_sage_const_1p0 ])):
                  self.i[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j]]=inval[i,j]
                 else:
                  self.i[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j]]=inval
     elif dim==_sage_const_3 :
         for i in range(len(indexlist[_sage_const_0 ])):
             for j in range(len(indexlist[_sage_const_1 ])):   
                 for k in range(len(indexlist[_sage_const_2 ])):
                     self.v[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j],indexlist[_sage_const_2 ][k]]=var(name+'_'+indexlist[_sage_const_0 ][i]+'_'+indexlist[_sage_const_1 ][j]+'_'+indexlist[_sage_const_2 ][k], domain=RR)
                     self.b[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j],indexlist[_sage_const_2 ][k]]=bounds
                     if type(inval)==type(np.array([_sage_const_1p0 ])):
                      self.i[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j],indexlist[_sage_const_2 ][k]]=inval[i,j,k]
                     else:
                      self.i[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j],indexlist[_sage_const_2 ][k]]=inval

