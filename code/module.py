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


class Equation:
     
   def __init__(self,descr=None,dim=_sage_const_0 ,indexlist=None,eqstr=None):
     "Indexed set of equations (each equation is a string, to be evaluated) es. e=Equation(descr,dim,indexlist,eqstr)\n§i,§j,§k replace index elements"  
     self.description=descr
     self.dimension=dim
     self.v={}
     if dim==_sage_const_0 :
         self.v[_sage_const_0 ]=eqstr
     if dim==_sage_const_1 :
         for i in range(len(indexlist)):
             self.v[indexlist[i]]=eqstr.replace('§i','"'+indexlist[i]+'"')
     elif dim==_sage_const_2 :
         for i in range(len(indexlist[_sage_const_0 ])):
             for j in range(len(indexlist[_sage_const_1 ])):
                 self.v[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j]]=eqstr.replace('§i','"'+indexlist[_sage_const_0 ][i]+'"').replace('§j','"'+indexlist[_sage_const_1 ][j]+'"')
     elif dim==_sage_const_3 :
         for i in range(len(indexlist[_sage_const_0 ])):
             for j in range(len(indexlist[_sage_const_1 ])):   
                 for k in range(len(indexlist[_sage_const_2 ])):
                     self.v[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j],indexlist[_sage_const_2 ][k]]=eqstr.replace('§i','"'+indexlist[_sage_const_0 ][i]+'"').replace('§j','"'+indexlist[_sage_const_1 ][j]+'"').replace('§k','"'+indexlist[_sage_const_2 ][k]+'"')


class Parameter: # indexed parameters
    
  def __init__(self,descr=None,dim=_sage_const_0 ,indexlist=None,val=_sage_const_0 ):
     "Indexed parameter es. p=Parameter(descr,dim,indexlist,mat)"
     self.description=descr
     self.dimension=dim
     self.v={}
     if dim==_sage_const_0 :
         self.v[_sage_const_0 ]=val
     if dim==_sage_const_1 :
         for i in range(len(indexlist)):
             self.v[indexlist[i]]=val[i]
     elif dim==_sage_const_2 :
         for i in range(len(indexlist[_sage_const_0 ])):
             for j in range(len(indexlist[_sage_const_1 ])):
                 self.v[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j]]=val[i,j]
     elif dim==_sage_const_3 :
         for i in range(len(indexlist[_sage_const_0 ])):
             for j in range(len(indexlist[_sage_const_1 ])):   
                 for k in range(len(indexlist[_sage_const_2 ])):
                     self.v[indexlist[_sage_const_0 ][i],indexlist[_sage_const_1 ][j],indexlist[_sage_const_2 ][k]]=val[i,j,k]


class Model: # define models
    
    def __init__(self,eqlist,varlist):
      "Initialize equations and variables list es. m=Model(eqlist,varlist)\nlists only include names (without .v)"
      self.elist=[]
      self.vlist=[]
      self.initialvalues=[]
      self.ranges=[]
      self.solution={}
      for i in range(len(eqlist)):
        for j in range(len(eqlist[i].v.values())):
            self.elist=self.elist+[eval(list(eqlist[i].v.values())[j])]
      for i in range(len(varlist)):
          for j in range(len(varlist[i].v.values())):
              self.vlist=self.vlist+[list(varlist[i].v.values())[j]]
              self.ranges=self.ranges+[list(varlist[i].b.values())[j]]
              self.initialvalues=self.initialvalues+[list(varlist[i].i.values())[j]]
                 
    def plf(self,v):
          self.vlist=map(float,v)
          return self.slf(*v)         
                    
    def equations(self):
        "Show model equations"
        for i in range(len(self.elist)):
          print (self.elist[i])
        
    def variables(self):
        "Show model variables"
        for i in range(len(self.vlist)):
           print (self.vlist[i])
    
    def solve(self):
        "Solve the model, returning a dictionary (self.solution)"
        solout=solve(self.elist,self.vlist,solution_dict=True)
        if len(solout)==_sage_const_0 :
          print('No solution found! Try nsolve')
        elif len(solout)==_sage_const_1 :  
          self.solution=solout[_sage_const_0 ]
          for i in range(len(solout[_sage_const_0 ].values())):
            try:
              print(str(solout[_sage_const_0 ].keys()[i])+' -> '+str(solout[_sage_const_0 ].values()[i].N(digits=_sage_const_10 )))
            except TypeError:
              print(str(solout[_sage_const_0 ].keys()[i])+' -> '+str(solout[_sage_const_0 ].values()[i]))
        else:
          print('Multiple solutions:')
          print(solout)
          
    def nopt(self):
         "Numerical optimization (bundled), returning a dictionary (self.solution)"
         self.sf=_sage_const_0 
         self.max=False
         self.cl=[]
         self.conl=[]
         for i in range(len(self.elist)):
          if str(self.elist[i].rhs())=='MIN':
           self.sf=self.sf+self.elist[i].lhs()
          elif str(self.elist[i].rhs())=='MAX':
           self.sf=self.sf-self.elist[i].lhs()
           self.max=True
          elif str(self.elist[i]).find(">")>-_sage_const_1 :
           self.cl.append([True,self.elist[i].lhs()-self.elist[i].rhs()])
           print (self.cl)
          elif str(self.elist[i]).find("<")>-_sage_const_1 :
           self.cl.append([True,-self.elist[i].lhs()+self.elist[i].rhs()])
          else:
           self.cl.append([False,self.elist[i].lhs()-self.elist[i].rhs()])
         def pf(v): 
           vdict = dict(zip(self.vlist, v)) 
           return self.sf.subs(vdict)
         def sos(i, v):
           vdict = dict(zip(self.vlist, v)) 
           return self.cl[i][_sage_const_1 ].subs(vdict) 
         for i in range(len(self.cl)):
          if self.cl[i][_sage_const_0 ]:
            ks='ineq'
          else:
            ks='eq'
          self.conl.append({'type':ks,'fun':functools.partial(sos, i)})
         self.solout=op.minimize(pf,self.initialvalues,bounds=self.ranges,constraints=self.conl,options={'disp':True})
         if self.max:
           self.of=-self.solout.fun 
         else:
           self.of=self.solout.fun 
         print('-------------------------------------------------------------------')
         for i in range(len(self.vlist)):
           self.solution[self.vlist[i]]=self.solout.x[i]
           print(str(self.vlist[i])+' = '+str(round(self.solout.x[i])))
         print('-------------------------------------------------------------------')
         print('Objective Function:')
         print(self.of)
         print('-------------------------------------------------------------------')

    def fix(self,v,val):
        "Fix the value of the variable v to value val"
        pos=self.vlist.index(eval(v))
        self.initialvalues[pos]=val
        self.ranges[pos]=(val,val)
        
    def nsolve(self,method='hybr'):
         """
         Numerical solve (root finding), returning a dictionary (self.solution)
         Available methods: 
         hybr(def) lm broyden1 broyden2 anderson linearmixing 
         diagbroyden excitingmixing krylov
         """ 
         self.slf=[]
         for i in range(len(self.elist)):
            self.slf.append(self.elist[i].lhs()-self.elist[i].rhs())
         def pf(v):
          vdict = dict(zip(self.vlist, v))
          r=[]
          for i in range(len(self.slf)):
            a=self.slf[i].subs(vdict)
            r.append(a)
          return r 
         self.solout=op.root(pf,self.initialvalues,jac=False,method=method)
         print(self.solout)
         print('-------------------------------------------------------------------')
         for i in range(len(self.vlist)):
           self.solution[self.vlist[i]]=self.solout.x[i]
           print(str(self.vlist[i])+' = '+str(round(self.solout.x[i],_sage_const_2 )))
         print('-------------------------------------------------------------------')