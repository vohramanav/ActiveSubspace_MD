import active_subspaces as ac
import numpy as np

from borehole_functions import *

# Draw M random samples
M = 1000

#Sample the input space according to the distributions in the table above
rw = np.random.normal(.1, .0161812, (M, 1))
r = np.exp(np.random.normal(7.71, 1.0056, (M, 1)))
Tu = np.random.uniform(63070, 115600, (M, 1))
Hu = np.random.uniform(990, 1110, (M, 1))
Tl = np.random.uniform(63.1, 116, (M, 1))
Hl = np.random.uniform(700, 820, (M, 1))
L = np.random.uniform(1120, 1680, (M, 1))
Kw = np.random.uniform(9855, 12045, (M, 1))

x = np.hstack((rw, r, Tu, Hu, Tl, Hl, L, Kw))

#Upper and lower limits for uniform-bounded inputs
xl = np.array([63070, 990, 63.1, 700, 1120, 9855])
xu = np.array([115600, 1110, 116, 820, 1680, 12045])

#XX = normalized input matrix
XX = ac.utils.misc.BoundedNormalizer(xl, xu).normalize(x[:, 2:])

#normalize non-uniform inputs
rw_norm = ((rw - .1)/.0161812).reshape(M, 1)
r_norm = np.log(r); r_norm = ((r_norm - 7.71)/1.0056).reshape(M, 1)

XX = np.hstack((rw_norm, r_norm, XX))

#output values (f) and gradients (df)
f = borehole(XX)
df = borehole_grad(XX)

#Set up our subspace using the gradient samples
ss = ac.subspaces.Subspaces()
ss.compute(df=df, nboot=500)

#Component labels
in_labels = ['rw', 'r', 'Tu', 'Hu', 'Tl', 'Hl', 'L', 'Kw']

#plot eigenvalues, subspace errors
#ac.utils.plotters.eigenvalues(ss.eigenvals, ss.e_br)
#ac.utils.plotters.subspace_errors(ss.sub_br)

#manually make the subspace 2D for the eigenvector and 2D summary plots
ss.partition(2)
#Compute the active variable values
y = XX.dot(ss.W1)

#ac.utils.plotters.eigenvectors(ss.W1, in_labels=in_labels)
ac.utils.plotters.sufficient_summary(y, f)
