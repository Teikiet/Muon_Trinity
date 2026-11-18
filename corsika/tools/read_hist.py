import numpy as np
import matplotlib.pyplot as plt
import boost_histogram as bh

def read_hist(filename):
    """
    read numpy file produced with CORSIKA 8's save_hist() function into
    boost-histogram object.
    """

    d = np.load(filename)
    axistypes = d['axistypes'].view('c')
    overflow = d['overflow']
    underflow = d['underflow']

    axes = []
    for i, (at, has_overflow, has_underflow) in enumerate(zip(axistypes, overflow, underflow)):
        if at == b'c':
            axes.append(bh.axis.Variable(d[f'binedges_{i}'], overflow=has_overflow, underflow=has_underflow))
        elif at == b'd':
            axes.append(bh.axis.IntCategory(d[f'bins_{i}'], growth=(not has_overflow)))
 
    h = bh.Histogram(*axes)
    h.view(flow=True)[:] = d['data']
    
    return h

