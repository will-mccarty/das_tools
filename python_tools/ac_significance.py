import numpy as np

class ac_significance:

   def __init__(self): 
      self.confidence_int = 0.9



   def get_ztran_diff(self, expcor1, expcor2):
      from gmaopy.stats.critval import critval


      diff = expcor1 - expcor2
      ztdiff = 0.5 * np.log( (1.0 + 0.5*diff)  / (1.0 - 0.5*diff) )
   
      ztmn  = np.mean(ztdiff)
      ztvar = np.var(ztdiff)

      dof   = diff.size - 1
      crit  = critval(self.confidence_int, dof)
      zcrit = crit * ( (ztvar / dof) ** (0.5) )

      cordiff = 2 * ( (np.exp(2 * ztmn) - 1)  / (np.exp(2 * ztmn) + 1) )
      corup   = 2 * ( (np.exp(2 * zcrit) - 1)  / (np.exp(2 * zcrit) + 1) )
      corlow  = 2 * ( (np.exp(-2 * zcrit) - 1)  / (np.exp(-2 * zcrit) + 1) )

      corsig = False

      if (cordiff > corup or cordiff < corlow):
         corsig = True

      print( cordiff, corup, corlow, corsig)
      

      return cordiff, [corup, corlow], corsig
   
      print( np.mean(diff), ztmn, ztvar, crit)
      print(cordiff, corup, corlow, corsig)



   def oplot_sig_hatch(self,ax,sigarr,x,y):
      import matplotlib.patches as mpatches
      from matplotlib.collections import PatchCollection

      sz = sigarr.shape
      print(sz)
      patches = []

      for i in np.arange(sz[0]):
         for j in np.arange(sz[1]):
#            if True:
            if (sigarr[i,j]):
               cx = x[i]
               cy = y[j]
               wd = x[i+1] - x[i]
               ht = y[j+1] - y[j]
#               cpatch = mpatches.Rectangle([cx,cy],wd,ht,color='White',alpha=0.7, linewidth=0)
               cpatch = mpatches.Rectangle([cx,cy],wd,ht,hatch='..',color='Gray',linewidth=0, fill=None )
               ax.add_patch(cpatch)
#               patches.append(cpatch)

#      collection = PatchCollection(patches)

#      ax.add_collection(collection)


