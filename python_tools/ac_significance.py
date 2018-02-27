import numpy as np

class ac_significance:

   def __init__(self,confidence=0.9): 
#      print('setting confidence to ',confidence)
      self.confidence_int = confidence



   def get_ztran_diff(self, expcor1, expcor2):
      from gmaopy.stats.critval import critval

      tiny = 1e-6

      diff = expcor1 - expcor2
      zexpcor1 = 0.5 * np.log( (1.0 + expcor1)  / (1.0 - expcor1 + tiny) )
      zexpcor2 = 0.5 * np.log( (1.0 + expcor2)  / (1.0 - expcor2 + tiny) )

#      ztdiff = 0.5 * np.log( (1.0 + diff)  / (1.0 - diff) )
#      ztdiff = zexpcor1 - zexpcor2
   
      ztmn1  = np.mean(zexpcor1)
      ztmn2  = np.mean(zexpcor2)

      expcor1mn = (np.exp(2*ztmn1) - 1) / (np.exp(2*ztmn1) + 1) 
      expcor2mn = (np.exp(2*ztmn2) - 1) / (np.exp(2*ztmn2) + 1)

      ztdiff = 0.5 * np.log( (1.0 + 0.5*diff)  / (1.0 - 0.5*diff) )

      ztmn  = np.mean(ztdiff)
      ztvar = np.var(ztdiff)

      dof   = diff.size #- 1
#      print('setting confidence to ',self.confidence_int)
      crit  = critval(self.confidence_int, dof)
#      print(crit)
      zcrit = crit * ( np.sqrt(ztvar / dof) )
#      cordiff =  2 * ( (np.exp(2 * ztmn) - 1)  / (np.exp(2 * ztmn) + 1) )
      cordiff = expcor1mn-expcor2mn
      corup   =  2 * ( (np.exp(2 * zcrit) - 1)  / (np.exp(2 * zcrit) + 1) )
      corlow  =  2 * ( (np.exp(-2 * zcrit) - 1)  / (np.exp(-2 * zcrit) + 1) )

#      print(crit,cordiff,corup,corlow)
      corsig = False

      if (cordiff > corup or cordiff < corlow):
         corsig = True

      

      return cordiff, [corup, corlow], corsig, expcor1mn, expcor2mn
   
#      print( np.mean(diff), ztmn, ztvar, crit)
#      print(cordiff, corup, corlow, corsig)



   def oplot_sig_hatch(self,ax,sigarr,x,y):
      import matplotlib.patches as mpatches
      from matplotlib.collections import PatchCollection

      sz = sigarr.shape
#      print(sz)
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


