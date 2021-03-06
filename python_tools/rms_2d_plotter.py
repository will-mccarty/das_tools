import psycopg2
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import colormaps
import matplotlib.cm as mcm
from gmaopy.modules.scoreplot import *
from gmaopy.stats.statistics import Statistics, statisticList
from string import Template

class rms_2d_plotter:
#   import psycopg2
#   import numpy as np
#   import numpy.ma as ma
#   import matplotlib.pyplot as plt
#   import colormaps
#   import matplotlib.cm as mcm
#   import ac_significance as ac_sig

#   def __init__(self):

   
   def rms_2d_plotter(self,startdate=2015080200,enddate=2015093000,expid1='x0016_cld.21z',expid2='x0016_ctl.21z',variables=['h','t','u','v','q'],domains=['n.hem','s.hem'],sigplot=False,verbose=False,straightmean=False,confidence=0.95,title=None,verif1='gmao',verif2='gmao',dblevs=None,maxlev=1000.,minlev=100.,levs=None,range=None,normed=False):

      from scipy.stats import ttest_ind   
      import ac_significance as ac_sig

#      import psycopg2
#      import numpy as np
#      import numpy.ma as ma
#      import matplotlib.pyplot as plt
#      import colormaps 
#      import matplotlib.cm as mcm
#      from gmaopy.modules.scoreplot import *
#      from gmaopy.stats.statistics import Statistics, statisticList


      ac = ac_sig.ac_significance(confidence=confidence)
      if levs is None:
          levs=np.array([1000.0,850.0,700.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0])
      else:
          levs = np.array(levs)

      
      conn = psycopg2.connect("dbname='semper' user='gmao_user' host='edb1' ")
      cur = conn.cursor()
      cm = colormaps.colormaps()
      
      
#      startdate=2015080200
#      enddate=2015093000
      statistic='rms'
#      #statistic='rms'
#      #expid1='e572p2_fp'
#      #expid2='g5ncep'
#      expid1='x0016_cld.21z'
#      expid2='x0016_ctl.21z'
#      variables=['h','t','u','v']
#      domains=['n.hem','s.hem']
      
      if (dblevs):
          print("Determining levels from database")
          query = " SELECT distinct level from fc_exp.v_view where date = {} and step = 0 and expver = '{}' and verify = '{}' ORDER BY level ".format(startdate, expid1,verif1)
          cur.execute(query)
          val = cur.fetchall()
          levs = np.array(val)[:,0]
          levs = levs[::-1]
          msk = (levs <= maxlev) & (levs >= minlev)
          levs = levs[msk]
#          levs = val
          
#      else:
#          levs=np.array([1000.0,850.0,700.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0])

      print("Masking levels between {} and {}".format(minlev,maxlev))
      msk = (levs <= maxlev) & (levs >= minlev)
      levs = levs[msk]

      print(levs.size,'Levels=',levs)

      #levs=np.array([1000.0,850.0,700.0,500.0,400.0,300.0,200.0,100.0])
      
      steps=np.array([0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114, 120])
      steps=np.array([0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120])
      
      #steps=np.array([0, 24, 48, 72, 96, 120])
      
      corarr=np.zeros((steps.size,levs.size))
      pctarr=np.zeros((steps.size,levs.size))
      sigarr=np.zeros((steps.size,levs.size))
      exp1arr = np.zeros((steps.size,levs.size))
      exp2arr = np.zeros((steps.size,levs.size))

      levs2d=np.zeros((steps.size,levs.size))
      stps2d=np.zeros((steps.size,levs.size))
      
#      pltlevs  = np.array([1000.0,850.0,700.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0])
#      pltlevs  = np.array([1075.0,925.0,750,600.,450.,350.,275.,225.,175.,125.,75.])
      pltlevs = (levs + np.roll(levs,1))/2.
      pltlevs[0] = levs[0] - (pltlevs[1]-levs[0])
      pltlevs = np.append(pltlevs,levs[-1] - (levs[-2]-levs[-1])*0.5)
      pltsteps = np.array([0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132]) - 6
      
      
      levs2d=np.zeros((pltsteps.size,pltlevs.size))
      stps2d=np.zeros((pltsteps.size,pltlevs.size))
      
      
      #500.0,400.0,300.0,250.0,200.0,150.0,100.0])
      
      for j in np.arange(pltlevs.size):
         levs2d[:,j] = pltlevs[j]
      
      for i in np.arange(pltsteps.size):
         stps2d[i,:] = pltsteps[i]
      
      for var in variables:
         for dom in domains:
            for i in np.arange(steps.size):
              for j in np.arange(levs.size):
            #     j = 1
                 cstp = steps[i]
                 clev = levs[j]
            
                 query = " SELECT value from fc_exp.v_view where date >= {} and date <= {} and level = {} and step = {} and domain_name = '{}' and statistic = '{}' and expver = '{}' and variable = '{}' and verify = '{}' ORDER BY date ".format(startdate, enddate, clev, cstp, dom, statistic, expid1,var,verif1)
#                 print(query)
                 cur.execute(query)
                 val = cur.fetchall()
                 exp1 = np.array(val)
            
                 query = " SELECT value from fc_exp.v_view where date >= {} and date <= {} and level = {} and step = {} and domain_name = '{}' and statistic = '{}' and expver = '{}' and variable = '{}' and verify = '{}' ORDER BY date ".format(startdate, enddate, clev, cstp, dom, statistic, expid2,var,verif2)
#                 print(query)
                 cur.execute(query)
                 val = cur.fetchall()
                 exp2 = np.array(val)

                 n = exp2.size
                 print('Starting...',cstp,clev,dom,statistic,expid1,expid2,var,startdate,enddate,'size=',n,'confidence=',confidence)
                 if (verbose):
                     print('  i,j=',i,j)
                     print(exp1.T)	
                     print(exp2.T)
                 #print(exp2)
#                 ztdiff, sigrange, sig, expcor1, expcor2 = ac.get_ztran_diff(exp1, exp2)
                 if (verbose): print('  NP Mean exp1, exp2, Diff:',np.mean(exp1),np.mean(exp2),np.mean(exp1) - np.mean(exp2))
                 t, prob = ttest_ind(exp1,exp2,equal_var=False,nan_policy='omit')
                 conf = 1 - prob
                 if (verbose): print('  confidence(exp,thresh):',conf,confidence)
 
                 exp1arr[i,j] = np.mean(exp1)
                 exp2arr[i,j] = np.mean(exp2)

                 if (normed):
                     corarr[i,j] = (exp1arr[i,j] - exp2arr[i,j])/exp2arr[i,j]*100 # mndiff
                 else:
                     corarr[i,j] = exp1arr[i,j] - exp2arr[i,j] # mndiff
                 sigarr[i,j] = False if conf < confidence else True
#                 sigarr[i,j] = False if mndiff < upper and mndiff > lower else True
#                 pctarr[i,j] = mndiff /  np.mean(exp2)
            #     corarr[i,j] = val
#                 print(cstp,clev,i,j,corarr[i,j],sigarr[i,j])


            fig= plt.figure(figsize=(6,5)) 
            ax = fig.add_subplot(111)
            if (verbose): print('min/max diff:',np.min(corarr),np.max(corarr))
            plt.subplots_adjust(top=0.85, left=0.18)
#            plt.tight_layout() 
            #ac.oplot_sig_hatch(ax,sigarr,pltsteps,pltlevs)
            mx=np.max(corarr)
            mn=np.min(corarr)
            rng=max((np.abs(mx),np.abs(mn)))
            
            pctarr = pctarr * 100
            pmx=np.max(pctarr)
            pmn=np.min(pctarr)
            prng=max((np.abs(pmx),np.abs(pmn)))
            
            if range is not None:
                rng = range

            fg = plt.pcolor(stps2d,levs2d,ma.masked_invalid(corarr),cmap=cm.newjet,vmin=rng * -1.0 ,vmax=rng)
            #plt.pcolormesh(stps2d,levs2d,pctarr ,cmap=cm.newjet,vmin=prng * -1.0 ,vmax=prng)
            plt.xlim([0 - 6,120 + 6])
#            plt.ylim([1075,75])
            plt.ylim([pltlevs[0],pltlevs[-1]])           
            plt.yscale('log')
            
            plt.yticks(levs, levs)
            plt.xticks(steps)
            
            plt.xlabel('Forecast Hour')
            plt.ylabel('Pressure (hPa)')

            if (title == None):
                outtitle = "RMS Difference\n{} minus {}\n{} - {}, n={}".format(expid1, expid2, startdate, enddate, n)
            else:
                if (dom == 's.hem'): longdom='Southern Hemisphere'
                if (dom == 'n.hem'): longdom='Northern Hemisphere'
                if (dom == 'tropics'): longdom='Tropics'
                if (dom == 'global'): longdom='Global'

                tmpl = Template(title)
#                print(dom)
                outtitle = tmpl.substitute(expid1=expid1, expid2=expid2, startdate=startdate, enddate=enddate, n=n, var=var, dom=dom, longdom=longdom)
            plt.title(outtitle)
            if (sigplot):
                ac.oplot_sig_hatch(ax,sigarr,pltsteps,pltlevs)
            
            cbar= plt.colorbar()
            if (normed):
                cbar_label = 'Difference (%)'
            else:
                cbar_label = 'Difference'
            cbar.set_label(cbar_label)

            
            #sigarr = ma.masked_where(corarr > 0.005, corarr)
            #sigarr = np.nan_to_num(sigarr)
            #print sigarr
            #plt.pcolormesh(stps2d,levs2d,sigarr,vmin=0.0,vmax=0.006,cmap=mcm.binary,alpha=0.5)
            
            fn = '{}_z.{}.{}-{}.{}-{}.{}.verif-{}.{}.png'.format(statistic, var, expid1, expid2, startdate, enddate, dom,verif1,verif2)
            print('Writing: ',fn) 
      
            fig.savefig(fn)
            
            
            plt.close(fig)

      if (verbose):print('returning ',dom,var)            
      return(corarr,sigarr,exp1arr,exp2arr)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     #plt.show()
