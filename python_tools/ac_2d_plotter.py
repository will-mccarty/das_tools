class ac_2d_plotter:
#   import psycopg2
#   import numpy as np
#   import numpy.ma as ma
#   import matplotlib.pyplot as plt
#   import colormaps
#   import matplotlib.cm as mcm
#   import ac_significance as ac_sig

#   def __init__(self):

   
   def ac_2d_plotter(self,startdate=2015080200,enddate=2015093000,expid1='x0016_cld.21z',expid2='x0016_ctl.21z',variables=['h','t','u','v','q'],domains=['n.hem','s.hem']):
   
      import psycopg2
      import numpy as np
      import numpy.ma as ma
      import matplotlib.pyplot as plt
      import colormaps 
      import matplotlib.cm as mcm
      import ac_significance as ac_sig
      
      conn = psycopg2.connect("dbname='semper' user='gmao_user' host='dpdb' ")
      cur = conn.cursor()
      cm = colormaps.colormaps()
      ac = ac_sig.ac_significance()
      
      
#      startdate=2015080200
#      enddate=2015093000
      statistic='cor'
#      #statistic='rms'
#      #expid1='e572p2_fp'
#      #expid2='g5ncep'
#      expid1='x0016_cld.21z'
#      expid2='x0016_ctl.21z'
#      variables=['h','t','u','v']
#      domains=['n.hem','s.hem']
      
      
      levs=np.array([1000.0,850.0,700.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0])
      #levs=np.array([1000.0,850.0,700.0,500.0,400.0,300.0,200.0,100.0])
      
      steps=np.array([0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114, 120])
      steps=np.array([0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120])
      
      #steps=np.array([0, 24, 48, 72, 96, 120])
      
      print (steps.size,levs.size)
      corarr=np.zeros((steps.size,levs.size))
      pctarr=np.zeros((steps.size,levs.size))
      sigarr=np.zeros((steps.size,levs.size))
      levs2d=np.zeros((steps.size,levs.size))
      stps2d=np.zeros((steps.size,levs.size))
      
      pltlevs  = np.array([1000.0,850.0,700.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0])
      pltlevs  = np.array([1075.0,925.0,750,600.,450.,350.,275.,225.,175.,125.,75.])
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
            ac = ac_sig.ac_significance()
            for i in np.arange(steps.size):
              for j in np.arange(levs.size):
            #     j = 1
                 cstp = steps[i]
                 clev = levs[j]
            
                 query = " SELECT value from fc_exp.v_view where date >= {} and date <= {} and level = {} and step = {} and domain_name = '{}' and statistic = '{}' and expver = '{}' and variable = '{}' ".format(startdate, enddate, clev, cstp, dom, statistic, expid1,var)
                 cur.execute(query)
                 val = cur.fetchall()
                 exp1 = np.array(val)
            
                 query = " SELECT value from fc_exp.v_view where date >= {} and date <= {} and level = {} and step = {} and domain_name = '{}' and statistic = '{}' and expver = '{}' and variable = '{}'  ".format(startdate, enddate, clev, cstp, dom, statistic, expid2,var)
                 cur.execute(query)
                 val = cur.fetchall()
                 exp2 = np.array(val)

                 n = exp2.size
#                 print('Starting...',cstp,clev,dom,statistic,expid1,expid2,var,startdate,enddate)
                 print(exp1.size,exp2.size)	
                 diff, sigrange, sig = ac.get_ztran_diff(exp1, exp2)
                 corarr[i,j] = diff
                 sigarr[i,j] = sig
                 pctarr[i,j] = diff /  np.mean(exp2)
            #     corarr[i,j] = val
                 print(cstp,clev,i,j,corarr[i,j],sig)


            fig= plt.figure(figsize=(6,5)) 
            ax = fig.add_subplot(111)

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
            
            fg = plt.pcolor(stps2d,levs2d,corarr,cmap=cm.newjet,vmin=rng * -1.0 ,vmax=rng)
            #plt.pcolormesh(stps2d,levs2d,pctarr ,cmap=cm.newjet,vmin=prng * -1.0 ,vmax=prng)
            plt.xlim([0 - 6,120 + 6])
            plt.ylim([1075,75])
            
            plt.yscale('log')
            
            plt.yticks(levs, levs)
            plt.xticks(steps)
            
            plt.xlabel('Forecast Hour')
            plt.ylabel('Pressure (hPa)')

            title = "Anomaly Correlation Difference\n{} minus {}\n{} - {}, n={}".format(expid1, expid2, startdate, enddate, n)
            plt.title(title)
            ac.oplot_sig_hatch(ax,sigarr,pltsteps,pltlevs)
            
            cbar= plt.colorbar()
            cbar.set_label('Difference')

            
            #sigarr = ma.masked_where(corarr > 0.005, corarr)
            #sigarr = np.nan_to_num(sigarr)
            #print sigarr
            #plt.pcolormesh(stps2d,levs2d,sigarr,vmin=0.0,vmax=0.006,cmap=mcm.binary,alpha=0.5)
            
            fn = '{}_z.{}.{}-{}.{}-{}.{}.png'.format(statistic, var, expid1, expid2, startdate, enddate, dom)
            print(fn) 
      
            fig.savefig(fn)
            
            
            plt.close(fig)
            
            
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     #plt.show()
