import obs_database_access as obdb
import pandas as pd
import numpy as np

def mask_mutual_invalids(a,b):
    a = np.ma.masked_invalid(a)
    b = np.ma.array(b,mask = a.mask)
    b = np.ma.masked_invalid(b)
    a = np.ma.array(a,mask = b.mask)

    return(a,b)


class scorecard_from_config():

    def __init__(self,cfg_file, startdate=None, enddate=None, exp=None, ctl=None, exp_database=None, ctl_database=None):
        self.load_config(cfg_file)

        self.startdate=startdate if startdate is not None else self.config['startdate']
        self.enddate=enddate if enddate is not None else self.config['enddate']
        self.exp=exp if exp is not None else self.config['exp']['name']
        self.ctl=ctl if ctl is not None else self.config['ctl']['name']
        self.exp_database=exp_database if exp_database is not None else self.config['exp']['database']
        self.ctl_database=ctl_database if ctl_database is not None else self.config['ctl']['database']

        
        self.sc_cfg = scorecard(self.exp, self.ctl,
                            exp_database=self.exp_database,
                            ctl_database=self.ctl_database,
                            startdate=self.startdate,
                            enddate=self.enddate)
        self.loop_obs_classes()

        self.sc_cfg.plot_scores(self.score_df)

        print(self.val_df)
        print(self.score_df)

    def load_config(self,cfg_file):
        import yaml
        
        print('Initializing from %s' % cfg_file)
        with open(cfg_file, 'r') as f:
            self.cfg_all    = yaml.load(f)
            self.config     = self.cfg_all['config']
            self.class_type = self.cfg_all['class_type']
            self.classes    = self.cfg_all['obs_classes']

    def loop_obs_classes(self):

        query = obdb.query_structure(variable='omf',usage='used',domain_name='global')

        self.val_df = None
        self.score_df = None

        for cls in self.classes:
            cur_obs = self.classes[cls]['types']
            cur_obs_list = sorted(self.classes[cls]['types'])
            print(cur_obs)
            for ob in cur_obs_list:
                query = obdb.query_structure(variable='omf',usage='used',domain_name='global',kt=self.class_type[cls]['kt'])
                cur_ob = cur_obs[ob]

                for key in cur_ob:
                   if key in query: query[key] = cur_ob[key]
              
                print('Retrieving %s' % cur_ob['name'])         
                idx = pd.MultiIndex.from_arrays([[cls],[cur_ob['name']]],names=['Class','Type'])
                self.val_df, self.score_df = self.sc_cfg.scorecard_get_diff_scores(query,
                                                                                   in_val_df=self.val_df,
                                                                                   in_score_df=self.score_df,
                                                                                   index_name=idx)
#                                                                                   index_name=cur_ob['name'])


                

class scorecard():
#    import obs_database_access as obdb
#    import pandas as pd
#    import numpy as np

    def __init__(self, exp, ctl, exp_database='ob_exp', ctl_database='ob_exp',
                 startdate=None, enddate=None):

        self.startdate=startdate
        self.enddate=enddate
        self.exp = obdb.experiment(exp,database=exp_database,startdate=self.startdate,enddate=self.enddate)
        self.exp_name=exp
        self.ctl = obdb.experiment(ctl,database=ctl_database,startdate=self.startdate,enddate=self.enddate)
        self.ctl_name=ctl


    def scorecard_get_difference(self,query,index_name=None,in_df=None):
        exp_bias = self.exp.get_bias_from_list(query,'domain_name',[('global'),('n.hem'), ('tropics'), ('s.hem')])
        ctl_bias = self.ctl.get_bias_from_list(query,'domain_name',[('global'),('n.hem'), ('tropics'), ('s.hem')])
        bias_abs_diff = np.abs(exp_bias)-np.abs(ctl_bias)
    
        exp_rms = self.exp.get_rms_from_list(query,'domain_name',[('global'),('n.hem'), ('tropics'), ('s.hem')])
        ctl_rms = self.ctl.get_rms_from_list(query,'domain_name',[('global'),('n.hem'), ('tropics'), ('s.hem')])
        rms_diff = exp_rms-ctl_rms
        arr = np.array([np.append(bias_abs_diff,rms_diff)])

        if (index_name is not None):
            index_name=[index_name]
        df = pd.DataFrame(arr,columns =['GL_Bias','NH_Bias','TR_Bias','SH_Bias','GL_RMS','NH_RMS','TR_RMS','SH_RMS'],index=index_name)
        if (in_df is not None):
            df = in_df.append(df)
        return(df)

    def scorecard_get_diff_scores(self,in_query,index_name=None,in_score_df=None, in_val_df=None):
        from copy import deepcopy
        from scipy.stats import ttest_ind

        domain = [('global'),('n.hem'), ('tropics'), ('s.hem')]
        query = deepcopy(in_query)

        ts = self.exp.get_ndate_timeseries()

        bias_abs_diff_arr = np.array([])
        rms_diff_arr = np.array([])
        bias_diff_score_arr = np.array([])
        rms_diff_score_arr  = np.array([])

        for dom in domain:
            query['domain_name'] = dom
            exp_bias = self.exp.get_bias_from_list(query,'date',ts)
            ctl_bias = self.ctl.get_bias_from_list(query,'date',ts)
            exp_total_bias = self.exp.get_bias(query)
            ctl_total_bias = self.ctl.get_bias(query)
#            exp_bias, ctl_bias = mask_mutual_invalids(exp_bias,ctl_bias)            
            bias_abs_diff = np.abs(exp_total_bias) - np.abs(ctl_total_bias)
#            for t, eb, cb in zip(ts, exp_bias, ctl_bias):
#               print(t,eb,cb)
#            print(exp_bias)
#            print(ctl_bias)
            t, prob = ttest_ind(exp_bias,ctl_bias,equal_var=False,nan_policy='omit')
            conf = 1.0 - prob
#            print(dom,bias_abs_diff,prob)
            bias_abs_diff_arr = np.append(bias_abs_diff_arr,bias_abs_diff)
            if (conf < 0.05):
                score = 0
            elif (conf > 0.05 and conf < 0.90):
                score = -1 if bias_abs_diff > 0.0 else 1
            elif (conf > 0.90 and conf < 0.99):
                score = -2 if bias_abs_diff > 0.0 else 2
            elif (conf > 0.99):
                score = -3 if bias_abs_diff > 0.0 else 3
            else:
                score = None
            bias_diff_score_arr = np.append(bias_diff_score_arr,score)

            exp_rms = self.exp.get_rms_from_list(query,'date',ts)
            ctl_rms = self.ctl.get_rms_from_list(query,'date',ts)
            exp_total_rms = self.exp.get_rms(query)
            ctl_total_rms = self.ctl.get_rms(query)

            rms_diff = exp_total_rms - ctl_total_rms

            t, prob = ttest_ind(exp_rms,ctl_rms,equal_var=False,nan_policy='omit')
            conf = 1.0 - prob
#            print(dom,rms_diff,prob)

            rms_diff_arr = np.append(rms_diff_arr,rms_diff)
            if (conf < 0.05):
                score = 0
            elif (conf > 0.05 and conf < 0.90):
                score = -1 if rms_diff > 0.0 else 1
            elif (conf > 0.90 and conf < 0.99):
                score = -2 if rms_diff > 0.0 else 2
            elif (conf > 0.99):
                score = -3 if rms_diff > 0.0 else 3

            rms_diff_score_arr = np.append(rms_diff_score_arr,score)


#            exp_rms = self.exp.get_rms_from_list(query,'date',ts)
#            ctl_rms = self.ctl.get_rms_from_list(query,'date',ts)
#            rms_diff = exp_rms-ctl_rms

        val_arr = np.array([np.append(bias_abs_diff_arr,rms_diff_arr)]) #  arr = np.array([np.append(bias_abs_diff,rms_diff)])        
        score_arr = np.array([np.append(bias_diff_score_arr,rms_diff_score_arr)])

        if (index_name is not None): 
            index_name=[index_name]
#        df = pd.DataFrame(arr,columns =['GL_Bias','NH_Bias','TR_Bias','SH_Bias','GL_RMS','NH_RMS','TR_RMS','SH_RMS'],index=index_name)
        c_val_df   = pd.DataFrame(val_arr   ,columns =['GL_Bias','NH_Bias','TR_Bias','SH_Bias','GL_RMS','NH_RMS','TR_RMS','SH_RMS'],index=index_name)
        c_score_df = pd.DataFrame(score_arr ,columns =['GL_Bias','NH_Bias','TR_Bias','SH_Bias','GL_RMS','NH_RMS','TR_RMS','SH_RMS'],index=index_name)
        
        if (in_score_df is not None):
            in_score_df = in_score_df.append(c_score_df)
        else:
            in_score_df = c_score_df

        if (in_val_df is not None):
#            in_val_df = in_val_df.append(c_val_df)
            in_val_df = c_val_df.append(in_val_df)
        else:
            in_val_df = c_val_df

        return(in_val_df,in_score_df)

    def plot_scores(self,score_df,filename=None,show=True,save=True):
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        import seaborn as sns

        cmap = colors.LinearSegmentedColormap.from_list('custom scores', 
                                             [(0,    '#FF0000'),
                                              (0.167,'#EC9D84'),
                                              (0.333,'#FBDED0'),
                                              (0.5,  '#FFFFFF'),
                                              (0.667,'#DFEEF5'),
                                              (0.833,'#79AED2'),
                                              (1.0  ,'#0000FF')], N=256)


        if (filename is None):
            filename='scorecard.{}.{}.{}-{}.png'.format(self.exp_name,self.ctl_name,self.startdate,self.enddate)
#            filename='scorecard-'+self.exp_name+'-'+self.ctl_name+'.png'

        titlestr='EXP: {} CTL: {}\n{}-{}'.format(self.exp_name,self.ctl_name,self.startdate,self.enddate)

        fig, ax = plt.subplots()
        height = score_df.shape[0] * 0.225 + 0.5
        print(height)
        fig.set_size_inches(5,height)
        bounds = np.array([-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5])
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

        print(score_df)
        score_df.sort_index(inplace=True)
        #score_df = pltdf
        print(score_df)
        pcm = ax.pcolor(np.ma.masked_invalid(score_df),edgecolor='White',linewidths=4,
                        cmap=cmap,norm=norm)
#        pcm = sns.heatmap(score_df, cmap=cmap, linewidths=0.5, annot=True)
        plt.gca().invert_yaxis()
        plt.autoscale(enable=True,axis='y',tight=True)
        ax.set_yticks(np.arange(score_df.shape[0]) + 0.5, minor=False)
        idx = (['{}'.format(i[1]) for i in score_df.index])
#        idx = score_df.index.get_level_values(0).astype(str).values + ' ' + score_df.index.get_level_values(1).astype(str).values 
        ax.set_yticklabels(idx, minor=False)
        ax.set_xticks(np.arange(score_df.shape[1]) + 0.5, minor=False)
        ax.set_xticklabels(score_df.columns, minor=False,rotation=90)

        plt.title(titlestr)

        for spine in plt.gca().spines.values():
            spine.set_visible(False)
    

        plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on', labelbottom='on')

        cbar = fig.colorbar(pcm, orientation='horizontal', ticks=[-2.5,-1.5,-0.5,0.5,1.5,2.5,-2,2])
        cbar.ax.set_xticklabels(['99%','90%','5%','5%','90%','99%','\nDegraded','\nImproved'])#,rotation=90)
        cbar.set_label('Confidence in Difference')
#        cbar.ax.yaxis.set_label_position('left')
        cbar.ax.tick_params(axis='both', which='both',length=0)
        plt.tight_layout()

        if (save): 
            plt.savefig(filename)
        if (show): 
            plt.show()
        return

#    def load_config(self,cfg_file):
#        import yaml
#
#        with open(cfg_file, 'r') as f:
#            self.cfg_all = yaml.load(f)
#            self.config=self.cfg_all['config']
#            self.class_type = self.cfg_all['class_type']
#            self.classes = self.cfg_all['obs_classes']

#    def scorecard_from_config(self, cfg_file):
#        self.load_config(cfg_file)
#
#        f517 = sc.scorecard(self.config['exp']['name'], self.config['ctl']['name'],
#                            exp_database=self.config['exp']['database'],
#                            ctl_database=self.config['ctl']['database'],
#                            startdate=self.config['startdate'],
#                            enddate=self.config['enddate']) 

