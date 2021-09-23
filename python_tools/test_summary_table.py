import gmao_tools as gt
import ioda_access as ia

fn = '/discover/nobackup/projects/gmao/obsdev/wrmccart/jcsda/pyioda_test/geos.a005.diag.PT6H.aircraft.2020-12-15T03:00:00Z.PT6H.nc4'

acft_file = ia.obs(fn)
acft_file.set_obs_var('air_temperature')

config_dict = gt.yaml_to_dict('summary_table.example.yaml')

df = gt.summary_table(acft_file,config_dict,in_msk='(kx > 0)')#,verbose=True)

print(df)
