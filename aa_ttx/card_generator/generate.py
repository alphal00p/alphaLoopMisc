import os

replace_dict = {}
replace_dict['batch_size'] = 50_000
replace_dict['n_start'] = 100_000_000
replace_dict['n_max'] = 500_000_000
replace_dict['n_workers'] = 500
#replace_dict['options'] = '--fresh '
replace_dict['options'] = ''

output_path = './run_scan_nnlo_all_singlets.gL'

with open('./header.gL','r') as f:
    header = f.read()
with open('./integrate.gL','r') as f:
    integrate = f.read()

#sqrt_s_values = [346.2,347.0,348.0,349.0,350.0,352.0,354.0,360.0,370.0,380.0,400.0,420.0,450.0,470.0,500.0,550.0,600.0,700.0,800.0,900.0,1000.0,1100.0,1200.0,1500.0,1700.0,2000.0,2500.0,3000.0,4000.0,5000.0,6000.0,8000.0,10000.0,12000.0,15000.0,20000.0,100000.0]

sqrt_s_values = [ 346.05, 346.1, 346.6, 347.6, 348.6, 349.6, 350.4, 351.0, 353.0, 357.0, 362.0, 376.0, 390.0, 410.0, 430.0, 460.0, 490.0, 540.0, 590.0, 650.0, 750.0, 850.0, 950.0, 1050.0, 1150.0, 1400.0, 1600.0, 1900.0, 2300.0, 2800.0, 3400.0, 4400.0, 5400.0, 9000.0, 13000.0, 30000.0, 80000.0 ]

offset = 0

res = [header,]
for i, sqrts in enumerate(sqrt_s_values):
    rep_dict = dict(replace_dict)
    rep_dict['run_id'] = offset+i+1
    suffix = f"{i+1}_sqrts_{str(sqrts).replace('.','_')}"
    rep_dict['run_description'] = f"nnlo_all_singlets_{suffix}"
    rep_dict['sqrts_over_2'] = sqrts/2.
    res.append(integrate.format(**rep_dict))

with open(output_path,'w') as f:
    f.write('\n'.join(res))
