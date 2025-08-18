import os
script_path = os.path.realpath(__file__)
script_dir  = os.path.dirname(script_path)
proj_dir    = os.path.join(script_dir, "..")
os.sys.path.append(proj_dir)

prov_dir = os.path.join(proj_dir, 'results', '01-model-tuning',
                        '03-provenance')

prov_f = os.path.join(prov_dir, 'prov-models.tsv')
base_mdl_dir = os.path.join(proj_dir, 'models', 'deephlaffy', 'both')

import pandas as pd
from nbslpy._utils import readPickle


datf = pd.read_csv(prov_f, sep = '\t')

res_lst = []
for _, x in datf.iterrows():
    mdl      = x['Model']
    mdl_key  = x['ModelKey']
    mdl_seed = x['Seed']
    mdl_dir  = os.path.join(base_mdl_dir, mdl_key, mdl_seed)
    hist_f   = os.path.join(mdl_dir, 'hist.pkl')
    hist     = readPickle(hist_f)
    tmp      = pd.DataFrame(hist.history)
    tmp.loc[:,"Epoch"] = [i + 1 for i in range(len(tmp))]
    tmp.loc[:,"Model"]    = mdl
    tmp.loc[:,"ModelKey"] = mdl_key
    tmp.loc[:,"Seed"]  = mdl_seed
    res_lst.append(tmp)

res = pd.concat(res_lst, axis = 0)
res_f = os.path.join(prov_dir, 'prov-models-data.tsv')
res.to_csv(res_f, sep = '\t', index = False)
