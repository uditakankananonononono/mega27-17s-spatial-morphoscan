import json, numpy as np
from math import comb
d=json.load(open('results/p2_multitask.json'))['folds']
rng=np.random.default_rng(20261008)
out={"prespec":"results/p5_pls_robust_prespec.json","n_folds":len(d)}
ok=True
for name,key in [("committed","ridge_m1_committed"),("refold","ridge_m1_refold")]:
    x=np.array([f['pls_median_r']-f[key] for f in d]); w=int((x>0).sum()); n=len(x)
    p=sum(comb(n,k) for k in range(w,n+1))/2**n
    bs=np.array([np.median(rng.choice(x,n,replace=True)) for _ in range(10000)])
    lo,hi=np.percentile(bs,[2.5,97.5])
    out[name]={"wins":w,"median_delta":float(np.median(x)),"sign_p_one_sided":p,"boot_ci95":[float(lo),float(hi)]}
    ok&= (w>=9 and lo>0)
out["verdict"]="SURVIVES" if ok else "FAILS"
json.dump(out,open('results/p5_pls_robust.json','w'),indent=1); print(json.dumps(out,indent=1))
