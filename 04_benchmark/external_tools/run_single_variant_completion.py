#!/usr/bin/env python3
"""Immutable same-model SV completion. Failed/null-nonconverged fits never publish NPZs."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, time
import numpy as np
import pandas as pd

ROOT=Path('${ADSP_ROOT}')
HERE=ROOT/'revision/sv_completion_20261015'
sys.path.insert(0,str(ROOT/'revision/scripts'))
import benchmark_null_harness as H

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--scheme',choices=['null','power'],required=True)
    p.add_argument('--replicate',type=int,required=True)
    p.add_argument('--tool',choices=['genesis','gmmat'],default='genesis')
    p.add_argument('--threads',type=int,default=8)
    p.add_argument('--save-null',action='store_true')
    a=p.parse_args()
    target=200 if a.scheme=='null' else 40
    assert 1<=a.replicate<=target
    native=ROOT/'outputs/plasmode/power'/('bignull' if a.scheme=='null' else 'bigpower')
    hits=sorted((native/'npz').glob(f'rep{a.replicate:03d}_seed*.npz'))
    assert len(hits)==1, hits
    native_file=hits[0]
    out=HERE/'results'/a.scheme
    covset=a.tool+'_project'
    dest=out/'npz'/f'regenie_{covset}_rep{a.replicate:03d}.npz'
    meta_dest=out/'meta'/f'{a.tool}_rep{a.replicate:03d}.json'
    if dest.exists() or meta_dest.exists():
        assert dest.exists() and meta_dest.exists(), 'Incomplete prior publication requires audit'
        m=json.loads(meta_dest.read_text())
        assert m['complete'] and m['npz_sha256']==sha(dest)
        print('Validated output already exists',dest)
        return 0
    attempt=time.strftime('%Y%m%dT%H%M%S')+'_'+str(os.getpid())
    work=HERE/'work'/a.scheme/f'{a.tool}_rep{a.replicate:03d}'/attempt
    work.mkdir(parents=True,exist_ok=False)
    with np.load(native_file,allow_pickle=False) as z:
        pool=pd.read_csv(ROOT/'outputs/plasmode/plasmode_pool.csv',usecols=['IID'],dtype={'IID':str})
        n=len(pool)
        y=z['ystar'].astype(np.int8)
        keep=z['keep'].astype(bool) if 'keep' in z.files else np.ones(n,dtype=bool)
        if len(y)==int(keep.sum()) and len(y)!=n:
            full=np.zeros(n,dtype=np.int8); full[keep]=y; y=full
        assert len(y)==n and len(keep)==n
        seed=int(z['seed'])
    ycsv=work/'phenotype.csv'
    pd.DataFrame({'IID':pool.IID,'Y':np.where(keep,y.astype(float),np.nan)}).to_csv(ycsv,index=False,na_rep='NA')
    prefix=work/'fit'
    cmd=['Rscript',str(HERE/'fit_single_variant_validated.R'),'--ystar',str(ycsv),'--out',str(prefix),'--which',a.tool,'--threads',str(a.threads)]
    if a.save_null: cmd+=['--save-null','TRUE']
    with open(work/'R.stdout.log','w') as f:
        rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
    if rc:
        print('R process failed',rc,'see',work/'R.stdout.log'); return rc
    meta=json.loads(Path(str(prefix)+'_meta.json').read_text())
    tmeta=meta.get(a.tool,{})
    if not tmeta.get('converged',False) or tmeta.get('fallback_glm',False):
        print('No valid original-model fit. Not publishing',work); return 3
    pv=pd.read_csv(str(prefix)+'_pvals.tsv.gz',sep='\t',dtype={'snp':str})
    bim=pd.read_csv(ROOT/'outputs/plasmode/pool_281k.bim',sep=r'\s+',header=None,
                   names=['chr','id','cm','pos','a1','a2'],dtype={'id':str})
    assert np.array_equal(pv.snp.to_numpy(),bim.id.to_numpy())
    values=pv['p_genesis_spa' if a.tool=='genesis' else 'p_gmmat'].to_numpy(float)
    tested=np.isfinite(values)&(values>=0)&(values<=1)
    if tested.sum()<0.99*len(values):
        print('Fewer than99%of benchmark sites have valid test outputs. Not publishing',work);return 4
    ann=H.snp_annotation(bim,H.PERSNP)
    for dd in [dest.parent,meta_dest.parent]:dd.mkdir(parents=True,exist_ok=True)
    tmp=dest.with_suffix('.pending.npz')
    np.savez_compressed(tmp,p=values.astype(np.float32),tested=tested,excluded=ann.excluded.to_numpy(),
                        covset=np.array([covset]),tool=np.array([a.tool]),seed=np.uint64(seed),
                        replicate=np.int64(a.replicate),ystar=y,keep=keep,
                        chr=ann.chr.to_numpy().astype(np.int16),pos=ann.pos.to_numpy().astype(np.int64))
    m=dict(complete=True,scheme=a.scheme,replicate=a.replicate,tool=a.tool,seed=seed,
           native_npz=str(native_file),native_sha256=sha(native_file),phenotype_sha256=sha(ycsv),
           driver_sha256=sha(HERE/'fit_single_variant_validated.R'),raw_metadata=str(prefix)+'_meta.json',
           npz_sha256=sha(tmp),n_tested=int(tested.sum()),r_metadata=meta)
    os.replace(tmp,dest)
    mt=meta_dest.with_suffix('.pending.json');mt.write_text(json.dumps(m,indent=2));os.replace(mt,meta_dest)
    print('Published converged original-model output',dest)
    return 0
if __name__=='__main__':sys.exit(main())
