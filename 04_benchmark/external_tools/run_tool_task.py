#!/usr/bin/env python3
"""Run one frozen null replicate/design/tool task without touching legacy results."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, time, traceback
import numpy as np
import pandas as pd

ROOT=Path('${ADSP_ROOT}')
HERE=Path(__file__).resolve().parent
PYTHON=Path('${CONDA_ROOT}/envs/graph/bin/python')
REGENIE='${HOME_DIR}/bin/regenie'
PLINK='${HOME_DIR}/bin/plink2'
ARMS=['retained_additive','retained_joint','full_joint_diagnostic']
TOOLS=['regenie','gmmat','genesis']

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def write_json(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.pending.json');tmp.write_text(json.dumps(obj,indent=2,default=str));os.replace(tmp,path)

def run(cmd,log):
    print('RUN',' '.join(map(str,cmd)),flush=True)
    with open(log,'w') as f:
        rc=subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT).returncode
    if rc:raise RuntimeError(f'Command failed with exit {rc}. See {log}')

def run_regenie(work,prep,bim,threads):
    d=pd.read_csv(str(prep)+'_phenocovar.tsv',sep='\t',dtype={'IID':str})
    fam=pd.read_csv(ROOT/'outputs/plasmode/pool_281k.fam',sep=r'\s+',header=None,dtype=str)
    byid=fam.set_index(1)
    assert not byid.index.duplicated().any() and d.IID.isin(byid.index).all()
    ph=pd.DataFrame({'FID':byid.loc[d.IID,0].to_numpy(),'IID':d.IID,'Y':d.Y})
    ph.to_csv(work/'pheno.tsv',sep='\t',index=False)
    ph[['FID','IID']].to_csv(work/'keep.txt',sep='\t',index=False,header=False)
    cv=pd.concat([ph[['FID','IID']],d.drop(columns=['IID','Y'])],axis=1)
    assert cv.iloc[:,2:].notna().all().all()
    cv.to_csv(work/'covar.tsv',sep='\t',index=False)
    bfile=ROOT/'outputs/plasmode/pool_281k'
    # Derive the step-1 marker list on the actual retained participants.
    # Both a sufficient MAC and positive genotype variance are required.
    run([PLINK,'--bfile',bfile,'--keep',work/'keep.txt','--geno-counts',
         'cols=homref,refalt1,homalt1','--threads',threads,'--memory','16000',
         '--out',work/'counts'],work/'counts.stdout.log')
    gc=pd.read_csv(work/'counts.gcount',sep=r'\s+',dtype={'ID':str})
    assert gc.shape[1]==4 and not gc.ID.duplicated().any(),gc.columns.tolist()
    g=gc.iloc[:,1:].to_numpy(float)
    n=g.sum(axis=1); ac=g[:,1]+2*g[:,2]
    mac=np.minimum(ac,2*n-ac)
    varsum=g[:,1]+4*g[:,2]-np.divide(ac*ac,n,out=np.zeros_like(ac),where=n>0)
    use=(mac>=5)&(varsum>1e-8)&(n>0)
    gc.loc[use,'ID'].to_csv(work/'step1_extract.txt',index=False,header=False)
    assert use.sum()>1000
    common=['--bed',bfile,'--keep',work/'keep.txt','--phenoFile',work/'pheno.tsv',
            '--covarFile',work/'covar.tsv','--covarColList',','.join(cv.columns[2:]),
            '--bt','--threads',threads]
    s1=work/'step1';s2=work/'step2'
    t1=time.time()
    run([REGENIE,'--step','1']+common+['--extract',work/'step1_extract.txt','--bsize','1000',
         '--loocv','--lowmem','--lowmem-prefix',work/'l0','--out',s1],work/'step1.stdout.log')
    step1_s=time.time()-t1
    pred=Path(str(s1)+'_pred.list');assert pred.exists()
    # Validate the numerical prediction files before any association scan.
    pl=pd.read_csv(pred,sep=r'\s+',header=None,dtype=str)
    pred_checks=[]
    for predfile in pl.iloc[:,-1]:
        predpath=Path(predfile)
        if not predpath.is_absolute():predpath=work/predpath
        x=pd.read_csv(predpath,sep=r'\s+')
        assert x.shape[1]==len(d)+1,(str(predpath),x.shape,len(d))
        assert np.isfinite(x.iloc[:,1:].to_numpy(float)).all()
        expected={f'{fid}_{iid}' for fid,iid in zip(ph.FID,ph.IID)}
        assert set(x.columns[1:])==expected,'LOCO prediction participant IDs differ'
        pred_checks.append({'file':str(predpath),'shape':list(x.shape),'sha256':sha(predpath)})
    t2=time.time()
    run([REGENIE,'--step','2']+common+['--bsize','400','--firth','--approx','--pThresh','0.05',
         '--pred',pred,'--out',s2],work/'step2.stdout.log')
    step2_s=time.time()-t2
    rr=pd.read_csv(work/'step2_Y.regenie',sep=r'\s+',dtype={'ID':str})
    if 'TEST' in rr:rr=rr[rr.TEST.isin(['ADD','ADD-FIRTH','ADD-SPA'])]
    assert not rr.ID.duplicated().any()
    unknown=set(rr.ID)-set(bim.id);assert not unknown
    rr=rr.set_index('ID').reindex(bim.id)
    logp=pd.to_numeric(rr.LOG10P,errors='coerce').to_numpy(float)
    assert not np.any(np.isfinite(logp)&(logp<0))
    p=np.power(10.,-logp)
    logs='\n'.join((work/f).read_text() for f in ['step1.stdout.log','step2.stdout.log'])
    warnings=[x for x in logs.splitlines() if any(t in x.lower() for t in ['warn','fail','converg'])]
    # Overall process completion is necessary; site-level failures remain
    # explicit through missing p values and the common-panel coverage gate.
    meta=dict(accepted=True,converged=True,fallback_glm=False,step1_seconds=step1_s,step2_seconds=step2_s,
              n_step1_markers=int(use.sum()),predictions=pred_checks,log_diagnostics=warnings,
              procedure='LOCO ridge step 1; approximate Firth correction at normal P <= 0.05',
              null_fit_validation='Shared finite fixed-effect preflight, completed native solver, finite identity-checked LOCO predictions')
    return p,{'log10p':logp},meta

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--replicate',type=int,required=True)
    ap.add_argument('--arm',choices=ARMS,required=True)
    ap.add_argument('--tool',choices=TOOLS,required=True)
    ap.add_argument('--threads',type=int,default=8)
    ap.add_argument('--preflight-only',action='store_true')
    a=ap.parse_args();assert 1<=a.replicate<=200
    key=f'{a.tool}_{a.arm}_rep{a.replicate:03d}'
    status=HERE/'status'/f'{key}.json'
    if status.exists():
        old=json.loads(status.read_text())
        if old.get('status') in ['accepted','structurally_separated']:
            if old['status']=='accepted':assert sha(old['result'])==old['result_sha256']
            print('Already finished',key,old['status']);return 0
        raise RuntimeError(f'Prior attempt requires audit before rerun: {status}')
    stamp=time.strftime('%Y%m%dT%H%M%S')+'_'+str(os.getpid())
    work=HERE/'work'/key/stamp;work.mkdir(parents=True,exist_ok=False)
    meta=dict(key=key,tool=a.tool,arm=a.arm,replicate=a.replicate,status='running',accepted=False,
              attempt_dir=str(work),slurm_job_id=os.getenv('SLURM_JOB_ID'),slurm_array_job_id=os.getenv('SLURM_ARRAY_JOB_ID'),
              slurm_array_task_id=os.getenv('SLURM_ARRAY_TASK_ID'),start=time.time(),
              driver_sha256={f:sha(HERE/f) for f in ['run_tool_task.py','prepare_design_for_tools.R','fit_tools_replicates.R']})
    write_json(status,meta)
    try:
        hits=list((ROOT/'outputs/plasmode/power/bignull/npz').glob(f'rep{a.replicate:03d}_seed*.npz'))
        assert len(hits)==1
        source=hits[0]
        pool=pd.read_csv(ROOT/'outputs/plasmode/plasmode_pool.csv',usecols=['IID'],dtype=str)
        bim=pd.read_csv(ROOT/'outputs/plasmode/pool_281k.bim',sep=r'\s+',header=None,
                        names=['chr','id','cm','pos','a1','a2'],dtype={'id':str})
        with np.load(source,allow_pickle=False) as z:
            y=z['ystar'].astype(np.int8);seed=int(z['seed']);native_rows=z['rows__strat_project'].astype(int)
            assert len(y)==len(pool)==27306 and np.all(z['keep'])
            excluded=z['excluded'].astype(bool);chrom=z['chr'].astype(int);pos=z['pos'].astype(int)
        assert np.array_equal(chrom,bim.chr.to_numpy(int)) and np.array_equal(pos,bim.pos.to_numpy(int))
        yc=work/'phenotype.csv';pd.DataFrame({'IID':pool.IID,'Y':y}).to_csv(yc,index=False)
        meta.update(native_source=str(source),native_sha256=sha(source),seed=seed,phenotype_sha256=sha(yc))
        prep=work/'design'
        run(['Rscript',HERE/'prepare_design_for_tools.R','--arm',a.arm,'--ystar',yc,'--out',prep],work/'design.stdout.log')
        dm=json.loads(Path(str(prep)+'_meta.json').read_text());meta['design']=dm
        sel=pd.read_csv(str(prep)+'_selection.tsv',sep='\t',dtype={'IID':str})
        assert np.array_equal(sel.IID,pool.IID)
        keep=sel.keep.to_numpy(bool)
        if a.arm!='full_joint_diagnostic':assert np.array_equal(np.flatnonzero(keep),native_rows)
        meta['n_retained']=int(keep.sum());meta['selection_sha256']=sha(str(prep)+'_selection.tsv')
        if dm['status']=='structurally_separated':
            meta.update(status='structurally_separated',reason=dm['reason']);return 0
        if not dm['accepted_design']:raise RuntimeError('Shared fixed-effects design failed preflight')
        if a.preflight_only:
            meta.update(status='preflight_passed');return 0
        extra={}
        if a.tool=='regenie':
            values,extra,tm=run_regenie(work,prep,bim,a.threads)
        else:
            prefix=work/'fit'
            run(['Rscript',HERE/'fit_tools_replicates.R','--ystar',yc,'--design',str(prep)+'.rds',
                 '--out',prefix,'--which',a.tool,'--threads',a.threads,'--save-null','TRUE'],work/'fit.stdout.log')
            fm=json.loads(Path(str(prefix)+'_meta.json').read_text());meta['fit_metadata']=fm
            tm=fm.get(a.tool,{})
            if not tm.get('converged') or tm.get('fallback_glm') or not tm.get('accepted'):
                raise RuntimeError('Native original-model null fit failed acceptance checks')
            pv=pd.read_csv(str(prefix)+'_pvals.tsv.gz',sep='\t',dtype={'snp':str})
            assert np.array_equal(pv.snp.to_numpy(),bim.id.to_numpy())
            values=pv['p_genesis_spa' if a.tool=='genesis' else 'p_gmmat'].to_numpy(float)
            if a.tool=='genesis':extra['p_normal']=pv.p_genesis_score.to_numpy(float)
        meta['tool_validation']=tm
        tested=np.isfinite(values)&(values>=0)&(values<=1)
        hardy=pd.read_csv(ROOT/'revision/tables/plasmode_power/hardy/all.hardy',sep='\t')
        bad=hardy[(hardy['O(HET_A1)']>hardy['E(HET_A1)'])&(hardy.P<1e-4)]
        badset=set(zip(bad['#CHROM'].astype(int),bad.ID.str.split('_').str[1].astype(int)))
        evaluation=~excluded & ~np.array([(c,p) in badset for c,p in zip(chrom,pos)])
        meta.update(n_variants=len(values),n_tested=int(tested.sum()),n_evaluation=int(evaluation.sum()),
                    n_evaluation_tested=int((evaluation&tested).sum()))
        assert (evaluation&tested).sum()>=0.99*evaluation.sum(),'Less than 99% valid coverage on fixed evaluation panel'
        dest=HERE/'results'/f'{key}.npz';dest.parent.mkdir(parents=True,exist_ok=True)
        assert not dest.exists()
        tmp=dest.with_suffix('.pending.npz')
        np.savez_compressed(tmp,p=values,tested=tested,evaluation=evaluation,excluded=excluded,chr=chrom,pos=pos,
                            snp=bim.id.to_numpy(str),ystar=y,keep=keep,seed=np.uint64(seed),replicate=a.replicate,
                            arm=np.array(a.arm),tool=np.array(a.tool),**extra)
        os.replace(tmp,dest)
        meta.update(status='accepted',accepted=True,result=str(dest),result_sha256=sha(dest))
        return 0
    except Exception as e:
        meta.update(status='failed',error=str(e),traceback=traceback.format_exc())
        print(meta['traceback'],flush=True);return 1
    finally:
        meta['elapsed_seconds']=time.time()-meta['start'];meta['finished']=time.time()
        write_json(status,meta);write_json(work/'final_status.json',meta)
        print('STATUS',key,meta['status'],flush=True)

if __name__=='__main__':sys.exit(main())
