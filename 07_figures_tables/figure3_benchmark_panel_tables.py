#!/usr/bin/env python
"""Build figure and table payloads from audited summaries. No manuscript or workbook edits."""
from pathlib import Path
import sys, json, hashlib, shutil, os
import numpy as np
import pandas as pd
import openpyxl
HERE=Path(__file__).resolve().parent
REV=HERE.parents[2]
SCRIPTS=REV/'scripts'
TABLES=REV/'tables/plasmode_power'
EXT=REV/'five_analysis_execution_20261016/external_tools'
WORK=REV/'manuscript/Supplementary_Tables.xlsx'
BACKUP=HERE/'Supplementary_Tables.before_20261027.xlsx'
if not BACKUP.exists(): shutil.copy2(WORK,BACKUP)
w=openpyxl.load_workbook(BACKUP)
base=[[c.value for c in row] for row in w['ST15a'].iter_rows(min_row=5)]
cols=[c.value for c in w['ST15a'][4]]
old=pd.DataFrame(base,columns=cols)
old.to_csv(HERE/'ST15a_original_rows.csv',index=False)
newcols=cols+['Participants','Replicates','Fit and analysis status']
rows=[]
for vals in base:
 rows.append(vals+[None,1000,'Primary benchmark, original normal or saddlepoint statistic'])
GM=EXT/'retained_1000_20261019/audit/20261020T083739_final_5728839/summary_by_arm.tsv'
RG=EXT/'regenie_expansion_20261019/audit/run_20261020T083739/counts_summary.csv'
gm=pd.read_csv(GM,sep='\t'); rg=pd.read_csv(RG)
external=[]
for r in gm[gm['set']=='primary'].to_dict('records'):
 external.append(dict(tool=r['tool'],arm=r['arm'],tier='primary',alpha=r['alpha'],n_replicates=r['n_replicates'],total_tests=r['total_tests'],observed=r['observed_raw'],expected=r['expected_raw'],ratio=r['ratio_raw'],lo=r['ratio_raw_boot_lo'],hi=r['ratio_raw_boot_hi'],N=20786,status='Accepted converged mixed model'))
for r in rg[rg.tool=='regenie'].to_dict('records'):
 if r['tier_set'] not in ('clean','accepted'): continue
 if r['tier_set']=='accepted' and r['arm']=='retained_legacy_additive': continue
 external.append(dict(tool='regenie',arm=r['arm'],tier=r['tier_set'],alpha=r['alpha'],n_replicates=r['n_replicates'],total_tests=r['total_tests'],observed=r['observed_raw'],expected=r['expected'],ratio=r['ratio_raw'],lo=r['ratio_raw_ci_lo'],hi=r['ratio_raw_ci_hi'],N=27306 if r['arm'].startswith('full') else 20786,status='Accepted without transient null retry' if r['tier_set']=='clean' else 'Sensitivity including one transient null retry'))
ext=pd.DataFrame(external)
ext.to_csv(HERE/'external_calibration_source.csv',index=False)
arm_names={'retained_joint':'Retained, joint project × center','retained_legacy_additive':'Retained, historical additive coding','full_joint':'Complete pool, joint project × center','full_joint_regenie':'Complete pool, joint project × center'}
for r in external:
 design=arm_names[r['arm']]
 test={'gmmat':'GMMAT normal score','genesis':'GENESIS Score.SPA','regenie':'REGENIE LOCO and approximate Firth'}[r['tool']]
 rows.append(['Heterozygote-excess filter and AD-window exclusions',design,test,r['alpha'],r['total_tests'],r['observed'],r['expected'],r['ratio'],r['lo'],r['hi'],r['N'],r['n_replicates'],r['status']])
for tool in ('GMMAT','GENESIS'):
 for arm in ('joint project × center','historical additive coding'):
  rows.append(['Same prespecified evaluation panel','Complete pool, '+arm,tool,None,0,None,None,None,None,None,27306,0,'No accepted fit in 3 attempted pilots. No false-positive ratio estimated.'])
new=pd.DataFrame(rows,columns=newcols)
new.to_csv(HERE/'ST15a_updated.csv',index=False)
rvsrc=TABLES/'fpr_rv_native_tests_1000.csv'
rv=pd.read_csv(rvsrc)
rv=rv[np.isclose(rv.alpha,2.5e-6,rtol=1e-12,atol=0)].copy()
assert len(rv)==24
name={'saturated_cell':'Oracle (generating cells)','strat_project':'Stratified (project x center)','pooled_project':'Pooled + center + project + PCs'}
tname={'SB':'SKAT–burden omnibus test','STAAR-O':'STAAR-O','burden':'burden','SKAT':'SKAT'}
rvrows=[]
for floor in (10,0):
 for design in name:
  for test in tname:
   z=rv[(rv.cmac_floor==floor)&(rv.design==design)&(rv.test==test)]
   assert len(z)==1
   r=z.iloc[0]
   rvrows.append(['cMAC >= 10' if floor else 'none',name[design],tname[test],int(r.n_rep),r.n_tests_mean,int(r.observed),r.ratio,r.boot_lo,r.boot_hi,r.share_any_fdr,r.share_lo,r.share_hi])
rvtable=pd.DataFrame(rvrows,columns=[c.value for c in w['ST15b'][4]])
rvtable.to_csv(HERE/'ST15b_updated.csv',index=False)
cap_a=('Observed null-test counts, expected counts and their ratios at the stated threshold. Primary benchmark rows retain the original 1,000-permutation normal and saddlepoint summaries, with and without removal of heterozygote-excess sites. The filtered external-tool evaluation panel contains 201,603 variants after AD-region exclusions. New GMMAT and GENESIS analyses use 1,000 shared null phenotypes on 20,786 retained participants. All 4,000 mixed-model fits across the two tools and two retained specifications converged with the sparse genetic relatedness matrix. REGENIE uses 200 phenotypes per specification. Its stringent summary excludes one retained joint-unit scan and one complete-pool joint-unit scan with a transient null-model retry. The all-accepted sensitivity rows include these scans. Additive project coding was defined on the complete pool before retention. Confidence intervals resample phenotype replicates with 2,000 resamples for the primary and GMMAT/GENESIS summaries and 4,000 for REGENIE. Complete-pool GMMAT and GENESIS rows report failed fitting diagnostics and have no calibration estimate. Historical failed or fallback external-tool results are not used as valid mixed-model comparisons.')
cap_b=('Null coding-gene scans under the oracle, stratified and pooled designs. The stratified design has 1,000 scans, pooled has 100 and oracle has 50. False-positive ratios at P < 2.5 × 10−6 have 95% percentile-bootstrap intervals over scans. The proportion of complete-null scans with any Benjamini–Hochberg discovery at FDR 5% estimates the realized FDR and has an exact Clopper–Pearson interval. The SKAT–burden omnibus test combines the SKAT and burden components and excludes ACAT-V. With cMAC ≥ 10, 51 of 1,000 stratified scans selected at least one gene, giving 5.1% with a 95% interval of 3.8% to 6.7%. The allele-count floor uses the frozen full-pool gene annotation. These simulations evaluate the benchmark statistic on unrelated participants and do not reproduce every feature of the deployed weighted mixed model or its annotation-based coding and noncoding families.')
payload={'ST15a':{'title':w['ST15a']['A1'].value,'caption':cap_a,'columns':newcols,'rows':rows},'ST15b':{'title':w['ST15b']['A1'].value,'caption':cap_b,'columns':list(rvtable.columns),'rows':rvrows}}
(HERE/'table15_payload.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False))

# Keep established visual style while separating new calibration from historical power.
os.environ.setdefault('MPLCONFIGDIR','/tmp/riga_integration_20261027_mpl')
sys.path.insert(0,str(SCRIPTS))
import fig_plasmode_benchmark_ng_20261014 as F
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as mticker
F.original.REV=REV;F.original.TABP=TABLES;F.original.LEVEL=1;F.original.ALPHA=5e-8
F.original.test_of=lambda arm:'firth' if arm.startswith('regenie') else 'chi2'
single=F.original.load('big',TABLES/'fpr_bignull_figure_table_hetfilt_chi2.csv')
keys=['oracle','strat_project','riga_projectkey','pooled_project']
single=single[single.key.isin(keys)].copy()
# Replace independent-test intervals with the published whole-replicate intervals.
labels={'oracle':'Oracle (generating cells)','strat_project':'Stratified (project x center)','riga_projectkey':'Matched 1:1 (project x center)','pooled_project':'Pooled + center + project + PCs'}
for idx,r in single.iterrows():
 match=old[(old['Panel filter'].str.startswith('heterozygote'))&(old.Design==labels[r.key])&(old.Test=='normal approximation')&np.isclose(old['α'],5e-8,atol=0,rtol=1e-12)]
 assert len(match)==1
 z=match.iloc[0]
 assert abs(float(z.Ratio)-r.ratio)<.001
 single.loc[idx,'lo']=z['Bootstrap lower'];single.loc[idx,'hi']=z['Bootstrap upper']
single.to_csv(HERE/'figure_panel_a_source.csv',index=False)
rv10=rv[rv.cmac_floor==10].copy();rv10.to_csv(HERE/'figure_panel_b_source.csv',index=False)
plotext=ext[(ext.tier.isin(['primary','clean']))&ext.arm.str.startswith('retained')].copy()
assert len(plotext)==18
plotext.to_csv(HERE/'figure_panel_c_source.csv',index=False)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],'font.size':8,'text.color':F.INK,'axes.labelcolor':F.INK,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False,'mathtext.fontset':'dejavusans','savefig.facecolor':'white'})
W,H=183,180
fig=plt.figure(figsize=(W/25.4,H/25.4),dpi=300,facecolor='white')
def axmm(x,y,w,h):return fig.add_axes([x/W,y/H,w/W,h/H])
def txt(x,y,s,**kw):return fig.text(x/W,y/H,s,**kw)
a= axmm(17,101,66,55)
F.base_style(a);a.set_yscale('log');a.set_ylim(.25,620);a.set_xlim(.1,.65)
a.set_xticks(np.arange(.1,.7,.1));a.set_yticks([.3,1,3,10,30,100,300]);a.set_yticklabels(['0.3','1','3','10','30','100','300']);a.yaxis.set_minor_locator(mticker.NullLocator());a.grid(axis='y',color=F.GRID,lw=.5);a.axhline(1,color='#888888',lw=.8,ls=(0,(3,2)))
for r in single.itertuples():
 color=F.STYLES[r.key][1]
 a.plot([r.power]*2,[r.lo,r.hi],color=color,lw=.95)
 a.plot([r.power-1.96*r.power_se,r.power+1.96*r.power_se],[r.ratio]*2,color=color,lw=.95)
 F.point(a,r.power,r.ratio,r.key)
a.set_xlabel('Power',fontsize=8.5,labelpad=4);a.set_ylabel('Type-I error / nominal',fontsize=8.5,labelpad=4)
a.text(.11,1.23,'Nominal',fontsize=7.5,color=F.SECONDARY)
b=axmm(130,101,49,55);F.panel_b(b,rv10)
for x,label,title in [(3,'a','Single-variant tests'),(95,'b','Rare-variant tests')]:
 txt(x,175,label,fontsize=11,fontweight='bold',va='top')
 txt(17 if label=='a' else 103,175,title,fontsize=9,fontweight='bold',va='top')
txt(17,169,r'Calibration and power at $\alpha = 5\times10^{-8}$',fontsize=7.8,va='top',color=F.SECONDARY)
txt(17,164,'Effect size equivalent to OR 1.20 at MAF 0.30',fontsize=7.1,va='top',color=F.SECONDARY)
txt(103,169,r'Coding sets with cMAC $\geq$ 10',fontsize=7.8,va='top',color=F.SECONDARY)
txt(103,164,r'Calibration at $\alpha = 2.5\times10^{-6}$',fontsize=7.8,va='top',color=F.SECONDARY)
handles=[]
for key in keys:
 label,color,marker,size,hollow=F.STYLES[key]
 if key=='riga_projectkey':label='Case–control matched'
 handles.append(Line2D([],[],ls='none',marker=marker,ms=size,mfc='white' if hollow else color,mec=color if hollow else 'white',color=color,label=label))
fig.legend(handles=handles,ncol=4,loc='center',frameon=False,bbox_to_anchor=(.51,88/H),fontsize=7.6,handletextpad=.4,columnspacing=1.4)
txt(3,78,'c',fontsize=11,fontweight='bold',va='top')
txt(17,78,'Calibration across association engines',fontsize=9,fontweight='bold',va='top')
txt(17,72,'Same 20,786 retained participants and shared null phenotypes',fontsize=7.8,color=F.SECONDARY,va='top')
armsty={'retained_joint':(F.BLUE,'s',False,.12),'retained_legacy_additive':('#009E73','o',True,-.12)}
for ix,(alpha,alabel) in enumerate([(1e-5,r'$\alpha = 10^{-5}$'),(1e-6,r'$\alpha = 10^{-6}$'),(5e-8,r'$\alpha = 5\times10^{-8}$')]):
 ax=axmm(32+52*ix,25,43,38);F.base_style(ax);ax.set_xscale('log');ax.set_xlim(.3,400);ax.set_ylim(-.5,2.5)
 ax.set_xticks([1,10,100]);ax.set_xticklabels(['1','10','100']);ax.xaxis.set_minor_locator(mticker.NullLocator());ax.axvline(1,color='#888888',lw=.8,ls=(0,(3,2)));ax.grid(axis='x',color=F.GRID,lw=.5)
 ax.spines['left'].set_visible(False);ax.set_yticks([2,1,0]);ax.set_yticklabels(['GMMAT','GENESIS','REGENIE'] if ix==0 else ['','','']);ax.tick_params(axis='y',length=0)
 ax.set_title(alabel,fontsize=8,pad=6);ax.set_xlabel('Type-I error / nominal',fontsize=7.6,labelpad=4)
 for ti,tool in enumerate(['gmmat','genesis','regenie']):
  for arm,(color,marker,hollow,dy) in armsty.items():
   z=plotext[(plotext.tool==tool)&(plotext.arm==arm)&np.isclose(plotext.alpha,alpha,atol=0,rtol=1e-12)].iloc[0]
   y=2-ti+dy;ax.plot([z.lo,z.hi],[y,y],color=color,lw=1.05)
   ax.plot(z.ratio,y,ls='none',marker=marker,ms=4.8,mfc='white' if hollow else color,mec=color,mew=.9)
handles=[Line2D([],[],ls='none',marker=m,ms=4.8,mfc='white' if h else c,mec=c,label='Joint-unit intercepts' if ar=='retained_joint' else 'Historical additive coding') for ar,(c,m,h,dy) in armsty.items()]
fig.legend(handles=handles,ncol=2,loc='center',frameon=False,bbox_to_anchor=(.53,12/H),fontsize=7.6,handletextpad=.4,columnspacing=2)
txt(91.5,4,'GMMAT and GENESIS  1,000 scans per arm       REGENIE  199–200 scans per arm',ha='center',fontsize=7,color=F.SECONDARY)
out=REV/'manuscript/figures/fig2_benchmark_updated_20261027'
for suffix in ('.pdf','.svg','.png'):fig.savefig(out.with_suffix(suffix),dpi=300,facecolor='white')
plt.close(fig)
inputs=[BACKUP,GM,RG,rvsrc,TABLES/'fpr_bignull_figure_table_hetfilt_chi2.csv',TABLES/'big/power_snp.csv']
manifest={'figure':str(out),'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'notes':['Panel a retains original power and four primary design point estimates, with replicate-bootstrap null intervals from ST15a.','Panel b updates only stratified data to 1000 scans. Oracle50 and pooled100 unchanged.','Panel c contains only validated retained null summaries and no power estimate.','Full-pool diagnostics and REGENIE sensitivities are in ST15a, not plotted.']}
(HERE/'benchmark_provenance.json').write_text(json.dumps(manifest,indent=2))
print('Wrote',out)
print('ST15a rows',len(rows),'ST15b rows',len(rvrows),'external plotted',len(plotext))
