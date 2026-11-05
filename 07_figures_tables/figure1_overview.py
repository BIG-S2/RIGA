#!/usr/bin/env python3
"""Render Figure 1 from the existing analysis pool, with an aligned cohort display.

The original DA-EM BIC panel is embedded unchanged: numerical values for that
particular saved fit were not exported. No curve is reconstructed from pixels.
Run with the gwas conda environment. Writes PDF, SVG and 400-dpi PNG.
"""
from pathlib import Path
import re
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'revision/manuscript/figures'
BLUE, ORANGE, GREY = '#0072B2', '#D55E00', '#7A8793'
INK, LIGHT = '#243746', '#DCE3E8'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                     'axes.labelsize': 8, 'xtick.labelsize': 7,
                     'ytick.labelsize': 7.5, 'text.color': INK,
                     'axes.labelcolor': INK, 'axes.edgecolor': GREY,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'axes.spines.top': False, 'axes.spines.right': False})


def pool():
    p = pd.read_csv(ROOT / 'pheno/adsp_combined.csv').dropna(subset=['AD', 'Sex', 'Age'])
    p['Age'] = pd.to_numeric(p.Age.astype(str).str.replace('+', '', regex=False))
    c = pd.read_csv(ROOT / 'adsp_lpa_clusters_20PC.tsv', sep='\t', usecols=['IID', 'lpa_cluster'])
    p = p.merge(c, on='IID', validate='one_to_one')
    # The original R cut(Age, seq(30,100,3), right=TRUE) ends at 99.
    p = p[(p.Age > 30) & (p.Age <= 99)].copy()
    q = pd.read_csv(ROOT / 'pheno/gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv',
                    usecols=['Sample_name', 'Study', 'SeqCenter']).drop_duplicates('Sample_name')
    p = p.merge(q, left_on='IID', right_on='Sample_name', how='left', validate='one_to_one')
    assert (len(p), int(p.AD.sum()), p.Cohort.nunique(), p.Study.nunique(), p.lpa_cluster.nunique()) == (38129, 12756, 50, 66, 13)
    manual = {'Arizona APOE Cohort Study':'AZ_APOE', 'Cache County Study':'CacheCounty',
              'Case Western Reserve University (CWRU) Autopsy Cohort':'CWRU_Autopsy',
              'Case Western Reserve University (CWRU) Rapid Decline Cohort':'CWRU_RDC',
              'IndyAA - Indianapolis African Americans':'IndyAA', 'Mount Sinai Brain Bank':'MSBB',
              'Mayo Clinic':'Mayo', 'University of Miami':'UM',
              'University of Pittsburgh ADRC- Kamboh':'UPitt_Kamboh',
              'University of Washington Families':'UW_Families', 'Vanderbilt University':'Vanderbilt'}
    def short(x):
        m = re.search(r'\(([^()]+)\)\s*$', x)
        return manual.get(x, m.group(1) if m else x)
    p['cohort_label'] = p.Cohort.map(short)
    return p


def main():
    p = pool()
    cohorts = sorted(p.cohort_label.unique(), key=str.casefold)
    # Retain alphabetical cohort ordering, common to every panel.
    n = len(cohorts)
    y = np.arange(n)
    fig = plt.figure(figsize=(7.2, 8.65), facecolor='white')
    bottom, height = .350, .575
    specs = [(0.197,.092),(.305,.080),(.407,.092),(.520,.108),(.653,.317)]
    axes = [fig.add_axes([x,bottom,w,height]) for x,w in specs]
    for ax in axes:
        ax.set_ylim(n-.5,-.5)
        ax.set_yticks(y)
        ax.tick_params(axis='y',length=0,labelleft=False)
        ax.spines[['left','right','top']].set_visible(False)
        ax.spines['bottom'].set_linewidth(.6)
        ax.tick_params(axis='x',length=2,width=.6,pad=3)
        for yy in range(0,n,2):
            ax.axhspan(yy-.5,yy+.5,color='#F4F6F8',zorder=0,lw=0)
    axes[0].set_yticklabels(cohorts, fontsize=7.2)
    axes[0].tick_params(axis='y',labelleft=True,pad=5)
    names = ['AD status','Sex','Age (years)','Race/ethnicity','Mean genetic PCs']
    for ax,title in zip(axes,names): ax.set_title(title,loc='center',fontweight='bold',fontsize=8,pad=9)
    g = p.groupby('cohort_label')
    ad = g.AD.mean().reindex(cohorts).values
    female = g.Sex.mean().reindex(cohorts).values
    for ax,vals,colors in [(axes[0],ad,[GREY,BLUE]),(axes[1],female,['#B7C1C9','#4C526B'])]:
        ax.barh(y,1-vals,height=.78,color=colors[0],edgecolor='none')
        ax.barh(y,vals,left=1-vals,height=.78,color=colors[1],edgecolor='none')
        ax.set_xlim(0,1);ax.set_xticks([0,1],['0','100%'])
        ax.get_xticklabels()[0].set_ha('left');ax.get_xticklabels()[-1].set_ha('right')
    ag = g.Age.agg(['min','median','max']).reindex(cohorts)
    axes[2].hlines(y,ag['min'],ag['max'],color='#A4B5C1',linewidth=.7)
    axes[2].scatter(ag['median'],y,s=7,c=BLUE,zorder=3)
    axes[2].set_xlim(28,102); axes[2].set_xticks([40,60,80,100],['40','60','80','100'],fontsize=6.3)
    race_colors = ['#D55E00','#E6AB02','#5A7D3A','#009E73','#A9B6C1','#8064A2']
    race_labels = ['American Indian / Alaska Native','Asian','Native Hawaiian / Pacific Islander',
                   'Black or African American','White','Hispanic / Latino']
    ra = pd.crosstab(p.cohort_label,p.Race,normalize='index').reindex(index=cohorts,columns=range(1,7),fill_value=0)
    left=np.zeros(n)
    for col,color in zip(ra.columns,race_colors):
        vals=ra[col].values
        axes[3].barh(y,vals,left=left,height=.78,color=color,edgecolor='none')
        left+=vals
    axes[3].set_xlim(0,1);axes[3].set_xticks([0,1],['0','100%'])
    axes[3].get_xticklabels()[0].set_ha('left');axes[3].get_xticklabels()[-1].set_ha('right')
    means=g[[f'PC{i}' for i in range(1,21)]].mean().reindex(cohorts)
    lim=np.abs(means.values).max()
    cmap=LinearSegmentedColormap.from_list('blue_white_orange',[BLUE,'#FFFFFF',ORANGE])
    hm=axes[4].imshow(means.values,aspect='auto',cmap=cmap,vmin=-lim,vmax=lim,interpolation='nearest',zorder=2)
    axes[4].set_xticks([0,4,9,14,19],['1','5','10','15','20'],fontsize=7)
    axes[4].set_xlabel('Principal component',labelpad=4)
    for yy in np.arange(.5,n,1): axes[4].axhline(yy,color='white',lw=.35,zorder=3)
    cax=fig.add_axes([.74,.297,.17,.008])
    cb=fig.colorbar(hm,cax=cax,orientation='horizontal',ticks=[-.04,0,.04])
    cb.ax.tick_params(labelsize=6.5,length=2,pad=2);cb.outline.set_visible(False)
    cb.set_label('Mean PC value',fontsize=7,labelpad=2)
    fig.text(.018,.951,'a',fontsize=12,fontweight='bold')
    fig.text(.653,.951,'b',fontsize=12,fontweight='bold')
    fig.text(.197,.979,'38,129 participants across 50 cohorts',fontsize=9,color=INK)
    handles=[Patch(color=GREY,label='Control'),Patch(color=BLUE,label='AD'),
             Patch(color='#B7C1C9',label='Male'),Patch(color='#4C526B',label='Female')]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.185,.337),ncol=4,
               frameon=False,fontsize=7,handlelength=1,handletextpad=.4,columnspacing=1)
    fig.legend(handles=[Patch(color=c,label=l) for c,l in zip(race_colors,race_labels)],
               loc='upper left',bbox_to_anchor=(.185,.307),ncol=2,frameon=False,
               fontsize=6.8,handlelength=1,handletextpad=.5,columnspacing=1.3)
    # Preserve the complete existing BIC panel and all labels verbatim.
    bic=Image.open(ROOT/'bic_curve_20PC_unrelated.png')
    ac=fig.add_axes([.015,.005,.510,.245]);ac.imshow(bic);ac.axis('off')
    fig.text(.018,.247,'c',fontsize=12,fontweight='bold')
    axd=fig.add_axes([.600,.050,.270,.165])
    pal=['#0072B2','#D55E00','#009E73','#CC79A7','#6F4C9B','#A6761D','#56B4E9',
         '#7A8793','#B79F00','#00A6A6','#3B4B87','#B56576','#E69F00']
    rng=np.random.default_rng(1)
    order=rng.permutation(len(p))
    z=p.iloc[order]
    axd.scatter(z.PC1,z.PC2,c=[pal[int(k)-1] for k in z.lpa_cluster],s=.30,
                alpha=.55,linewidths=0,rasterized=True)
    axd.set_xlabel('PC1',labelpad=2);axd.set_ylabel('PC2',labelpad=2)
    axd.set_title('Latent genetic classes',fontsize=8,fontweight='bold',pad=5)
    axd.tick_params(labelsize=6.5,length=2,pad=2)
    axd.grid(color=LIGHT,lw=.4);axd.set_axisbelow(True)
    axd.legend(handles=[Line2D([],[],marker='o',ls='',ms=3,color=c,label=str(i)) for i,c in enumerate(pal,1)],
               loc='center left',bbox_to_anchor=(1.01,.5),ncol=2,frameon=False,
               fontsize=6,handlelength=.6,handletextpad=.3,columnspacing=.5,labelspacing=.3,title='Class',title_fontsize=6.5)
    fig.text(.600,.247,'d',fontsize=12,fontweight='bold')
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'fig1_ng_20261014.{ext}',dpi=400,facecolor='white')
    plt.close(fig)
    audit={'participants':len(p),'cases':int(p.AD.sum()),'cohorts':n,'projects':int(p.Study.nunique()),
           'classes':int(p.lpa_cluster.nunique()),'cohort_order':cohorts,
           'sex_codes':'0=Male, 1=Female; ADSP phenotype dictionaries',
           'bic_panel':'Original unmodified bic_curve_20PC_unrelated.png embedded; no numerical reconstruction',
           'pca_points':len(p),'pca_rasterization':'All participants retained; dense scatter only rasterized'}
    (OUT/'fig1_ng_20261014_provenance.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:audit[k] for k in ['participants','cases','cohorts','projects','classes']}))


if __name__=='__main__': main()
