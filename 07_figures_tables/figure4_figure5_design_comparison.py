#!/usr/bin/env python3
"""Complete QC_label7 comparisons in STAARpipelineSummary visual conventions.

No association reruns and no intersection of native eligible universes.
Dense point fields map every observation to a documented raster pixel.
Axes, symbols, highlighted tails and annotations remain editable vectors.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/riga_qc7_mpl')
import json
import hashlib
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
from matplotlib.ticker import MaxNLocator
from scipy.ndimage import maximum_filter
from scipy.special import ndtri

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
OUT = BASE.parent / 'figures'
NEW = ROOT / 'revision/pooled_qc7_update_20261015'
OLD = BASE / 'sv_cache' / 'riga'          # 2026-10-19: persistent copies of the former /tmp caches
CACHE = BASE / 'sv_cache' / 'cmp'
N = {'riga': 30894, 'matched': 12724, 'pooled': 38138}
NAMES = {'riga': 'RIGA', 'matched': 'Matched, project key', 'pooled': 'Pooled'}
COLORS = {'riga': '#000000', 'matched': '#0072B2', 'pooled': '#D55E00'}
TAIL_MARKERS = {'riga': 'o', 'matched': '^', 'pooled': 's'}
NAVY, GOLD, RED, GRAY = '#00008B', '#CD8500', '#FF0000', '#606060'
SEL = '#BDBDBD'   # 2026-10-19: BH-selected associations in light gray (threshold line stays red)
LENGTH = np.array([248956422,242193529,198295559,190214555,181538259,
                  170805979,159345973,145138636,138394717,133797422,135086622,
                  133275309,114364328,107043718,101991189,90338345,83257441,
                  80373285,58617616,64444167,46709983,50818468])
OFFSET = np.r_[0, 0, np.cumsum(LENGTH)[:-1]]
XLIM = (-10, LENGTH.sum()/1e6+10)
TICKS = (OFFSET[1:] + LENGTH/2)/1e6
MASKS = {
 'coding': [('plof','pLoF','s'),('plof_ds','pLoF + disruptive missense','o'),
            ('missense','Missense','^'),('disruptive_missense','Disruptive missense','+'),
            ('synonymous','Synonymous','x')],
 'noncoding': [('UTR','UTR','s'),('upstream','Upstream','o'),('downstream','Downstream','^'),
               ('promoter_CAGE','Promoter CAGE','+'),('promoter_DHS','Promoter DHS','x'),
               ('enhancer_CAGE','Enhancer CAGE','D'),('enhancer_DHS','Enhancer DHS','v')],
 'ncRNA': [('ncRNA','ncRNA','o')]}
TITLES = {'sv': 'Single variants', 'coding': 'Coding variant sets',
          'noncoding': 'Noncoding variant sets', 'ncRNA': 'Noncoding RNA genes'}
ANNOTATIONS = []
ANNOTATION_DATA = []

plt.rcParams.update({'font.family':'Liberation Sans', 'font.size':7,
                     'pdf.fonttype':42, 'ps.fonttype':42, 'svg.fonttype':'none',
                     'axes.linewidth':.6, 'savefig.facecolor':'white'})


def compress(value):
    a = np.asarray(value)
    return np.where(a <= 15, a, 15+(a-15)/16)


def svpath(arm, part):
    return OLD / f'sv_{part}' if arm == 'riga' else CACHE / f'{arm}_sv_{part}'


def raster(ax, x, y, color, xlim, ylim, width=1750, height=1250, radius=1, alpha=1):
    """Map every finite observation without subsampling or capping."""
    assert len(x) == len(y)
    occupied = np.zeros((height, width), dtype=bool)
    clipped = 0
    for start in range(0, len(x), 1000000):
        xx, yy = np.asarray(x[start:start+1000000]), np.asarray(y[start:start+1000000])
        good = (np.isfinite(xx) & np.isfinite(yy) & (xx >= xlim[0]) & (xx <= xlim[1])
                & (yy >= ylim[0]) & (yy <= ylim[1]))
        clipped += int((~good).sum())
        ix = np.rint((xx[good]-xlim[0])/(xlim[1]-xlim[0])*(width-1)).astype(int)
        iy = np.rint((yy[good]-ylim[0])/(ylim[1]-ylim[0])*(height-1)).astype(int)
        occupied[iy, ix] = True
    pixels = int(occupied.sum())
    if radius:
        occupied = maximum_filter(occupied, size=2*radius+1)
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    rgba[:, :, :3] = np.array(to_rgba(color)[:3])*255
    rgba[:, :, 3] = occupied.astype(np.uint8)*round(alpha*255)
    ax.imshow(rgba, origin='lower', extent=(*xlim, *ylim), aspect='auto',
              interpolation='nearest', zorder=1)
    assert clipped == 0, clipped
    return {'input_count': len(x), 'rendered_count': len(x)-clipped,
            'clipped_count': clipped, 'occupied_pixels_before_dilation': pixels,
            'raster_size': [width, height], 'radius_pixels': radius}


def load_gene(family, arm, coords):
    if arm == 'riga':
        src = (BASE/'figure_inputs'/f'{family}.csv' if family != 'ncRNA' else
               ROOT/'revision/tables/review/genesets_primary_ncRNA.csv.gz')
    elif family == 'coding':
        src = BASE/'figure_inputs'/f'{arm}_coding.csv'
    else:
        src = NEW/f'{arm}_{family}.csv'
    d = pd.read_csv(src).rename(columns={'gene':'Gene name'})
    if 'tested' in d:
        d = d[d.tested.astype(str).str.lower().isin(['true','1'])].copy()
    d = d[d.cMAC >= 10].copy()
    assert np.isfinite(d.SB).all() and d.SB.between(0,1,inclusive='right').all()
    assert not d.duplicated(['Gene name','Chr','Category']).any()
    d = d.merge(coords[['hgnc_symbol','chromosome_name','pos']],
                left_on=['Gene name','Chr'], right_on=['hgnc_symbol','chromosome_name'],
                how='left',validate='many_to_one')
    assert d.pos.notna().all(), d[d.pos.isna()][['Gene name','Chr']].to_dict('records')
    d['x'] = (d.pos + OFFSET[d.Chr.to_numpy(dtype=int)])/1e6
    d['lp'] = -np.log10(d.SB)
    ps = np.sort(d.SB.to_numpy()); nn=len(ps)
    ranks=np.arange(1,nn+1); passing=ranks[ps <= .05*ranks/nn]
    k=int(passing[-1]) if len(passing) else 0
    critical=.05*k/nn if k else None
    d['selected']=d.SB <= critical if k else False
    if 'q_SB' in d:
        assert int((d.q_SB < .05).sum()) == k
    stats={'n':nn,'k':k,'critical':critical,'lambda':float((ndtri(np.median(ps)/2)/ndtri(.25))**2),
           'max_logp':float(d.lp.max()),'source':str(src),
           'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
           'selected_genes':int(d.loc[d.selected,'Gene name'].nunique())}
    return d,stats


def load_inputs():
    coordinates={
        'protein_coding':pd.read_csv(BASE/'figure_inputs/gene_coordinates.csv'),
        'ncRNA':pd.read_csv(NEW/'ncRNA_gene_coordinates.csv')}
    for key,coords in coordinates.items():
        coords['chromosome_name']=pd.to_numeric(coords.chromosome_name,errors='coerce')
        coords=coords[coords.chromosome_name.between(1,22)].copy()
        coords['pos']=(coords.start_position+coords.end_position)/2
        assert not coords.duplicated(['hgnc_symbol','chromosome_name']).any()
        coordinates[key]=coords
    sv={a:pd.read_csv(svpath(a,'stats.csv')).iloc[0].to_dict() for a in N}
    for s in sv.values():
        s['n'],s['k']=int(s['n']),int(s['k'])
    genes={}; stats={}
    for family in MASKS:
        genes[family],stats[family]={},{}
        for arm in (N if family in ('coding','noncoding') else ['riga','pooled']):   # 2026-10-19: matched noncoding sets now exist (NEW/matched_noncoding.csv)
            coords=coordinates['ncRNA' if family=='ncRNA' else 'protein_coding']
            genes[family][arm],stats[family][arm]=load_gene(family,arm,coords)
    assert (sv['riga']['n'],sv['riga']['k']) == (37900595,296)
    assert (stats['coding']['riga']['n'],stats['coding']['riga']['k']) == (69555,3)
    assert (stats['noncoding']['riga']['n'],stats['noncoding']['riga']['k']) == (101641,10)
    assert (stats['ncRNA']['riga']['n'],stats['ncRNA']['riga']['k']) == (7704,1)
    return sv,genes,stats


def annotate(ax, x, y, name, target, arm, family, record, ha='center', size=6.5):
    col = COLORS[arm] if record.get('plot') == 'qq' else '#161616'
    art=ax.annotate(name,xy=(x,y),xytext=target,textcoords='axes fraction',
                    fontsize=size,ha=ha,va='bottom',color=col,
                    fontstyle='normal' if 'region' in name else 'italic',
                    bbox={'facecolor':'white','edgecolor':'none','pad':.2},
                    arrowprops={'arrowstyle':'-','color':col,'lw':.45,
                                'shrinkA':2,'shrinkB':1.8},
                    zorder=10,annotation_clip=False)
    ANNOTATIONS.append((ax,name,art))
    ANNOTATION_DATA.append(dict(name=name,arm=arm,family=family,x=float(x),
                                y_display=float(y),**record))


def axes_style(ax, xlim, ylim, compressed=False, genomic=False):
    ax.set_xlim(xlim);ax.set_ylim(ylim);ax.grid(False);ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_visible(True);spine.set_linewidth(.6);spine.set_color('black')
    ax.tick_params(axis='both',labelsize=6.5,length=2.3,width=.55,pad=2.4)
    if genomic:
        ch=np.array([1,3,5,7,9,11,13,15,17,19,22])
        ax.set_xticks(TICKS[ch-1],[str(c) for c in ch]);ax.tick_params(axis='x',length=0)
    else:
        ax.set_xticks(np.arange(0,xlim[1]+.1,2 if xlim[1]>6 else 1))
        ax.set_xlabel(r'Expected $-\log_{10} P$',fontsize=7,fontweight='bold',labelpad=4)
        ax.set_ylabel(r'Observed $-\log_{10} P$',fontsize=7,fontweight='bold',labelpad=4)
    if compressed:
        ticks=np.array([0,5,10,15,100,200,300,400,500,600,800,1000])
        ticks=ticks[compress(ticks)<=ylim[1]]
        ax.set_yticks(compress(ticks),[str(t) for t in ticks])
        y=15/ylim[1]
        for side in [0,1]:
            for shift in [-.013,.013]:
                ax.plot([side-.012,side+.012],[y+shift-.009,y+shift+.009],
                        transform=ax.transAxes,lw=.6,color='black',clip_on=False,zorder=12)
    else:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5,integer=True))


def manhattan_sv(ax,arm,stats,ymax):
    cc=np.memmap(svpath(arm,'chr.i32'),dtype='<i4',mode='r')
    pp=np.memmap(svpath(arm,'pos.i32'),dtype='<i4',mode='r')
    lp=np.memmap(svpath(arm,'logp.f32'),dtype='<f4',mode='r')
    assert len(cc)==len(pp)==len(lp)==stats['n']
    audits=[]
    for c in range(1,23):
        keep=cc==c
        xx=(np.asarray(pp[keep],dtype=float)+OFFSET[c])/1e6
        audits.append(raster(ax,xx,compress(lp[keep]),NAVY if c%2 else GOLD,XLIM,(0,ymax)))
    if arm=='riga':
        folder=ROOT/'output/matched_lpa_pc_pair33_projcell/individual_analysis_qc7'
        sig=pd.read_csv(folder/'manhattan_MAC_20_sig_annotated_v4.csv');sig['logp']=-np.log10(sig.pvalue)
        labs=pd.read_csv(folder/'manhattan_MAC_20_labels_v4.csv')
        targets={'BIN1':(.11,.39),'DGKB':(.29,.52),'BAIAP2L1':(.27,.83),'KCNV1':(.29,.68),
                 'SHARPIN':(.56,.83),'TUBB8':(.56,.52),'CRTAC1':(.59,.68),'MS4A4A':(.82,.52),
                 'MMP17':(.81,.83),'KLHL1':(.87,.68),'CELF4':(.85,.37)}
        for row in labs.itertuples():
            name='APOE region' if row.label.startswith('APOE') else row.label
            target=(.98,.93) if name=='APOE region' else targets[name]
            annotate(ax,(row.POS+OFFSET[row.CHR])/1e6,compress(-np.log10(row.pvalue)),name,
                     target,arm,'sv',{'plot':'manhattan','chr':int(row.CHR),'position':int(row.POS),
                                      'p':row.pvalue},ha='right' if name=='APOE region' else 'center')
    else:
        sig=pd.read_csv(CACHE/f'{arm}_sv_sig.csv')
        apoe=sig[(sig.CHR==19)&sig.POS.between(44.4e6,46.4e6)]
        row=apoe.loc[apoe.logp.idxmax()]
        annotate(ax,(row.POS+OFFSET[19])/1e6,compress(row.logp),'APOE region',(.98,.93),arm,'sv',
                 {'plot':'manhattan','chr':19,'position':int(row.POS),'logp':row.logp},ha='right')
    assert len(sig)==stats['k']
    ax.scatter((sig.POS+OFFSET[sig.CHR.to_numpy(dtype=int)])/1e6,compress(sig.logp),
               s=3.5,c=SEL,linewidths=0,rasterized=len(sig)>5000,zorder=4)
    ax.axhline(compress(-np.log10(stats['critical'])),c=RED,lw=.6,ls=(0,(1,2)),zorder=3)
    axes_style(ax,XLIM,(0,ymax),True,True)
    return {'input_count':stats['n'],'rendered_count':sum(x['rendered_count'] for x in audits),
            'clipped_count':sum(x['clipped_count'] for x in audits)}


def gene_hit(data,gene):
    hit=data[data['Gene name']==gene]
    return None if hit.empty else hit.loc[hit.lp.idxmax()]


def manhattan_gene(ax,data,stats,family,arm,ymax):
    count=0
    for category,_,marker in MASKS[family]:
        for parity,col in [(1,NAVY),(0,GOLD)]:
            s=data[(data.Category==category)&(data.Chr%2==parity)]
            kw={'color':col} if marker in ('+','x') else {'facecolors':'none','edgecolors':col}
            ax.scatter(s.x,s.lp,s=3.7,marker=marker,linewidths=.32,rasterized=True,zorder=1,**kw)
            count+=len(s)
        s=data[(data.Category==category)&data.selected]
        kw={'color':SEL} if marker in ('+','x') else {'facecolors':'none','edgecolors':SEL}
        ax.scatter(s.x,s.lp,s=9,marker=marker,linewidths=.55,zorder=4,
                   rasterized=len(s)>5000,**kw)
    assert count==stats['n'] and data.lp.between(0,ymax).all()
    if stats['k']:
        ax.axhline(-np.log10(stats['critical']),c=RED,lw=.6,ls=(0,(1,2)),zorder=3)
    targets=({'coding':{'KAZN':(.13,.36),'TREM2':(.48,.49),'NECTIN2':(.97,.36)},
              'noncoding':{'TREM2':(.32,.30),'PGAM5':(.60,.40),'MIPOL1':(.73,.51),'TOMM40':(.97,.35)},
              'ncRNA':{'RN7SL293P':(.60,.39)}}[family] if arm=='riga' else {})
    for gene,target in targets.items():
        r=gene_hit(data[data.selected],gene)
        if r is not None:
            annotate(ax,r.x,r.lp,gene,target,arm,family,
                     {'plot':'manhattan','chr':int(r.Chr),'position':float(r.pos),'p':r.SB,'category':r.Category},
                     ha='right' if gene in ('NECTIN2','TOMM40') else 'center')
    if arm=='pooled' and family in ('coding','noncoding'):
        d=data[(data.Chr==19)&data.pos.between(44.4e6,46.4e6)]
        r=d.loc[d.lp.idxmax()]
        annotate(ax,r.x,r.lp,'APOE region',(.97,.91),arm,family,
                 {'plot':'manhattan','chr':19,'position':float(r.pos),'p':r.SB,
                  'target_gene':r['Gene name'],'category':r.Category},ha='right')
    axes_style(ax,XLIM,(0,ymax),False,True)
    return {'input_count':stats['n'],'rendered_count':count,'clipped_count':0}


def mask_legend(fig,family,y):
    handles=[Line2D([],[],marker=m,ls='none',color=GRAY,markerfacecolor='none',
                    markeredgewidth=.5,markersize=3.1,label=t) for _,t,m in MASKS[family]]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.073,y),ncol=len(handles),
               frameon=False,fontsize=5.9,columnspacing=.85,handletextpad=.25,borderaxespad=0)


def check_annotations(fig):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();overlaps=[];outside=[]
    for i,(ax,name,art) in enumerate(ANNOTATIONS):
        box=art.get_bbox_patch().get_window_extent(renderer)
        if not fig.bbox.contains(box.x0,box.y0) or not fig.bbox.contains(box.x1,box.y1):
            outside.append(name)
        for other_ax,other,other_art in ANNOTATIONS[i+1:]:
            if ax is other_ax and box.overlaps(other_art.get_bbox_patch().get_window_extent(renderer)):
                overlaps.append([name,other])
    assert not overlaps,overlaps
    assert not outside,outside
    return {'label_overlap_pairs':overlaps,'labels_outside_canvas':outside}


def save_figure(fig,name,record):
    record['annotation_checks']=check_annotations(fig)
    record['annotations']=list(ANNOTATION_DATA)
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{name}.{ext}',dpi=600)
    fig.savefig(BASE/f'{name}_preview.png',dpi=160)
    (OUT/f'{name}.provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    plt.close(fig)


def make_fig4(sv,genes,stats):
    ANNOTATIONS.clear();ANNOTATION_DATA.clear()
    fig=plt.figure(figsize=(183/25.4,218/25.4))
    lefts=[.075,.580];bottoms=[.705,.415,.125];width=.405;height=.195
    families=['sv','coding','noncoding'];audits={};ranges={}
    for col,arm in enumerate(['riga','pooled']):
        fig.text(lefts[col]+width/2,.987,NAMES[arm],fontsize=7,fontweight='bold',ha='center',va='top')
        fig.text(lefts[col]+width/2,.969,f'N = {N[arm]:,}',fontsize=7,ha='center',va='top',color=GRAY)
    for row,family in enumerate(families):
        fig.text(.075,bottoms[row]+height+.021,TITLES[family],fontsize=7,fontweight='bold',va='bottom')
        st=sv if family=='sv' else stats[family]
        maximum=max(st[a]['max_logp'] for a in ['riga','pooled'])
        ymax=float(np.ceil((compress(maximum) if family=='sv' else maximum)*1.08))
        ranges[family]=[0,ymax]
        for col,arm in enumerate(['riga','pooled']):
            ax=fig.add_axes([lefts[col],bottoms[row],width,height])
            ax.text(-.005,1.035,'abcdefgh'[row*2+col],transform=ax.transAxes,
                    fontweight='bold',fontsize=8,va='bottom')
            if family=='sv':audits[f'{arm}_{family}']=manhattan_sv(ax,arm,st[arm],ymax)
            else:audits[f'{arm}_{family}']=manhattan_gene(ax,genes[family][arm],st[arm],family,arm,ymax)
            if col==0:ax.set_ylabel(r'$-\log_{10} P$',fontsize=7,labelpad=2)
            if row==2:ax.set_xlabel('Chromosome',fontsize=7,labelpad=3)
    fig.text(.075,.673,r'Single-variant axes use 16-fold compression above $-\log_{10} P = 15$.',fontsize=6.3,color=GRAY)
    mask_legend(fig,'coding',.377);mask_legend(fig,'noncoding',.063)
    h=[Line2D([],[],c=RED,lw=.65,ls=(0,(1,2))),Line2D([],[],c=SEL,marker='o',ls='none',markersize=3)]
    fig.legend(h,['Arm-specific BH threshold at nominal FDR 5%','BH-selected associations'],
               loc='lower left',bbox_to_anchor=(.073,.030),frameon=False,ncol=2,fontsize=6.5,
               handlelength=2.2,borderaxespad=0)
    save_figure(fig,'fig4_qc7_complete_20261015',{
        'size_mm':[183,218],'panel_mapping':dict(zip('abcdef',
         [f'{f} {a}' for f in families for a in ['riga','pooled']])),
        'sample_sizes':N,'single_variant':sv,'gene_sets':stats,'rendering':audits,
        'shared_y_ranges':ranges,'style':'STAARpipelineSummary navy/gold chromosomes, red BH symbols and dotted lines',
        'native_eligible_universes_preserved':True,'capped_values':0,
        'single_variant_transform':{'threshold':15,'compression_factor':16}})


def qq_points(ax,observed,arm,xlim,ylim,compressed):
    if not isinstance(observed,np.memmap):observed=np.sort(np.asarray(observed))[::-1]
    n=len(observed)
    # A sorted cached memmap or the complete sorted gene-family array supplies every rank.
    for start in range(0,n-1,1000000):
        chunk=observed[start:min(n,start+1000001)]
        assert (chunk[:-1]>=chunk[1:]).all()
    expected=-np.log10(np.arange(1,n+1,dtype=float)/(n+1))
    yy=compress(observed) if compressed else observed
    audit=raster(ax,expected,yy,COLORS[arm],xlim,ylim,width=1700,height=1450,radius=1,alpha=.75)
    tail=expected>=3.4
    ax.scatter(expected[tail],yy[tail],s=3.1,c=COLORS[arm],marker=TAIL_MARKERS[arm],linewidths=0,zorder=3)
    audit.update(vector_tail_count=int(tail.sum()),rank_definition='i/(N+1)')
    return observed,audit


def qq_gene_labels(ax,data,family,sorted_lp,compact=False):
    targets={'coding':{'TREM2':(.20,.44),'NECTIN2':(.42,.30),'KAZN':(.76,.39)},
             'noncoding':{'TOMM40':(.30,.34),'TREM2':(.60,.24)},
             'ncRNA':{'RN7SL293P':(.50,.35)}}[family]
    for gene,target in targets.items():
        r=gene_hit(data[data.selected],gene)
        if r is None:continue
        rank=int(np.searchsorted(-sorted_lp,-r.lp,side='left'))+1
        x=-np.log10(rank/(len(data)+1))
        annotate(ax,x,r.lp,gene,target,'riga',family,{'plot':'qq','rank':rank,'p':r.SB,
                 'observed_logp':r.lp,'category':r.Category,'chr':int(r.Chr),
                 'position':float(r.pos)},size=5.8 if compact else 7)


def qq_sv_labels(ax,sorted_lp):
    folder=ROOT/'output/matched_lpa_pc_pair33_projcell/individual_analysis_qc7'
    labels=pd.read_csv(folder/'manhattan_MAC_20_labels_v4.csv')
    cc=np.memmap(svpath('riga','chr.i32'),dtype='<i4',mode='r')
    pp=np.memmap(svpath('riga','pos.i32'),dtype='<i4',mode='r')
    lp=np.memmap(svpath('riga','logp.f32'),dtype='<f4',mode='r')
    for name,target in {'APOE region':(.40,.86),'BIN1':(.38,.34)}.items():
        hit=labels[labels.label.str.startswith('APOE')] if name=='APOE region' else labels[labels.label==name]
        r=hit.iloc[0];keep=(cc==r.CHR)&(pp==r.POS);actual=lp[keep]
        assert len(actual)>0
        value=float(np.max(actual))
        # Float32 ranks are determined from the exact stored plotting value.
        rank=int(np.searchsorted(-sorted_lp,-value,side='left'))+1
        x=-np.log10(rank/(len(sorted_lp)+1))
        annotate(ax,x,compress(value),name,target,'riga','sv',{'plot':'qq','rank':rank,
                 'observed_logp':value,'chr':int(r.CHR),'position':int(r.POS)},size=5.8)


def make_fig5(sv,genes,stats):
    ANNOTATIONS.clear();ANNOTATION_DATA.clear()
    fig=plt.figure(figsize=(183/25.4,90/25.4))
    families=['sv','coding','noncoding'];audits={};ranges={};mapping={}
    positions=[(.075,.25),(.405,.25),(.735,.25)]
    for i,(family,(left,bottom)) in enumerate(zip(families,positions)):
        ax=fig.add_axes([left,bottom,.245,.53])
        st=sv if family=='sv' else stats[family]
        arms=list(st);compressed=family=='sv'
        maxlp=max(st[a]['max_logp'] for a in arms)
        ymax=float(np.ceil((compress(maxlp) if compressed else maxlp)*1.08))
        xmax=float(np.ceil(np.log10(max(st[a]['n'] for a in arms)+1)))
        xlim,ylim=(0,xmax),(0,ymax);ranges[family]={'x':xlim,'y_display':ylim}
        audits[family]={};riga_lp=None
        for arm in ['pooled','matched','riga']:
            if arm not in arms:continue
            lp=(np.memmap(svpath(arm,'qq_logp.f32'),dtype='<f4',mode='r') if compressed
                else genes[family][arm].lp.to_numpy())
            sorted_lp,audits[family][arm]=qq_points(ax,lp,arm,xlim,ylim,compressed)
            if arm=='riga':riga_lp=sorted_lp
        diag=np.linspace(0,xmax,200)
        ax.plot(diag,compress(diag) if compressed else diag,c=RED,lw=.8,zorder=5)
        axes_style(ax,xlim,ylim,compressed,False)
        ax.tick_params(labelsize=5.8)
        ax.xaxis.label.set_size(6.3);ax.yaxis.label.set_size(6.3)
        ax.text(-.16,1.045,'abcd'[i],transform=ax.transAxes,fontsize=8,fontweight='bold',va='bottom')
        ax.set_title(TITLES[family],fontsize=7,fontweight='bold',loc='left',pad=12)
        mapping['abcd'[i]]={'family':family,'methods':arms}
        if family=='sv':qq_sv_labels(ax,riga_lp)
        else:qq_gene_labels(ax,genes[family]['riga'],family,riga_lp,compact=True)
        if 'matched' not in arms:
            ax.text(.035,.96,'RIGA and pooling',transform=ax.transAxes,fontsize=5.8,color=GRAY,va='top')
    handles=[Line2D([],[],marker=TAIL_MARKERS[a],ls='none',c=COLORS[a],markersize=4,label=NAMES[a]) for a in N]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.994),ncol=3,
               frameon=False,fontsize=7,columnspacing=1.6,handletextpad=.4)
    fig.text(.075,.075,r'Single-variant axes use 16-fold compression above $-\log_{10} P = 15$.',fontsize=6.3,color=GRAY)
    fig.text(.075,.027,'Labels identify selected RIGA signals. Red lines show the uniform-P reference.',fontsize=6.3,color=GRAY)
    save_figure(fig,'fig5_qc7_overlay_20261015',{
        'size_mm':[183,90],'panel_mapping':mapping,'sample_sizes':N,'single_variant':sv,
        'gene_sets':stats,'rendering':audits,'shared_ranges':ranges,
        'style':'STAARpipelineSummary boxed axes, filled points and solid red reference, with method colors',
        'colors':COLORS,'tail_markers':TAIL_MARKERS,'lambda_annotations':False,'test_count_annotations':False,
        'label_arm':'riga','native_eligible_universes_preserved':True,
        'incomplete_matched_regulatory_or_ncRNA_used':False,'capped_values':0,
        'single_variant_transform':{'threshold':15,'compression_factor':16}})


def make_supp8(sv,genes,stats):
    ANNOTATIONS.clear();ANNOTATION_DATA.clear()
    fig=plt.figure(figsize=(183/25.4,150/25.4))
    family='ncRNA';st=stats[family]
    ymax=float(np.ceil(max(s['max_logp'] for s in st.values())*1.08))
    audits={}
    for i,(arm,left) in enumerate([('riga',.075),('pooled',.580)]):
        ax=fig.add_axes([left,.63,.405,.27])
        audits[f'manhattan_{arm}']=manhattan_gene(ax,genes[family][arm],st[arm],family,arm,ymax)
        ax.set_xlabel('Chromosome',fontsize=7,labelpad=3)
        if i==0:ax.set_ylabel(r'$-\log_{10} P$',fontsize=7,labelpad=2)
        ax.text(-.005,1.045,'ab'[i],transform=ax.transAxes,fontsize=8,fontweight='bold',va='bottom')
        ax.set_title(NAMES[arm],fontsize=7,fontweight='bold',pad=13)
    ax=fig.add_axes([.32,.14,.405,.33])
    xmax=float(np.ceil(np.log10(max(s['n'] for s in st.values())+1)))
    for arm in ['pooled','riga']:
        sorted_lp,audits[f'qq_{arm}']=qq_points(ax,genes[family][arm].lp.to_numpy(),arm,(0,xmax),(0,ymax),False)
        if arm=='riga':riga_lp=sorted_lp
    diag=np.linspace(0,xmax,100);ax.plot(diag,diag,c=RED,lw=.8,zorder=5)
    axes_style(ax,(0,xmax),(0,ymax),False,False)
    ax.text(-.16,1.045,'c',transform=ax.transAxes,fontsize=8,fontweight='bold',va='bottom')
    qq_gene_labels(ax,genes[family]['riga'],family,riga_lp)
    handles=[Line2D([],[],marker=TAIL_MARKERS[a],ls='none',c=COLORS[a],markersize=3.5,label=NAMES[a])
             for a in ['riga','pooled']]
    ax.legend(handles=handles,loc='upper left',frameon=False,fontsize=6.5,handletextpad=.4)
    fig.text(.075,.985,'Noncoding RNA association scans',fontsize=7,fontweight='bold',va='top')
    handles=[Line2D([],[],c=RED,lw=.65,ls=(0,(1,2))),Line2D([],[],c=SEL,marker='o',ls='none',markersize=3)]
    fig.legend(handles,['Arm-specific BH threshold at nominal FDR 5%','BH-selected associations'],
               loc='lower left',bbox_to_anchor=(.073,.525),frameon=False,ncol=2,fontsize=6.3,
               handlelength=2.2,borderaxespad=0)
    fig.text(.075,.045,'Labels identify the RIGA ncRNA signal. The solid red line shows the uniform-P reference.',fontsize=6.3,color=GRAY)
    save_figure(fig,'supp_fig8_qc7_comparison_20261015',{
        'size_mm':[183,150],'panel_mapping':{'a':'RIGA Manhattan','b':'Pooled Manhattan','c':'RIGA and pooled QQ'},
        'sample_sizes':{'riga':N['riga'],'pooled':N['pooled']},'gene_sets':{'ncRNA':st},
        'rendering':audits,'shared_y_range':[0,ymax],'colors':COLORS,'tail_markers':TAIL_MARKERS,
        'native_eligible_universes_preserved':True,'capped_values':0,
        'lambda_annotations':False,'test_count_annotations':False,'label_arm':'riga',
        'coordinate_source':str(NEW/'ncRNA_gene_coordinates.csv')})


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--figures',nargs='+',choices=['4','5','8'],default=['4','5','8']);args=ap.parse_args()
    sv,genes,stats=load_inputs()
    for f in args.figures:
        {'4':make_fig4,'5':make_fig5,'8':make_supp8}[f](sv,genes,stats)
        print('Completed Figure '+f,flush=True)
    summary=[]
    for arm,st in sv.items():summary.append({'family':'sv','arm':arm,**st})
    for fam,arms in stats.items():
        for arm,st in arms.items():summary.append({'family':fam,'arm':arm,**st})
    pd.DataFrame(summary).to_csv(BASE/'figure_comparison_summary.csv',index=False)


if __name__=='__main__':main()
