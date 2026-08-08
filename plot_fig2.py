#!/usr/bin/env python3
"""
covalent5 · Fig 2 합성 — (a) DFT 최적화 3D 구조(ball-and-stick, 수소제외) + 원자 범례(구조 바로 아래)
                          (b) 쌍극자  (c) P 전하
실행: cd "$COV5" && python plot_fig2.py
"""
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams.update({
 'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,
 'legend.fontsize':6.5,'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.8,'axes.axisbelow':True,
 'axes.grid':True,'grid.linewidth':0.5,'grid.alpha':0.4,'legend.frameon':False,
 'savefig.dpi':300,'savefig.bbox':'tight','savefig.pad_inches':0.06})
import matplotlib.pyplot as plt, numpy as np, os
from matplotlib.lines import Line2D
from PIL import Image

HERE=os.path.dirname(os.path.abspath(__file__))
S2=os.path.join(HERE,"figure","mol3d")
struct=["A-230","A-232","A-234","A-242"]
ELEM=[("C","#999999"),("N","#0000FF"),("O","#FF0000"),("P","#FF8000"),("F","#00FFFF")]

dip={"A-230":[11.56,8.76,9.74,10.67],"A-232":[7.83,8.07,11.06,12.80],"A-234":[9.17,10.60,10.61,9.30],
"A-242":[8.11,8.46],"VX":[4.72,3.84],"Sarin":[4.10,4.12],"Soman":[4.20,4.20,4.45,4.21],"Tabun":[5.81,5.70],
"Paraoxon":[2.99],"Parathion":[2.59],"Chlorpyrifos":[3.38],"Malathion":[6.15,7.48],"Dichlorvos":[3.94]}
qP={"A-230":[.773,.801,.787,.773],"A-232":[.895,.923,.845,.834],"A-234":[.870,.868,.867,.855],
"A-242":[.805,.815],"VX":[.698,.716],"Sarin":[.873,.871],"Soman":[.860,.865,.868,.862],"Tabun":[.801,.802],
"Paraoxon":[.850],"Parathion":[.675],"Chlorpyrifos":[.733],"Malathion":[.632,.617],"Dichlorvos":[.936]}
cls={"A-230":"Novichok","A-232":"Novichok","A-234":"Novichok","A-242":"Novichok","Sarin":"G","Soman":"G",
 "Tabun":"G","VX":"V","Parathion":"OP","Paraoxon":"OP","Chlorpyrifos":"OP","Malathion":"OP","Dichlorvos":"OP"}
col={"Novichok":"#DC267F","G":"#648FFF","V":"#FE6100","OP":"#009E73"}
order=["A-230","A-232","A-234","A-242","VX","Sarin","Soman","Tabun","Paraoxon","Parathion","Chlorpyrifos","Malathion","Dichlorvos"]

def autocrop(p):
    im=Image.open(p).convert("RGBA"); bb=im.getbbox()
    if bb:
        pad=20; bb=(max(0,bb[0]-pad),max(0,bb[1]-pad),min(im.width,bb[2]+pad),min(im.height,bb[3]+pad)); im=im.crop(bb)
    return im

fig=plt.figure(figsize=(7.2,5.4),constrained_layout=True)
outer=fig.add_gridspec(2,1,height_ratios=[1.12,1.05],hspace=0.10)
top=outer[0].subgridspec(1,4,wspace=0.06)
bot=outer[1].subgridspec(1,2,wspace=0.12,width_ratios=[1.0,1.0])

saxes=[]
for j,n in enumerate(struct):
    ax=fig.add_subplot(top[0,j]); ax.imshow(autocrop(f"{S2}/{n}.png")); ax.axis('off')
    saxes.append(ax)

def strip(ax,d,ylab):
    for i,name in enumerate(order):
        ys=d[name]; xs=np.full(len(ys),i)+(np.linspace(-0.12,0.12,len(ys)) if len(ys)>1 else 0)
        ax.scatter(xs,ys,s=20,color=col[cls[name]],edgecolor='white',linewidth=0.4,zorder=3)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order,rotation=45,ha='right'); ax.set_ylabel(ylab)

axB=fig.add_subplot(bot[0,0]); strip(axB,dip,"Dipole moment (Debye)")
leg=[Line2D([0],[0],marker='o',color='w',markerfacecolor=col[c],markersize=6,label=c) for c in ["Novichok","G","V","OP"]]
axB.legend(handles=leg,loc='upper right',ncol=2,title="Class",columnspacing=0.8,handletextpad=0.3)
axC=fig.add_subplot(bot[0,1]); strip(axC,qP,"P Mulliken charge (e)")

fig.canvas.draw()
# 원자 범례: 구조 축 하단 바로 아래에 부착
sbot=min(a.get_position().y0 for a in saxes)
for ax,n in zip(saxes,struct):
    p=ax.get_position(); fig.text((p.x0+p.x1)/2, sbot-0.004, n, fontsize=8.5, fontweight='bold', ha='center', va='top')
h=[Line2D([0],[0],marker='o',color='w',markerfacecolor=c,markeredgecolor='gray',markersize=8,label=e) for e,c in ELEM]
fig.legend(handles=h,loc='upper center',bbox_to_anchor=(0.5,sbot-0.05),ncol=5,frameon=False,handletextpad=0.3,columnspacing=1.3)

pB=axB.get_position(); pC=axC.get_position(); ly=max(pB.y1,pC.y1)+0.012
fig.text(0.008,0.997,"(a) Novichok agent structures (DFT-optimised, "+r"$\omega$B97X-D3)",fontsize=8.5,fontweight='bold',va='top',ha='left')
fig.text(0.008,ly,"(b) Molecular polarity",fontsize=8.5,fontweight='bold',va='bottom',ha='left')
fig.text(pC.x0,ly,"(c) Electrophilicity at P",fontsize=8.5,fontweight='bold',va='bottom',ha='left')
fig.savefig(os.path.join(HERE,"figure","fig2_descriptors.png")); print("saved")
