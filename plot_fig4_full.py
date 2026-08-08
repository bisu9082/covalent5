#!/usr/bin/env python3
"""covalent5 · Fig 4 full (a,b,c data + d misranked poses). 실행: cd "$COV5" && python plot_fig4_full.py"""
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,
 'ytick.labelsize':7,'legend.fontsize':6.5,'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.8,'axes.axisbelow':True,
 'axes.grid':True,'grid.linewidth':0.5,'grid.alpha':0.4,'legend.frameon':False,
 'savefig.dpi':300,'savefig.bbox':'tight','savefig.pad_inches':0.06})
import matplotlib.pyplot as plt, numpy as np, os
from matplotlib.lines import Line2D
from PIL import Image
HERE=os.path.dirname(os.path.abspath(__file__)); POSE=os.path.join(HERE,"figure","poses")
def crop(p):
    im=Image.open(p).convert("RGBA"); bb=im.getbbox()
    if bb:
        pad=15; bb=(max(0,bb[0]-pad),max(0,bb[1]-pad),min(im.width,bb[2]+pad),min(im.height,bb[3]+pad)); im=im.crop(bb)
    return im
D=[("BSAR",-5.08,3.425,8.79,"G"),("GF",-5.54,3.776,8.69,"G"),("VR",-5.98,4.333,8.64,"V"),
("VX-Sp",-6.02,5.089,8.15,"V"),("GD",-6.23,3.895,7.96,"G"),("GB",-3.75,3.04,7.43,"G"),
("GA",-4.67,3.525,6.87,"G"),("PXNE",-6.3,3.987,6.34,"OP"),("VX-Rp",-5.97,4.96,6.08,"V"),
("MPXN",-7.01,3.865,6.08,"OP"),("DETAB",-5.76,4.228,5.94,"G"),("DFP",-5.53,3.958,5.11,"OP"),
("METH",-4.13,3.966,3.28,"OP"),("FEN",-6.59,5.426,1.30,"OP")]
col={"G":"#648FFF","V":"#FE6100","OP":"#009E73"}
aff=np.array([d[1] for d in D]); cnn=np.array([d[2] for d in D]); ki=np.array([d[3] for d in D]); cc=[col[d[4]] for d in D]
fig=plt.figure(figsize=(7.2,5.9),constrained_layout=True)
outer=fig.add_gridspec(3,1,height_ratios=[1.0,0.12,0.92])
top=outer[0].subgridspec(1,3,width_ratios=[1,1,1.18]); bot=outer[2].subgridspec(1,2,wspace=0.05)
ax=fig.add_subplot(top[0,0]); ax.scatter(aff,ki,s=30,c=cc,edgecolor='white',linewidth=0.5,zorder=3)
ax.set_xlabel("Docking affinity (kcal/mol)"); ax.set_ylabel("Experimental log $k_i$")
ax.set_title("(a) Vina affinity vs potency",loc='left',fontweight='bold',fontsize=8.5)
ax.text(0.04,0.06,r"$\rho=+0.12$ (n.s.)",transform=ax.transAxes,fontsize=8,fontweight='bold',color='#444',va='bottom')
leg=[Line2D([0],[0],marker='o',color='w',markerfacecolor=col[c],markersize=6,label=c) for c in ["G","V","OP"]]
ax.legend(handles=leg,loc='lower right',handletextpad=0.3)
ax2=fig.add_subplot(top[0,1]); ax2.scatter(cnn,ki,s=30,c=cc,edgecolor='white',linewidth=0.5,zorder=3)
for d in D:
    if d[0] in ("FEN","BSAR"): ax2.annotate(d[0],(d[2],d[3]),fontsize=6.5,fontweight='bold',xytext=(6,-2 if d[0]=="FEN" else 4),textcoords='offset points')
ax2.set_xlabel("GNINA CNN affinity"); ax2.set_ylabel("Experimental log $k_i$")
ax2.set_title("(b) CNN affinity vs potency",loc='left',fontweight='bold',fontsize=8.5)
ax2.text(0.04,0.06,r"$\rho=-0.40$ (all, n=14)",transform=ax2.transAxes,fontsize=8,fontweight='bold',color='#B00',va='bottom')
ax3=fig.add_subplot(top[0,2]); meth=["Vina\naffinity","GNINA\nCNN","GNINA\ncovalent","Physics\ndescriptors"]; rho=[0.12,-0.40,0.23,0.55]
ci_lo=[-0.47,-0.82,-0.57,0.04]; ci_hi=[0.64,0.21,0.87,0.74]
colb=["#999999","#DC267F","#FE6100","#009E73"]; y=np.arange(4)
xerr=np.array([[r-lo for r,lo in zip(rho,ci_lo)],[hi-r for r,hi in zip(rho,ci_hi)]])
ax3.barh(y,rho,color=colb,edgecolor='white',height=0.6,zorder=3)
ax3.errorbar(rho,y,xerr=xerr,fmt='none',ecolor='#444',elinewidth=1.0,capsize=3,capthick=1.0,zorder=4)
ax3.axvline(0,color='k',lw=0.9,zorder=2)
ax3.set_yticks(y); ax3.set_yticklabels(meth,fontsize=7); ax3.set_xlabel(r"Spearman $\rho$ (rank vs $k_i$)"); ax3.set_xlim(-1.0,1.08)
for yi,v,hi in zip(y,rho,ci_hi): ax3.text(hi+0.05,yi,f"{v:+.2f}",va='center',ha='left',fontsize=7,fontweight='bold',color='#333')
ax3.set_title("(c) Method comparison (95% CI)",loc='left',fontweight='bold',fontsize=8.5)
poses=[("FEN_00_P3R","Fenamiphos","strong docking / weak inhibitor","#009E73"),
       ("BSAR_00_P5S","Butylsarin","weak docking / strong inhibitor","#648FFF")]
for j,(sid,nm,sub,c) in enumerate(poses):
    axp=fig.add_subplot(bot[0,j]); axp.imshow(crop(f"{POSE}/{sid}.png")); axp.axis('off')
    axp.set_title(f"{nm}\n{sub}",fontsize=7.5,color=c,fontweight='bold')
fig.text(0.008,0.455,"(d) Docking misranks the two extremes (same gorge, opposite potency)",fontsize=8.5,fontweight='bold',va='bottom',ha='left')
fig.savefig("figure/fig4_rq1.png"); print("saved figure/fig4_rq1.png")
