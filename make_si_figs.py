#!/usr/bin/env python3
"""
covalent5 · R2 Supporting Information 그림 2종
  S2  도킹 박스 방향·크기 + 대표 top-ranked 포즈       (Reviewer 1 Minor 8, Major 2)
  S3  14개 kinetic anchor 구조 격자 + 측정 k_i          (Reviewer 1 Minor 1)
출처: docking/4EY7.pdb, docking/out_noncovalent/*.sdf.gz,
      docking/run_gnina_noncovalent.sh, data/anchor_ki.csv, data/structures/*.sdf
"""
import csv, gzip, os, warnings
import numpy as np
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams.update({
    'font.size': 8, 'axes.titlesize': 8.5, 'axes.labelsize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 6.8,
    'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.8,
    'axes.axisbelow': True, 'grid.linewidth': 0.5, 'grid.alpha': 0.35,
    'legend.frameon': False, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.06})
import matplotlib.pyplot as plt
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
warnings.filterwarnings("ignore")

C_BOX="#6F5BD5"; C_NEAR="#0F9D8F"; C_FAR="#F4611A"; C_PROT="#B8BDC4"
C_TRIAD="#E05C5C"; INK="#1f2328"; MUTED="#6B7280"

CENTER=np.array([-5.337,-41.735,29.543]); SIZE=22.5      # run_gnina_noncovalent.sh
OG=CENTER                                                 # = Ser203 OG (chain A)
NEAR=4.5

anchors=list(csv.DictReader(open('data/anchor_ki.csv')))
PRETTY={'VX_00_P3S':'VX ($S_\\mathrm{P}$)','VX_01_P3R':'VX ($R_\\mathrm{P}$)',
 'GB_00_P4S':'Sarin','GD_00_C1S_P3S':'Soman','GA_00_P3R':'Tabun',
 'PXNE_00_achiral':'Paraoxon-ethyl','BSAR_00_P5S':'Butylsarin','GF_00_P1S':'Cyclosarin',
 'VR_00_P8S':'VR','MPXN_00_achiral':'Paraoxon-methyl','DETAB_00_P3R':'Diethyltabun',
 'DFP_00_achiral':'DFP','METH_00_P2R':'Methamidophos','FEN_00_P3R':'Fenamiphos'}

# ── 4EY7 chain A Cα + 촉매 삼중쌍 ────────────────────────────────────────
ca=[]; triad={}
for ln in open('docking/4EY7.pdb'):
    if not ln.startswith('ATOM'): continue
    if ln[21]!='A': continue
    name=ln[12:16].strip(); res=ln[17:20].strip(); num=int(ln[22:26])
    xyz=np.array([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])
    if name=='CA': ca.append(xyz)
    if (res,num) in [('SER',203),('HIS',447),('GLU',334)] and name=='CA':
        triad[f"{res.title()}{num}"]=xyz
ca=np.array(ca)

# ── 각 앵커 top pose 의 P···OG 거리 ──────────────────────────────────────
rows=[]
for a in anchors:
    sid=a['stereo_id']
    p=f'docking/out_noncovalent/{sid}_nc.sdf.gz'
    m=[x for x in Chem.ForwardSDMolSupplier(gzip.open(p)) if x][0]
    cf=m.GetConformer()
    Ps=[np.array(cf.GetAtomPosition(at.GetIdx())) for at in m.GetAtoms() if at.GetSymbol()=='P']
    d=min(np.linalg.norm(q-OG) for q in Ps)
    rows.append((sid, PRETTY[sid], d, Ps[int(np.argmin([np.linalg.norm(q-OG) for q in Ps]))]))
rows.sort(key=lambda r: r[2])
n_near=sum(1 for r in rows if r[2]<=NEAR)

# ══════════════════ Figure S2 ═══════════════════════════════════════════
fig=plt.figure(figsize=(7.0,3.3))
gs=fig.add_gridspec(1,2,width_ratios=[1.05,1.0],wspace=0.34)

ax=fig.add_subplot(gs[0,0],projection='3d')
near_prot=ca[np.linalg.norm(ca-CENTER,axis=1)<26]
ax.plot(near_prot[:,0],near_prot[:,1],near_prot[:,2],'-',color=C_PROT,lw=0.7,alpha=0.85,zorder=1)
h=SIZE/2; c=CENTER
corners=np.array([[c[0]+sx*h,c[1]+sy*h,c[2]+sz*h] for sx in(-1,1) for sy in(-1,1) for sz in(-1,1)])
edges=[(0,1),(0,2),(0,4),(1,3),(1,5),(2,3),(2,6),(3,7),(4,5),(4,6),(5,7),(6,7)]
for i,j in edges:
    ax.plot(*zip(corners[i],corners[j]),color=C_BOX,lw=1.0,alpha=0.9,zorder=3)
_off={'Ser203':np.array([1.0,1.0,2.6]),'His447':np.array([-6.5,1.0,-1.4]),
      'Glu334':np.array([1.0,1.0,-3.6])}
for nm,xyz in triad.items():
    ax.scatter(*xyz,s=26,color=C_TRIAD,edgecolor='white',linewidth=0.5,zorder=5)
    ax.text(*(xyz+_off.get(nm,np.array([0.6,0.6,1.4]))),nm,fontsize=5.8,color=C_TRIAD)
P=np.array([r[3] for r in rows]); dd=np.array([r[2] for r in rows])
ax.scatter(P[dd<=NEAR,0],P[dd<=NEAR,1],P[dd<=NEAR,2],s=16,color=C_NEAR,
           edgecolor='white',linewidth=0.4,zorder=6,label=f'P, $\\leq${NEAR} Å')
ax.scatter(P[dd>NEAR,0],P[dd>NEAR,1],P[dd>NEAR,2],s=16,color=C_FAR,
           edgecolor='white',linewidth=0.4,zorder=6,label=f'P, $>${NEAR} Å')
ax.set_xlabel('x (Å)',labelpad=-8,fontsize=6.5); ax.set_ylabel('y (Å)',labelpad=-8,fontsize=6.5)
ax.set_zlabel('z (Å)',labelpad=-8,fontsize=6.5)
ax.tick_params(labelsize=5.0,pad=-3)
ax.set_title('(a) Search box on hAChE (4EY7 chain A)',loc='left',fontweight='bold',pad=-2)
ax.view_init(elev=16,azim=-64)
ax.legend(loc='upper left',bbox_to_anchor=(-0.02,0.98),handletextpad=0.3,labelspacing=0.25)
ax.text2D(0.01,0.76,f'box {SIZE}$\\times${SIZE}$\\times${SIZE} Å\ncentre = Ser203 O$_\\gamma$\n'
          f'({c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f})',transform=ax.transAxes,ha='left',
          va='top',fontsize=5.5,color=MUTED,linespacing=1.35)

ax=fig.add_subplot(gs[0,1])
y=np.arange(len(rows))[::-1]
cols=[C_NEAR if r[2]<=NEAR else C_FAR for r in rows]
ax.barh(y,[r[2] for r in rows],height=0.62,color=cols,zorder=3)
ax.axvline(NEAR,color=INK,lw=0.9,ls=(0,(4,3)),zorder=4)
for yi,r in zip(y,rows):
    ax.text(r[2]+0.12,yi,f'{r[2]:.2f}',va='center',fontsize=6.0,color=INK)
ax.set_yticks(y); ax.set_yticklabels([r[1] for r in rows],fontsize=6.4)
ax.set_xlim(0,10.6); ax.set_xlabel('P$\\cdots$Ser203 O$_\\gamma$, top-ranked pose (Å)')
ax.set_title('(b) Near-attack distance, 14 kinetic anchors',loc='left',fontweight='bold')
ax.grid(axis='x'); ax.grid(axis='y',visible=False)
ax.text(0.985,0.96,f'{n_near} of {len(rows)} within {NEAR} Å',transform=ax.transAxes,
        ha='right',va='top',fontsize=6.6,color=INK,fontweight='bold')
fig.savefig('figS2_docking_box.png'); plt.close(fig)
print(f'figS2_docking_box.png  ({n_near}/{len(rows)} within {NEAR} Å)')

# ══════════════════ Figure S3 ═══════════════════════════════════════════
from rdkit.Chem import Draw, AllChem
mols=[]; legends=[]
for a in anchors:
    sid=a['stereo_id']
    path=None
    for root,_,fs in os.walk('data/structures'):
        if sid+'.sdf' in fs: path=os.path.join(root,sid+'.sdf'); break
    m=Chem.MolFromMolFile(path) if path else None
    if m is None:
        m=[x for x in Chem.ForwardSDMolSupplier(gzip.open(f'docking/out_noncovalent/{sid}_nc.sdf.gz')) if x][0]
    m=Chem.RemoveHs(m); AllChem.Compute2DCoords(m)
    mols.append(m)
    nm=PRETTY[sid].replace('$\\mathrm{P}$','P').replace('$S_\\mathrm{P}$','Sp').replace('$R_\\mathrm{P}$','Rp')
    legends.append(f"{nm}\nlog k_i = {float(a['pKi']):.2f}")
img=Draw.MolsToGridImage(mols,molsPerRow=4,subImgSize=(330,300),legends=legends,
                         useSVG=False,returnPNG=False)
img.save('figS3_anchor_structures.png')
print(f'figS3_anchor_structures.png  ({len(mols)} structures)')

# ══ 검증표 CSV (SI 표 + 응답레터 근거용) ═════════════════════════════════
with open('near_attack_distances.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['stereo_id','name','P_OG_top_pose_A','within_4.5A'])
    for r in rows: w.writerow([r[0],r[1].replace('$','').replace('\\mathrm{P}','P'),
                               f'{r[2]:.2f}', int(r[2]<=NEAR)])
print('near_attack_distances.csv')
