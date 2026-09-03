import sys, itertools, numpy as np
from collections import Counter, deque
R={'H':0.31,'C':0.76,'N':0.71,'O':0.66,'P':1.07,'S':1.05}
def load(p):
    L=open(p).read().split('\n'); n=int(L[0].split()[0]); el=[];xyz=[]
    for ln in L[2:2+n]:
        w=ln.split(); el.append(w[0]); xyz.append([float(x) for x in w[1:4]])
    return el,np.array(xyz)
for path in sys.argv[1:]:
    el,X=load(path); n=len(el); adj={i:set() for i in range(n)}
    for i,j in itertools.combinations(range(n),2):
        if np.linalg.norm(X[i]-X[j]) < R.get(el[i],.8)+R.get(el[j],.8)+0.40:
            adj[i].add(j); adj[j].add(i)
    seen=set(); frags=[]
    for s in range(n):
        if s in seen: continue
        q=deque([s]); comp=set()
        while q:
            u=q.popleft()
            if u in comp: continue
            comp.add(u); q.extend(adj[u]-comp)
        seen|=comp; frags.append(sorted(comp))
    print(f"\n=== {path}  {n} atoms, {len(frags)} fragments ===")
    for f in sorted(frags,key=len,reverse=True):
        c=Counter(el[i] for i in f)
        formula=''.join(f"{e}{c[e]}" for e in ['C','H','N','O','P','S'] if c[e])
        ring = any(el[i]=='N' for i in f) and len([i for i in f if el[i]=='N'])==2 and len(f)<16
        note=''
        if 'P' in c: note='  <- VX'
        elif ring: note='  <- imidazole (His447 side chain)'
        elif c['O']==1 and c['C']<=2: note='  <- methoxide/methanol (Ser203 side chain)'
        elif c['O']==2: note='  <- carboxyl/carboxylate (Glu334 side chain)'
        elif c['N']==1 and c['O']==1: note='  <- amide (oxyanion hole)'
        print(f"  {len(f):3d} atoms  {formula:<14}{note}")
