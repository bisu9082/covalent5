import sys, itertools, numpy as np
R={'H':0.31,'C':0.76,'N':0.71,'O':0.66,'P':1.07,'S':1.05}
def load(p):
    L=open(p).read().split('\n'); n=int(L[0].split()[0])
    el=[];xyz=[]
    for ln in L[2:2+n]:
        w=ln.split(); el.append(w[0]); xyz.append([float(x) for x in w[1:4]])
    return el, np.array(xyz)
for path in sys.argv[1:]:
    el,X=load(path); n=len(el)
    adj={i:set() for i in range(n)}
    for i,j in itertools.combinations(range(n),2):
        d=np.linalg.norm(X[i]-X[j]); cut=R.get(el[i],.8)+R.get(el[j],.8)+0.40
        if d<cut: adj[i].add(j); adj[j].add(i)
    # 5원 고리 탐색
    rings=set()
    for a in range(n):
        for b in adj[a]:
            for c in adj[b]-{a}:
                for d in adj[c]-{a,b}:
                    for e in adj[d]-{a,b,c}:
                        if a in adj[e]: rings.add(frozenset((a,b,c,d,e)))
    from collections import Counter
    print(f"\n=== {path}  ({n} atoms, {Counter(el)}) ===")
    print(f"5원 고리 {len(rings)}개")
    hit=0
    for r in rings:
        f=Counter(el[i] for i in r)
        tag=""
        if f['N']==2 and f['C']==3: tag="  <<< 이미다졸 (His 곁사슬)"; hit+=1
        print("   ", dict(f), tag)
    # 아마이드 개수 (C=O 에 N 이 붙은 탄소)
    am=sum(1 for i in range(n) if el[i]=='C'
           and any(el[j]=='O' and np.linalg.norm(X[i]-X[j])<1.35 for j in adj[i])
           and any(el[j]=='N' for j in adj[i]))
    print(f"아마이드 카보닐(C(=O)-N) {am}개")
    print(">>> 판정:", "His447 포함" if hit else "His447 없음")
