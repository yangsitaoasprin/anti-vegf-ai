# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# =============================================================================
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#5 接触残基级结合模式分析 (复用 analyze_binding_mode_4.py 的位点映射逻辑)
输入: designs/candidates_6/structures/vegf_len100_3_T0.1_s2.complex.pdb
      (链A=VEGF 1-98, 链B=binder 1-100, 原始 Boltz 无穿模, ΔG=-38.3)
输出: results/binding_mode/contacts_vegf_len100_3_T0.1_s2.csv
      results/binding_mode/binding_mode_report.md (追加 #5 章节)
"""
import os, csv, math
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, ".."))
DATA = os.path.join(PROJ, "data")
RES  = os.path.join(PROJ, "results", "binding_mode")
os.makedirs(RES, exist_ok=True)

AA = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
      'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
      'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
CUT = 5.0
REF = open(os.path.join(DATA, "vegfA_W.seq")).read().strip()
assert len(REF) == 98

def load_chains(path):
    chains={}; cur=None; last=None
    for line in open(path):
        if not line.startswith(("ATOM","HETATM")): continue
        ch=line[21]
        try:
            resnum=int(line[22:26]); resname=line[17:20].strip()
            atom=line[12:16].strip(); x,y,z=float(line[30:38]),float(line[38:46]),float(line[46:54])
        except ValueError: continue
        if atom.startswith("H"): continue
        key=(ch,resnum,line[26])
        if cur is None or key!=last:
            cur={"resnum":resnum,"resname":resname,"aa":AA.get(resname,'?'),"atoms":[]}
            chains.setdefault(ch,[]).append(cur); last=key
        cur["atoms"].append((x,y,z))
    return chains

def min_dist(ra,rb):
    best=None
    for ax,ay,az in ra["atoms"]:
        for bx,by,bz in rb["atoms"]:
            d=math.sqrt((ax-bx)**2+(ay-by)**2+(az-bz)**2)
            if best is None or d<best: best=d
    return best

def interface(vseg,rseg,cutoff=CUT):
    out=[]
    for i,v in enumerate(vseg):
        best=None; bpart=[]
        for j,r in enumerate(rseg):
            d=min_dist(v,r)
            if d is not None and (best is None or d<best): best=d
            if d is not None and d<=cutoff: bpart.append(j)
        if best is not None and best<=cutoff: out.append((i,best,bpart))
    return out

def nw(a,b):
    n,m=len(a),len(b)
    dp=[[0]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): dp[i][0]=-i
    for j in range(1,m+1): dp[0][j]=-j
    for i in range(1,n+1):
        for j in range(1,m+1):
            s=1 if a[i-1]==b[j-1] else -1
            dp[i][j]=max(dp[i-1][j-1]+s,dp[i-1][j]-1,dp[i][j-1]-1)
    i,j=n,m; ali=[]
    while i>0 or j>0:
        if i>0 and j>0 and dp[i][j]==dp[i-1][j-1]+(1 if a[i-1]==b[j-1] else -1):
            ali.append((i-1,j-1)); i-=1; j-=1
        elif i>0 and dp[i][j]==dp[i-1][j]-1: ali.append((i-1,None)); i-=1
        else: ali.append((None,j-1)); j-=1
    ali.reverse(); return ali

def build_map(vseg):
    vseq="".join(r["aa"] for r in vseg)
    ali=nw(vseq,REF); vmap=[None]*len(vseg)
    for vi,ri in ali:
        if vi is not None and ri is not None: vmap[vi]=ri+1
    return vmap

def site_residues(pdb,veg_ch,rec_chs):
    chains=load_chains(os.path.join(DATA,pdb)); vseg=chains[veg_ch]; vmap=build_map(vseg)
    idxs=set()
    for rc in rec_chs:
        for vi,d,_ in interface(vseg,chains[rc]):
            oi=vmap[vi]
            if oi: idxs.add(oi)
    return idxs

kdr_idx=site_residues("3V2A.pdb","A",["R"])
flt_idx=site_residues("1FLT.pdb","W",["X","Y"])
dim_idx=site_residues("1FLT.pdb","W",["V"])
print("KDR 位点:",sorted(kdr_idx)); print("FLT1 位点:",sorted(flt_idx)); print("二聚面:",sorted(dim_idx))

# ---- #5 ----
cid="vegf_len100_3_T0.1_s2"
pdb=os.path.join(PROJ,"designs","candidates_6","structures",f"{cid}.complex.pdb")
chains=load_chains(pdb); A=chains["A"]; B=chains["B"]
vmap=build_map(A); bseq="".join(r["aa"] for r in B)
iface=interface(A,B)
rows=[]; kdr_h=[]; flt_h=[]; dim_h=[]; other_h=[]
for vi,d,bpart in iface:
    our=vmap[vi]; bparts=sorted({b+1 for b in bpart}); bstr=",".join(str(x) for x in bparts)
    rows.append({"idx":our,"veg":f"{A[vi]['resname']}{our}","dist":round(d,2),
                 "binder":bstr,"nb":len(bparts),
                 "KDR":"Y" if our in kdr_idx else "","FLT1":"Y" if our in flt_idx else "",
                 "DIM":"Y" if our in dim_idx else ""})
    if our in kdr_idx: kdr_h.append(our)
    elif our in flt_idx: flt_h.append(our)
    elif our in dim_idx: dim_h.append(our)
    else: other_h.append(our)
with open(os.path.join(RES,f"contacts_{cid}.csv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["idx","veg","dist","binder","nb","KDR","FLT1","DIM"])
    w.writeheader(); [w.writerow(r) for r in rows]
kdr_names=",".join(f"{A[o-1]['resname']}{o}" for o in kdr_h)
flt_names=",".join(f"{A[o-1]['resname']}{o}" for o in flt_h)
dim_names=",".join(f"{A[o-1]['resname']}{o}" for o in dim_h)
oth_names=",".join(f"{A[o-1]['resname']}{o}" for o in other_h[:60])
print(f"[{cid}] 接触={len(rows)} KDR重叠={kdr_h} FLT1重叠={flt_h} 二聚面重叠={dim_h}")

# 主导位点判定
def dominant():
    scores={"KDR":len(kdr_h),"Flt-1":len(flt_h),"二聚面":len(dim_h)}
    top=max(scores,key=scores.get)
    return top, scores

top,scores=dominant()

# ---- 追加到主报告 ----
rep=os.path.join(RES,"binding_mode_report.md")
with open(rep,"a") as f:
    f.write(f"\n---\n\n## 补充：#5 ({cid}) 结合模式\n\n")
    f.write(f"- 使用结构：`designs/candidates_6/structures/{cid}.complex.pdb`（原始 Boltz 无穿模结构，MM-GBSA ΔG=-38.3 kcal/mol）\n")
    f.write(f"- VEGF 接触残基数：**{len(rows)}**\n")
    f.write(f"- 重叠 **KDR(VEGFR2)** 位点 {len(kdr_h)} 个： `{kdr_names}`\n")
    f.write(f"- 重叠 **Flt-1(VEGFR1)** 位点 {len(flt_h)} 个： `{flt_names}`\n")
    f.write(f"- 重叠 **VEGF 二聚面** {len(dim_h)} 个： `{dim_names}`\n")
    f.write(f"- 其余 VEGF 表面接触 {len(other_h)} 个： `{oth_names}`\n\n")
    f.write(f"**解读**：#5 的主导结合区为 **{top}**（KDR {scores['KDR']} / Flt-1 {scores['Flt-1']} / 二聚面 {scores['二聚面']} 个重叠残基）。至此 6 个候选表位全部分析完毕，无未标注候选。\n")
print("[DONE] #5 binding mode ->", RES)
print("dominant:",top,scores)
