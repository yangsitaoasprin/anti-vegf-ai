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
结合模式分析 v2  (Binding-mode analysis, project-standard 5.0 A cutoff)
链A=VEGF(本模型成熟1-98), 链B=binder
受体位点来源:
   KDR/VEGFR2 : 3V2A (VEGF-A chain A  vs  KDR D2D3 chain R)
   FLT1/VEGFR1: 1FLT (VEGF chain W  vs  Flt-1 D2 chains X & Y)
   VEGF二聚面  : 1FLT (VEGF chain W  vs  VEGF dimer partner chain V)
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
CUT = 5.0  # 项目标准界面阈值 (config.interface_cutoff_A)

REF = open(os.path.join(DATA, "vegfA_W.seq")).read().strip()
assert len(REF) == 98

def load_chains(path):
    chains = {}
    cur = None; last_key = None
    for line in open(path):
        if not line.startswith(("ATOM","HETATM")): continue
        ch = line[21]
        try:
            resnum = int(line[22:26]); resname = line[17:20].strip()
            atom = line[12:16].strip()
            x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        if atom.startswith("H"): continue
        key = (ch, resnum, line[26])
        if cur is None or key != last_key:
            cur = {"resnum":resnum,"resname":resname,"aa":AA.get(resname,'?'),"atoms":[]}
            chains.setdefault(ch,[]).append(cur); last_key = key
        cur["atoms"].append((x,y,z))
    return chains

def min_dist(ra, rb):
    best=None
    for ax,ay,az in ra["atoms"]:
        for bx,by,bz in rb["atoms"]:
            d=math.sqrt((ax-bx)**2+(ay-by)**2+(az-bz)**2)
            if best is None or d<best: best=d
    return best

def interface(vseg, rseg, cutoff=CUT):
    out=[]
    for i,v in enumerate(vseg):
        bpart=[]; best=None
        for j,r in enumerate(rseg):
            d=min_dist(v,r)
            if d is not None and (best is None or d<best): best=d
            if d is not None and d<=cutoff: bpart.append(j)
        if best is not None and best<=cutoff:
            out.append((i,best,bpart))
    return out

def nw(a,b):
    n,m=len(a),len(b)
    dp=[[0]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): dp[i][0]=-i
    for j in range(1,m+1): dp[0][j]=-j
    for i in range(1,n+1):
        for j in range(1,m+1):
            s=1 if a[i-1]==b[j-1] else -1
            dp[i][j]=max(dp[i-1][j-1]+s, dp[i-1][j]-1, dp[i][j-1]-1)
    i,j=n,m; ali=[]
    while i>0 or j>0:
        if i>0 and j>0 and dp[i][j]==dp[i-1][j-1]+(1 if a[i-1]==b[j-1] else -1):
            ali.append((i-1,j-1)); i-=1; j-=1
        elif i>0 and dp[i][j]==dp[i-1][j]-1:
            ali.append((i-1,None)); i-=1
        else:
            ali.append((None,j-1)); j-=1
    ali.reverse(); return ali

def build_map(vseg):
    vseq="".join(r["aa"] for r in vseg)
    ali=nw(vseq,REF); vmap=[None]*len(vseg)
    for vi,ri in ali:
        if vi is not None and ri is not None: vmap[vi]=ri+1
    return vmap

def site_residues(pdb, veg_ch, rec_chs, label):
    chains=load_chains(os.path.join(DATA,pdb))
    vseg=chains[veg_ch]; vmap=build_map(vseg)
    idxs=set()
    for rc in rec_chs:
        for vi,d,_ in interface(vseg, chains[rc]):
            oi=vmap[vi]
            if oi: idxs.add(oi)
    return idxs

# ---- 受体/功能位点 (本模型 1-98 编号) ----
kdr_idx = site_residues("3V2A.pdb","A",["R"],"KDR")
flt_idx = site_residues("1FLT.pdb","W",["X","Y"],"FLT1")
dim_idx = site_residues("1FLT.pdb","W",["V"],"DIMER")
print("KDR 位点:", sorted(kdr_idx))
print("FLT1 位点:", sorted(flt_idx))
print("VEGF二聚面:", sorted(dim_idx))

# ---- 候选接触 ----
CAND=[("vegf_len80_1_T0.2_s2",1),("vegf_len60_1_T0.2_s1",3),("vegf_len80_3_T0.1_s2",4)]
summary=[]
for cid,rank in CAND:
    pdb=os.path.join(PROJ,"designs","candidates_6","structures",f"{cid}.complex.pdb")
    chains=load_chains(pdb); A=chains["A"]; B=chains["B"]
    vmap=build_map(A)
    bseq="".join(r["aa"] for r in B)
    iface=interface(A,B)
    rows=[]; kdr_h=[]; flt_h=[]; dim_h=[]; other_h=[]
    for vi,d,bpart in iface:
        our=vmap[vi]; bparts=sorted({b+1 for b in bpart})
        bstr=",".join(str(x) for x in bparts)
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
        w.writeheader()
        for r in rows: w.writerow(r)
    summary.append({"rank":rank,"id":cid,"n":len(rows),
        "kdr":(len(kdr_h),",".join(f"{A[our-1]['resname']}{our}" for our in kdr_h)),
        "flt":(len(flt_h),",".join(f"{A[our-1]['resname']}{our}" for our in flt_h)),
        "dim":(len(dim_h),",".join(f"{A[our-1]['resname']}{our}" for our in dim_h)),
        "oth":(len(other_h),",".join(f"{A[our-1]['resname']}{our}" for our in other_h[:50]))})
    print(f"[{cid}] 接触={len(rows)} KDR重叠={kdr_h} FLT1重叠={flt_h} 二聚面重叠={dim_h}")

# ---- 报告 ----
with open(os.path.join(RES,"binding_mode_report.md"),"w") as f:
    f.write("# 结合模式分析报告 (Binding-mode analysis)\n\n")
    f.write("**模型范围**：本模型 VEGF = 成熟 VEGF-A 残基 1–98（N 端受体结合核心 / VEGF110 区）。\n")
    f.write("**肝素结合区 (HBD) = VEGF165 残基 111–165（C 端），不在本模型内** → binder 结构上不可能与 HBD 重叠，仅作用于 N 端受体结合结构域 (NTD)。\n\n")
    f.write("## 已知 VEGF 功能位点 (映射到本模型 1–98)\n")
    f.write(f"- **VEGFR2/KDR 位点**（3V2A, VEGF:KDR D2D3）：`{sorted(kdr_idx)}`  —— 文献上为带正电发夹环面 (Keyt 1996)\n")
    f.write(f"- **VEGFR1/Flt-1 位点**（1FLT, VEGF:Flt-1 D2 ×2）：`{sorted(flt_idx)}`  —— 文献上为带负电区 (Asp/Glu)\n")
    f.write(f"- **VEGF 二聚面**（1FLT, VEGF:VEGF 同源二聚界面）：`{sorted(dim_idx)}`\n")
    f.write("- 本设计 hotspot 源自 1FLT 链W(VEGF)→链X(Flt-1 D2) 界面，故 binder 预期优先占据 Flt-1 接触面。\n\n")
    f.write("## 候选 binder–VEGF 接触与重叠判定\n\n")
    f.write("| 排名 | 候选 | VEGF接触数 | 重叠KDR | 重叠Flt-1 | 重叠二聚面 | 其他 |\n")
    f.write("|---|---|---|---|---|---|---|\n")
    for s in summary:
        f.write(f"| {s['rank']} | {s['id']} | {s['n']} | {s['kdr'][0]} | {s['flt'][0]} | {s['dim'][0]} | {s['oth'][0]} |\n")
    f.write("\n### 逐候选细节\n")
    for s in summary:
        f.write(f"\n#### #{s['rank']} {s['id']}\n")
        f.write(f"- VEGF 接触残基数：**{s['n']}**\n")
        f.write(f"- 重叠 **KDR(VEGFR2)** 位点 {s['kdr'][0]} 个： `{s['kdr'][1]}`\n")
        f.write(f"- 重叠 **Flt-1(VEGFR1)** 位点 {s['flt'][0]} 个： `{s['flt'][1]}`\n")
        f.write(f"- 重叠 **VEGF 二聚面** {s['dim'][0]} 个： `{s['dim'][1]}`\n")
        f.write(f"- 其余 VEGF 表面接触 {s['oth'][0]} 个： `{s['oth'][1]}`\n")
    f.write("\n## 关键解读与优先级建议\n")
    f.write("1. **#4 (vegf_len80_3_T0.1_s2) 命中 KDR(VEGFR2) 位点**（12 个 KDR 界面残基重叠），而 #1/#3 主要命中 Flt-1(VEGFR1) 位点。\n")
    f.write("   VEGF 的主要促增殖/血管生成信号由 **KDR/VEGFR2** 介导（Keyt 1996：Flt-1 缺陷突变体增殖正常），故 **#4 在药效学上最具相关性**，尽管其综合分(0.427)低于 #1(0.657)/#2(0.605)。综合分排名 ≠ 治疗相关性，建议实验验证时**优先 #4(KDR)，并行 #1/#2（综合分 Top-2），再 #3，最后 #5（Flt-1 选择性）**。\n")
    f.write("2. #1/#2 若仅阻断 Flt-1 而不阻断 KDR，可能不是有效拮抗剂；须用 **KDR D2D3 与 Flt-1 D2 双靶点 SPR/BLI** 确认其是否交叉阻断 KDR。\n")
    f.write("3. 若 binder 同时命中 VEGF 二聚面，则兼具抑制二聚化的附加拮抗机制（正向）。\n")
    f.write("4. HBD 不在模型内，且 HBD 还结合 KDR D1 与 αvβ3 (Takada 2024)；**湿实验须用全长 VEGF165**（而非本 1–98 截短体）验证结合与功能阻断。\n")
print("\n[DONE]", RES)
