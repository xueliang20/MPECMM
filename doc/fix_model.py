import scipy.io
import cobra
from cobra import Model, Reaction, Metabolite
import re
import numpy as np

def clean_id(name):
    if name is None: return None
    if isinstance(name, bytes): name = name.decode('utf-8', errors='ignore')
    elif not isinstance(name, str): name = str(name)
    name = name.strip()
    if not name: return None
    name = re.sub(r'[^a-zA-Z0-9_]', '_', name)
    if name[0].isdigit(): name = 'M_' + name
    return name

def get_cell_item(arr, idx):
    try:
        flat = np.asarray(arr).flatten()
        val = flat[idx]
        if isinstance(val, bytes): return val.decode('utf-8', errors='ignore').strip()
        if isinstance(val, np.ndarray):
            if val.dtype == 'object' and val.size > 0:
                v = val.flatten()[0]
                return v.decode('utf-8', errors='ignore').strip() if isinstance(v, bytes) else str(v).strip()
            return str(val).strip()
        return str(val).strip()
    except:
        return ""

def dedup(ids_list):
    seen = {}
    result = []
    for old_id in ids_list:
        nid = clean_id(old_id)
        if not nid: 
            result.append(None)
            continue
        if nid in seen:
            seen[nid] += 1
            nid = f"{nid}_{seen[nid]}"
        else:
            seen[nid] = 0
        result.append(nid)
    return result

# ================= 主流程 =================
print("正在加载 iCG875.mat ...")
data = scipy.io.loadmat('D:/MS-mSystems/iCG875.mat', struct_as_record=False, squeeze_me=True)
model_data = data['model']

mets = [get_cell_item(model_data.mets, i) for i in range(model_data.mets.size)]
rxns = [get_cell_item(model_data.rxns, i) for i in range(model_data.rxns.size)]

# 兼容处理：如果 S 已经是普通 numpy 数组，就直接用；如果是稀疏矩阵，才调用 todense()
S_raw = model_data.S
if hasattr(S_raw, 'todense'):
    S = np.asarray(S_raw.todense())
else:
    S = np.asarray(S_raw)

lbs = np.array(model_data.lb).flatten()
ubs = np.array(model_data.ub).flatten()
rs = np.array(model_data.rs).flatten() if hasattr(model_data, 'rs') else np.zeros(len(rxns))

clean_mets = dedup(mets)
clean_rxns = dedup(rxns)

valid_mets = [m for m in clean_mets if m]
valid_rxns = [r for r in clean_rxns if r]

print(f"有效代谢物: {len(valid_mets)}, 有效反应: {len(valid_rxns)}")

model = Model("iCG875_fixed")

# 创建代谢物
met_dict = {}
for mid in valid_mets:
    m = Metabolite(id=mid)
    met_dict[mid] = m
model.add_metabolites(list(met_dict.values()))

# 创建反应
rxn_list = []
for i, rid in enumerate(valid_rxns):
    if rid is None: continue
    rxn = Reaction(id=rid)
    rxn.lower_bound = float(lbs[i]) if i < len(lbs) else -1000.0
    rxn.upper_bound = float(ubs[i]) if i < len(ubs) else 1000.0
    
    # 添加代谢物计量关系
    col = S[:, i]
    met_stoich = {}
    for j, coeff in enumerate(col):
        if coeff != 0 and j < len(clean_mets) and clean_mets[j] in met_dict:
            met_stoich[met_dict[clean_mets[j]]] = float(coeff)
    rxn.add_metabolites(met_stoich)
    
    # 目标系数（生物质等）
    if i < len(rs) and rs[i] != 0:
        rxn.objective_coefficient = float(rs[i])
    
    rxn_list.append(rxn)

print("正在将反应加入模型...")
model.add_reactions(rxn_list)
print(f"成功创建并添加了 {len(rxn_list)} 个反应！")

# 兜底：确保所有 metabolite 都有合法的 compartment
for met in model.metabolites:
    if not met.compartment or not str(met.compartment).strip():
        met.compartment = 'c'
        
# 导出 SBML
out_path = 'D:/MS-mSystems/iCG875_fixed.sbml'
print(f"正在导出 SBML 到: {out_path} ...")
cobra.io.write_sbml_model(model, out_path)
print("完成！")