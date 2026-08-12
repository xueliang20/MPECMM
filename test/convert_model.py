import cobra
import scipy.io
import numpy as np
import re
import libsbml

mat_path = 'D:/MS-mSystems/修稿/上传/iCG875.mat'

print("Step 1: Loading .mat file and building model in memory...")
mat_data = scipy.io.loadmat(mat_path)

model_var = [k for k in mat_data if k not in ['__header__', '__version__', '__globals__']][0]
model_struct = mat_data[model_var]

# 修复代谢物 ID
original_mets = model_struct['mets'][0, 0]
mets = [str(m[0]) if isinstance(m, np.ndarray) else str(m) for m in original_mets]
seen = set()
unique_mets = []
for i, met_id in enumerate(mets):
    if met_id in seen:
        unique_mets.append(f'{met_id}_dup{i}')
    else:
        seen.add(met_id)
        unique_mets.append(met_id)

# 【核心修改】强制指定正确的区室，不再盲目推断
valid_comps = {'c': 'Cytoplasm', 'e': 'Extracellular', 'p': 'Periplasmic space'}

# 构建 COBRA 模型
model = cobra.Model("iCG875")
for comp_id, comp_name in valid_comps.items():
    model.compartments[comp_id] = comp_name

# 添加代谢物，并强制将未知区室归为细胞质 'c'
met_objects = []
for m in unique_mets:
    if '[' in m and ']' in m:
        comp = m.split('[')[1].split(']')[0]
        if comp not in valid_comps:
            comp = 'c'  # 如果区室不在标准列表中，默认归为细胞质
    else:
        comp = 'c'
    met = cobra.Metabolite(m, compartment=comp)
    met_objects.append(met)

# 提取并清洗反应 ID
raw_rxns = model_struct['rxns'][0, 0]
rxns = [re.sub(r"[\[\]\'\"]", "", str(r[0]) if isinstance(r, np.ndarray) else str(r)).strip().replace(" ", "_") for r in raw_rxns]

# 提取化学计量矩阵 S
S = model_struct['S'][0, 0].toarray() if hasattr(model_struct['S'][0, 0], 'toarray') else model_struct['S'][0, 0]

# 添加反应
for i, rxn_id in enumerate(rxns):
    reaction = cobra.Reaction(rxn_id)
    met_dict = {met_objects[j]: float(coeff) for j, coeff in enumerate(S[:, i]) if coeff != 0}
    reaction.add_metabolites(met_dict)
    model.add_reactions([reaction])

print(f"Model built: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites.")

# 直接在内存中转换为 libSBML 对象并评估
print("Step 2: Converting to SBML in memory and evaluating...")
sbml_document = libsbml.SBMLDocument(2, 1)
sbml_model = sbml_document.createModel()
sbml_model.setId("iCG875")

# 添加正确的区室
for comp_id, comp_name in valid_comps.items():
    comp = sbml_model.createCompartment()
    comp.setId(comp_id)
    comp.setName(comp_name)
    comp.setConstant(True)

# 添加代谢物
for met in model.metabolites:
    species = sbml_model.createSpecies()
    species.setId(met.id)
    species.setCompartment(met.compartment)
    species.setConstant(False)
    species.setBoundaryCondition(False)

# 添加反应
for rxn in model.reactions:
    reaction = sbml_model.createReaction()
    reaction.setId(rxn.id)
    reaction.setReversible(rxn.reversibility)
    
    for met, coeff in rxn.metabolites.items():
        if coeff < 0:
            reactant = reaction.createReactant()
            reactant.setSpecies(met.id)
            reactant.setStoichiometry(-coeff)
        else:
            product = reaction.createProduct()
            product.setSpecies(met.id)
            product.setStoichiometry(coeff)

# --- 直接在内存中生成 Memote 报告 ---

print("\nStep 3: Generating Memote report directly in memory...")

try:
    # 1. 导入 memote 的核心测试模块
    from memote.suite.api import test_model
    from memote.suite.results import MemoteResult
    
    # 2. 直接在内存中对 COBRA 模型对象运行测试
    # 这会跳过所有文件读写，直接评估模型质量
    results = test_model(model)
    
    # 3. 将测试结果保存为 JSON 文件（Memote 的核心数据格式）
    result_path = 'D:/MS-mSystems/修稿/上传/memote_result.json'
    with open(result_path, 'w', encoding='utf-8') as f:
        import json
        json.dump(results, f, indent=2)
        
    print(f"SUCCESS! Memote result saved to: {result_path}")
    print("\nYou can view the detailed results by opening the JSON file,")
    print("or by running: memote report snapshot --filename 'report.html' 'D:/MS-mSystems/修稿/上传/memote_result.json'")
    
except Exception as e:
    print(f"ERROR: Failed to generate Memote report. Details: {e}")
    print("\nThis usually means memote is not installed correctly.")
    print("Please run: pip install memote")