import cobra
import pandas as pd

# === 路径配置 ===
MODEL_PATH = r'D:\MS-mSystems\iCG875_mass_fixed_v6.sbml'
OUTPUT_PATH = r'D:\MS-mSystems\iCG875_with_gpr.sbml'
GPR_TABLE = r'D:\MS-mSystems\generxn.txt'

# 加载模型
model = cobra.io.read_sbml_model(MODEL_PATH)
print(f'模型加载成功：{len(model.metabolites)} 个代谢物，{len(model.reactions)} 个反应。')

# 读取 txt 文件，跳过表头
df = pd.read_csv(GPR_TABLE, sep='\t', dtype=str, header=0)
print(f'读取到 {len(df)} 条反应-基因对应关系。')

# 遍历写入 GPR
success, skipped, not_found = 0, 0, 0
for idx, row in df.iterrows():
    rxn_id = str(row[0]).strip()
    gene_rule = str(row[1]).strip()
    
    if rxn_id == 'nan' or gene_rule == 'nan' or not gene_rule:
        skipped += 1
        continue
    
    try:
        reaction = model.reactions.get_by_id(rxn_id)
        reaction.gene_reaction_rule = gene_rule
        success += 1
    except KeyError:
        not_found += 1
        print(f'  ⚠️ 未找到反应: {rxn_id}')

# 保存模型
cobra.io.write_sbml_model(model, OUTPUT_PATH)
print(f'\n完成：成功 {success} | 跳过 {skipped} | 未找到 {not_found}')
print(f'已保存至：{OUTPUT_PATH}')