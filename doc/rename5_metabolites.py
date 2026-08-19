import cobra
import pandas as pd

MODEL_FILE = "iCG875v11met.xml"
EXCEL_FILE = "iCG875v1011metabolites_info.xlsx"
OUTPUT_FILE = "iCG875v12met.xml"

model = cobra.io.read_sbml_model(MODEL_FILE)
df = pd.read_excel(EXCEL_FILE)

updated_count = 0
not_found = []

for idx, row in df.iterrows():
    # 去掉 M_ 前缀，得到模型中的真实 ID
    met_id = str(row["id"]).replace("M_", "", 1)
    
    if met_id in model.metabolites:
        met = model.metabolites.get_by_id(met_id)
        
        # 写入名称
        if pd.notna(row["name"]):
            met.name = row["name"]
        
        # 写入 compartment
        if "compartment" in row and pd.notna(row["compartment"]):
            met.compartment = str(row["compartment"])
        
        # 写入 SBO term
        if "sboTerm" in row and pd.notna(row["sboTerm"]):
            met.annotation["sbo"] = str(row["sboTerm"])
        
        updated_count += 1
    else:
        not_found.append(met_id)

# 保存模型
cobra.io.write_sbml_model(model, OUTPUT_FILE)

print(f"✅ 成功更新 {updated_count} 个代谢物")
print(f"💾 已保存至：{OUTPUT_FILE}")

if not_found:
    print(f"\n⚠️ 以下 {len(not_found)} 个代谢物在模型中未找到：")
    print(not_found)