import cobra
import pandas as pd

# ================= 1. 路径配置区 =================
model_path = "iCG875v16G.xml"          # 你的原始 SBML 模型文件
excel_path = "reaction_list.xlsx"         # 你提供的 Excel 表格
output_path = "iCG875v17R.xml"            # 修改后保存的新模型文件
# ==================================================

# 1. 加载模型和表格
print("正在加载模型与表格，请稍候...")
model = cobra.io.read_sbml_model(model_path)
df = pd.read_excel(excel_path, engine='openpyxl')

# 统计修改次数
count_name = 0
count_sbo = 0
count_type = 0

print(f"✅ 模型加载成功，共 {len(model.reactions)} 个反应。")
print(f"✅ 表格读取成功，共 {len(df)} 行记录。\n")

# 2. 遍历表格，依照反应ID修改模型
for index, row in df.iterrows():
    rxn_id = row["反应ID"]
    
    # 检查模型中是否存在该反应
    try:
        rxn = model.reactions.get_by_id(rxn_id)
    except KeyError:
        print(f"⚠️  警告：模型中找不到反应 ID: {rxn_id}，已跳过。")
        continue

    # (1) 修改反应名称 (Name)
    new_name = str(row["反应名称"]) if pd.notna(row["反应名称"]) else ""
    if rxn.name != new_name:
        rxn.name = new_name
        count_name += 1

    # (2) 修改 SBO 编号
    new_sbo = str(row["SBO编号"]) if pd.notna(row["SBO编号"]) else ""
    # 确保 SBO 格式统一（有些表里可能没写 SBO: 前缀）
    if new_sbo and not new_sbo.startswith("SBO:"):
        new_sbo = f"SBO:{new_sbo}"
        
    if rxn.annotation.get("sbo") != new_sbo:
        rxn.annotation["sbo"] = new_sbo
        count_sbo += 1

    # (3) 修改反应类型 (写入到注释/子系统/SubSystem 中)
    if "反应类型" in row and pd.notna(row["反应类型"]):
        new_type = str(row["反应类型"])
        rxn.annotation["subsystem"] = new_type
        rxn.annotation["reaction_type"] = new_type # 备用字段
        count_type += 1

    # (4) 修改上下界与可逆性 (可选，如果需要严格同步)
    if pd.notna(row["下界(LB)"]):
        rxn.lower_bound = float(row["下界(LB)"])
    if pd.notna(row["上界(UB)"]):
        rxn.upper_bound = float(row["上界(UB)"])
    if pd.notna(row["可逆"]):
        # 根据表格的 TRUE/FALSE 自动调整边界和可逆属性
        if str(row["可逆"]).upper() == "TRUE":
            rxn.reversibility = True
        else:
            rxn.reversibility = False

# 3. 保存修改后的模型
print("正在保存修改后的模型...")
cobra.io.write_sbml_model(model, output_path)

# 4. 输出修正报告
print("\n" + "="*30)
print("🎉 模型修正完成！")
print(f"修改反应名称：{count_name} 个")
print(f"修改 SBO 编号：{count_sbo} 个")
print(f"修改反应类型：{count_type} 个")
print(f"模型已保存至：{output_path}")
print("="*30)