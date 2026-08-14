"""
基于本地 BiGG 数据库文件的高效注释脚本（终极修复版：含目标函数自动修复）
"""

import cobra
import csv
import re
from cobra.io import write_sbml_model

# ================= 配置区 =================
INPUT_MODEL  = "D:/MS-mSystems/iCG875_fixed.sbml"
OUTPUT_MODEL = "D:/MS-mSystems/iCG875_annotated_final.sbml"
BIGG_DB_FILE = "D:/MS-mSystems/bigg-met.txt"
# ==========================================

def parse_bigg_database(filepath):
    """解析 bigg-met.txt，返回 {bigg_id: {db: id, ...}}"""
    db = {}
    print(f"正在解析数据库文件: {filepath} ...")
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            bigg_id = row.get('bigg_id', '').strip()
            if not bigg_id:
                continue
            links = row.get('database_links', '')
            annotations = {}
            for link in links.split(', '):
                if ':' in link:
                    k, v = link.split(':', 1)
                    annotations[k.strip().lower()] = v.strip()
            db[bigg_id] = annotations
    print(f"数据库加载完成，共包含 {len(db)} 个代谢物条目。\n")
    return db

def clean_metabolite_id(met_id):
    """
    清理代谢物 ID：
    1. 去掉首尾的 '_'
    2. 去掉常见的区室后缀（如 _c, _e, _p, _m 等）
    3. 再次去掉尾部可能残留的 '_'
    """
    # 1. 先去掉首尾的下划线
    clean_id = met_id.strip('_')
    
    # 2. 去掉区室后缀（如 _c, _e, _p, _m, _b, _w, _x, _r, _l, _s）
    #    使用正则表达式匹配：以一个 '_' 结尾，前面是一个小写字母
    clean_id = re.sub(r'_[a-z]$', '', clean_id)
    
    # 3. 再次去掉可能残留的尾部下划线
    clean_id = clean_id.rstrip('_')
    
    return clean_id

def annotate_model(model, db):
    """使用 BiGG 数据库注释模型中的代谢物，并添加 SBO 术语"""
    metabolites = model.metabolites
    total = len(metabolites)
    success_count = 0
    
    print(f"开始为 {total} 个代谢物匹配注释 ...\n")
    
    for i, met in enumerate(metabolites, 1):
        original_id = met.id
        # 清理 ID
        clean_id = clean_metabolite_id(original_id)
        
        match_found = False
        
        # 策略 1：直接用清理后的 ID 去匹配
        if clean_id in db:
            match_found = True
            annotations = db[clean_id]
        else:
            # 策略 2：如果失败，尝试全小写匹配
            clean_id_lower = clean_id.lower()
            for db_key, db_annotations in db.items():
                if db_key.lower() == clean_id_lower:
                    match_found = True
                    annotations = db_annotations
                    break
        
        if match_found:
            success_count += 1
            # 写入注释
            for db_name, db_id in annotations.items():
                if db_name and db_id:
                    met.annotation[f"{db_name}"] = f"{db_id}"
            # 添加 SBO 术语 (SBO:0000247 代表 simple chemical)
            met.annotation["sbo"] = "SBO:0000247"
            print(f"[{i:04d}/{total}] {original_id:<20} -> √ 匹配 (clean: {clean_id})")
        else:
            print(f"[{i:04d}/{total}] {original_id:<20} -> 未找到 (clean: {clean_id})")
    
    print(f"\n匹配完成！成功: {success_count}/{total}，未找到: {total - success_count}/{total}")
    return model

if __name__ == '__main__':
    # 1. 加载模型
    print("正在加载模型 ...")
    model = cobra.io.read_sbml_model(INPUT_MODEL)
    print(f"模型加载完成，共 {len(model.metabolites)} 个代谢物。\n")

    # ================= 关键修复：自动设置目标函数 =================
    # 自动寻找含有 'BIOMASS' 的反应，并将其设为目标函数
    biomass_reactions = [r for r in model.reactions if 'BIOMASS' in r.id.upper()]
    if biomass_reactions:
        model.objective = biomass_reactions[0].id
        print(f"✅ 已成功设置目标函数: {biomass_reactions[0].id}\n")
    else:
        print("⚠️  警告：未在模型中找到 BIOMASS 反应，无法自动设置目标函数！\n")
    # ================================================================

    # 2. 加载数据库并注释
    bigg_db = parse_bigg_database(BIGG_DB_FILE)
    annotated_model = annotate_model(model, bigg_db)

    # 3. 导出模型
    print("\n正在导出注释后的模型 ...")
    write_sbml_model(annotated_model, OUTPUT_MODEL)
    print(f"✅ 完成！输出文件: {OUTPUT_MODEL}")