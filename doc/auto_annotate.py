"""
高级代谢物注释脚本
自动匹配 KEGG 和 ChEBI 数据库 ID 并填入 SBML 模型
"""

import cobra
import requests
import time
import re
from urllib.parse import quote

# ====== 配置区 ======
INPUT_MODEL  = "D:/MS-mSystems/iCG875_fixed.sbml"
OUTPUT_MODEL = "D:/MS-mSystems/iCG875_annotated.sbml"
DELAY        = 0.4       # KEGG API 限流：≥ 0.33s
MAX_RESULTS  = 1         # 每库只取最匹配的一个结果

def clean_name(met_id):
    """
    从模型内部ID推测代谢物名称
    例: M_glc__D_c → glc-D
        M_atp_c    → atp
    """
    name = re.sub(r'_[cepmnx]$', '', met_id)   # 去 compartment 后缀
    name = re.sub(r'^M_', '', name)             # 去 M_ 前缀
    name = name.replace('__', '-')              # glc__D → glc-D
    return name

def search_kegg(query):
    """KEGG Compound 搜索，返回最匹配的 KEGG ID (如 'C00031')"""
    url = f"https://rest.kegg.jp/find/compound/{quote(query)}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200 and r.text.strip():
            # 返回格式: "cpd:C00031\tD-Glucose"
            kegg_id = r.text.strip().split('\n')[0].split('\t')[0].replace('cpd:', '')
            return kegg_id
    except Exception:
        pass
    return None

def search_chebi(query):
    """ChEBI 搜索，返回最匹配的 ChEBI ID (如 'CHEBI:17234')"""
    url = f"https://www.ebi.ac.uk/chebi/webservices/rest/v1/search?query={quote(query)}&category=chebiName&searchCategory=ALL"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            import xml.etree.ElementTree as ET
            root = ET.fromstring(r.text)
            for hit in root.findall('.//listElement'):
                cid = hit.find('chebiId')
                if cid is not None and cid.text:
                    return cid.text
    except Exception:
        pass
    return None

def annotate_metabolite(met, name_query):
    """对一个代谢物执行 KEGG + ChEBI 搜索并写入 annotation"""
    annotated = 0

    # 1) KEGG
    kegg_id = search_kegg(name_query)
    if kegg_id:
        met.annotation['kegg.compound'] = f"urn:miriam:kegg.compound:{kegg_id}"
        annotated += 1

    # 2) ChEBI
    time.sleep(DELAY)
    chebi_id = search_chebi(name_query)
    if chebi_id:
        met.annotation['chebi'] = f"urn:miriam:chebi:{chebi_id}"
        annotated += 1

    return annotated

def main():
    print("正在加载模型 ...")
    model = cobra.io.read_sbml_model(INPUT_MODEL)

    total = len(model.metabolites)
    success = 0

    print(f"开始为 {total} 个代谢物搜索注释 ...\n")

    for i, met in enumerate(model.metabolites, 1):
        name_query = met.name if met.name else clean_name(met.id)
        n = annotate_metabolite(met, name_query)
        if n > 0:
            success += 1

        # 进度输出
        status = f"[{i:04d}/{total}] {met.id:20s} → "
        if n == 2:
            status += "KEGG ✓  ChEBI ✓"
        elif n == 1:
            status += "部分匹配"
        else:
            status += "未匹配"
        print(status)

    # 导出带注释的模型
    print(f"\n正在导出注释后的模型 ...")
    cobra.io.write_sbml_model(model, OUTPUT_MODEL)

    # 统计
    print(f"\n{'='*40}")
    print(f"注释完成！成功: {success}/{total} 个代谢物")
    print(f"输出文件: {OUTPUT_MODEL}")
    print(f"{'='*40}")

if __name__ == "__main__":
    main()