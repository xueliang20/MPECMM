import re

# 读取原始文件
with open('D:/MS-mSystems/iCG875.xml', 'r', encoding='utf-8') as f:
    content = f.read()

# 记录已经出现过的 compartment id
seen_ids = {}

def add_unique_id(match):
    attrs = match.group(1)
    # 尝试提取 compartment 名称
    nm = re.search(r'name="([^"]+)"', attrs)
    base = nm.group(1) if nm else 'compartment'
    # 替换非字母数字字符为下划线
    base = re.sub(r'[^a-zA-Z0-9]', '_', base)
    
    # 如果这个 id 已经存在，就在后面加数字（避免重复）
    if base in seen_ids:
        seen_ids[base] += 1
        uid = f"{base}_{seen_ids[base]}"
    else:
        seen_ids[base] = 0
        uid = base
        
    return f'<species id="{uid}" {attrs}>'

# 自动给没有 id 的 species 标签补上唯一的 id
new_content = re.sub(r'<species\s+(?![^>]*id=)([^>]+)>', add_unique_id, content)

# 写入修复后的新文件
with open('D:/MS-mSystems/iCG875_fixed.xml', 'w', encoding='utf-8') as f:
    f.write(new_content)

print('✅ XML 修复完成！')