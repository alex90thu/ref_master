import streamlit as st
import requests
import pandas as pd
import os
import time
import hashlib
import zipfile
from io import BytesIO

st.set_page_config(page_title="Reverse-RAG Manager", page_icon="🧬", layout="wide")
API_URL = "http://127.0.0.1:8020"

# 高对比度配色
COLOR_MAP = {"blue": "#E3F2FD", "green": "#F1F8E9", "orange": "#FFF3E0", "red": "#FCE4EC", "purple": "#F3E5F5", "teal": "#E0F2F1"}

def get_tag_style(tag):
    colors = list(COLOR_MAP.values())
    return colors[int(hashlib.md5(tag.encode()).hexdigest(), 16) % len(colors)]

st.title("👍 你说得对 Here are your refs")

tab1, tab2 = st.tabs(["🚀 提交新任务", "📋 任务管理大厅"])

# --- Tab 1: 提交 ---
with tab1:
    with st.form("task_submission"):
        task_tag = st.text_input("任务标签 (Tag)：", value="Default")
        content = st.text_area("待处理文本：", height=300)
        if st.form_submit_button("部署后台队列"):
            if content.strip():
                try:
                    requests.post(f"{API_URL}/submit-task", json={"content": content, "tag": task_tag})
                    st.success(f"任务 [{task_tag}] 提交成功！")
                    time.sleep(0.5)
                    st.rerun()
                except Exception as e: st.error(f"连接失败: {e}")

# --- Tab 2: 列表管理 ---
with tab2:
    with st.expander("📝 Overleaf 配置模板 (点击右上角图标复制)", expanded=False):
        st.markdown("""
        **使用步骤：**
        1. 编译器选 **XeLaTeX**。
        2. 新建 **refs.bib**，粘贴 Report 中的 BibTeX 内容。
        3. 复制下方代码到 **main.tex**。
        """)
        try:
            with open("main.tex", "r", encoding="utf-8") as f:
                # 使用 st.code 渲染，它自带官方的复制按钮
                st.code(f.read(), language="latex")
        except: st.warning("根目录下未找到 main.tex")

    try:
        all_tasks = requests.get(f"{API_URL}/tasks").json()
    except: all_tasks = {}

    if all_tasks:
        for tid, info in sorted(all_tasks.items(), key=lambda x: x[1]['create_time'], reverse=True):
            tag = info.get("tag", "Default")
            bg_color = get_tag_style(tag)
            
            with st.container():
                st.markdown(f"""
                    <div style="background-color:{bg_color}; padding:12px; border-radius:10px; border-left:8px solid #555; margin-bottom:5px;">
                        <h4 style="margin:0;">🏷️ {tag} | <small>ID: {tid}</small></h4>
                        <p style="margin:0; font-size:0.9rem;">状态: <b>{info['status']}</b> | 进度: {info['progress']} | 时间: {info['create_time']}</p>
                    </div>
                """, unsafe_allow_html=True)
                
                with st.expander("任务操作与统计"):
                    c1, c2 = st.columns([1, 1])
                    
                    with c1:
                        if info['status'] == 'completed' and info['result_files']:
                            # 创建内存中的 ZIP 文件
                            zip_buffer = BytesIO()
                            with zipfile.ZipFile(zip_buffer, "w") as zf:
                                for f_path in info['result_files']:
                                    if os.path.exists(f_path):
                                        zf.write(f_path, os.path.basename(f_path))
                            
                            st.download_button(
                                label="📦 一键下载结果压缩包 (ZIP)",
                                data=zip_buffer.getvalue(),
                                file_name=f"RAG_{tag}_{tid}.zip",
                                mime="application/zip",
                                key=f"dl_zip_{tid}"
                            )
                        else:
                            st.write("⏳ 任务排队中或正在处理...")

                    with c2:
                        if info['status'] == 'completed':
                            try:
                                sdf = pd.read_csv("output/summary.csv")
                                task_s = sdf[sdf['task_id'] == tid]
                                if not task_s.empty:
                                    st.write(f"📈 命中率: {task_s.iloc[0]['hit_rate']}")
                                    st.write(f"📚 引用总数: {task_s.iloc[0]['refs']}")
                            except: pass
            st.divider()

if any(t.get("status") in ["pending", "running"] for t in all_tasks.values()):
    time.sleep(5)
    st.rerun()