import streamlit as st
import requests
import pandas as pd
import os
import time
import re

# --- 配置区 ---
st.set_page_config(page_title="Reverse-RAG Manager", page_icon="🧬", layout="wide")

API_URL = "http://127.0.0.1:8020"
OUTPUT_DIR = "output"

# --- 样式逻辑 ---
st.markdown("""
    <style>
    .status-done { color: #28a745; font-weight: bold; }
    .status-running { color: #007bff; font-weight: bold; }
    .status-pending { color: #ffc107; font-weight: bold; }
    </style>
    """, unsafe_allow_html=True)

# --- 主界面 ---
st.title("🧬 Reverse-RAG 任务管理系统")

tab1, tab2 = st.tabs(["🚀 提交新任务", "📋 任务管理大厅"])

# --- Tab 1: 提交任务 ---
with tab1:
    with st.form("task_submission"):
        content = st.text_area("输入待处理文本：", height=300, placeholder="在此粘贴需要溯源的学术段落...")
        if st.form_submit_button("提交后台排队"):
            if content.strip():
                try:
                    r = requests.post(f"{API_URL}/submit-task", json={"content": content})
                    st.success(f"任务已提交！ID: {r.json()['task_id']}")
                    time.sleep(1)
                    st.rerun()
                except Exception as e:
                    st.error(f"连接失败: {e}")
            else:
                st.warning("内容不能为空")

# --- Tab 2: 列表式管理 ---
with tab2:
    st.subheader("所有任务状态")
    
    try:
        # 从后端获取所有任务数据
        # 注意：这里假设后端 main.py 已经增加了一个 GET /tasks 接口，如果没有，我们先尝试获取全局列表
        response = requests.get(f"{API_URL}/tasks") 
        if response.status_code == 200:
            all_tasks = response.json()
        else:
            all_tasks = {}
    except:
        st.error("无法获取任务列表，请检查后端是否运行。")
        all_tasks = {}

    if not all_tasks:
        st.info("暂无活跃任务。")
    else:
        # 将字典转换为 DataFrame 方便展示，倒序排列（最新在上）
        task_list = []
        for tid, info in all_tasks.items():
            task_list.append({
                "任务ID": tid,
                "创建时间": info.get("create_time", "-"),
                "当前状态": info.get("status", "unknown"),
                "进度": info.get("progress", "0%"),
                "文件": info.get("result_files", [])
            })
        
        df_tasks = pd.DataFrame(task_list).iloc[::-1]

        # 遍历展示
        for index, row in df_tasks.iterrows():
            with st.expander(f"ID: {row['任务ID']} | 状态: {row['当前状态']} | 时间: {row['创建时间']}", expanded=(row['当前状态'] == 'running')):
                c1, c2, c3 = st.columns([1, 2, 2])
                
                with c1:
                    st.write(f"**进度**: {row['进度']}")
                
                with c2:
                    if row['当前状态'] == 'completed' and row['文件']:
                        for f_path in row['文件']:
                            if os.path.exists(f_path):
                                with open(f_path, "rb") as f:
                                    label = "📥 下载 Output" if "output" in f_path else "📊 下载 Report"
                                    st.download_button(label, f, file_name=os.path.basename(f_path), key=f"{f_path}_{row['任务ID']}")
                    elif "failed" in row['当前状态']:
                        st.error("任务出错，请检查日志")
                    else:
                        st.write("⏳ 正在排队或处理中...")

                with c3:
                    if row['当前状态'] == 'completed':
                        # 快速预览摘要
                        try:
                            # 假设 summary.csv 里有对应记录
                            summary_df = pd.read_csv(f"{OUTPUT_DIR}/summary.csv")
                            task_summary = summary_df[summary_df['task_id'] == row['任务ID']]
                            if not task_summary.empty:
                                st.write(f"📈 命中率: {task_summary.iloc[0]['hit_rate']}")
                                st.write(f"引用数: {task_summary.iloc[0]['refs']}")
                        except:
                            pass

    if st.button("🔄 刷新列表"):
        st.rerun()

# 自动刷新逻辑
running_exists = any(t.get("status") in ["pending", "running"] for t in all_tasks.values())
if running_exists:
    time.sleep(5)
    st.rerun()