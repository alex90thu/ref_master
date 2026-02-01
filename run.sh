#!/bin/bash

# 1. 显式定义目标路径（同步 .env 里的设置）
REAL_OUTPUT="/data/guozehua/ref_master/output"
REAL_LOGS="/data/guozehua/ref_master/logs"

# 2. 确保 /data 下的真实物理目录存在
mkdir -p "$REAL_OUTPUT"
mkdir -p "$REAL_LOGS"

# 3. 在当前项目目录下创建指向那里的软链接（如果不存在）
# 这样 VSCode 侧边栏就能看到正确的“绿色”映射，且不会带引号
[ -L "output" ] || ln -s "$REAL_OUTPUT" "./output"
[ -L "logs" ] || ln -s "$REAL_LOGS" "./logs"

export PYTHONPATH=$PYTHONPATH:$(pwd)

echo "正在启动 Reverse-RAG 后端服务..."
nohup python main.py > "$REAL_LOGS/backend.log" 2>&1 &
BACKEND_PID=$!

echo "正在启动 Reverse-RAG 前端服务..."
nohup streamlit run app.py --server.port 8021 --server.address 0.0.0.0 > "$REAL_LOGS/frontend.log" 2>&1 &
FRONTEND_PID=$!

echo "------------------------------------------------"
echo "软链接已建立: ./output -> $REAL_OUTPUT"
echo "后端 PID: $BACKEND_PID | 前端 PID: $FRONTEND_PID"
echo "------------------------------------------------"