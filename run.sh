#!/bin/bash

# 加载 .env 环境变量（简单解析）
if [ -f .env ]; then
    export $(cat .env | grep -v '#' | awk '/=/ {print $1}')
fi

# 使用环境变量定义的目录，若未定义则使用默认值 [cite: 2026-01-14]
LOG_DIR=${LOG_DIR:-"logs"}
OUTPUT_DIR=${OUTPUT_DIR:-"output"}

# 自动创建目录
mkdir -p "$LOG_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p scripts

export PYTHONPATH=$PYTHONPATH:$(pwd)

echo "正在启动 Reverse-RAG 后端服务 (127.0.0.1:8020)..."
nohup python main.py > "$LOG_DIR/backend.log" 2>&1 &
BACKEND_PID=$!
echo "后端启动成功，PID: $BACKEND_PID"

echo "正在启动 Reverse-RAG 前端服务 (0.0.0.0:8021)..."
nohup streamlit run app.py --server.port 8021 --server.address 0.0.0.0 > "$LOG_DIR/frontend.log" 2>&1 &
FRONTEND_PID=$!
echo "前端启动成功，PID: $FRONTEND_PID"

echo "------------------------------------------------"
echo "日志目录: $LOG_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "后端日志: tail -f $LOG_DIR/backend.log"
echo "前端日志: tail -f $LOG_DIR/frontend.log"
echo "停止服务: kill $BACKEND_PID $FRONTEND_PID"
echo "------------------------------------------------"