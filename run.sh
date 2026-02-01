#!/bin/bash

# 确保在项目根目录下创建必要的文件夹
mkdir -p logs output scripts

# 将当前目录加入 Python 搜索路径，确保 main.py 能找到 scripts.processor
export PYTHONPATH=$PYTHONPATH:$(pwd)

echo "正在启动 Reverse-RAG 后端服务 (127.0.0.1:8020)..."
# 使用 nohup 运行后端，并将 PID 存入变量
nohup python main.py > logs/backend.log 2>&1 &
BACKEND_PID=$!

# 检查进程是否真的跑起来了（简单校验）
if ps -p $BACKEND_PID > /dev/null
then
    echo "后端启动成功，PID: $BACKEND_PID"
else
    echo "后端启动失败，请检查 logs/backend.log"
fi

echo "正在启动 Reverse-RAG 前端服务 (0.0.0.0:8021)..."
# 使用 nohup 运行前端，并将 PID 存入变量
nohup streamlit run app.py --server.port 8021 --server.address 0.0.0.0 > logs/frontend.log 2>&1 &
FRONTEND_PID=$!

if ps -p $FRONTEND_PID > /dev/null
then
    echo "前端启动成功，PID: $FRONTEND_PID"
else
    echo "前端启动失败，请检查 logs/frontend.log"
fi

echo "------------------------------------------------"
echo "服务已在后台运行："
echo "- 后端日志: tail -f logs/backend.log"
echo "- 前端日志: tail -f logs/frontend.log"
echo "停止服务请执行: kill $BACKEND_PID $FRONTEND_PID"
echo "------------------------------------------------"