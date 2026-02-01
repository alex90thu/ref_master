import os
import uuid
import asyncio
import time  # <--- 修复点：必须导入 time 库
from fastapi import FastAPI
from pydantic import BaseModel
from dotenv import load_dotenv
from scripts.processor import TaskProcessor

# 加载环境变量
load_dotenv()

app = FastAPI(title="Reverse-RAG Pro Backend")

# 初始化处理器
# 注意：确保你的 .env 中有 PUBMED_EMAIL
processor = TaskProcessor(
    ollama_url="http://127.0.0.1:11434/api/generate",
    model="deepseek-r1:70b",
    email=os.getenv("PUBMED_EMAIL", "your_email@sjtu.edu.cn"),
    api_key=os.getenv("PUBMED_API_KEY")
)

# 全局任务数据库和队列
tasks_db = {}
queue = asyncio.Queue()

class TaskRequest(BaseModel):
    content: str

async def worker():
    """
    后台消费者进程：从队列中获取任务 ID 并调用处理器执行
    """
    while True:
        # 获取一个任务 ID
        task_id = await queue.get()
        try:
            # 执行实际的重活
            await processor.run_task(
                task_id=task_id, 
                content=tasks_db[task_id]["raw_content"], 
                tasks_db=tasks_db
            )
        except Exception as e:
            tasks_db[task_id]["status"] = f"failed: {str(e)}"
            print(f"Error processing task {task_id}: {e}")
        finally:
            # 标记任务完成
            queue.task_done()

@app.on_event("startup")
async def startup_event():
    """
    应用启动时自动开启后台 Worker
    """
    asyncio.create_task(worker())

@app.post("/submit-task")
async def submit(req: TaskRequest):
    """
    接收任务并放入队列，立即返回任务 ID
    """
    tid = str(uuid.uuid4())[:8]
    # 初始化任务状态
    tasks_db[tid] = {
        "raw_content": req.content, 
        "status": "pending", 
        "progress": "0%", 
        "create_time": time.strftime("%Y-%m-%d %H:%M:%S"), # 这里之前报错，现已修复
        "result_files": []
    }
    
    # 放入队列排队
    await queue.put(tid)
    return {"task_id": tid}

@app.get("/task/{tid}")
async def get_task(tid: str):
    """
    供前端查询任务状态
    """
    if tid not in tasks_db:
        return {"error": "task not found"}
    return tasks_db.get(tid)

# 在 main.py 中添加：

@app.get("/tasks")
async def get_all_tasks():
    """返回所有任务的精简信息"""
    return tasks_db

if __name__ == "__main__":
    import uvicorn
    # 监听 8020 端口
    uvicorn.run(app, host="127.0.0.1", port=8020)