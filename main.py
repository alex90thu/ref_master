import os
import uuid
import asyncio
import time
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from dotenv import load_dotenv
from scripts.processor import TaskProcessor

# 加载配置
load_dotenv()

app = FastAPI(title="Reverse-RAG Pro Backend")

# 初始化处理器
processor = TaskProcessor(
    ollama_url="http://127.0.0.1:11434/api/generate",
    model="deepseek-r1:70b",
    email=os.getenv("PUBMED_EMAIL", "your_email@sjtu.edu.cn"),
    api_key=os.getenv("PUBMED_API_KEY")
)

# 全局存储与队列
tasks_db = {}
queue = asyncio.Queue()

class TaskRequest(BaseModel):
    content: str
    tag: str = "Default" # 确保 Pydantic 匹配字段名

async def worker():
    """异步执行器"""
    while True:
        task_id = await queue.get()
        try:
            await processor.run_task(
                task_id=task_id, 
                content=tasks_db[task_id]["raw_content"], 
                tasks_db=tasks_db
            )
        except Exception as e:
            tasks_db[task_id]["status"] = f"failed: {str(e)}"
            print(f"Task {task_id} Error: {e}")
        finally:
            queue.task_done()

@app.on_event("startup")
async def startup_event():
    asyncio.create_task(worker())

@app.post("/submit-task")
async def submit(req: TaskRequest):
    tid = str(uuid.uuid4())[:8]
    # 显式提取 tag，防止空字符串导致 Default 失效
    actual_tag = req.tag.strip() if req.tag and req.tag.strip() else "Default"
    
    tasks_db[tid] = {
        "raw_content": req.content,
        "tag": actual_tag,
        "status": "pending",
        "progress": "0%",
        "create_time": time.strftime("%Y-%m-%d %H:%M:%S"),
        "result_files": []
    }
    await queue.put(tid)
    return {"task_id": tid}

@app.get("/tasks")
async def get_all_tasks():
    return tasks_db

@app.get("/task/{tid}")
async def get_single_task(tid: str):
    if tid not in tasks_db:
        raise HTTPException(status_code=404, detail="Task not found")
    return tasks_db.get(tid)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8020)