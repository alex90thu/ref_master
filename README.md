# Reverse-RAG

[![中文](https://img.shields.io/badge/中文-切换-blue)](#中文) [![English](https://img.shields.io/badge/English-Switch-brightgreen)](#english)

---

## 中文

**Reverse-RAG** 是一个“学术段落溯源”小系统：
- 前端：Streamlit 任务管理台
- 后端：FastAPI 异步队列
- 核心：
  1) 使用 Ollama LLM 从句子中抽取 PubMed 检索关键词
  2) PubMed 检索并生成引用
  3) 产出 Output 与 Report 文件，同时写入 summary.csv

### 功能
- 提交文本并排队处理
- 实时查看任务状态/进度
- 下载 Output / Report 文件
- 快速预览命中率与引用数

### 目录结构
```
.
├── app.py                 # Streamlit 前端
├── main.py                # FastAPI 后端
├── scripts/processor.py   # 任务处理逻辑
├── output/                # 产出文件与 summary.csv
├── logs/                  # 运行日志
├── run.sh                 # 一键启动脚本
├── requirement.txt        # 依赖
└── environment.yml        # Conda 依赖
```

### 环境准备
#### 方式一：Conda
```bash
conda env create -f environment.yml
conda activate reverse_rag
```

#### 方式二：pip
```bash
pip install -r requirement.txt
```

### 运行
```bash
bash ./run.sh
```

### 访问
前端应用可通过 [http://127.0.0.1:8021](http://127.0.0.1:8021) 访问。
后端 API 可通过 [http://127.0.0.1:8020](http://127.0.0.1:8020) 访问。
- 后端：127.0.0.1:8020
- 前端：0.0.0.0:8021

### 环境变量
复制 [/.env.example](/.env.example) 为 [/.env](/.env) 并按需修改：
```
# 路径配置（请确保运行用户对这些目录有读写权限）
OUTPUT_DIR="/mnt/data/rag_project/output"
LOG_DIR="/mnt/data/rag_project/logs"

OPENAI_API_BASE="http://127.0.0.1:11434/v1"
OPENAI_API_KEY=""
PUBMED_EMAIL="your_email@domain.com"
# 如果有 PubMed API Key 速度更快，没有可留空
PUBMED_API_KEY=""
```

### API 简述
- `POST /submit-task` 提交任务
- `GET /task/{tid}` 查询单任务
- `GET /tasks` 获取全部任务

---

## English

**Reverse-RAG** is a small “academic paragraph tracing” system:
- Frontend: Streamlit task dashboard
- Backend: FastAPI async queue
- Core flow:
  1) Use an Ollama LLM to extract PubMed search keywords per sentence
  2) Query PubMed and build citations
  3) Generate Output / Report files and append to summary.csv

### Features
- Submit text and process in queue
- Real-time task status/progress
- Download Output / Report files
- Quick preview of hit rate and refs count

### Structure
```
.
├── app.py                 # Streamlit frontend
├── main.py                # FastAPI backend
├── scripts/processor.py   # Task pipeline
├── output/                # Generated files + summary.csv
├── logs/                  # Runtime logs
├── run.sh                 # One-click start
├── requirement.txt        # Dependencies
└── environment.yml        # Conda env
```

### Setup
#### Option A: Conda
```bash
conda env create -f environment.yml
conda activate reverse_rag
```

#### Option B: pip
```bash
pip install -r requirement.txt
```

### Run
```bash
bash ./run.sh
```
- Backend: 127.0.0.1:8020
- Frontend: 0.0.0.0:8021

### Environment Variables
Copy [/.env.example](/.env.example) to [/.env](/.env) and adjust as needed:
```
# Path settings (ensure the runtime user has read/write permission)
OUTPUT_DIR="/mnt/data/rag_project/output"
LOG_DIR="/mnt/data/rag_project/logs"

OPENAI_API_BASE="http://127.0.0.1:11434/v1"
OPENAI_API_KEY=""
PUBMED_EMAIL="your_email@domain.com"
# PubMed API Key is optional (faster if provided)
PUBMED_API_KEY=""
```

### API
- `POST /submit-task` submit task
- `GET /task/{tid}` query one task
- `GET /tasks` list all tasks
