import os
import re
import time
import asyncio
import pandas as pd
from datetime import datetime
import httpx
from Bio import Entrez

class TaskProcessor:
    def __init__(self, ollama_url, model, email, api_key=None):
        self.ollama_url = ollama_url
        self.model = model
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        
        # 确保输出目录存在
        os.makedirs("output", exist_ok=True)
        self.summary_file = "output/summary.csv"

    async def get_keywords(self, sentence: str):
        """调用 Ollama 提取关键词并统计 Token"""
        prompt = (f"Extract 3-5 specific English search keywords for PubMed from this statement: '{sentence}'. "
                  f"Return only keywords separated by commas, no preamble.")
        
        payload = {"model": self.model, "prompt": prompt, "stream": False}
        async with httpx.AsyncClient(timeout=150.0) as client:
            try:
                resp = await client.post(self.ollama_url, json=payload)
                data = resp.json()
                tks = data.get("prompt_eval_count", 0) + data.get("eval_count", 0)
                content = re.sub(r'<think>.*?</think>', '', data.get("response", ""), flags=re.DOTALL)
                return content.strip().strip('"'), tks
            except Exception as e:
                print(f"Ollama Error: {e}")
                return "", 0

    def search_pubmed_multi(self, keywords: str, max_results=3):
        """同步检索 PubMed"""
        if not keywords: return []
        query = f"({keywords}) AND (2020:2026[pdat])"
        results = []
        try:
            with Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance") as h:
                id_list = Entrez.read(h)["IdList"]
            for pmid in id_list:
                with Entrez.esummary(db="pubmed", id=pmid) as h:
                    s = Entrez.read(h)[0]
                    results.append({
                        "title": s["Title"],
                        "author": s["LastAuthor"],
                        "year": s["PubDate"].split(' ')[0],
                        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        "bib": f"@article{{pmid{pmid},\n  title={{{s['Title']}}},\n  author={{{s['LastAuthor']} et al.}},\n  year={{{s['PubDate'].split(' ')[0]}}},\n  journal={{PubMed}}\n}}"
                    })
        except Exception as e:
            print(f"PubMed Error: {e}")
        return results

    async def run_task(self, task_id, content, tasks_db):
        """执行完整任务流"""
        tasks_db[task_id]["status"] = "running"
        start_time_stamp = time.time()
        
        sentences = [s.strip() for s in re.split(r'(?<=[。！？.!?;])', content) if s.strip()]
        total_tokens, total_refs, hit_sentences = 0, 0, 0
        results_data = []
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        for i, sent in enumerate(sentences):
            tasks_db[task_id]["progress"] = f"Sent {i+1}/{len(sentences)}"
            kw, tks = await self.get_keywords(sent)
            total_tokens += tks
            refs = self.search_pubmed_multi(kw)
            
            if refs:
                hit_sentences += 1
                total_refs += len(refs)
            
            results_data.append({"sentence": sent, "refs": refs, "keywords": kw})
            await asyncio.sleep(0.5) # 频率限制保护

        # 文件持久化
        out_path = f"output/{ts}_output.md"
        rep_path = f"output/{ts}_report.md"
        
        # 1. 生成 Output.md
        with open(out_path, "w", encoding="utf-8") as f:
            for item in results_data:
                marks = "".join([f"[{r['author']} et al., {r['year']}]" for r in item['refs']])
                f.write(f"{item['sentence']}{marks} ")

        # 2. 生成 Report.md
        duration = round(time.time() - start_time_stamp, 2)
        with open(rep_path, "w", encoding="utf-8") as f:
            f.write(f"# Reverse-RAG Report\n\n- Created: {tasks_db[task_id]['create_time']}\n")
            f.write(f"- Duration: {duration}s\n- Tokens: {total_tokens}\n")
            f.write(f"- Hit Rate: {hit_sentences}/{len(sentences)}\n\n")
            f.write("## BibTeX\n```bibtex\n" + "\n".join([r['bib'] for item in results_data for r in item['refs']]) + "\n```\n")

        # 3. 更新 Summary CSV
        new_summary = {
            "task_id": task_id, "time": ts, "duration": duration,
            "tokens": total_tokens, "hit_rate": f"{hit_sentences}/{len(sentences)}", "refs": total_refs
        }
        df = pd.DataFrame([new_summary])
        df.to_csv(self.summary_file, mode='a', index=False, header=not os.path.exists(self.summary_file))

        tasks_db[task_id].update({
            "status": "completed", "progress": "100%", "result_files": [out_path, rep_path]
        })