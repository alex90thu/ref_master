import os
import re
import time
import asyncio
import hashlib
import json
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
        
        # 核心修复：直接读取 .env 绝对路径，不依赖相对路径
        self.output_dir = os.getenv("OUTPUT_DIR", "/data/guozehua/ref_master/output")
        os.makedirs(self.output_dir, exist_ok=True)
        
        self.pubmed_semaphore = asyncio.Semaphore(3) 
        
        # 缓存文件强制锁定在绝对路径下的 output 目录
        self.cache_file = os.path.join(self.output_dir, "pubmed_cache.json")
        self.cache = self._load_cache()
        print(f"DEBUG: Cache file path is {self.cache_file}")

    def _load_cache(self):
        # 增加容错检查
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    return data if isinstance(data, dict) else {}
            except Exception as e:
                print(f"Cache load error: {e}")
                return {}
        return {}

    def _save_cache(self):
        try:
            with open(self.cache_file, 'w', encoding='utf-8') as f:
                json.dump(self.cache, f, ensure_ascii=False, indent=2)
        except Exception as e:
            print(f"Cache save error: {e}")

    def latex_escape(self, text):
        conv = {
            '&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#', '_': r'\_',
            '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}', '^': r'\^{}',
        }
        regex = re.compile('|'.join(re.escape(str(key)) for key in sorted(conv.keys(), key=lambda item: -len(item))))
        return regex.sub(lambda mo: conv[mo.group()], text)

    async def get_keywords(self, sentence: str):
        prompt = (f"Extract 3-5 specific English search keywords for PubMed from: '{sentence}'. "
                  f"Return only keywords separated by commas.")
        payload = {"model": self.model, "prompt": prompt, "stream": False}
        async with httpx.AsyncClient(timeout=150.0) as client:
            try:
                resp = await client.post(self.ollama_url, json=payload)
                data = resp.json()
                content = re.sub(r'<think>.*?</think>', '', data.get("response", ""), flags=re.DOTALL)
                return content.strip().strip('"'), (data.get("prompt_eval_count", 0) + data.get("eval_count", 0))
            except: return "", 0

    async def search_pubmed_safe(self, keywords: str, max_results=3):
        if not keywords: return []
        kw_hash = hashlib.md5(keywords.encode()).hexdigest()
        if kw_hash in self.cache: return self.cache[kw_hash]

        async with self.pubmed_semaphore:
            query = f"({keywords}) AND (2020:2026[pdat])"
            results = []
            try:
                await asyncio.sleep(0.3)
                with Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance") as h:
                    id_list = Entrez.read(h)["IdList"]
                for pmid in id_list:
                    with Entrez.esummary(db="pubmed", id=pmid) as h:
                        s = Entrez.read(h)[0]
                        authors = s.get("AuthorList", [])
                        full_authors = " and ".join(authors) if authors else s.get("LastAuthor", "Anon")
                        journal = s.get("Source", "Unknown Journal")
                        doi = ""
                        eloc = s.get("elocationid", "")
                        if "doi:" in eloc.lower(): doi = re.sub(r'(?i)doi:\s*', '', eloc).strip()
                        elif "ArticleIds" in s and isinstance(s["ArticleIds"], dict): doi = s["ArticleIds"].get("doi", "")
                        
                        year = s.get("PubDate", "2026").split(' ')[0]
                        bib_key = f"pmid{pmid}"
                        results.append({
                            "id": bib_key, "title": s.get("Title", "No Title"), "author": full_authors,
                            "year": year, "journal": journal, "doi": doi,
                            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                            "bib": (f"@article{{{bib_key},\n  title={{{s.get('Title', 'No Title')}}},\n"
                                    f"  author={{{full_authors}}},\n  year={{{year}}},\n"
                                    f"  journal={{{journal}}},\n  doi={{{doi}}}\n}}")
                        })
                self.cache[kw_hash] = results
                self._save_cache()
                return results
            except: return []

    async def run_task(self, task_id, content, tasks_db):
        tasks_db[task_id]["status"] = "running"
        start_ts = time.time()
        sentences = [s.strip() for s in re.split(r'(?<=[。！？.!?;])', content) if s.strip()]
        total_tokens, total_refs, results_data = 0, 0, []

        for i, sent in enumerate(sentences):
            tasks_db[task_id]["progress"] = f"Sent {i+1}/{len(sentences)}"
            kw, tks = await self.get_keywords(sent)
            total_tokens += tks
            refs = await self.search_pubmed_safe(kw)
            total_refs += len(refs)
            results_data.append({"sentence": sent, "refs": refs, "keywords": kw})

        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_path = os.path.join(self.output_dir, f"{ts}_output.md")
        rep_path = os.path.join(self.output_dir, f"{ts}_report.md")

        with open(out_path, "w", encoding="utf-8") as f:
            for item in results_data:
                # 对正文进行 LaTeX 转义
                safe_sent = self.latex_escape(item['sentence'])
                
                # 寻找句末标点符号（。！？.!?;）
                punc_match = re.search(r'([。！？.!?;])$', safe_sent)
                
                if item['refs']:
                    cite_keys = ",".join([r['id'] for r in item['refs']])
                    cite_str = f"\\cite{{{cite_keys}}}"
                    
                    if punc_match:
                        # 关键修复：将引用插在标点符号前面
                        punc = punc_match.group(1)
                        sentence_body = safe_sent[:-1]
                        f.write(f"{sentence_body}{cite_str}{punc} ")
                    else:
                        f.write(f"{safe_sent}{cite_str} ")
                else:
                    f.write(f"{safe_sent} ")

        with open(rep_path, "w", encoding="utf-8") as f:
            f.write(f"# Reverse-RAG Report\n- Tag: {tasks_db[task_id].get('tag')}\n\n## BibTeX\n```bibtex\n")
            seen = set()
            for item in results_data:
                for r in item['refs']:
                    if r['id'] not in seen:
                        f.write(r['bib'] + "\n")
                        seen.add(r['id'])
            f.write("```\n")

        # 统计数据存入迁移后的路径
        summary_csv = os.path.join(self.output_dir, "summary.csv")
        hit_rate = f"{sum(1 for x in results_data if x['refs'])}/{len(sentences)}"
        pd.DataFrame([{
            "task_id": task_id, "tag": tasks_db[task_id].get('tag'), "time": ts,
            "duration": round(time.time()-start_ts, 2), "tokens": total_tokens,
            "hit_rate": hit_rate, "refs": total_refs
        }]).to_csv(summary_csv, mode='a', index=False, header=not os.path.exists(summary_csv))

        tasks_db[task_id].update({"status": "completed", "progress": "100%", "result_files": [out_path, rep_path]})