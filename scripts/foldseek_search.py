#!/usr/bin/env python3
"""
Foldseek 结构相似性搜索脚本
在 MASTER 节点运行（需要联网）
用法: python foldseek_search.py --pdb <receptor.pdb> --output <output_dir> [--databases pdb100,afdb50,...]
"""
import requests
import os
import argparse
from time import sleep

ALL_DATABASES = [
    "afdb50", "BFVD", "afdb-swissprot", "afdb-proteome",
    "bfmd", "mgnify_esm30", "pdb100", "cath50", "gmgcl_id",
]

def run_foldseek(pdb_path, output_dir, databases=None):
    """Submit a PDB file to Foldseek API and download results."""
    if databases is None:
        databases = ALL_DATABASES

    url = "https://search.foldseek.com/api/ticket"

    print(f"[INFO] Submitting {pdb_path} to Foldseek...")
    print(f"[INFO] Databases: {databases}")

    with open(pdb_path, "rb") as f:
        files = {"q": f}
        data = [("mode", "3diaa")]
        for db in databases:
            data.append(("database[]", db))

        response = requests.post(url, files=files, data=data)
        ticket = response.json()
        print(f"[INFO] Ticket: {ticket}")

    # Poll until complete
    repeat = True
    while repeat:
        status = requests.get(f'https://search.foldseek.com/api/ticket/{ticket["id"]}').json()
        if status["status"] == "ERROR":
            print(f"[ERROR] Foldseek job failed for {pdb_path}")
            return False
        sleep(2)
        repeat = status["status"] != "COMPLETE"
        if repeat:
            print(f"[INFO] Waiting... status: {status['status']}")

    ticket_id = ticket["id"]
    print(f"[INFO] Job completed! Downloading results...")

    os.makedirs(output_dir, exist_ok=True)
    download_url = f'https://search.foldseek.com/api/result/download/{ticket_id}'
    tar_path = f'{output_dir}/foldseek_results.tar.gz'

    download = requests.get(download_url, stream=True)
    with open(tar_path, 'wb') as fd:
        for chunk in download.iter_content(chunk_size=128):
            fd.write(chunk)

    # Extract
    print(f"[INFO] Extracting to {output_dir}...")
    os.system(f'tar -xzf {tar_path} -C {output_dir}')
    print(f"[INFO] Done! Results in {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Foldseek structure similarity search")
    parser.add_argument("--pdb", required=True, help="Path to query PDB file")
    parser.add_argument("--output", required=True, help="Output directory for results")
    parser.add_argument("--databases", default="pdb100,afdb50,afdb-swissprot,afdb-proteome",
                        help="Comma-separated list of databases")
    args = parser.parse_args()

    databases = [db.strip() for db in args.databases.split(",")]
    run_foldseek(args.pdb, args.output, databases)
