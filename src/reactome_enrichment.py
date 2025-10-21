"""
Reactome 通路富集（AnalysisService API）并下载叠加分析的通路图
================================================================

功能
- 读取基因列表（CSV 含 gene 列或纯文本一行一个基因）
- 调用 Reactome AnalysisService 提交标识符，获取 token
- 拉取通路富集结果（与官网一致的统计口径，包含 P 值/FDR 等）
- 将显著通路保存为 CSV
- 下载通路图（叠加本次分析结果的着色），与官网导出一致

用法
  python src/reactome_enrichment.py --input_file results/brm_analysis_csv/significant_DEGs.csv \
      --analysis_name SMARCA2_analysis --output_dir results/reactome_enrichment --top_n 20

依赖
- requests, pandas, numpy

说明
- Reactome AnalysisService 文档：reactome.org（REST API）。本脚本采用公开端点：
  * POST  /AnalysisService/identifiers  提交基因，返回 token
  * GET   /AnalysisService/token/{token}            获取整体摘要（含映射统计）
  * GET   /AnalysisService/token/{token}/pathways   获取通路结果表
  * GET   /ContentService/exporter/diagram/{stId}.png?token={token}&resource=TOTAL  下载叠加通路图

"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import List, Dict, Any

import requests
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]


ANALYSIS_BASE = "https://reactome.org/AnalysisService"
CONTENT_BASE = "https://reactome.org/ContentService"


def setup_logging(output_dir: Path) -> logging.Logger:
    logger = logging.getLogger("ReactomeEnrichment")
    logger.setLevel(logging.INFO)
    # 清空旧 handler（避免在多次运行时重复）
    for h in logger.handlers[:]:
        logger.removeHandler(h)
    output_dir.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(output_dir / "log.txt", mode="w", encoding="utf-8")
    ch = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', '%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def sanitize_filename(name: str) -> str:
    name = re.sub(r"[\\/:*?\"<>|]", "_", name)
    name = re.sub(r"\s+", " ", name).strip()
    return name[:180]


def load_gene_list(input_file: Path, logger: logging.Logger) -> List[str]:
    if not input_file.exists():
        raise FileNotFoundError(f"输入文件不存在: {input_file}")

    genes: List[str] = []
    if input_file.suffix.lower() in {".csv", ".tsv"}:
        sep = "," if input_file.suffix.lower() == ".csv" else "\t"
        df = pd.read_csv(input_file, sep=sep)
        col = None
        for cand in ["gene", "Gene", "symbol", "Symbol", "GENE"]:
            if cand in df.columns:
                col = cand
                break
        if col is None:
            raise ValueError("未找到 gene 列。请确保文件包含列名 'gene'.")
        genes = df[col].astype(str).tolist()
    else:
        # 纯文本：每行一个标识符
        genes = [line.strip() for line in input_file.read_text(encoding="utf-8").splitlines() if line.strip()]

    # 清洗：处理 GENE|ID、?|ID 等格式，统一大写
    cleaned: List[str] = []
    for g in genes:
        if g.startswith("?|"):
            # 保留 ID 或直接丢弃前缀：Reactome 支持多类型标识符，这里取 ID
            parts = g.split("|")
            cleaned.append(parts[1].strip().upper())
        elif "|" in g:
            cleaned.append(g.split("|")[0].strip().upper())
        else:
            cleaned.append(g.strip().upper())

    # 去重
    unique_genes = sorted(list({x for x in cleaned if x}))
    logger.info(f"载入基因：原始 {len(genes)}，清洗去重后 {len(unique_genes)}")
    return unique_genes


def submit_identifiers(genes: List[str], species: str, logger: logging.Logger) -> Dict[str, Any]:
    url = f"{ANALYSIS_BASE}/identifiers/?interactors=false&species={requests.utils.quote(species)}"
    headers = {"Content-Type": "text/plain"}
    payload = "\n".join(genes)
    logger.info("提交标识符至 Reactome AnalysisService…")
    r = requests.post(url, headers=headers, data=payload.encode("utf-8"), timeout=120)
    r.raise_for_status()
    data = r.json()
    token = data.get("summary", {}).get("token") or data.get("token")
    if not token:
        raise RuntimeError("未从 Reactome 获取 token。响应中缺少 token 字段。")
    logger.info(f"获取 token: {token}")
    return {"token": token, "raw": data}


def get_analysis_summary(token: str, logger: logging.Logger) -> Dict[str, Any]:
    url = f"{ANALYSIS_BASE}/token/{token}"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    data = r.json()
    # 尝试读取识别统计
    not_found = data.get("summary", {}).get("notFound", [])
    found = data.get("summary", {}).get("found", [])
    # 有时返回的是 numbers
    found_size = data.get("summary", {}).get("found", {}).get("size") if isinstance(found, dict) else len(found)
    not_found_size = data.get("summary", {}).get("notFound", {}).get("size") if isinstance(not_found, dict) else len(not_found)
    logger.info("标识符映射统计（Reactome 口径）：")
    logger.info(f"  - 成功映射: {found_size}")
    logger.info(f"  - 未映射:   {not_found_size}")
    if isinstance(not_found, list) and 0 < len(not_found) <= 20:
        logger.info("  - 未映射列表: " + ", ".join(map(str, not_found)))
    return data


def get_pathway_results(token: str, logger: logging.Logger, page_size: int = 2000) -> pd.DataFrame:
    url = f"{ANALYSIS_BASE}/token/{token}/pathways?resource=TOTAL&order=ASC&sortBy=ENTITIES_FDR&pageSize={page_size}"
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    items = r.json() or []
    logger.info(f"通路结果条目数: {len(items)}")
    # 解析为 DataFrame
    rows = []
    for it in items:
        st_id = it.get("stId")
        name = it.get("name")
        species = it.get("speciesName") or it.get("species")
        entities = it.get("entities", {})
        p_value = entities.get("pValue")
        fdr = entities.get("fdr") or entities.get("FDR")
        found = entities.get("found")
        total = entities.get("total")
        rows.append({
            "stId": st_id,
            "name": name,
            "species": species,
            "entities_found": found,
            "entities_total": total,
            "entities_pValue": p_value,
            "entities_fdr": fdr,
        })
    df = pd.DataFrame(rows)
    return df


def get_pathways_via_reactome2py(genes: List[str], species: str, logger: logging.Logger):
    """
    备用方案：通过 reactome2py 直接获取 pathways 列表（无需手动拼 REST 子端点）。
    返回与 get_pathway_results 相同列。
    """
    try:
        from reactome2py import analysis as rx_analysis
    except Exception as e:
        raise RuntimeError("reactome2py 未安装，无法使用备用方案。请先 pip install reactome2py") from e

    logger.info("使用 reactome2py 备用方案获取通路结果…")
    gene_str = ",".join(genes)
    res = rx_analysis.identifiers(gene_str, species=species)
    # reactome2py 新版本直接返回 dict，旧版本返回 JSON 字符串
    try:
        import json as _json
        res_obj = res if isinstance(res, dict) else _json.loads(res)
    except Exception:
        res_obj = res
    items = res_obj.get("pathways", [])

    rows = []
    for it in items:
        st_id = it.get("stId")
        name = it.get("name")
        species_name = it.get("speciesName") or it.get("species")
        entities = it.get("entities", {})
        p_value = entities.get("pValue")
        fdr = entities.get("fdr") or entities.get("FDR")
        found = entities.get("found")
        total = entities.get("total")
        rows.append({
            "stId": st_id,
            "name": name,
            "species": species_name,
            "entities_found": found,
            "entities_total": total,
            "entities_pValue": p_value,
            "entities_fdr": fdr,
        })
    df = pd.DataFrame(rows)
    token = None
    try:
        token = res_obj.get("summary", {}).get("token")
    except Exception:
        token = None
    return df, token


def download_pathway_diagram(token: str, st_id: str, out_png: Path, logger: logging.Logger) -> None:
    # 叠加分析结果的通路图导出端点
    url = (
        f"{CONTENT_BASE}/exporter/diagram/{st_id}.png?"
        f"quality=7&diagramProfile=Modern&analysisProfile=Standard&"
        f"token={token}&resource=TOTAL"
    )
    r = requests.get(url, stream=True, timeout=120)
    r.raise_for_status()
    with open(out_png, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    logger.info(f"✓ 保存通路图: {out_png.name}")


def main() -> None:
    parser = argparse.ArgumentParser(description="使用 Reactome AnalysisService 进行通路富集并下载通路图")
    parser.add_argument("--input_file", type=Path, default=PROJECT_ROOT / "results" / "brm_analysis_csv" / "significant_DEGs.csv",
                        help="输入基因列表文件：CSV需包含'gene'列，或纯文本每行一个基因")
    parser.add_argument("--output_dir", type=Path, default=PROJECT_ROOT / "results" / "reactome_enrichment",
                        help="输出目录")
    parser.add_argument("--analysis_name", type=str, default="SMARCA2_analysis",
                        help="用于区分本次运行的子目录名")
    parser.add_argument("--species", type=str, default="Homo sapiens", help="物种，例如 'Homo sapiens'")
    parser.add_argument("--fdr", type=float, default=0.05, help="显著性阈值（entities_fdr）")
    parser.add_argument("--top_n", type=int, default=20, help="下载前多少条显著通路的叠加通路图")
    args = parser.parse_args()

    out_dir = args.output_dir / args.analysis_name
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(out_dir)

    try:
        logger.info("=== Reactome 富集分析（API）===")
        genes = load_gene_list(args.input_file, logger)
        # 首选 REST token + pathways 接口；若失败则退回 reactome2py
        try:
            submit = submit_identifiers(genes, args.species, logger)
            token = submit["token"]
            _ = get_analysis_summary(token, logger)
            df = get_pathway_results(token, logger)
        except Exception as e_rest:
            logger.warning(f"通过 REST /token/.../pathways 获取失败，将退回 reactome2py：{e_rest}")
            df, token = get_pathways_via_reactome2py(genes, args.species, logger)
        if df.empty:
            logger.warning("未获取到任何通路结果。")
        else:
            # 保存完整结果
            df_all_path = out_dir / "reactome_pathways_all.csv"
            df.to_csv(df_all_path, index=False)
            logger.info(f"保存完整通路结果: {df_all_path.name} ({df.shape[0]} 行)")

            # 筛选显著（优先按 FDR，否则 pValue）
            if "entities_fdr" in df.columns and df["entities_fdr"].notna().any():
                sig_df = df[df["entities_fdr"].astype(float) < args.fdr].copy()
                sig_df = sig_df.sort_values(["entities_fdr", "entities_pValue"], na_position="last")
            else:
                sig_df = df[df["entities_pValue"].astype(float) < args.fdr].copy()
                sig_df = sig_df.sort_values(["entities_pValue"], na_position="last")

            sig_csv = out_dir / "reactome_pathways_significant.csv"
            sig_df.to_csv(sig_csv, index=False)
            logger.info(f"显著通路（阈值 {args.fdr}）: {len(sig_df)} 条 → {sig_csv.name}")

            # 下载前 top_n 通路图
            to_download = sig_df.head(args.top_n)
            if to_download.empty:
                logger.warning("无显著通路可供下载通路图。")
            else:
                logger.info(f"开始下载前 {len(to_download)} 条显著通路的叠加通路图…")
                for i, row in enumerate(to_download.itertuples(index=False), 1):
                    st_id: str = getattr(row, "stId")
                    name: str = getattr(row, "name") or st_id
                    rank = f"{i:02d}"
                    fname = sanitize_filename(f"{rank}_{name}_{st_id}.png")
                    # 尝试使用 token 叠加图；若没有 token（走了 reactome2py）或失败，则回退导出静态图
                    try:
                        token_val = locals().get('token')
                        if token_val:
                            download_pathway_diagram(token_val, st_id, out_dir / fname, logger)
                        else:
                            raise RuntimeError("无 token")
                    except Exception:
                        # 回退：使用 ContentService 的事件导出（不含分析叠加，但可用）
                        from reactome2py import content as rx_content
                        png_bytes = rx_content.export_event(st_id, format='png')
                        (out_dir / fname).write_bytes(png_bytes)
                        logger.info(f"✓ 保存通路图(静态): {fname}")

        logger.info("分析完成。输出目录：" + str(out_dir.relative_to(PROJECT_ROOT)))

    except Exception as e:
        logger.exception(f"运行失败: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()


