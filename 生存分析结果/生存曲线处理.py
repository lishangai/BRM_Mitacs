#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量为生存曲线PDF添加右上角标注（基因名 + Logrank p-value）
使用方法：
    直接在包含PDF文件的目录运行本脚本（或 python 生存曲线处理.py）
依赖：PyPDF2, reportlab
    pip install PyPDF2 reportlab
"""

from pathlib import Path
import re
from io import BytesIO

from PyPDF2 import PdfReader, PdfWriter
from reportlab.pdfgen import canvas
from reportlab.lib.colors import Color, black, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# ===== 可调样式 =====
FONT_SIZE = 36       # 标注字体大小（进一步增大）
MARGIN = 32          # 边距适当增大
CLEAR_WIDTH = 560    # 覆盖旧标注的白底矩形宽度（points）
CLEAR_HEIGHT = 110   # 覆盖旧标注的白底矩形高度（points）
TEXT_COLOR = black
BG_COLOR = white
# ====================

# 尝试注册中文字体（可选）
try:
    # 如果系统有 SimHei 或者其他中文字体，可替换为对应路径
    pdfmetrics.registerFont(TTFont('SimHei', 'SimHei.ttf'))
    DEFAULT_FONT = 'SimHei'
except Exception:
    DEFAULT_FONT = 'Helvetica'

# 解析文件名中的基因名与p值
# 示例：OV_BRM_Logrank p-value=0.5.pdf 或 OV_VLDLR _Logrank p-value=0.5.pdf
FILENAME_PATTERN = re.compile(r"OV[_\s]+(?P<gene>[^_\s]+)[_\s]+Logrank\s*p-value\s*=\s*(?P<p>[0-9]*\.?[0-9]+)", re.IGNORECASE)


def parse_gene_pvalue(filename: str):
    """从文件名中解析基因名和p值，解析失败返回(None, None)。"""
    name = Path(filename).stem
    m = FILENAME_PATTERN.search(name)
    if not m:
        return None, None
    gene = m.group('gene')
    pval = m.group('p')
    return gene, pval


def build_overlay(width: float, height: float, text: str) -> BytesIO:
    """构建与页面同尺寸的叠加PDF，右上角绘制白底矩形+文本。返回BytesIO。"""
    buf = BytesIO()
    c = canvas.Canvas(buf, pagesize=(width, height))

    # 计算文本宽度
    text_width = pdfmetrics.stringWidth(text, DEFAULT_FONT, FONT_SIZE)
    x_text = max(MARGIN, width - MARGIN - text_width)
    y_text = height - MARGIN - FONT_SIZE

    # 在右上角放置白底，尽量覆盖此前的标注区域
    rect_right = width - MARGIN
    rect_top = height - MARGIN
    rect_left = max(MARGIN, rect_right - CLEAR_WIDTH)
    rect_bottom = max(0, rect_top - CLEAR_HEIGHT)

    c.setFillColor(BG_COLOR)
    c.setStrokeColor(BG_COLOR)
    c.rect(rect_left, rect_bottom, rect_right - rect_left, rect_top - rect_bottom, fill=1, stroke=0)

    # 绘制文字
    c.setFillColor(TEXT_COLOR)
    c.setFont(DEFAULT_FONT, FONT_SIZE)
    c.drawString(x_text, y_text, text)

    c.showPage()
    c.save()
    buf.seek(0)
    return buf


def annotate_pdf(pdf_path: Path, gene: str, pval: str):
    """为PDF的每一页右上角添加标注（Gene + p值），覆盖原文件。"""
    try:
        reader = PdfReader(str(pdf_path))
    except Exception as e:
        print(f"✗ 无法读取PDF: {pdf_path.name} -> {e}")
        return

    writer = PdfWriter()
    label = f"{gene} | Logrank p-value={pval}"

    for i, page in enumerate(reader.pages):
        # 获取页面尺寸（单位：points）
        try:
            width = float(page.mediabox.width)
            height = float(page.mediabox.height)
        except Exception:
            # 常见A4默认尺寸（595 x 842 points）
            width, height = 595.0, 842.0

        overlay_buf = build_overlay(width, height, label)
        overlay_reader = PdfReader(overlay_buf)
        overlay_page = overlay_reader.pages[0]

        # 合并叠加层
        try:
            page.merge_page(overlay_page)  # PyPDF2 >= 2.0支持
        except Exception:
            try:
                page.mergePage(overlay_page)  # 旧版本
            except Exception as e:
                print(f"  - 第{i+1}页合并失败: {e}")
        writer.add_page(page)

    # 覆盖写回（先写到临时文件，再替换）
    tmp_path = pdf_path.with_suffix('.tmp.pdf')
    with open(tmp_path, 'wb') as f:
        writer.write(f)

    tmp_path.replace(pdf_path)
    print(f"✓ 已标注: {pdf_path.name} -> {label}")


def main():
    base_dir = Path(__file__).resolve().parent
    print(f"工作目录: {base_dir}")

    pdf_files = sorted([p for p in base_dir.glob('*.pdf') if p.is_file()])
    if not pdf_files:
        print("未在当前目录找到PDF文件。")
        return

    total = len(pdf_files)
    print(f"共发现 {total} 个PDF文件，开始处理……")

    for pdf in pdf_files:
        gene, pval = parse_gene_pvalue(pdf.name)
        if not gene or not pval:
            print(f"! 跳过（文件名未匹配基因与p值）：{pdf.name}")
            continue
        annotate_pdf(pdf, gene, pval)

    print("全部处理完成。")


if __name__ == '__main__':
    main()
