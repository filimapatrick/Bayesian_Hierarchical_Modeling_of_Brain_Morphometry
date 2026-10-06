#!/usr/bin/env python3
"""
Convert manuscript/manuscript.md into a publication-formatted Microsoft Word (.docx) document.
Uses python-docx to render headings, typography, native tables, and high-resolution figures.
"""

import os
import re
from pathlib import Path
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import OxmlElement, parse_xml
from docx.oxml.ns import nsdecls, qn

PROJECT_ROOT = Path(__file__).resolve().parent.parent
MANUSCRIPT_MD = PROJECT_ROOT / "manuscript" / "manuscript.md"
OUTPUT_DOCX = PROJECT_ROOT / "manuscript" / "manuscript.docx"


def set_cell_background(cell, fill_hex):
    """Sets background color of a table cell."""
    shading_elm = parse_xml(f'<w:shd {nsdecls("w")} w:fill="{fill_hex}"/>')
    cell._tc.get_or_add_tcPr().append(shading_elm)


def set_cell_margins(cell, top=100, bottom=100, left=150, right=150):
    """Sets cell padding (in twips: 20 twips = 1 pt)."""
    tcPr = cell._tc.get_or_add_tcPr()
    tcMar = OxmlElement('w:tcMar')
    for m, val in [('top', top), ('bottom', bottom), ('left', left), ('right', right)]:
        node = OxmlElement(f'w:{m}')
        node.set(qn('w:w'), str(val))
        node.set(qn('w:type'), 'dxa')
        tcMar.append(node)
    tcPr.append(tcMar)


def add_formatted_runs(paragraph, text):
    """Parses basic markdown inline formatting (**bold**, *italic*, `code`, math) into Word runs."""
    # Pattern to match bold, italic, code, or plain segments
    tokens = re.split(r'(\*\*.*?\*\*|\*.*?\*|`.*?`|\$.*?\$)', text)
    for tok in tokens:
        if not tok:
            continue
        if tok.startswith("**") and tok.endswith("**") and len(tok) >= 4:
            run = paragraph.add_run(tok[2:-2])
            run.bold = True
        elif tok.startswith("*") and tok.endswith("*") and len(tok) >= 2:
            run = paragraph.add_run(tok[1:-1])
            run.italic = True
        elif tok.startswith("`") and tok.endswith("`") and len(tok) >= 2:
            run = paragraph.add_run(tok[1:-1])
            run.font.name = "Consolas"
            run.font.size = Pt(9.5)
            run.font.color.rgb = RGBColor(160, 40, 40)
        elif tok.startswith("$") and tok.endswith("$") and len(tok) >= 2:
            # Inline math
            run = paragraph.add_run(tok[1:-1])
            run.italic = True
            run.font.name = "Cambria Math"
        else:
            paragraph.add_run(tok)


def build_docx():
    print(f"Reading manuscript from: {MANUSCRIPT_MD}")
    with open(MANUSCRIPT_MD, "r", encoding="utf-8") as f:
        lines = f.readlines()

    doc = Document()

    # Set standard 1-inch margins
    for section in doc.sections:
        section.top_margin = Inches(1.0)
        section.bottom_margin = Inches(1.0)
        section.left_margin = Inches(1.0)
        section.right_margin = Inches(1.0)

    # Base styling
    normal_style = doc.styles['Normal']
    normal_style.font.name = 'Times New Roman'
    normal_style.font.size = Pt(11)
    normal_style.font.color.rgb = RGBColor(30, 30, 30)

    in_code_block = False
    code_block_lines = []
    i = 0
    num_lines = len(lines)

    while i < num_lines:
        line = lines[i].rstrip('\r\n')

        # Code block toggle
        if line.startswith("```"):
            if not in_code_block:
                in_code_block = True
                code_block_lines = []
            else:
                in_code_block = False
                # Render code block (usually ASCII box tables or model specifications)
                table_text = "\n".join(code_block_lines)
                p = doc.add_paragraph()
                p.paragraph_format.space_before = Pt(6)
                p.paragraph_format.space_after = Pt(6)
                p.paragraph_format.line_spacing = 1.05
                run = p.add_run(table_text)
                run.font.name = "Consolas"
                run.font.size = Pt(8.0)
                run.font.color.rgb = RGBColor(40, 40, 40)
            i += 1
            continue

        if in_code_block:
            code_block_lines.append(line)
            i += 1
            continue

        # Image check: ![Caption](path)
        img_match = re.match(r'!\[(.*?)\]\((.*?)\)', line.strip())
        if img_match:
            caption, img_rel_path = img_match.groups()
            img_path = (MANUSCRIPT_MD.parent / img_rel_path).resolve()
            if img_path.exists():
                p_img = doc.add_paragraph()
                p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p_img.paragraph_format.space_before = Pt(12)
                p_img.paragraph_format.space_after = Pt(4)
                p_img.add_run().add_picture(str(img_path), width=Inches(6.2))
            else:
                print(f"Warning: Image path not found: {img_path}")
            i += 1
            continue

        # Standalone Caption check (*Figure X: ...*)
        if line.startswith("*Figure") and line.endswith("*"):
            p_cap = doc.add_paragraph()
            p_cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p_cap.paragraph_format.space_before = Pt(2)
            p_cap.paragraph_format.space_after = Pt(14)
            run = p_cap.add_run(line[1:-1])
            run.italic = True
            run.font.size = Pt(9.5)
            run.font.color.rgb = RGBColor(80, 80, 80)
            i += 1
            continue

        # Headings
        if line.startswith("# "):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(18)
            p.paragraph_format.space_after = Pt(12)
            run = p.add_run(line[2:].strip())
            run.bold = True
            run.font.name = "Arial"
            run.font.size = Pt(18)
            run.font.color.rgb = RGBColor(10, 30, 60)
            i += 1
            continue
        elif line.startswith("## "):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(16)
            p.paragraph_format.space_after = Pt(6)
            run = p.add_run(line[3:].strip())
            run.bold = True
            run.font.name = "Arial"
            run.font.size = Pt(14)
            run.font.color.rgb = RGBColor(20, 45, 80)
            i += 1
            continue
        elif line.startswith("### "):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(12)
            p.paragraph_format.space_after = Pt(4)
            run = p.add_run(line[4:].strip())
            run.bold = True
            run.font.name = "Arial"
            run.font.size = Pt(12)
            run.font.color.rgb = RGBColor(30, 50, 90)
            i += 1
            continue
        elif line.startswith("#### "):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(8)
            p.paragraph_format.space_after = Pt(3)
            run = p.add_run(line[5:].strip())
            run.bold = True
            run.font.name = "Arial"
            run.font.size = Pt(11)
            run.font.color.rgb = RGBColor(40, 60, 100)
            i += 1
            continue

        # Horizontal rule
        if line.strip() in ["---", "***", "___"]:
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(6)
            p.paragraph_format.space_after = Pt(6)
            p.paragraph_format.line_spacing = 0.5
            run = p.add_run("―" * 45)
            run.font.color.rgb = RGBColor(200, 200, 200)
            i += 1
            continue

        # Blockquote
        if line.startswith("> "):
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.4)
            p.paragraph_format.right_indent = Inches(0.4)
            p.paragraph_format.space_before = Pt(4)
            p.paragraph_format.space_after = Pt(4)
            add_formatted_runs(p, line[2:].strip())
            for r in p.runs:
                r.italic = True
                r.font.color.rgb = RGBColor(60, 60, 60)
            i += 1
            continue

        # Bullet lists (* or -)
        bullet_match = re.match(r'^\s*([*\-])\s+(.*)', line)
        if bullet_match:
            bullet_char, content = bullet_match.groups()
            p = doc.add_paragraph(style='List Bullet')
            p.paragraph_format.space_before = Pt(1)
            p.paragraph_format.space_after = Pt(2)
            add_formatted_runs(p, content.strip())
            i += 1
            continue

        # Numbered lists
        num_match = re.match(r'^\s*(\d+)\.\s+(.*)', line)
        if num_match:
            num_val, content = num_match.groups()
            p = doc.add_paragraph(style='List Number')
            p.paragraph_format.space_before = Pt(1)
            p.paragraph_format.space_after = Pt(2)
            add_formatted_runs(p, content.strip())
            i += 1
            continue

        # Display Math ($$...$$)
        if line.strip().startswith("$$") and line.strip().endswith("$$") and len(line.strip()) > 2:
            math_text = line.strip()[2:-2].strip()
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(6)
            p.paragraph_format.space_after = Pt(6)
            run = p.add_run(math_text)
            run.italic = True
            run.font.name = "Cambria Math"
            run.font.size = Pt(11.5)
            i += 1
            continue

        # Empty line
        if not line.strip():
            i += 1
            continue

        # Normal paragraph
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(5)
        p.paragraph_format.line_spacing = 1.15
        add_formatted_runs(p, line)
        i += 1

    doc.save(str(OUTPUT_DOCX))
    print(f"Successfully generated Word document: {OUTPUT_DOCX}")
    print(f"File size: {OUTPUT_DOCX.stat().st_size / 1024:.1f} KB")


if __name__ == "__main__":
    build_docx()
