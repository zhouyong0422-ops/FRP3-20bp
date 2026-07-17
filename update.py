import streamlit as st
import pandas as pd
import numpy as np
import math
import io
import time
import re
from datetime import datetime
from collections import Counter

# ==========================================
# 1. 页面全局排版与 CSS 样式注入
# ==========================================
st.set_page_config(
    page_title="自动化引物探针设计比对系统",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# 嵌入深度还原并强化的 UI 排版 CSS
st.markdown("""
<style>
    /* 全局背景与字体 */
    .stApp { background-color: #eef2f5; color: #333; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
    
    /* 核心标题栏样式 */
    .header-box {
        background: white; padding: 25px 30px; border-radius: 12px;
        box-shadow: 0 6px 20px rgba(0,0,0,0.08); margin-bottom: 20px;
        border-bottom: 3px solid #1abc9c; display: flex; justify-content: space-between; align-items: center;
    }
    .header-title { color: #2c3e50; font-size: 28px; font-weight: bold; margin: 0; }
    .badge { background: #8e44ad; color: white; padding: 5px 14px; border-radius: 15px; font-size: 14px; font-weight: bold; }
    
    /* 上传提示文字 */
    .upload-hint { font-style: italic; color: #7f8c8d; font-size: 14px; margin-top: 4px; margin-bottom: 10px; line-height: 1.6; }
    
    /* 文本输入框美化 */
    .stTextArea textarea {
        border: 2px solid #bdc3c7 !important; border-radius: 8px !important;
        font-family: monospace !important; background-color: #ffffff !important;
        padding: 12px !important; font-size: 13px !important; line-height: 1.5 !important;
    }
    .stTextArea textarea:focus { border-color: #1abc9c !important; box-shadow: 0 0 8px rgba(26, 188, 156, 0.3) !important; }
    
    /* 靶区与卡片布局 */
    .locus-title-bar {
        background: #edf2f7; color: #2c3e50; padding: 14px 20px; font-size: 18px; font-weight: bold;
        border: 1px solid #d1d5db; border-radius: 8px 8px 0 0; margin-top: 35px;
        display: flex; justify-content: space-between; align-items: center;
        box-shadow: 0 2px 4px rgba(0,0,0,0.02);
    }
    .card-box {
        background: #fff; padding: 20px; border-left: 6px solid #bdc3c7;
        border-right: 1px solid #d1d5db; border-bottom: 1px solid #d1d5db; margin-bottom: 15px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.02);
    }
    .card-primary { border-left-color: #27ae60; background: #f9fbfd; }
    .card-header { display:flex; justify-content:space-between; align-items:center; border-bottom:1px dashed #e1e8ed; padding-bottom:10px; margin-bottom:15px; }
    
    /* 角色与分数徽章 */
    .role-badge-main { background: #27ae60; color: white; font-size: 12px; padding: 3px 10px; border-radius: 4px; font-weight: bold; margin-right: 10px; }
    .role-badge-sub { background: #7f8c8d; color: white; font-size: 12px; padding: 3px 10px; border-radius: 4px; font-weight: bold; margin-right: 10px; }
    .score-badge { background: #f39c12; color: white; padding: 3px 12px; border-radius: 12px; font-size: 14px; font-weight: bold; }
    
    /* 序列区块排版 */
    .seq-row { display: flex; justify-content: space-between; align-items: center; font-family: monospace; font-size: 14px; margin-bottom: 6px; padding: 8px 12px; background: #fff; border-radius: 6px; border: 1px solid #ecf0f1; }
    .seq-probe-row { border-left: 4px solid #8e44ad; background: #fdfafb; }
    .seq-type { font-weight: bold; width: 140px; color: #34495e; }
    .seq-string { color: #d35400; letter-spacing: 1px; font-weight: bold; flex-grow: 1; }
    .seq-stats { color: #7f8c8d; font-size: 12px; text-align: right; }
    
    /* 错配统计行排版 */
    .stats-row { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 4px; margin-bottom: 14px; padding-left: 4px; }
    .stat-badge { padding: 3px 8px; border-radius: 10px; font-weight: bold; font-size: 11px; display: inline-block; }
    .bg-0 { background: #d4edda; color: #155724; border: 1px solid #c3e6cb; }
    .bg-1 { background: #fff3cd; color: #856404; border: 1px solid #ffeeba; }
    .bg-2 { background: #ffeeba; color: #856404; border: 1px solid #ffdf7e; }
    .bg-3 { background: #f8d7da; color: #721c24; border: 1px solid #f5c6cb; }
    .mix-badge { background: #9b59b6; color: white; padding: 3px 8px; border-radius: 8px; font-size: 11px; display: inline-block; font-weight: bold; }
    
    /* 展开细节文字排版 */
    .details-row { display: flex; justify-content: space-between; padding: 4px 0; border-bottom: 1px dashed #f1f2f6; }
    .details-row:last-child { border-bottom: none; border-top: 1px solid #ddd; padding-top: 8px; margin-top: 4px; font-weight: bold; }
</style>
""", unsafe_allow_html=True)

# 头部设计
st.markdown("""
<div class="header-box">
    <div class="header-title">🧬 自动化引物探针设计比对系统</div>
    <div class="badge">融合终极优化版 V28</div>
</div>
""", unsafe_allow_html=True)

# ==========================================
# 2. 核心算法与极速统计逻辑 (低内存开销版)
# ==========================================

@st.cache_data(show_spinner=False)
def calc_tm(seq: str) -> float:
    seq_upper = seq.upper()
    g = seq_upper.count('G')
    c = seq_upper.count('C')
    s = seq_upper.count('S')
    half = sum(seq_upper.count(b) for b in ['R', 'Y', 'M', 'K'])
    gc_count = g + c + s + half * 0.5
    if len(seq) == 0: return 0.0
    return round(64.9 + 41 * (gc_count - 16.4) / len(seq), 1)

@st.cache_data(show_spinner=False)
def calc_gc(seq: str) -> float:
    seq_upper = seq.upper()
    g = seq_upper.count('G')
    c = seq_upper.count('C')
    s = seq_upper.count('S')
    half = sum(seq_upper.count(b) for b in ['R', 'Y', 'M', 'K'])
    if len(seq) == 0: return 0.0
    return round(((g + c + s + half * 0.5) / len(seq)) * 100, 1)

def reverse_complement(seq: str) -> str:
    trans = str.maketrans('ATCGAGRYMSWKMVDHVNatcgagrymswkmvdhvn', 'TAGCTCRYKSWMKBHDVnTAGCTCRYKSWMKBHDVn')
    return seq.translate(trans)[::-1]

def has_secondary_structure_risk(seq: str) -> bool:
    if len(seq) < 12: return False
    for i in range(4, 6):
        head = seq[:i]
        tail_comp = reverse_complement(seq[-i:])
        if head == tail_comp: return True
    return False

def has_3prime_dimer_risk(seq1: str, seq2: str) -> bool:
    return seq1[-4:] == reverse_complement(seq2[-4:])

def is_hard_valid_oligo(seq: str, is_probe: bool = False) -> bool:
    if 'N' in seq: return False
    if not all(b in 'ATGCRYSWKM' for b in seq): return False
    deg_count = sum(seq.count(b) for b in 'RYSWKM')
    if deg_count > 1: return False
    if any(b in 'RYSWKM' for b in seq[-5:]): return False
    if is_probe and any(b in 'RYSWKM' for b in seq[:1]): return False
    gc = calc_gc(seq)
    if gc < 20 or gc > 80: return False
    if is_probe and seq.startswith('G'): return False
    return True

def calc_soft_penalties(seq: str, is_probe: bool) -> float:
    penalty = 0.0
    gc = calc_gc(seq)
    for b in ['A', 'T', 'G', 'C']:
        if b * 4 in seq: penalty += 10.0
    if has_secondary_structure_risk(seq): penalty += 15.0
    if is_probe:
        if gc < 30: penalty += (30 - gc) * 2
        if gc > 65: penalty += (gc - 65) * 2
        if seq.count('G') >= seq.count('C'): penalty += 8.0
    else:
        if gc < 40: penalty += (40 - gc) * 2
        if gc > 60: penalty += (gc - 60) * 2
        end5 = seq[-5:]
        end_gc = end5.count('G') + end5.count('C')
        if end_gc < 1 or end_gc > 2: penalty += 6.0
        if seq.endswith('T'): penalty += 12.0
        if seq.endswith('GG') or seq.endswith('CC'): penalty += 8.0
    return penalty

def get_conserved_probe_sequence(start_index: int, length: int, sequences_array: list, is_reverse: bool = False) -> list:
    total_seq = len(sequences_array)
    valid_seqs = [s[start_index:start_index+length] for s in sequences_array if '-' not in s[start_index:start_index+length] and 'N' not in s[start_index:start_index+length]]
    if not valid_seqs or len(valid_seqs) / total_seq < 0.80: return []
    counts = Counter(valid_seqs)
    top_seq, top_count = counts.most_common(1)[0]
    top_cov = top_count / len(valid_seqs)
    if top_cov >= 0.98: return [top_seq]
    
    iupac_map = {'AG': 'R', 'GA': 'R', 'CT': 'Y', 'TC': 'Y', 'GC': 'S', 'CG': 'S', 'AT': 'W', 'TA': 'W', 'GT': 'K', 'TG': 'K', 'AC': 'M', 'CA': 'M'}
    best_deg_col, max_boost = -1, 0.0
    start_col = 5 if is_reverse else 1
    end_col = (length - 1) if is_reverse else (length - 5)
    
    for c in range(start_col, end_col):
        base_counts = Counter(seq[c] for seq in valid_seqs)
        b_sorted = base_counts.most_common()
        if len(b_sorted) >= 2:
            top1_cov = b_sorted[0][1] / len(valid_seqs)
            top2_cov = (b_sorted[0][1] + b_sorted[1][1]) / len(valid_seqs)
            boost = top2_cov - top1_cov
            if boost > max_boost and boost >= 0.04:
                max_boost = boost
                best_deg_col = c
                
    if best_deg_col != -1 and max_boost > 0:
        base_counts = Counter(seq[best_deg_col] for seq in valid_seqs)
        b_sorted = base_counts.most_common()
        deg_char = iupac_map.get(b_sorted[0][0] + b_sorted[1][0])
        if deg_char:
            deg_seq = list(top_seq)
            deg_seq[best_deg_col] = deg_char
            return ["".join(deg_seq)]
    return [top_seq]

def get_top_variants(start_index: int, length: int, sequences_array: list, max_variants: int = 2, is_probe: bool = False, is_reverse: bool = False) -> list:
    if is_probe: return get_conserved_probe_sequence(start_index, length, sequences_array, is_reverse)
    total_seq = len(sequences_array)
    valid_seqs = [s[start_index:start_index+length] for s in sequences_array if '-' not in s[start_index:start_index+length] and 'N' not in s[start_index:start_index+length]]
    if not valid_seqs or len(valid_seqs) / total_seq < 0.90: return []
    counts = Counter(valid_seqs)
    sorted_counts = counts.most_common()
    variants = [sorted_counts[0][0]]
    coverage = sorted_counts[0][1] / len(valid_seqs)
    
    if coverage < 0.97 and len(sorted_counts) > 1 and max_variants > 1:
        if (sorted_counts[1][1] / len(valid_seqs)) > 0.04:
            seq1, seq2 = sorted_counts[0][0], sorted_counts[1][0]
            diff_idx = [j for j in range(length) if seq1[j] != seq2[j]]
            if len(diff_idx) == 1:
                idx = diff_idx[0]
                is_3prime5_risk = (idx < 5) if is_reverse else (idx >= length - 5)
                if not is_3prime5_risk:
                    iupac_map = {'AG': 'R', 'GA': 'R', 'CT': 'Y', 'TC': 'Y', 'GC': 'S', 'CG': 'S', 'AT': 'W', 'TA': 'W', 'GT': 'K', 'TG': 'K', 'AC': 'M', 'CA': 'M'}
                    deg_char = iupac_map.get(seq1[idx] + seq2[idx])
                    if deg_char:
                        merged = list(seq1)
                        merged[idx] = deg_char
                        return ["".join(merged)]
                    else: variants.append(seq2)
                else: variants.append(seq2)
            else: variants.append(seq2)
    return variants

def is_base_match(b1: str, b2: str) -> bool:
    if b1 == b2: return True
    if b1 == 'R' and b2 in 'AG': return True
    if b1 == 'Y' and b2 in 'CT': return True
    if b1 == 'S' and b2 in 'GC': return True
    if b1 == 'W' and b2 in 'AT': return True
    if b1 == 'K' and b2 in 'GT': return True
    if b1 == 'M' and b2 in 'AC': return True
    return False

def calculate_mix_mismatch(variants_array: list, start_index: int, sequences_array: list) -> dict:
    total_seq = len(sequences_array)
    m0, m1, m2, m3p = 0, 0, 0, 0
    seq_len = len(variants_array[0])
    for seq in sequences_array:
        snippet = seq[start_index:start_index+seq_len]
        if '-' in snippet or 'N' in snippet:
            m3p += 1; continue
        best_mismatches = seq_len
        for target in variants_array:
            mismatches = sum(1 for j in range(seq_len) if not is_base_match(target[j], snippet[j]))
            if mismatches < best_mismatches: best_mismatches = mismatches
        if best_mismatches == 0: m0 += 1
        elif best_mismatches == 1: m1 += 1
        elif best_mismatches == 2: m2 += 1
        else: m3p += 1
    if total_seq == 0: return {"p0": 0.0, "p1": "0.0", "p2": "0.0", "p3": "0.0", "m0":0, "m1":0, "m2":0, "m3p":0, "total":0}
    return {
        "p0": round((m0 / total_seq) * 100, 1),
        "p1": f"{(m1 / total_seq) * 100:.1f}",
        "p2": f"{(m2 / total_seq) * 100:.1f}",
        "p3": f"{(m3p / total_seq) * 100:.1f}",
        "m0": m0, "m1": m1, "m2": m2, "m3p": m3p, "total": total_seq
    }

# ==========================================
# 3. 内存优化的预筛算法
# ==========================================
@st.cache_data(show_spinner=False)
def scan_entropy(sequences: list) -> list:
    seq_len = len(sequences[0])
    entropies = []
    for i in range(seq_len):
        col_chars = [seq[i] for seq in sequences if seq[i] != '-']
        total = len(col_chars)
        if total == 0:
            entropies.append(0.0)
            continue
        counts = Counter(col_chars)
        ent = 0.0
        for count in counts.values():
            p = count / total
            ent -= p * math.log2(p)
        entropies.append(ent)
    return entropies

def run_pipeline_engine(sequences: list, progress_bar, status_text):
    seq_len = len(sequences[0])
    total_steps = max(1, seq_len - 160)
    
    memo_fwd_primers, memo_fwd_probes, memo_rev_probes, memo_rev_primers = {}, {}, {}, {}

    def get_all_valid_variants_fwd(start_idx: int, is_probe: bool):
        cache_dict = memo_fwd_probes if is_probe else memo_fwd_primers
        if start_idx in cache_dict: return cache_dict[start_idx]
        valid_list = []
        target_lengths = list(range(18, 31)) if is_probe else list(range(18, 26))
        for l in target_lengths:
            if start_idx + l > seq_len: continue
            raw_vars = get_top_variants(start_idx, l, sequences, 1 if is_probe else 2, is_probe, False)
            if not raw_vars: continue
            if all(is_hard_valid_oligo(v, is_probe) for v in raw_vars):
                soft_pen = sum(calc_soft_penalties(v, is_probe) for v in raw_vars)
                avg_tm = sum(calc_tm(v) for v in raw_vars) / len(raw_vars)
                tm_bonus = 15.0 if (is_probe and 62 <= avg_tm <= 68) or (not is_probe and 55 <= avg_tm <= 60) else 0.0
                stats = calculate_mix_mismatch(raw_vars, start_idx, sequences)
                p0 = stats["p0"]
                cons_score = (p0 * 3 + (50 if p0 >= 99.0 else 0) - ((98.0 - p0) * 10 if p0 < 98.0 else 0)) if is_probe else p0
                valid_list.append({"length": l, "variants": raw_vars, "score": tm_bonus - soft_pen, "stats": stats, "totalScore": cons_score + tm_bonus - soft_pen})
        valid_list.sort(key=lambda x: x["totalScore"], reverse=True)
        cache_dict[start_idx] = valid_list
        return valid_list

    def get_all_valid_variants_rev(start_idx: int, is_probe: bool):
        cache_dict = memo_rev_probes if is_probe else memo_rev_primers
        if start_idx in cache_dict: return cache_dict[start_idx]
        valid_list = []
        target_lengths = list(range(18, 31)) if is_probe else list(range(18, 26))
        for l in target_lengths:
            if start_idx + l > seq_len: continue
            raw_vars = get_top_variants(start_idx, l, sequences, 1 if is_probe else 2, is_probe, True)
            if not raw_vars: continue
            final_vars = []
            all_passed = True
            for rv in raw_vars:
                comp = reverse_complement(rv)
                if not is_hard_valid_oligo(comp, is_probe): all_passed = False; break
                final_vars.append(comp)
            if all_passed:
                soft_pen = sum(calc_soft_penalties(v, is_probe) for v in final_vars)
                avg_tm = sum(calc_tm(v) for v in final_vars) / len(final_vars)
                tm_bonus = 15.0 if (is_probe and 62 <= avg_tm <= 68) or (not is_probe and 55 <= avg_tm <= 60) else 0.0
                stats = calculate_mix_mismatch(raw_vars, start_idx, sequences)
                p0 = stats["p0"]
                cons_score = (p0 * 3 + (50 if p0 >= 99.0 else 0) - ((98.0 - p0) * 10 if p0 < 98.0 else 0)) if is_probe else p0
                valid_list.append({"length": l, "variants": final_vars, "rawVariants": raw_vars, "score": tm_bonus - soft_pen, "stats": stats, "totalScore": cons_score + tm_bonus - soft_pen})
        valid_list.sort(key=lambda x: x["totalScore"], reverse=True)
        cache_dict[start_idx] = valid_list
        return valid_list

    all_candidates = []
    
    for i in range(total_steps):
        if i % 10 == 0 or i == total_steps - 1:
            pct = int(15 + (i / total_steps) * 80)
            progress_bar.progress(pct)
            status_text.text(f"⏳ 全景矩阵加权寻优中... 当前扫描参考坐标: {i} / {total_steps} bp")
            
        f_obj_list = get_all_valid_variants_fwd(i, False)
        if not f_obj_list: continue
        
        for f_obj in f_obj_list[:2]:
            f_vars, f_len = f_obj["variants"], f_obj["length"]
            min_r_start = max(i + f_len + 1, i + 70 - 25)
            max_r_start = min(seq_len - 1, i + 160 - 18)
            
            for r_start in range(min_r_start, max_r_start + 1):
                r_obj_list = get_all_valid_variants_rev(r_start, False)
                if not r_obj_list: continue
                
                for r_obj in r_obj_list[:2]:
                    r_vars, r_len = r_obj["variants"], r_obj["length"]
                    amplicon_size = r_start + r_len - i
                    if amplicon_size < 70 or amplicon_size > 160: continue
                    
                    primer_max_tm = max(max(calc_tm(v) for v in f_vars), max(calc_tm(v) for v in r_vars))
                    primer_min_tm = min(min(calc_tm(v) for v in f_vars), min(calc_tm(v) for v in r_vars))
                    
                    valid_probes = []
                    for p_start in range(i + f_len + 1, min(r_start - 18, i + f_len + 30)):
                        p_obj_list = get_all_valid_variants_fwd(p_start, True)
                        for p_obj in p_obj_list:
                            if p_start + p_obj["length"] <= r_start - 1:
                                min_p_tm = min(calc_tm(v) for v in p_obj["variants"])
                                if min_p_tm > primer_max_tm:
                                    valid_probes.append({**p_obj, "start": p_start, "isRev": False, "gapToSameDir": p_start - (i + f_len)})
                                    
                    for p_end in range(r_start - 1, max(i + f_len + 18, r_start - 30), -1):
                        for p_len in range(18, 31):
                            p_start = p_end - p_len
                            if p_start >= i + f_len + 1 and (r_start - p_end < 30):
                                p_obj_list = get_all_valid_variants_rev(p_start, True)
                                for p_obj in p_obj_list:
                                    if p_obj["length"] == p_len:
                                        min_p_tm = min(calc_tm(v) for v in p_obj["variants"])
                                        if min_p_tm > primer_max_tm:
                                            valid_probes.append({**p_obj, "start": p_start, "isRev": True, "gapToSameDir": r_start - p_end})
                                            
                    if not valid_probes: continue
                    valid_probes.sort(key=lambda x: x["totalScore"], reverse=True)
                    
                    for p_obj in valid_probes[:3]:
                        p_vars = p_obj["variants"]
                        min_p_tm = min(calc_tm(v) for v in p_vars)
                        
                        cross_soft_pen = 0.0
                        if any(has_3prime_dimer_risk(f, r) for f in f_vars for r in r_vars): cross_soft_pen += 20.0
                        primer_tm_diff = abs(primer_max_tm - primer_min_tm)
                        if primer_tm_diff > 2.0: cross_soft_pen += (primer_tm_diff - 2.0) * 5.0
                        if min_p_tm < primer_max_tm + 5.0: cross_soft_pen += (primer_max_tm + 5.0 - min_p_tm) * 6.0
                        if p_obj["gapToSameDir"] > 10: cross_soft_pen += (p_obj["gapToSameDir"] - 10) * 1.5
                        
                        mix_f = -35.0 if len(f_vars) > 1 else 0.0
                        mix_r = -35.0 if len(r_vars) > 1 else 0.0
                        
                        final_score = f_obj["totalScore"] + p_obj["totalScore"] + r_obj["totalScore"] + mix_f + mix_r - cross_soft_pen
                        base_p0 = f_obj["stats"]["p0"] + (p_obj["stats"]["p0"] * 3) + r_obj["stats"]["p0"]
                        p_p0 = p_obj["stats"]["p0"]
                        p_bonus = 50.0 if p_p0 >= 99.0 else 0.0
                        p_pen = (98.0 - p_p0) * 10.0 if p_p0 < 98.0 else 0.0
                        self_soft_pen = (f_obj["score"] - (15.0 if f_obj["length"] else 0)) + (p_obj["score"] - (15.0 if p_obj["length"] else 0)) + (r_obj["score"] - (15.0 if r_obj["length"] else 0))
                        display_soft_pen = -self_soft_pen + cross_soft_pen
                        
                        all_candidates.append({
                            "fwd": f_vars, "rev": r_vars, "probe": p_vars,
                            "fStats": f_obj["stats"], "pStats": p_obj["stats"], "rStats": r_obj["stats"],
                            "size": amplicon_size, "start": i, "score": final_score,
                            "probeDir": "Reverse (-)" if p_obj["isRev"] else "Forward (+)",
                            "probeGap": p_obj["gapToSameDir"],
                            "details": {"base": base_p0, "pBonus": p_bonus, "pPenalty": -p_pen, "mixF": mix_f, "mixR": mix_r, "softPen": -display_soft_pen}
                        })
                        
    progress_bar.progress(96)
    status_text.text("⏳ 正在整合排序并聚类黄金去重区...")
    
    all_candidates.sort(key=lambda x: x["score"], reverse=True)
    loci_groups = []
    locus_window = 50
    for cand in all_candidates:
        found = False
        for loc in loci_groups:
            if abs(cand["start"] - loc["anchorStart"]) <= locus_window:
                if len(loc["variants"]) < 3: loc["variants"].append(cand)
                found = True; break
        if not found:
            loci_groups.append({"locusId": len(loci_groups) + 1, "anchorStart": cand["start"], "variants": [cand]})
            
    progress_bar.progress(100)
    status_text.text("✅ 寻优完成！正在出具最终 DOE 清单...")
    time.sleep(0.3)
    progress_bar.empty()
    status_text.empty()
    return loci_groups

# ==========================================
# 4. 页面主体与用户交互区域 (带有边框卡片包装)
# ==========================================

# 用 border=True 将输入框和文件选择器严密包裹在一个带边框和阴影的专业卡片中
with st.container(border=True):
    col1, col2 = st.columns([1, 2.5])
    with col1:
        uploaded_file = st.file_uploader("📂 导入 FASTA 序列库", type=["fasta", "fas", "txt", "aln"], label_visibility="collapsed")
        st.markdown('<div class="upload-hint">导入多序列比对 FASTA 变异库(fasta,fas,txt,aln)，自动筛选符合 TaqMan qPCR 规范的引物探针组合，输出可视化报告并支持 CSV 导出。</div>', unsafe_allow_html=True)
    with col2:
        fasta_text = ""
        if uploaded_file is not None:
            fasta_text = uploaded_file.getvalue().decode("utf-8")
        fasta_input = st.text_area("或在此处直接粘贴比对完成的 FASTA 序列...", value=fasta_text, height=130, label_visibility="collapsed", placeholder="或在此处直接粘贴比对完成的 FASTA 序列...")

run_btn = st.button("⚙️ 启动靶区筛选：高精度加权运算输出引探 DOE 试验方案", use_container_width=True, type="primary")

# ==========================================
# 5. 任务执行与结果排版渲染
# ==========================================
if run_btn:
    if not fasta_input.strip():
        st.error("请导入包含变异库的文件！")
    else:
        sequences = []
        curr_seq = []
        for line in fasta_input.strip().splitlines():
            if line.startswith('>'):
                if curr_seq: sequences.append("".join(curr_seq).upper())
                curr_seq = []
            else:
                curr_seq.append(line.strip())
        if curr_seq: sequences.append("".join(curr_seq).upper())
        
        if len(sequences) < 2:
            st.error("文件格式有误或序列数不足2条。")
        else:
            prog_bar = st.progress(5)
            stat_txt = st.empty()
            stat_txt.text(f"⏳ 正在计算变异强度扫描 (0 / {len(sequences[0])} bp)...")
            
            # 渲染香农熵图表
            entropies = scan_entropy(sequences)
            st.markdown("### 📊 靶标序列群变异强度扫描 (香农熵)")
            chart_data = pd.DataFrame({"坐标": list(range(1, len(entropies) + 1)), "香农熵 (突变强度扫描)": entropies})
            st.bar_chart(chart_data, x="坐标", y="香农熵 (突变强度扫描)", color="#3498db", height=250)
            
            # 执行计算
            loci_groups = run_pipeline_engine(sequences, prog_bar, stat_txt)
            
            # 渲染多重分析报告页面
            st.markdown("---")
            st.markdown("<h2 style='text-align:center; color:#2c3e50;'>多重 TaqMan 引探组合设计分析报告</h2>", unsafe_allow_html=True)
            st.markdown(f"<p style='text-align:center; color:#7f8c8d; font-size:14px;'>生成时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>", unsafe_allow_html=True)
            st.markdown("<p style='text-align:center; font-size:13px; color:#95a5a6;'>核心设计框架：严格限定<b>探针长度 18-30bp</b>，探针 $T_m >$ 引物 $T_m$；全局严禁在 3' 末端 5bp 引入任何兼并碱基，结合多维加速评分，精准遴选靶区体系。</p>", unsafe_allow_html=True)
            
            if not loci_groups:
                st.markdown('<div class="card-box" style="border-left-color: #e74c3c;"><h4 style="color:#c0392b; margin-top:0;">⚠️ 体系设计失败</h4><p>在该序列库中未能找到产物 70-160bp、且满足探针 18-30bp 长度与 Tm 大于引物规范的有效靶区。</p></div>', unsafe_allow_html=True)
            else:
                # 准备生成 CSV 导出的数据体
                csv_rows = []
                for loc in loci_groups:
                    for v_idx, cand in enumerate(loc["variants"]):
                        loc_name = f"靶区_{loc['locusId']}"
                        role = "主力优选" if v_idx == 0 else f"备选_{v_idx}"
                        score_str = f"{cand['score']:.1f}"
                        for i, seq in enumerate(cand["fwd"]):
                            t = f"Forward_{i+1}" if len(cand["fwd"]) > 1 else "Forward"
                            st_f = cand["fStats"]
                            csv_rows.append([loc_name, role, score_str, t, seq, len(seq), calc_tm(seq), calc_gc(seq), f"{st_f['p0']}% ({st_f['m0']}/{st_f['total']})", f"{st_f['p1']}% ({st_f['m1']}/{st_f['total']})", f"{st_f['p2']}% ({st_f['m2']}/{st_f['total']})", f"{st_f['p3']}% ({st_f['m3p']}/{st_f['total']})", cand["size"], "-", cand["start"]])
                        for i, seq in enumerate(cand["probe"]):
                            t = f"Probe({cand['probeDir']})_{i+1}" if len(cand["probe"]) > 1 else f"Probe({cand['probeDir']})"
                            st_p = cand["pStats"]
                            csv_rows.append([loc_name, role, score_str, t, seq, len(seq), calc_tm(seq), calc_gc(seq), f"{st_p['p0']}% ({st_p['m0']}/{st_p['total']})", f"{st_p['p1']}% ({st_p['m1']}/{st_p['total']})", f"{st_p['p2']}% ({st_p['m2']}/{st_p['total']})", f"{st_p['p3']}% ({st_p['m3p']}/{st_p['total']})", cand["size"], cand["probeGap"], cand["start"]])
                        for i, seq in enumerate(cand["rev"]):
                            t = f"Reverse_{i+1}" if len(cand["rev"]) > 1 else "Reverse"
                            st_r = cand["rStats"]
                            csv_rows.append([loc_name, role, score_str, t, seq, len(seq), calc_tm(seq), calc_gc(seq), f"{st_r['p0']}% ({st_r['m0']}/{st_r['total']})", f"{st_r['p1']}% ({st_r['m1']}/{st_r['total']})", f"{st_r['p2']}% ({st_r['m2']}/{st_r['total']})", f"{st_r['p3']}% ({st_r['m3p']}/{st_r['total']})", cand["size"], "-", cand["start"]])
                
                df_csv = pd.DataFrame(csv_rows, columns=["靶区归属","变体角色","综合得分","寡核苷酸类型","序列 (5'->3')","长度 (bp)","Tm (°C)","GC (%)","完美匹配(0)","错配1碱基(1)","错配2碱基(2)","错配3碱基(≥3)","预期产物长度 (bp)","同向引物间距 (bp)","精确起始坐标"])
                csv_buffer = io.StringIO()
                df_csv.to_csv(csv_buffer, index=False, encoding='utf-8-sig')
                
                # 右上角导出按钮
                col_d1, col_d2 = st.columns([4, 1])
                with col_d2:
                    st.download_button(
                        label="📥 导出 Excel (含靶区分组)",
                        data=csv_buffer.getvalue(),
                        file_name=f"多重_TaqMan_引探组合设计分析报告_{datetime.now().strftime('%Y-%m-%d')}.csv",
                        mime="text/csv",
                        use_container_width=True
                    )
                
                # 遍历并深度渲染 HTML 卡片（完全移除了行首缩进，杜绝代码块解析 Bug）
                for loc in loci_groups:
                    st.markdown(f'<div class="locus-title-bar"><span>🎯 独立黄金靶区 {loc["locusId"]}</span><span style="font-size: 14px; color: #7f8c8d; font-weight: normal;">(参考起始坐标: {loc["anchorStart"]})</span></div>', unsafe_allow_html=True)
                    
                    for v_idx, cand in enumerate(loc["variants"]):
                        is_primary = (v_idx == 0)
                        card_class = "card-primary" if is_primary else ""
                        role_badge = '<span class="role-badge-main">主力优选</span>' if is_primary else f'<span class="role-badge-sub">微调备选 {v_idx}</span>'
                        
                        def render_rows(title, var_list, is_probe=False):
                            html_str = ""
                            for idx, seq in enumerate(var_list):
                                lbl = f"{title} {idx+1}" if len(var_list) > 1 else title
                                row_style = "seq-probe-row" if is_probe else ""
                                html_str += f'<div class="seq-row {row_style}"><span class="seq-type">{lbl}:</span><span class="seq-string">5\'- {seq} -3\'</span><span class="seq-stats">Len: {len(seq)}bp | Tm: {calc_tm(seq)}°C | GC: {calc_gc(seq)}%</span></div>'
                            return html_str

                        def render_stats(stats, is_mix):
                            mix_str = '<span class="mix-badge">混合套数扣分(-35)</span>' if is_mix else ''
                            return f'<div class="stats-row"><span class="stat-badge bg-0">完全匹配(0): {stats["p0"]}%, {stats["m0"]}/{stats["total"]}</span><span class="stat-badge bg-1">错配1碱基(1): {stats["p1"]}%, {stats["m1"]}/{stats["total"]}</span><span class="stat-badge bg-2">错配2碱基(2): {stats["p2"]}%, {stats["m2"]}/{stats["total"]}</span><span class="stat-badge bg-3">错配3碱基(≥3): {stats["p3"]}%, {stats["m3p"]}/{stats["total"]}</span>{mix_str}</div>'

                        # 核心卡片渲染：保证无缩进、整段解析
                        card_html = (
                            f'<div class="card-box {card_class}">'
                            f'<div class="card-header">'
                            f'<div>{role_badge} <span style="font-size:16px; font-weight:bold; color:#2c3e50;">综合得分: </span><span class="score-badge">{cand["score"]:.1f} 分</span></div>'
                            f'<div style="font-size: 13px; color: #95a5a6;">精确定位: {cand["start"]} | 产物长度: {cand["size"]} bp | 探针间距: {cand["probeGap"]} bp</div>'
                            f'</div>'
                            f'{render_rows("Forward", cand["fwd"])}'
                            f'{render_stats(cand["fStats"], len(cand["fwd"]) > 1)}'
                            f'{render_rows(f"Probe [{cand[\'probeDir\']}]", cand["probe"], True)}'
                            f'{render_stats(cand["pStats"], len(cand["probe"]) > 1)}'
                            f'{render_rows("Reverse", cand["rev"])}'
                            f'{render_stats(cand["rStats"], len(cand["rev"]) > 1)}'
                            f'</div>'
                        )
                        st.markdown(card_html, unsafe_allow_html=True)
                        
                        # 原生折叠面板对应详情项
                        with st.expander("🔍 展开查看当前候选体系全维度评估分"):
                            det = cand["details"]
                            st.markdown(f'<div style="font-size:13px; color:#555; line-height:2.0;"><div class="details-row"><span>基础匹配分 (探针 3 倍权重)</span><span style="color:#27ae60; font-weight:bold;">+{det["base"]:.1f}</span></div><div class="details-row"><span>探针匹配度奖励 (完美率 ≥99%)</span><span style="color:#27ae60; font-weight:bold;">+{det["pBonus"]:.1f}</span></div><div class="details-row"><span>探针错配惩罚 (低于98%十倍扣除)</span><span style="color:{"#c0392b" if det["pPenalty"]!=0 else "#333"}; font-weight:bold;">{det["pPenalty"]:.1f}</span></div><div class="details-row"><span>F/R 混合套数扣分 (-35分/次，极力优选单套)</span><span style="color:{"#c0392b" if (det["mixF"]+det["mixR"])!=0 else "#333"}; font-weight:bold;">{det["mixF"] + det["mixR"]}</span></div><div class="details-row"><span>偏离规则体系总扣分 (GC/温差/同向间距过长等)</span><span style="color:{"#c0392b" if det["softPen"]!=0 else "#333"}; font-weight:bold;">{det["softPen"]:.1f}</span></div></div>', unsafe_allow_html=True)
