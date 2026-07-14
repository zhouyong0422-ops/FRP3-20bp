import streamlit as st
import pandas as pd
import numpy as np
import math
import re
import io
import time
from datetime import datetime
from collections import Counter
import altair as alt

# ==========================================
# 0. 页面全局配置与原版样式还原 (CSS Inject)
# ==========================================
st.set_page_config(
    page_title="V02 原料齐套 - 工业级全景长度竞争寻优系统 V17",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# 注入与原版 HTML 完全一致的视觉样式
st.markdown("""
<style>
    body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; color: #333; }
    .main-title { color: #2c3e50; border-bottom: 3px solid #1abc9c; padding-bottom: 10px; margin-top: 0; display: flex; justify-content: space-between; align-items: flex-end; font-size: 28px; font-weight: bold;}
    .badge { font-size: 14px; background: #8e44ad; color: white; padding: 5px 12px; border-radius: 15px; font-weight: bold; vertical-align: middle; margin-left: 10px;}
    
    .report-box { margin-top: 30px; background: #fff; padding: 30px; border: 1px solid #e0e0e0; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.05); position: relative;}
    .report-header { text-align: center; border-bottom: 2px dashed #ccc; padding-bottom: 20px; margin-bottom: 20px; position: relative; }
    .report-header h2 { margin: 0; color: #2c3e50; font-size: 24px;}
    .report-header p { margin: 5px 0 0 0; color: #7f8c8d; }
    
    .locus-group { margin-bottom: 40px; background: #fdfefe; border: 1px solid #d1d5db; border-radius: 8px; overflow: hidden; box-shadow: 0 4px 6px rgba(0,0,0,0.02);}
    .locus-title { background: #edf2f7; color: #2c3e50; margin: 0; padding: 15px 20px; font-size: 18px; border-bottom: 1px solid #d1d5db; display: flex; justify-content: space-between; align-items: center; font-weight: bold;}
    
    .candidate-card { padding: 20px; border-bottom: 1px solid #ecf0f1; position: relative; }
    .candidate-card:last-child { border-bottom: none; }
    .cand-primary { border-left: 6px solid #27ae60; background: #f9fbfd;} 
    .cand-variant { border-left: 6px solid #bdc3c7; background: #fff;}    
    
    .candidate-card h4 { margin-top: 0; color: #2c3e50; font-size: 16px; border-bottom: 1px dashed #e1e8ed; padding-bottom: 8px; display: flex; justify-content: space-between; align-items: center;}
    .score-badge { background: #f39c12; color: white; padding: 3px 10px; border-radius: 12px; font-size: 13px; font-weight: bold;}
    .mix-badge { background: #9b59b6; color: white; padding: 2px 8px; border-radius: 8px; font-size: 11px; margin-left: 8px; vertical-align: middle;}
    .role-badge { font-size: 12px; padding: 2px 8px; border-radius: 4px; margin-right: 10px; color: white; font-weight: bold;}
    .role-main { background: #27ae60; }
    .role-sub { background: #7f8c8d; }
    
    .seq-block { margin-bottom: 12px; background: #fff; padding: 8px 12px; border-radius: 6px; border: 1px solid #ecf0f1; }
    .probe-block { border-left: 4px solid #8e44ad; background: #fdfafb; }
    
    .seq-row { display: flex; justify-content: space-between; font-family: monospace; font-size: 14px; margin-bottom: 4px;}
    .seq-type { font-weight: bold; width: 110px; display: inline-block; color: #34495e;}
    .seq-string { color: #d35400; letter-spacing: 1px; font-weight: bold; flex-grow: 1;}
    .seq-stats { color: #7f8c8d; font-size: 12px; text-align: right;}
    
    .mismatch-stats { display: flex; flex-wrap: wrap; gap: 8px; font-size: 12px; margin-top: 6px; padding-top: 6px; border-top: 1px dashed #ecf0f1; }
    .stat-badge { padding: 3px 8px; border-radius: 12px; font-weight: bold; display: flex; align-items: center; gap: 4px;}
    .bg-0 { background: #d4edda; color: #155724; border: 1px solid #c3e6cb;}
    .bg-1 { background: #fff3cd; color: #856404; border: 1px solid #ffeeba;}
    .bg-2 { background: #ffeeba; color: #856404; border: 1px solid #ffdf7e;}
    .bg-3 { background: #f8d7da; color: #721c24; border: 1px solid #f5c6cb;}

    details { margin-top: 10px; background: #fff; border: 1px solid #dcdde1; border-radius: 6px; padding: 6px 12px; }
    summary { cursor: pointer; font-size: 13px; color: #2980b9; font-weight: bold; outline: none; }
    summary:hover { color: #1f6391; }
    .details-content { font-size: 12px; color: #555; margin-top: 8px; line-height: 1.8; border-top: 1px dashed #eee; padding-top: 6px;}
    .score-item { display: flex; justify-content: space-between; }
    .score-plus { color: #27ae60; font-weight: bold;}
    .score-minus { color: #c0392b; font-weight: bold;}
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-title"><span>🧬 自动化引物探针柔性优先寻优系统</span><span class="badge">全量长度竞争 V17</span></div>', unsafe_allow_html=True)

# ==========================================
# 1. 核心热力学与生物信息学底板算法
# ==========================================
def calc_tm(seq: str) -> str:
    g = seq.count('G')
    c = seq.count('C')
    return f"{(64.9 + 41 * (g + c - 16.4) / len(seq)):.1f}"

def calc_gc(seq: str) -> str:
    g = seq.count('G')
    c = seq.count('C')
    return f"{(((g + c) / len(seq)) * 100):.1f}"

def reverse_complement(seq: str) -> str:
    trans_map = str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'})
    return seq[::-1].translate(trans_map)

def has_secondary_structure_risk(seq: str) -> bool:
    length = len(seq)
    if length < 12:
        return False
    for i in range(4, 6):
        head = seq[:i]
        tail_comp = reverse_complement(seq[-i:])
        if head == tail_comp:
            return True
    return False

def has_3prime_dimer_risk(seq1: str, seq2: str) -> bool:
    end1 = seq1[-4:]
    end2_comp = reverse_complement(seq2[-4:])
    return end1 == end2_comp

def is_hard_valid_oligo(seq: str, is_probe: bool = False) -> bool:
    if 'N' in seq or re.search(r'[^ATGC]', seq):
        return False
    gc = float(calc_gc(seq))
    if gc < 20 or gc > 80:
        return False
    if is_probe and seq.startswith('G'):
        return False
    return True

def calc_soft_penalties(seq: str, is_probe: bool) -> float:
    penalty = 0.0
    gc = float(calc_gc(seq))
    if re.search(r'([ATGC])\1{3,}', seq):
        penalty += 10.0
    if has_secondary_structure_risk(seq):
        penalty += 15.0

    if is_probe:
        if gc < 30: penalty += (30.0 - gc) * 2.0
        if gc > 65: penalty += (gc - 65.0) * 2.0
        if seq.count('G') >= seq.count('C'):
            penalty += 8.0
    else:
        if gc < 40: penalty += (40.0 - gc) * 2.0
        if gc > 60: penalty += (gc - 60.0) * 2.0
        end5 = seq[-5:]
        end_gc = end5.count('G') + end5.count('C')
        if end_gc < 1 or end_gc > 2:
            penalty += 6.0
        if seq.endswith('T'):
            penalty += 12.0
        if seq.endswith('GG') or seq.endswith('CC'):
            penalty += 8.0
    return penalty

def get_top_variants(start_index: int, length: int, sequences_array: list, max_variants: int = 2) -> list:
    total_seq = len(sequences_array)
    valid_seqs = []
    for seq in sequences_array:
        s = seq[start_index : start_index + length]
        if '-' not in s and 'N' not in s:
            valid_seqs.append(s)
            
    total_valid = len(valid_seqs)
    if total_seq == 0 or (total_valid / total_seq) < 0.90:
        return []
        
    counts = Counter(valid_seqs)
    if not counts:
        return []

    sorted_counts = counts.most_common()
    variants = [sorted_counts[0][0]]
    coverage = sorted_counts[0][1] / total_valid

    if coverage < 0.97 and len(sorted_counts) > 1 and max_variants > 1:
        second_cov = sorted_counts[1][1] / total_valid
        if second_cov > 0.04:
            variants.append(sorted_counts[1][0])
            
    return variants

def calculate_mix_mismatch(variants_array: list, start_index: int, sequences_array: list) -> dict:
    total_seq = len(sequences_array)
    stats = {'m0': 0, 'm1': 0, 'm2': 0, 'm3p': 0, 'total': total_seq}
    seq_len = len(variants_array[0])

    for i in range(total_seq):
        lib_seq_snippet = sequences_array[i][start_index : start_index + seq_len]
        if '-' in lib_seq_snippet or 'N' in lib_seq_snippet:
            stats['m3p'] += 1
            continue
            
        best_mismatches = seq_len
        for target_seq in variants_array:
            # 高性能矢量化对比
            mismatches = sum(1 for a, b in zip(target_seq, lib_seq_snippet) if a != b)
            if mismatches < best_mismatches:
                best_mismatches = mismatches
                
        if best_mismatches == 0: stats['m0'] += 1
        elif best_mismatches == 1: stats['m1'] += 1
        elif best_mismatches == 2: stats['m2'] += 1
        else: stats['m3p'] += 1

    if stats['total'] == 0:
        return {'p0': '0.0', 'p1': '0.0', 'p2': '0.0', 'p3': '0.0', 'm0': 0, 'm1': 0, 'm2': 0, 'm3p': 0, 'total': 0}
        
    return {
        'p0': round((stats['m0'] / stats['total']) * 100, 1),
        'p1': f"{((stats['m1'] / stats['total']) * 100):.1f}",
        'p2': f"{((stats['m2'] / stats['total']) * 100):.1f}",
        'p3': f"{((stats['m3p'] / stats['total']) * 100):.1f}",
        'm0': stats['m0'], 'm1': stats['m1'], 'm2': stats['m2'], 'm3p': stats['m3p'], 'total': stats['total']
    }

# =========================================================================
# 2. 内存优化版计算引擎 (加入局部 Memoization 缓存防止内存/CPU暴涨)
# =========================================================================
@st.cache_data(show_spinner=False)
def run_pipeline_engine(sequences: tuple):
    seq_len = len(sequences[0])
    
    # --- 步骤 A: 计算香农熵 ---
    entropies = []
    for i in range(seq_len):
        col_bases = [seq[i] for seq in sequences if seq[i] != '-']
        total = len(col_bases)
        if total == 0:
            entropies.append(0.0)
            continue
        counts = Counter(col_bases)
        entropy = 0.0
        for count in counts.values():
            p = count / total
            entropy -= p * math.log2(p)
        entropies.append(entropy)

    # --- 步骤 B: 构建高频重复查询的局部缓存 (极致降低内存和循环耗时) ---
    memo_fwd_probes = {}
    memo_rev_probes = {}
    memo_rev_primers = {}

    def get_all_valid_variants_memo(start_idx: int, is_probe: bool):
        cache_dict = memo_fwd_probes if is_probe else None
        if cache_dict is not None and start_idx in cache_dict:
            return cache_dict[start_idx]

        valid_list = []
        target_lengths = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30] if is_probe else [18, 19, 20, 21, 22, 23, 24, 25]
        for l in target_lengths:
            if start_idx + l > seq_len: continue
            raw_variants = get_top_variants(start_idx, l, sequences, 1 if is_probe else 2)
            if not raw_variants: continue

            final_variants = []
            all_passed = True
            for v in raw_variants:
                if not is_hard_valid_oligo(v, is_probe):
                    all_passed = False; break
                final_variants.append(v)

            if all_passed:
                soft_pen = sum(calc_soft_penalties(v, is_probe) for v in final_variants)
                avg_tm = sum(float(calc_tm(v)) for v in final_variants) / len(final_variants)
                tm_bonus = 0.0
                if is_probe:
                    if 62 <= avg_tm <= 68: tm_bonus = 15.0
                else:
                    if 55 <= avg_tm <= 60: tm_bonus = 15.0
                valid_list.append({'length': l, 'variants': final_variants, 'score': tm_bonus - soft_pen})

        valid_list.sort(key=lambda x: x['score'], reverse=True)
        if cache_dict is not None: cache_dict[start_idx] = valid_list
        return valid_list

    def get_reverse_valid_variants_memo(start_idx: int, is_probe: bool = False):
        cache_dict = memo_rev_probes if is_probe else memo_rev_primers
        if start_idx in cache_dict:
            return cache_dict[start_idx]

        valid_list = []
        target_lengths = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30] if is_probe else [18, 19, 20, 21, 22, 23, 24, 25]
        for l in target_lengths:
            if start_idx + l > seq_len: continue
            raw_variants = get_top_variants(start_idx, l, sequences, 1 if is_probe else 2)
            if not raw_variants: continue

            final_variants = []
            all_passed = True
            for rv in raw_variants:
                comp = reverse_complement(rv)
                if not is_hard_valid_oligo(comp, is_probe):
                    all_passed = False; break
                final_variants.append(comp)

            if all_passed:
                soft_pen = sum(calc_soft_penalties(v, is_probe) for v in final_variants)
                avg_tm = sum(float(calc_tm(v)) for v in final_variants) / len(final_variants)
                tm_bonus = 0.0
                if is_probe:
                    if 62 <= avg_tm <= 68: tm_bonus = 15.0
                else:
                    if 55 <= avg_tm <= 60: tm_bonus = 15.0
                valid_list.append({'length': l, 'variants': final_variants, 'rawVariants': raw_variants, 'score': tm_bonus - soft_pen})

        valid_list.sort(key=lambda x: x['score'], reverse=True)
        cache_dict[start_idx] = valid_list
        return valid_list

    # --- 步骤 C: 严格执行全量长度矩阵竞争与单向间距寻优 ---
    all_candidates = []
    for i in range(seq_len - 70):
        f_obj_list = get_all_valid_variants_memo(i, False)
        if not f_obj_list: continue

        for f_obj in f_obj_list[:2]:
            f_variants = f_obj['variants']
            f_len = f_obj['length']

            min_r_start = max(i + 70 - 25, i + f_len + 1)
            max_r_start = min(i + 160 - 18, seq_len - 1)

            for r_start in range(min_r_start, max_r_start + 1):
                r_obj_list = get_reverse_valid_variants_memo(r_start, False)
                if not r_obj_list: continue

                for r_obj in r_obj_list[:2]:
                    r_variants_raw = r_obj['rawVariants']
                    r_variants = r_obj['variants']
                    r_len = r_obj['length']

                    amplicon_size = r_start + r_len - i
                    if amplicon_size < 70 or amplicon_size > 160: continue

                    valid_probes = []
                    # 3.1 寻找正向探针 (紧跟着同向的 Forward 引物，间距 < 30 bp)
                    for p_start in range(i + f_len + 1, min(i + f_len + 30, r_start - 18)):
                        p_obj_list = get_all_valid_variants_memo(p_start, True)
                        for p_obj in p_obj_list[:2]:
                            if p_start + p_obj['length'] <= r_start:
                                valid_probes.append({**p_obj, 'start': p_start, 'isRev': False, 'gapToSameDir': p_start - (i + f_len)})

                    # 3.2 寻找反向互补探针 (紧贴着同向的 Reverse 引物延伸，距离 < 30bp)
                    for p_end in range(r_start - 1, max(r_start - 30, i + f_len + 18), -1):
                        for p_len in [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]:
                            p_start = p_end - p_len
                            if p_start >= i + f_len and (r_start - p_end < 30):
                                p_obj_list = get_reverse_valid_variants_memo(p_start, True)
                                for p_obj in p_obj_list[:1]:
                                    if p_obj['length'] == p_len:
                                        valid_probes.append({**p_obj, 'start': p_start, 'isRev': True, 'gapToSameDir': r_start - p_end})

                    if not valid_probes: continue
                    valid_probes.sort(key=lambda x: x['score'], reverse=True)

                    for p_obj in valid_probes[:3]:
                        p_variants = p_obj['variants']
                        p_start = p_obj['start']

                        min_f_tm = min(float(calc_tm(v)) for v in f_variants)
                        max_f_tm = max(float(calc_tm(v)) for v in f_variants)
                        min_r_tm = min(float(calc_tm(v)) for v in r_variants)
                        max_r_tm = max(float(calc_tm(v)) for v in r_variants)
                        min_p_tm = min(float(calc_tm(v)) for v in p_variants)

                        primer_max_tm = max(max_f_tm, max_r_tm)
                        primer_min_tm = min(min_f_tm, min_r_tm)

                        soft_penalty = sum(calc_soft_penalties(f, False) for f in f_variants)
                        soft_penalty += sum(calc_soft_penalties(p, True) for p in p_variants)
                        soft_penalty += sum(calc_soft_penalties(r, False) for r in r_variants)

                        dimer_risk = False
                        for f in f_variants:
                            for r in r_variants:
                                if has_3prime_dimer_risk(f, r): dimer_risk = true_val = True

                        if dimer_risk: soft_penalty += 20.0
                        primer_tm_diff = abs(primer_max_tm - primer_min_tm)
                        if primer_tm_diff > 2.0: soft_penalty += (primer_tm_diff - 2.0) * 5.0
                        if min_p_tm < primer_max_tm + 5.0: soft_penalty += (primer_max_tm + 5.0 - min_p_tm) * 6.0
                        if p_obj['gapToSameDir'] > 10: soft_penalty += (p_obj['gapToSameDir'] - 10) * 1.5

                        f_stats = calculate_mix_mismatch(f_variants, i, sequences)
                        p_stats = calculate_mix_mismatch(p_variants, p_start, sequences)
                        r_stats = calculate_mix_mismatch(r_variants_raw, r_start, sequences)

                        f_p0 = float(f_stats['p0'])
                        p_p0 = float(p_stats['p0'])
                        r_p0 = float(r_stats['p0'])

                        base_score = f_p0 + (p_p0 * 3.0) + r_p0
                        probe_bonus = 50.0 if p_p0 >= 99.0 else 0.0
                        probe_penalty = (98.0 - p_p0) * 10.0 if p_p0 < 98.0 else 0.0
                        mix_f = -35.0 if len(f_variants) > 1 else 0.0
                        mix_r = -35.0 if len(r_variants) > 1 else 0.0

                        total_score = base_score + probe_bonus - probe_penalty + mix_f + mix_r - soft_penalty

                        all_candidates.append({
                            'fwd': f_variants, 'rev': r_variants, 'probe': p_variants,
                            'fStats': f_stats, 'pStats': p_stats, 'rStats': r_stats,
                            'size': amplicon_size, 'start': i, 'score': total_score,
                            'probeDir': "Reverse (-)" if p_obj['isRev'] else "Forward (+)",
                            'probeGap': p_obj['gapToSameDir'],
                            'details': {'base': base_score, 'pBonus': probe_bonus, 'pPenalty': -probe_penalty, 'mixF': mix_f, 'mixR': mix_r, 'softPen': -soft_penalty}
                        })

    all_candidates.sort(key=lambda x: x['score'], reverse=True)
    
    # --- 步骤 D: 归类分组 (Locus 聚类) ---
    global_loci_groups = []
    locus_window = 50
    for cand in all_candidates:
        found_locus = False
        for locus in global_loci_groups:
            if abs(cand['start'] - locus['anchorStart']) <= locus_window:
                if len(locus['variants']) < 3:
                    locus['variants'].append(cand)
                found_locus = True; break
        if not found_locus:
            global_loci_groups.append({'locusId': len(global_loci_groups) + 1, 'anchorStart': cand['start'], 'variants': [cand]})

    return entropies, global_loci_groups

# ==========================================
# 3. Streamlit 页面交互与数据生成模块
# ==========================================
st.markdown('<div style="background:#f8f9fa; padding:15px; border-radius:8px; border:2px dashed #bdc3c7; margin-bottom:15px;">📂 <b>V17 修复版</b>：探针与同向引物间距严格 &lt;30bp；与不同向引物不限距离；产物 70-160bp。</div>', unsafe_allow_html=True)

uploaded_file = st.file_uploader("导入 FASTA 序列库", type=["fasta", "fas", "txt", "aln"], label_visibility="collapsed")
fasta_text_input = st.text_area("或在此处直接粘贴比对完成的 FASTA 序列...", height=120)

input_text = ""
if uploaded_file is not None:
    input_text = uploaded_file.getvalue().decode("utf-8")
elif fasta_text_input.strip():
    input_text = fasta_text_input.strip()

if st.button("⚙️ 启动寻优：基于单向间距限制与长度竞争出具 DOE", use_container_width=True, type="primary"):
    if not input_text:
        st.error("请导入包含变异库的文件或在文本框中粘贴序列！")
    else:
        with st.spinner("⏳ 正在进行全量长度矩阵竞争与热力学加权计算，请耐心等待..."):
            lines = input_text.split('\n')
            sequences = []
            current_seq = ""
            for line in lines:
                if line.startswith('>'):
                    if current_seq: sequences.append(current_seq.upper())
                    current_seq = ""
                else:
                    current_seq += line.strip()
            if current_seq: sequences.append(current_seq.upper())

            if len(sequences) < 2:
                st.error("文件格式有误或序列数不足2条。")
            else:
                # 转换成 tuple 以便满足 Streamlit @st.cache_data 的可哈希化缓存要求
                entropies, global_loci_groups = run_pipeline_engine(tuple(sequences))
                
                st.session_state['entropies'] = entropies
                st.session_state['loci_groups'] = global_loci_groups
                st.session_state['report_time'] = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
                st.session_state['has_results'] = True

# ==========================================
# 4. 结果展示与报告渲染 (完全保证 1:1 视觉内容)
# ==========================================
if st.session_state.get('has_results', False):
    entropies = st.session_state['entropies']
    global_loci_groups = st.session_state['loci_groups']
    
    st.markdown("### 📊 靶标序列群变异强度扫描 (香农熵)")
    df_chart = pd.DataFrame({'位置': range(1, len(entropies) + 1), '香农熵': entropies})
    df_chart['变异强度'] = np.where(df_chart['香农熵'] > 0.05, '高变异 (>0.05)', '低变异 (≤0.05)')
    
    chart = alt.Chart(df_chart).mark_bar().encode(
        x=alt.X('位置:Q', title='序列物理参考坐标 (bp)'),
        y=alt.Y('香农熵:Q', title='Shannon Entropy'),
        color=alt.Color('变异强度:N', scale=alt.Scale(domain=['高变异 (>0.05)', '低变异 (≤0.05)'], range=['#e74c3c', '#3498db'])),
        tooltip=['位置', '香农熵']
    ).properties(height=250)
    st.altair_chart(chart, use_container_width=True)

    # 生成与原版完全契合的报告 HTML 数据流
    report_html = f"""
    <div class="report-box">
        <div class="report-header">
            <h2>V02 原料齐套：工业级多重引探 DOE 报告</h2>
            <p>生成时间：{st.session_state['report_time']}</p>
            <p style="font-size: 13px; color: #95a5a6;">规范内核：产物 70–160 bp；探针跟同向引物 &lt;30bp 且跟不同向不限制；支持正反向探针全量竞争。</p>
        </div>
    """
    
    if not global_loci_groups:
        report_html += '<div class="candidate-card" style="border-left-color: #e74c3c;"><h4 style="color:#c0392b;">⚠️ 体系设计失败</h4><p>在该序列库中未能找到产物 70-160bp、且满足探针与同向引物 &lt;30bp 约束的有效靶区。</p></div>'
    else:
        def generate_seq_rows(title, variants_array):
            html = ''
            for idx, v in enumerate(variants_array):
                label = f"{title} {idx + 1}" if len(variants_array) > 1 else title
                html += f"""
                <div class="seq-row">
                    <span class="seq-type">{label}:</span>
                    <span class="seq-string">5'- {v} -3'</span>
                    <span class="seq-stats">Len: {len(v)}bp | Tm: {calc_tm(v)}°C | GC: {calc_gc(v)}%</span>
                </div>"""
            return html

        def get_stat_html(stats, is_mix):
            mix_badge = '<span class="mix-badge">混合套数扣分(-35)</span>' if is_mix else ''
            return f"""
                <div class="mismatch-stats">
                    <span class="stat-badge bg-0">完全匹配(0): {stats['p0']}%, {stats['m0']}/{stats['total']}</span>
                    <span class="stat-badge bg-1">错配1碱基(1): {stats['p1']}%, {stats['m1']}/{stats['total']}</span>
                    <span class="stat-badge bg-2">错配2碱基(2): {stats['p2']}%, {stats['m2']}/{stats['total']}</span>
                    <span class="stat-badge bg-3">错配3碱基(≥3): {stats['p3']}%, {stats['m3p']}/{stats['total']}</span>
                    {mix_badge}
                </div>"""

        for locus in global_loci_groups:
            report_html += f"""
            <div class="locus-group">
                <h3 class="locus-title">
                    <span>🎯 独立黄金靶区 {locus['locusId']}</span>
                    <span style="font-size: 14px; color: #7f8c8d; font-weight: normal;">(参考起始坐标: {locus['anchorStart']})</span>
                </h3>"""
                
            for v_idx, cand in enumerate(locus['variants']):
                is_primary = (v_idx == 0)
                card_class = 'cand-primary' if is_primary else 'cand-variant'
                role_badge = '<span class="role-badge role-main">主力优选</span>' if is_primary else f'<span class="role-badge role-sub">微调备选 {v_idx}</span>'

                report_html += f"""
                <div class="candidate-card {card_class}">
                    <h4>
                        <span>{role_badge} 综合得分: <span class="score-badge">{cand['score']:.1f} 分</span></span>
                        <span style="font-size: 13px; font-weight: normal; color: #95a5a6;">精确定位: {cand['start']} | 产物长度: {cand['size']} bp | 探针间距: {cand['probeGap']} bp</span>
                    </h4>
                    
                    <div class="seq-block">
                        {generate_seq_rows('Forward', cand['fwd'])}
                        {get_stat_html(cand['fStats'], len(cand['fwd']) > 1)}
                    </div>
                    <div class="seq-block probe-block">
                        {generate_seq_rows(f"Probe [{cand['probeDir']}]", cand['probe'])}
                        {get_stat_html(cand['pStats'], len(cand['probe']) > 1)}
                    </div>
                    <div class="seq-block">
                        {generate_seq_rows('Reverse', cand['rev'])}
                        {get_stat_html(cand['rStats'], len(cand['rev']) > 1)}
                    </div>

                    <details>
                        <summary>🔍 查看当前变体柔性评分明细 (单向间距优化版)</summary>
                        <div class="details-content">
                            <div class="score-item"><span>基础匹配分 (探针 3 倍权重)</span><span class="score-plus">+{cand['details']['base']:.1f}</span></div>
                            <div class="score-item"><span>探针卓越奖励 (完美率 ≥99%)</span><span class="score-plus">+{cand['details']['pBonus']:.1f}</span></div>
                            <div class="score-item"><span>探针错配惩罚 (低于98%十倍扣除)</span><span class="{'' if cand['details']['pPenalty'] == 0 else 'score-minus'}">{cand['details']['pPenalty']:.1f}</span></div>
                            <div class="score-item"><span>F/R 混合套数重罚 (-35分/次，极力优选单套)</span><span class="{'' if (cand['details']['mixF'] + cand['details']['mixR']) == 0 else 'score-minus'}">{cand['details']['mixF'] + cand['details']['mixR']}</span></div>
                            <div class="score-item" style="border-top: 1px solid #ddd; padding-top: 4px; font-weight: bold;"><span>柔性规则偏离总扣分 (GC/温差/同向间距过长等)</span><span class="{'' if cand['details']['softPen'] == 0 else 'score-minus'}">{cand['details']['softPen']:.1f}</span></div>
                        </div>
                    </details>
                </div>"""
            report_html += "</div>"
    report_html += "</div>"
    
    # 渲染结果 HTML
    st.markdown(report_html, unsafe_allow_html=True)
    
    # --- 准备 CSV 导出流 ---
    csv_content = "\uFEFF靶区归属,变体角色,综合得分,寡核苷酸类型,序列 (5'->3'),长度 (bp),Tm (°C),GC (%),完美匹配(0),错配1碱基(1),错配2碱基(2),错配3碱基(≥3),预期产物长度 (bp),同向引物间距 (bp),精确起始坐标\n"
    for locus in global_loci_groups:
        for v_idx, cand in enumerate(locus['variants']):
            locus_name = f"靶区_{locus['locusId']}"
            role = "主力优选" if v_idx == 0 else f"备选_{v_idx}"
            score = f"{cand['score']:.1f}"
            size = cand['size']
            start = cand['start']
            gap = cand['probeGap']

            for i, seq in enumerate(cand['fwd']):
                type_name = f"Forward_{i+1}" if len(cand['fwd']) > 1 else "Forward"
                stats = cand['fStats']
                csv_content += f"{locus_name},{role},{score},{type_name},{seq},{len(seq)},{calc_tm(seq)},{calc_gc(seq)},{stats['p0']}% ({stats['m0']}/{stats['total']}),{stats['p1']}% ({stats['m1']}/{stats['total']}),{stats['p2']}% ({stats['m2']}/{stats['total']}),{stats['p3']}% ({stats['m3p']}/{stats['total']}),{size},-,{start}\n"

            for i, seq in enumerate(cand['probe']):
                type_base = f"Probe({cand['probeDir']})"
                type_name = f"{type_base}_{i+1}" if len(cand['probe']) > 1 else type_base
                stats = cand['pStats']
                csv_content += f"{locus_name},{role},{score},{type_name},{seq},{len(seq)},{calc_tm(seq)},{calc_gc(seq)},{stats['p0']}% ({stats['m0']}/{stats['total']}),{stats['p1']}% ({stats['m1']}/{stats['total']}),{stats['p2']}% ({stats['m2']}/{stats['total']}),{stats['p3']}% ({stats['m3p']}/{stats['total']}),{size},{gap},{start}\n"

            for i, seq in enumerate(cand['rev']):
                type_name = f"Reverse_{i+1}" if len(cand['rev']) > 1 else "Reverse"
                stats = cand['rStats']
                csv_content += f"{locus_name},{role},{score},{type_name},{seq},{len(seq)},{calc_tm(seq)},{calc_gc(seq)},{stats['p0']}% ({stats['m0']}/{stats['total']}),{stats['p1']}% ({stats['m1']}/{stats['total']}),{stats['p2']}% ({stats['m2']}/{stats['total']}),{stats['p3']}% ({stats['m3p']}/{stats['total']}),{size},-,{start}\n"

    # 在顶部右侧布局专属下载按钮
    date_str = datetime.now().strftime("%Y-%m-%d")
    st.download_button(
        label="📥 导出 Excel (含靶区分组)",
        data=csv_content.encode('utf-8'),
        file_name=f"V02_单向间距竞争优选_DOE清单_{date_str}.csv",
        mime="text/csv",
        type="primary"
    )
