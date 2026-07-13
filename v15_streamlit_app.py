import streamlit as st
import pandas as pd
import numpy as np
from io import StringIO
import plotly.express as px
import plotly.graph_objects as go
import time
from datetime import datetime

st.set_page_config(page_title="V02 原料齐套 - 工业级柔性优先寻优系统 V15", layout="wide")

st.title("🧬 自动化引物探针柔性优先寻优系统")
st.markdown('<span style="font-size: 14px; background: #2c3e50; color: #f1c40f; padding: 5px 12px; border-radius: 15px; font-weight: bold; border: 1px solid #f1c40f;">纯净单套优选 V15</span>', unsafe_allow_html=True)

st.markdown("""
**规范内核**：绝对死守 70–150 bp 产物（无偏向）、探针 5'端排 G 及长度边界；严重扣分遏制 F/R 双套混合，保证生产极简。
""")

# Sidebar for inputs
with st.sidebar:
    st.header("输入序列")
    uploaded_file = st.file_uploader("📂 导入 FASTA 序列库", type=["fasta", "fas", "txt", "aln"])
    
    fasta_text = st.text_area("或在此处直接粘贴比对完成的 FASTA 序列...", height=300)

if uploaded_file is not None:
    fasta_text = uploaded_file.getvalue().decode("utf-8")
    st.success(f"已加载序列库: {uploaded_file.name}")

# Main run button
if st.button("⚙️ 启动寻优：基于纯净单套准则出具靶区 DOE", type="primary", use_container_width=True):
    if not fasta_text.strip():
        st.error("请导入包含变异库的文件或粘贴序列！")
        st.stop()
    
    with st.spinner("⏳ 正在进行全景矩阵与单套优先加权计算，请耐心等待..."):
        # Execute the logic (ported from JS)
        result = execute_logic(fasta_text)
    
    if result:
        display_results(result)

def calc_tm(seq):
    g = seq.count('G')
    c = seq.count('C')
    return round(64.9 + 41 * (g + c - 16.4) / len(seq), 1)

def calc_gc(seq):
    g = seq.count('G')
    c = seq.count('C')
    return round(((g + c) / len(seq)) * 100, 1)

def reverse_complement(seq):
    map_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(map_dict.get(base, base) for base in reversed(seq))

def has_secondary_structure_risk(seq):
    if len(seq) < 12:
        return False
    for i in range(4, 6):
        head = seq[:i]
        tail_comp = reverse_complement(seq[-i:])
        if head == tail_comp:
            return True
    return False

def has_3prime_dimer_risk(seq1, seq2):
    end1 = seq1[-4:]
    end2_comp = reverse_complement(seq2[-4:])
    return end1 == end2_comp

def is_hard_valid_oligo(seq, is_probe=False):
    if 'N' in seq or not all(base in 'ATGC' for base in seq):
        return False
    gc = calc_gc(seq)
    if gc < 20 or gc > 80:
        return False
    if is_probe and seq[0] == 'G':
        return False
    return True

def calc_soft_penalties(seq, is_probe):
    penalty = 0
    gc = calc_gc(seq)
    
    if any(seq.count(base * 4) > 0 for base in 'ATGC'):  # homopolymer
        penalty += 10
    if has_secondary_structure_risk(seq):
        penalty += 15

    if is_probe:
        if gc < 30:
            penalty += (30 - gc) * 2
        if gc > 65:
            penalty += (gc - 65) * 2
        g_count = seq.count('G')
        c_count = seq.count('C')
        if g_count >= c_count:
            penalty += 8
    else:
        if gc < 40:
            penalty += (40 - gc) * 2
        if gc > 60:
            penalty += (gc - 60) * 2
        end5 = seq[-5:]
        end_gc = sum(1 for b in end5 if b in 'GC')
        if end_gc < 1 or end_gc > 2:
            penalty += 6
        if seq[-1] == 'T':
            penalty += 12
        if seq.endswith('GG') or seq.endswith('CC'):
            penalty += 8
    return penalty

def get_top_variants(start_index, length, sequences, max_variants=2):
    counts = {}
    total_valid = 0
    total_seq = len(sequences)
    
    for seq in sequences:
        snippet = seq[start_index:start_index + length]
        if '-' in snippet or 'N' in snippet:
            continue
        counts[snippet] = counts.get(snippet, 0) + 1
        total_valid += 1
    
    if total_valid / total_seq < 0.90 or not counts:
        return []
    
    sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    variants = [sorted_counts[0][0]]
    coverage = sorted_counts[0][1] / total_valid
    
    if coverage < 0.97 and len(sorted_counts) > 1 and max_variants > 1:
        second_cov = sorted_counts[1][1] / total_valid
        if second_cov > 0.04:
            variants.append(sorted_counts[1][0])
    return variants

def get_all_valid_variants(start_idx, sequences, is_probe):
    valid_list = []
    target_lengths = [20, 22, 24, 18, 19, 21, 23, 25, 26, 27, 28, 29, 30] if is_probe else [20, 21, 22, 19, 18, 23, 24, 25]
    
    for length in target_lengths:
        if start_idx + length > len(sequences[0]):
            continue
        raw_variants = get_top_variants(start_idx, length, sequences, 1 if is_probe else 2)
        if not raw_variants:
            continue
        
        final_variants = []
        all_passed = True
        for v in raw_variants:
            if not is_hard_valid_oligo(v, is_probe):
                all_passed = False
                break
            final_variants.append(v)
        
        if all_passed:
            valid_list.append({'length': length, 'variants': final_variants})
            break
    return valid_list

def calculate_mix_mismatch(variants, start_index, sequences):
    total_seq = len(sequences)
    stats = {'m0': 0, 'm1': 0, 'm2': 0, 'm3p': 0, 'total': total_seq}
    seq_len = len(variants[0])
    
    for seq in sequences:
        lib_snippet = seq[start_index:start_index + seq_len]
        if '-' in lib_snippet or 'N' in lib_snippet:
            stats['m3p'] += 1
            continue
        
        best_mismatches = seq_len
        for target in variants:
            mismatches = sum(a != b for a, b in zip(target, lib_snippet))
            if mismatches < best_mismatches:
                best_mismatches = mismatches
        
        if best_mismatches == 0:
            stats['m0'] += 1
        elif best_mismatches == 1:
            stats['m1'] += 1
        elif best_mismatches == 2:
            stats['m2'] += 1
        else:
            stats['m3p'] += 1
    
    if stats['total'] == 0:
        return {'p0': 0.0, 'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'm0':0, 'm1':0, 'm2':0, 'm3p':0, 'total':0}
    
    return {
        'p0': round((stats['m0'] / stats['total']) * 100, 1),
        'p1': round((stats['m1'] / stats['total']) * 100, 1),
        'p2': round((stats['m2'] / stats['total']) * 100, 1),
        'p3': round((stats['m3p'] / stats['total']) * 100, 1),
        'm0': stats['m0'], 'm1': stats['m1'], 'm2': stats['m2'], 'm3p': stats['m3p'], 'total': stats['total']
    }

def execute_logic(text):
    lines = text.split('\n')
    sequences = []
    current_seq = ""
    for line in lines:
        if line.startswith('>'):
            if current_seq:
                sequences.append(current_seq.upper())
            current_seq = ""
        else:
            current_seq += line.strip()
    if current_seq:
        sequences.append(current_seq.upper())
    
    if len(sequences) < 2:
        st.error("文件格式有误或序列数不足2条。")
        return None
    
    seq_len = len(sequences[0])
    
    # Entropy calculation for chart
    entropies = []
    for i in range(seq_len):
        column = {}
        total = 0
        for seq in sequences:
            base = seq[i] if i < len(seq) else '-'
            if base and base != '-':
                column[base] = column.get(base, 0) + 1
                total += 1
        entropy = 0
        for count in column.values():
            p = count / total
            entropy -= p * np.log2(p)
        entropies.append(entropy)
    
    # Store for display
    result = {
        'entropies': entropies,
        'sequences': sequences,
        'seq_len': seq_len,
        'candidates': []
    }
    
    min_gap = 1
    max_gap = 15
    
    all_candidates = []
    for i in range(seq_len - 150):
        f_obj_list = get_all_valid_variants(i, sequences, False)
        if not f_obj_list:
            continue
        f_obj = f_obj_list[0]
        f_variants = f_obj['variants']
        f_len = f_obj['length']
        
        for gap1 in range(min_gap, max_gap + 1):
            p_start = i + f_len + gap1
            p_obj_list = get_all_valid_variants(p_start, sequences, True)
            if not p_obj_list:
                continue
            p_obj = p_obj_list[0]
            p_variants = p_obj['variants']
            p_len = p_obj['length']
            
            for gap2 in range(min_gap, max_gap + 1):
                r_start = p_start + p_len + gap2
                r_found = False
                r_variants_raw = []
                r_variants = []
                r_len = 0
                
                target_lengths = [20, 21, 22, 19, 18, 23, 24, 25]
                for test_len in target_lengths:
                    if r_start + test_len > seq_len:
                        continue
                    raw_vars = get_top_variants(r_start, test_len, sequences, 2)
                    if not raw_vars:
                        continue
                    
                    all_passed = True
                    processed_rev = []
                    for rv in raw_vars:
                        comp = reverse_complement(rv)
                        if not is_hard_valid_oligo(comp, False):
                            all_passed = False
                            break
                        processed_rev.append(comp)
                    
                    if all_passed:
                        r_found = True
                        r_variants_raw = raw_vars
                        r_variants = processed_rev
                        r_len = test_len
                        break
                
                if not r_found:
                    continue
                
                amplicon_size = r_start + r_len - i
                if amplicon_size < 70 or amplicon_size > 150:
                    continue
                
                # Tm calculations
                min_f_tm = min(calc_tm(v) for v in f_variants)
                max_f_tm = max(calc_tm(v) for v in f_variants)
                min_r_tm = min(calc_tm(v) for v in r_variants)
                max_r_tm = max(calc_tm(v) for v in r_variants)
                min_p_tm = min(calc_tm(v) for v in p_variants)
                
                primer_max_tm = max(max_f_tm, max_r_tm)
                primer_min_tm = min(min_f_tm, min_r_tm)
                
                # Soft penalties
                soft_penalty = 0
                for f in f_variants:
                    soft_penalty += calc_soft_penalties(f, False)
                for p in p_variants:
                    soft_penalty += calc_soft_penalties(p, True)
                for r in r_variants:
                    soft_penalty += calc_soft_penalties(r, False)
                
                dimer_risk = any(has_3prime_dimer_risk(f, r) for f in f_variants for r in r_variants)
                if dimer_risk:
                    soft_penalty += 20
                
                primer_tm_diff = abs(primer_max_tm - primer_min_tm)
                if primer_tm_diff > 2.0:
                    soft_penalty += (primer_tm_diff - 2.0) * 5
                
                if min_p_tm < primer_max_tm + 5.0:
                    soft_penalty += (primer_max_tm + 5.0 - min_p_tm) * 6
                
                if gap1 > 5:
                    soft_penalty += (gap1 - 5) * 2
                if gap2 > 5:
                    soft_penalty += (gap2 - 5) * 2
                
                f_stats = calculate_mix_mismatch(f_variants, i, sequences)
                p_stats = calculate_mix_mismatch(p_variants, p_start, sequences)
                r_stats = calculate_mix_mismatch(r_variants_raw, r_start, sequences)
                
                f_p0 = f_stats['p0']
                p_p0 = p_stats['p0']
                r_p0 = r_stats['p0']
                
                base_score = f_p0 + (p_p0 * 3) + r_p0
                probe_bonus = 50 if p_p0 >= 99.0 else 0
                probe_penalty = (98.0 - p_p0) * 10 if p_p0 < 98.0 else 0
                
                mix_f = -35 if len(f_variants) > 1 else 0
                mix_r = -35 if len(r_variants) > 1 else 0
                
                total_score = base_score + probe_bonus - probe_penalty + mix_f + mix_r - soft_penalty
                
                all_candidates.append({
                    'fwd': f_variants,
                    'rev': r_variants,
                    'probe': p_variants,
                    'f_stats': f_stats,
                    'p_stats': p_stats,
                    'r_stats': r_stats,
                    'size': amplicon_size,
                    'start': i,
                    'score': total_score,
                    'details': {
                        'base': base_score,
                        'p_bonus': probe_bonus,
                        'p_penalty': -probe_penalty,
                        'mix_f': mix_f,
                        'mix_r': mix_r,
                        'soft_pen': -soft_penalty
                    }
                })
    
    all_candidates.sort(key=lambda x: x['score'], reverse=True)
    
    # Group into loci
    global_loci_groups = []
    locus_window = 50
    for cand in all_candidates:
        found_locus = False
        for locus in global_loci_groups:
            if abs(cand['start'] - locus['anchor_start']) <= locus_window:
                if len(locus['variants']) < 3:
                    locus['variants'].append(cand)
                found_locus = True
                break
        if not found_locus:
            global_loci_groups.append({
                'locus_id': len(global_loci_groups) + 1,
                'anchor_start': cand['start'],
                'variants': [cand]
            })
    
    result['loci_groups'] = global_loci_groups
    return result

def display_results(result):
    st.subheader("📊 靶标序列群变异强度扫描 (香农熵)")
    df_entropy = pd.DataFrame({
        'Position': list(range(1, len(result['entropies']) + 1)),
        'Entropy': result['entropies']
    })
    fig = px.bar(df_entropy, x='Position', y='Entropy', 
                 title='香农熵 (突变强度扫描)',
                 color_discrete_sequence=['#3498db' if e <= 0.05 else '#e74c3c' for e in result['entropies']])
    st.plotly_chart(fig, use_container_width=True)
    
    st.subheader("V02 原料齐套：工业级多重引探 DOE 报告")
    st.caption(f"出具时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if not result['loci_groups']:
        st.error("⚠️ 体系设计失败：在绝对硬底线限制下，该序列库未能找到完整无缺失(Gap<10%)的组合区域。")
        return
    
    for locus in result['loci_groups']:
        st.markdown(f"### 🎯 独立黄金靶区 {locus['locus_id']} (参考起始坐标: {locus['anchor_start']})")
        
        for v_idx, cand in enumerate(locus['variants']):
            is_primary = v_idx == 0
            role = "主力优选" if is_primary else f"微调备选 {v_idx}"
            card_color = "#27ae60" if is_primary else "#bdc3c7"
            
            with st.expander(f"{role} - 综合得分: {cand['score']:.1f} 分 | 产物长度: {cand['size']} bp | 起始: {cand['start']}", expanded=is_primary):
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.markdown("**Forward**")
                    for idx, seq in enumerate(cand['fwd']):
                        st.code(f"5'- {seq} -3'\nLen: {len(seq)}bp | Tm: {calc_tm(seq)}°C | GC: {calc_gc(seq)}%", language="text")
                    display_stats(cand['f_stats'], len(cand['fwd']) > 1)
                
                with col2:
                    st.markdown("**Probe**")
                    for idx, seq in enumerate(cand['probe']):
                        st.code(f"5'- {seq} -3'\nLen: {len(seq)}bp | Tm: {calc_tm(seq)}°C | GC: {calc_gc(seq)}%", language="text")
                    display_stats(cand['p_stats'], len(cand['probe']) > 1)
                
                with col3:
                    st.markdown("**Reverse**")
                    for idx, seq in enumerate(cand['rev']):
                        st.code(f"5'- {seq} -3'\nLen: {len(seq)}bp | Tm: {calc_tm(seq)}°C | GC: {calc_gc(seq)}%", language="text")
                    display_stats(cand['r_stats'], len(cand['rev']) > 1)
                
                # Details
                st.markdown("**🔍 柔性评分明细**")
                details = cand['details']
                cols = st.columns(2)
                cols[0].metric("基础匹配分", f"+{details['base']:.1f}")
                cols[1].metric("探针卓越奖励", f"+{details['p_bonus']:.1f}")
                cols[0].metric("探针错配惩罚", f"{details['p_penalty']:.1f}")
                cols[1].metric("F/R 混合套数重罚", f"{details['mix_f'] + details['mix_r']}")
                st.metric("柔性规则偏离总扣分", f"{details['soft_pen']:.1f}", delta=None)
    
    # Download
    if st.button("📥 导出 Excel (含靶区分组)"):
        csv_data = generate_csv(result)
        st.download_button(
            label="下载 CSV",
            data=csv_data,
            file_name=f"V02_纯净单套柔性加权_DOE清单_{datetime.now().strftime('%Y-%m-%d')}.csv",
            mime="text/csv"
        )

def display_stats(stats, is_mix):
    st.markdown(f"""
    **匹配统计**  
    完全匹配(0): {stats['p0']}% ({stats['m0']}/{stats['total']})  
    错配1: {stats['p1']}%  
    错配2: {stats['p2']}%  
    错配≥3: {stats['p3']}%
    """)
    if is_mix:
        st.warning("混合套数扣分(-35)")

def generate_csv(result):
    output = StringIO()
    output.write("\uFEFF")
    output.write("靶区归属,变体角色,综合得分,寡核苷酸类型,序列 (5'->3'),长度 (bp),Tm (°C),GC (%),完美匹配(0),错配1碱基(1),错配2碱基(2),错配3碱基(≥3),预期产物长度 (bp),精确起始坐标\n")
    
    for locus in result['loci_groups']:
        for v_idx, cand in enumerate(locus['variants']):
            locus_name = f"靶区_{locus['locus_id']}"
            role = "主力优选" if v_idx == 0 else f"备选_{v_idx}"
            score = f"{cand['score']:.1f}"
            size = cand['size']
            start = cand['start']
            
            for seq in cand['fwd']:
                typ = f"Forward_{len(cand['fwd'])}" if len(cand['fwd']) > 1 else "Forward"
                stats = cand['f_stats']
                output.write(f"{locus_name},{role},{score},{typ},{seq},{len(seq)},{calc_tm(seq)},{calc_gc(seq)},{stats['p0']}% ({stats['m0']}/{stats['total']}),{stats['p1']}% ({stats['m1']}/{stats['total']}),{stats['p2']}% ({stats['m2']}/{stats['total']}),{stats['p3']}% ({stats['m3p']}/{stats['total']}),{size},{start}\n")
            
            for seq in cand['probe']:
                typ = f"Probe_{len(cand['probe'])}" if len(cand['probe']) > 1 else "Probe"
                stats = cand['p_stats']
                output.write(f"{locus_name},{role},{score},{typ},{seq},{len(seq)},{calc_tm(seq)},{calc_gc(seq)},{stats['p0']}% ({stats['m0']}/{stats['total']}),{stats['p1']}% ({stats['m1']}/{stats['total']}),{stats['p2']}% ({stats['m2']}/{stats['total']}),{stats['p3']}% ({stats['m3p']}/{stats['total']}),{size},{start}\n")
            
            for seq in cand['rev']:
                typ = f"Reverse_{len(cand['rev'])}" if len(cand['rev']) > 1 else "Reverse"
                stats = cand['r_stats']
                output.write(f"{locus_name},{role},{score},{typ},{seq},{len(seq)},{calc_tm(seq)},{calc_gc(seq)},{stats['p0']}% ({stats['m0']}/{stats['total']}),{stats['p1']}% ({stats['m1']}/{stats['total']}),{stats['p2']}% ({stats['m2']}/{stats['total']}),{stats['p3']}% ({stats['m3p']}/{stats['total']}),{size},{start}\n")
    
    return output.getvalue()

st.info("完整复刻版：所有计算逻辑、评分规则、硬底线、柔性惩罚均与原 HTML/JS 版本一致。")
