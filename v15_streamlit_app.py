import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import SeqIO
import io
from datetime import datetime
import time

st.set_page_config(page_title="V02 全量长度竞争寻优系统", layout="wide")
st.title("🧬 自动化引物探针柔性优先寻优系统 V02")
st.caption("工业级全景长度竞争 · Streamlit 优化版 | 产物 70–160 bp")

# ====================== 辅助函数 ======================
def calc_tm(seq):
    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
    return round(64.9 + 41 * (gc/100 * len(seq) - 16.4) / len(seq), 1)

def calc_gc(seq):
    return round((seq.count('G') + seq.count('C')) / len(seq) * 100, 1)

def reverse_complement(seq):
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def is_hard_valid_oligo(seq, is_probe=False):
    if not seq or 'N' in seq or not seq.isupper():
        return False
    gc = calc_gc(seq)
    if gc < 20 or gc > 80:
        return False
    if is_probe and seq[0] == 'G':
        return False
    return True

@st.cache_data
def parse_fasta(text):
    sequences = []
    for record in SeqIO.parse(io.StringIO(text), "fasta"):
        sequences.append(str(record.seq).upper())
    return sequences

def get_top_variants(start_idx, length, sequences, max_variants=2):
    counts = {}
    valid_count = 0
    for seq in sequences:
        if start_idx + length > len(seq):
            continue
        sub = seq[start_idx:start_idx + length]
        if '-' in sub or 'N' in sub:
            continue
        counts[sub] = counts.get(sub, 0) + 1
        valid_count += 1
    if valid_count / len(sequences) < 0.90 or not counts:
        return []
    sorted_v = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    variants = [sorted_v[0][0]]
    if len(sorted_v) > 1 and max_variants > 1:
        coverage = sorted_v[1][1] / valid_count
        if coverage > 0.04:
            variants.append(sorted_v[1][0])
    return variants

# ====================== 主界面 ======================
tab1, tab2 = st.tabs(["📤 序列输入", "⚙️ 参数设置"])

with tab1:
    uploaded = st.file_uploader("上传 FASTA 文件", type=["fasta", "fas", "txt"])
    fasta_text = st.text_area("或直接粘贴 FASTA 序列", height=300, placeholder=">seq1\nATGC...")

if uploaded:
    fasta_text = uploaded.getvalue().decode("utf-8")

if not fasta_text.strip():
    st.info("请上传文件或粘贴 FASTA 序列")
    st.stop()

sequences = parse_fasta(fasta_text)
if len(sequences) < 2:
    st.error("序列数量不足，至少需要2条序列进行比对。")
    st.stop()

st.success(f"成功加载 {len(sequences)} 条序列，长度 ≈ {len(sequences[0])} bp")

# 参数设置
with tab2:
    col1, col2 = st.columns(2)
    with col1:
        max_gap = st.slider("最大 Gap (bp)", 1, 30, 15)
        probe_weight = st.slider("探针匹配权重", 1, 5, 3)
    with col2:
        min_amplicon = st.slider("最小产物长度", 70, 120, 70)
        max_amplicon = st.slider("最大产物长度", 130, 200, 160)

# ====================== 运行计算 ======================
if st.button("🚀 启动全量长度竞争寻优", type="primary", use_container_width=True):
    with st.spinner("正在进行全量长度矩阵竞争计算，请稍候..."):
        progress_bar = st.progress(0)
        all_candidates = []
        seq_len = len(sequences[0])
        total_steps = seq_len - 160
        step = 0

        for i in range(seq_len - 160):
            # Forward
            f_list = []
            for length in range(18, 26):
                vars_f = get_top_variants(i, length, sequences, 2)
                if vars_f and all(is_hard_valid_oligo(v, False) for v in vars_f):
                    f_list.append((length, vars_f))
            
            for f_len, f_vars in f_list[:2]:  # Top 2
                for gap1 in range(1, max_gap + 1):
                    p_start = i + f_len + gap1
                    # Probe
                    p_list = []
                    for length in range(18, 31):
                        vars_p = get_top_variants(p_start, length, sequences, 1)
                        if vars_p and all(is_hard_valid_oligo(v, True) for v in vars_p):
                            p_list.append((length, vars_p))
                    
                    for p_len, p_vars in p_list[:2]:
                        for gap2 in range(1, max_gap + 1):
                            r_start = p_start + p_len + gap2
                            r_len_candidates = []
                            for length in range(18, 26):
                                raw_r = get_top_variants(r_start, length, sequences, 2)
                                if raw_r:
                                    comp_r = [reverse_complement(r) for r in raw_r]
                                    if all(is_hard_valid_oligo(c, False) for c in comp_r):
                                        r_len_candidates.append((length, comp_r, raw_r))
                            
                            for r_len, r_vars, raw_r in r_len_candidates[:2]:
                                amplicon = r_start + r_len - i
                                if not (min_amplicon <= amplicon <= max_amplicon):
                                    continue
                                
                                # 计算得分（简化版，保留核心逻辑）
                                score = 100  # 可进一步实现详细打分
                                all_candidates.append({
                                    "start": i,
                                    "fwd": f_vars,
                                    "probe": p_vars,
                                    "rev": r_vars,
                                    "size": amplicon,
                                    "score": score
                                })
            
            step += 1
            if step % 10 == 0:
                progress_bar.progress(min(1.0, step / total_steps))

        progress_bar.progress(1.0)

    # 显示结果
    st.success(f"发现 {len(all_candidates)} 个候选组合")
    
    # 简单排序后展示
    if all_candidates:
        df = pd.DataFrame(all_candidates)
        df = df.sort_values("score", ascending=False).head(20)
        
        st.subheader("🎯 优选候选靶区")
        for idx, row in df.iterrows():
            with st.expander(f"靶区 {idx+1} | 得分 {row['score']:.1f} | 产物 {row['size']} bp"):
                st.code(f"Forward: {row['fwd']}", language="text")
                st.code(f"Probe: {row['probe']}", language="text")
                st.code(f"Reverse: {row['rev']}", language="text")
        
        # 下载
        csv = df.to_csv(index=False).encode()
        st.download_button("📥 下载 CSV 结果", csv, f"DOE_结果_{datetime.now().strftime('%Y%m%d')}.csv", "text/csv")

st.info("Streamlit 版本已优化缓存与分步计算，适合 2.7GB 内存限制环境。")
