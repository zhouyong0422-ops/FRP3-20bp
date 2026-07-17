import streamlit as st
import math
import re
from collections import Counter
import plotly.graph_objects as go
import datetime

# ==========================================
# 页面与全局 CSS 样式设置 (保留原版 UI 风格)
# ==========================================
st.set_page_config(page_title="自动化引物探针设计比对系统", layout="wide")

st.markdown("""
<style>
    /* 核心报告卡片样式复刻 */
    .report-header { text-align: center; border-bottom: 2px dashed #ccc; padding-bottom: 20px; margin-bottom: 20px; }
    .report-header h2 { margin: 0; color: #2c3e50; font-size: 24px;}
    .report-header p { margin: 5px 0 0 0; color: #7f8c8d; font-size: 14px;}
    .locus-group { margin-bottom: 40px; background: #fdfefe; border: 1px solid #d1d5db; border-radius: 8px; overflow: hidden; box-shadow: 0 4px 6px rgba(0,0,0,0.02);}
    .locus-title { background: #edf2f7; color: #2c3e50; margin: 0; padding: 15px 20px; font-size: 18px; border-bottom: 1px solid #d1d5db; display: flex; justify-content: space-between; align-items: center;}
    .candidate-card { padding: 20px; border-bottom: 1px solid #ecf0f1; position: relative; }
    .candidate-card:last-child { border-bottom: none; }
    .cand-primary { border-left: 6px solid #27ae60; background: #f9fbfd;} 
    .cand-variant { border-left: 6px solid #bdc3c7; background: #fff;}    
    .candidate-card h4 { margin-top: 0; color: #2c3e50; font-size: 16px; border-bottom: 1px dashed #e1e8ed; padding-bottom: 8px; display: flex; justify-content: space-between; align-items: center;}
    .score-badge { background: #f39c12; color: white; padding: 3px 10px; border-radius: 12px; font-size: 13px; font-weight: bold;}
    .mix-badge { background: #9b59b6; color: white; padding: 2px 8px; border-radius: 8px; font-size: 11px; margin-left: 8px; vertical-align: middle;}
    .role-badge { font-size: 12px; padding: 2px 8px; border-radius: 4px; margin-right: 10px; color: white;}
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
    .details-content { font-size: 12px; color: #555; margin-top: 8px; line-height: 1.8; border-top: 1px dashed #eee; padding-top: 6px;}
    .score-item { display: flex; justify-content: space-between; }
    .score-plus { color: #27ae60; font-weight: bold;}
    .score-minus { color: #c0392b; font-weight: bold;}
</style>
""", unsafe_allow_html=True)

st.markdown("""
<h1 style="display: flex; justify-content: space-between; align-items: flex-end; border-bottom: 3px solid #1abc9c; padding-bottom: 10px;">
    <span style="color: #2c3e50;">🧬 自动化引物探针设计比对系统</span>
    <span style="font-size: 14px; background: #8e44ad; color: white; padding: 5px 12px; border-radius: 15px; font-weight: bold;">全局安全兼并 V26</span>
</h1>
""", unsafe_allow_html=True)

# ==========================================
# 核心算法辅助函数
# ==========================================
def calcTm(seq):
    g = len(re.findall(r'[gG]', seq))
    c = len(re.findall(r'[cC]', seq))
    s = len(re.findall(r'[sS]', seq))
    half = len(re.findall(r'[rRyYmMkK]', seq))
    gcCount = g + c + s + half * 0.5
    return round(64.9 + 41 * (gcCount - 16.4) / len(seq), 1)

def calcGC(seq):
    g = len(re.findall(r'[gG]', seq))
    c = len(re.findall(r'[cC]', seq))
    s = len(re.findall(r'[sS]', seq))
    half = len(re.findall(r'[rRyYmMkK]', seq))
    return round(((g + c + s + half * 0.5) / len(seq)) * 100, 1)

def reverseComplement(seq):
    mapping = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
        'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'
    }
    return ''.join([mapping.get(b.upper(), b) for b in reversed(seq)])

def hasSecondaryStructureRisk(seq):
    if len(seq) < 12: return False
    for i in range(4, 6):
        head = seq[:i]
        tail_comp = reverseComplement(seq[-i:])
        if head.upper() == tail_comp.upper(): return True
    return False

def has3PrimeDimerRisk(seq1, seq2):
    return seq1[-4:].upper() == reverseComplement(seq2[-4:]).upper()

def isHardValidOligo(seq, isProbe=False):
    if 'N' in seq.upper(): return False
    if re.search(r'[^ATGCRYSWKM]', seq, re.I): return False
    degCount = len(re.findall(r'[RYSWKM]', seq, re.I))
    if degCount > 1: return False

    if re.search(r'[RYSWKM]', seq[-5:], re.I): return False
    if isProbe and re.search(r'[RYSWKM]', seq[:1], re.I): return False
    
    gc = float(calcGC(seq))
    if gc < 20 or gc > 80: return False
    if isProbe and seq.upper().startswith('G'): return False
    return True

def calcSoftPenalties(seq, isProbe):
    penalty = 0
    gc = float(calcGC(seq))
    if re.search(r'([ATGC])\1{3,}', seq, re.I): penalty += 10
    if hasSecondaryStructureRisk(seq): penalty += 15

    if isProbe:
        if gc < 30: penalty += (30 - gc) * 2
        if gc > 65: penalty += (gc - 65) * 2
        gCount = len(re.findall(r'[gG]', seq))
        cCount = len(re.findall(r'[cC]', seq))
        if gCount >= cCount: penalty += 8
    else:
        if gc < 40: penalty += (40 - gc) * 2
        if gc > 60: penalty += (gc - 60) * 2
        end5 = seq[-5:]
        endGC = len(re.findall(r'[gGcC]', end5))
        if endGC < 1 or endGC > 2: penalty += 6
        if seq.upper().endswith('T'): penalty += 12
        if re.search(r'GG$|CC$', seq, re.I): penalty += 8
    return penalty

def isBaseMatch(b1, b2):
    b1, b2 = b1.upper(), b2.upper()
    if b1 == b2: return True
    if b1 == 'R' and b2 in ['A', 'G']: return True
    if b1 == 'Y' and b2 in ['C', 'T']: return True
    if b1 == 'S' and b2 in ['G', 'C']: return True
    if b1 == 'W' and b2 in ['A', 'T']: return True
    if b1 == 'K' and b2 in ['G', 'T']: return True
    if b1 == 'M' and b2 in ['A', 'C']: return True
    return False

def calculateMixMismatch(variantsArray, startIndex, sequencesArray):
    totalSeq = len(sequencesArray)
    stats = {'m0': 0, 'm1': 0, 'm2': 0, 'm3p': 0, 'total': totalSeq}
    seqLen = len(variantsArray[0])

    for seq_full in sequencesArray:
        libSeqSnippet = seq_full[startIndex : startIndex + seqLen]
        if '-' in libSeqSnippet or 'N' in libSeqSnippet.upper():
            stats['m3p'] += 1
            continue
        
        bestMismatches = seqLen
        for targetSeq in variantsArray:
            mismatches = 0
            for j in range(seqLen):
                if not isBaseMatch(targetSeq[j], libSeqSnippet[j]):
                    mismatches += 1
            if mismatches < bestMismatches:
                bestMismatches = mismatches
                
        if bestMismatches == 0: stats['m0'] += 1
        elif bestMismatches == 1: stats['m1'] += 1
        elif bestMismatches == 2: stats['m2'] += 1
        else: stats['m3p'] += 1

    if stats['total'] == 0:
        return {'p0': '0.0', 'p1': '0.0', 'p2': '0.0', 'p3': '0.0', 'm0': 0, 'm1': 0, 'm2': 0, 'm3p': 0, 'total': 0}
    
    return {
        'p0': round((stats['m0'] / stats['total']) * 100, 1),
        'p1': round((stats['m1'] / stats['total']) * 100, 1),
        'p2': round((stats['m2'] / stats['total']) * 100, 1),
        'p3': round((stats['m3p'] / stats['total']) * 100, 1),
        'm0': stats['m0'], 'm1': stats['m1'], 'm2': stats['m2'], 'm3p': stats['m3p'], 'total': stats['total']
    }

def getConservedProbeSequence(startIndex, length, sequencesArray, isReverse=False):
    validSeqs = []
    for seq_full in sequencesArray:
        snippet = seq_full[startIndex : startIndex + length]
        if '-' not in snippet and 'N' not in snippet.upper():
            validSeqs.append(snippet)
            
    if not validSeqs or (len(validSeqs) / len(sequencesArray)) < 0.80:
        return []

    counts = Counter(validSeqs)
    sorted_counts = counts.most_common()
    topSeq = sorted_counts[0][0]
    topCov = sorted_counts[0][1] / len(validSeqs)

    if topCov >= 0.98: return [topSeq]

    iupacMap = {
        'AG': 'R', 'GA': 'R', 'CT': 'Y', 'TC': 'Y',
        'GC': 'S', 'CG': 'S', 'AT': 'W', 'TA': 'W',
        'GT': 'K', 'TG': 'K', 'AC': 'M', 'CA': 'M'
    }

    bestDegCol = -1
    maxBoost = 0
    startCol = 5 if isReverse else 1
    endCol = (length - 1) if isReverse else (length - 5)

    for c in range(startCol, endCol):
        baseCounts = Counter([s[c] for s in validSeqs]).most_common()
        if len(baseCounts) >= 2:
            top1Cov = baseCounts[0][1] / len(validSeqs)
            top2Cov = (baseCounts[0][1] + baseCounts[1][1]) / len(validSeqs)
            boost = top2Cov - top1Cov
            if boost > maxBoost and boost >= 0.04:
                maxBoost = boost
                bestDegCol = c

    if bestDegCol != -1 and maxBoost > 0:
        baseCounts = Counter([s[bestDegCol] for s in validSeqs]).most_common()
        b1, b2 = baseCounts[0][0], baseCounts[1][0]
        degChar = iupacMap.get(b1 + b2)
        if degChar:
            degSeqChars = list(topSeq)
            degSeqChars[bestDegCol] = degChar
            return [''.join(degSeqChars)]

    return [topSeq]

def getTopVariants(startIndex, length, sequencesArray, maxVariants=2, isProbe=False, isReverse=False):
    if isProbe:
        return getConservedProbeSequence(startIndex, length, sequencesArray, isReverse)
        
    counts = Counter()
    totalValid = 0
    for seq_full in sequencesArray:
        snippet = seq_full[startIndex : startIndex + length]
        if '-' in snippet or 'N' in snippet.upper(): continue
        counts[snippet] += 1
        totalValid += 1

    if totalValid / len(sequencesArray) < 0.90 or not counts: return []

    sorted_counts = counts.most_common()
    variants = [sorted_counts[0][0]]
    coverage = sorted_counts[0][1] / totalValid

    if coverage < 0.97 and len(sorted_counts) > 1 and maxVariants > 1:
        if sorted_counts[1][1] / totalValid > 0.04:
            seq1 = sorted_counts[0][0]
            seq2 = sorted_counts[1][0]
            
            diffCount = sum(1 for a, b in zip(seq1, seq2) if a != b)
            diffIdx = next((i for i, (a, b) in enumerate(zip(seq1, seq2)) if a != b), -1)

            if diffCount == 1:
                is3Prime5Risk = (diffIdx < 5) if isReverse else (diffIdx >= length - 5)
                if not is3Prime5Risk:
                    iupacMap = {
                        'AG': 'R', 'GA': 'R', 'CT': 'Y', 'TC': 'Y',
                        'GC': 'S', 'CG': 'S', 'AT': 'W', 'TA': 'W',
                        'GT': 'K', 'TG': 'K', 'AC': 'M', 'CA': 'M'
                    }
                    b1, b2 = seq1[diffIdx], seq2[diffIdx]
                    degChar = iupacMap.get(b1 + b2)
                    if degChar:
                        merged = list(seq1)
                        merged[diffIdx] = degChar
                        return [''.join(merged)]
                    else:
                        variants.append(seq2)
                else:
                    variants.append(seq2)
            else:
                variants.append(seq2)
                
    return variants

# 包装后的获取函数 (带缓存)
def getAllValidVariants(startIdx, sequences, isProbe, memoDict):
    if startIdx in memoDict: return memoDict[startIdx]

    validList = []
    targetLengths = list(range(18, 37)) if isProbe else list(range(18, 26))
    
    for length in targetLengths:
        if startIdx + length > len(sequences[0]): continue
        rawVariants = getTopVariants(startIdx, length, sequences, 1 if isProbe else 2, isProbe, False)
        if not rawVariants: continue

        finalVariants = []
        allPassed = True
        for v in rawVariants:
            if not isHardValidOligo(v, isProbe):
                allPassed = False
                break
            finalVariants.append(v)
            
        if allPassed:
            softPen = sum(calcSoftPenalties(v, isProbe) for v in finalVariants)
            avgTm = sum(float(calcTm(v)) for v in finalVariants) / len(finalVariants)
            tmBonus = 0
            if isProbe:
                if 62 <= avgTm <= 68: tmBonus = 15
            else:
                if 55 <= avgTm <= 60: tmBonus = 15
                
            stats = calculateMixMismatch(finalVariants, startIdx, sequences)
            p0 = float(stats['p0'])
            conservationScore = (p0 * 3 + (50 if p0 >= 99.0 else 0) - ((98.0 - p0) * 10 if p0 < 98.0 else 0)) if isProbe else p0

            validList.append({
                'length': length, 'variants': finalVariants,
                'score': tmBonus - softPen, 'stats': stats,
                'totalScore': conservationScore + tmBonus - softPen
            })
            
    validList.sort(key=lambda x: x['totalScore'], reverse=True)
    memoDict[startIdx] = validList
    return validList

def getReverseValidVariants(startIdx, sequences, isProbe, memoDict):
    if startIdx in memoDict: return memoDict[startIdx]

    validList = []
    targetLengths = list(range(18, 37)) if isProbe else list(range(18, 26))

    for length in targetLengths:
        if startIdx + length > len(sequences[0]): continue
        rawVariants = getTopVariants(startIdx, length, sequences, 1 if isProbe else 2, isProbe, True)
        if not rawVariants: continue

        finalVariants = []
        allPassed = True
        for rv in rawVariants:
            comp = reverseComplement(rv)
            if not isHardValidOligo(comp, isProbe):
                allPassed = False
                break
            finalVariants.append(comp)
            
        if allPassed:
            softPen = sum(calcSoftPenalties(v, isProbe) for v in finalVariants)
            avgTm = sum(float(calcTm(v)) for v in finalVariants) / len(finalVariants)
            tmBonus = 0
            if isProbe:
                if 62 <= avgTm <= 68: tmBonus = 15
            else:
                if 55 <= avgTm <= 60: tmBonus = 15

            stats = calculateMixMismatch(rawVariants, startIdx, sequences)
            p0 = float(stats['p0'])
            conservationScore = (p0 * 3 + (50 if p0 >= 99.0 else 0) - ((98.0 - p0) * 10 if p0 < 98.0 else 0)) if isProbe else p0

            validList.append({
                'length': length, 'variants': finalVariants, 'rawVariants': rawVariants,
                'score': tmBonus - softPen, 'stats': stats,
                'totalScore': conservationScore + tmBonus - softPen
            })
            
    validList.sort(key=lambda x: x['totalScore'], reverse=True)
    memoDict[startIdx] = validList
    return validList

# ==========================================
# Streamlit UI 构建与逻辑触发
# ==========================================

st.markdown("""<p style="color: #7f8c8d; font-style: italic;">导入多序列比对 FASTA 变异库，自动筛选符合 TaqMan qPCR 规范的引物探针组合，输出可视化报告并支持 CSV 导出。</p>""", unsafe_allow_html=True)

uploaded_file = st.file_uploader("📂 导入 FASTA 序列库 (优先使用)", type=["fasta", "fas", "txt", "aln"])
fasta_text = st.text_area("或在此处直接粘贴比对完成的 FASTA 序列...", height=150)

if st.button("⚙️ 启动靶区筛选：高精度加权运算输出引探 DOE 试验方案", use_container_width=True):
    content = ""
    if uploaded_file is not None:
        content = uploaded_file.getvalue().decode("utf-8")
    elif fasta_text.strip():
        content = fasta_text.strip()
    
    if not content:
        st.warning("请导入或粘贴包含变异库的 FASTA 文件！")
    else:
        # 解析 FASTA
        lines = content.split('\n')
        sequences = []
        currentSeq = []
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if currentSeq:
                    sequences.append(''.join(currentSeq).upper())
                    currentSeq = []
            else:
                currentSeq.append(line)
        if currentSeq:
            sequences.append(''.join(currentSeq).upper())
            
        if len(sequences) < 2:
            st.error("文件格式有误或序列数不足2条。")
        else:
            seqLen = len(sequences[0])
            
            # 计算香农熵
            progress_bar = st.progress(0)
            status_text = st.empty()
            status_text.markdown(f"**⏳ 正在计算变异强度扫描 (0 / {seqLen} bp)...**")
            
            entropies = []
            for i in range(seqLen):
                column_chars = [seq[i] for seq in sequences if seq[i] != '-']
                total = len(column_chars)
                entropy = 0
                if total > 0:
                    counts = Counter(column_chars)
                    for count in counts.values():
                        p = count / total
                        entropy -= p * math.log2(p)
                entropies.append(entropy)
                
                if i % max(1, seqLen // 10) == 0:
                    progress_bar.progress(int(5 + (i / seqLen) * 10))
                    
            status_text.markdown("**✅ 变异分布扫描完成，准备执行靶区全景竞争...**")
            progress_bar.progress(15)

            # 显示图表
            st.markdown("<h3 style='margin-top: 30px;'>📊 靶标序列群变异强度扫描 (香农熵)</h3>", unsafe_allow_html=True)
            colors = ['#e74c3c' if v > 0.05 else '#3498db' for v in entropies]
            fig = go.Figure(data=[go.Bar(x=list(range(1, len(entropies)+1)), y=entropies, marker_color=colors)])
            fig.update_layout(margin=dict(l=0, r=0, t=10, b=0), height=250, plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
            st.plotly_chart(fig, use_container_width=True)

            # 核心变量
            allCandidates = []
            memoFwdPrimers, memoFwdProbes, memoRevProbes, memoRevPrimers = {}, {}, {}, {}
            totalSteps = max(1, seqLen - 160)
            
            # 主寻优循环
            for i in range(totalSteps):
                if i % 10 == 0 or i == totalSteps - 1:
                    pct = 15 + (i / totalSteps) * 80
                    progress_bar.progress(int(pct))
                    status_text.markdown(f"**⏳ 全景矩阵加权寻优中... 当前扫描参考坐标: {i} / {totalSteps} bp**")

                fObjList = getAllValidVariants(i, sequences, False, memoFwdPrimers)
                if not fObjList: continue

                for fObj in fObjList[:2]:
                    fVariants, fLen = fObj['variants'], fObj['length']
                    
                    minRStart = max(i + 70 - 25, i + fLen + 1)
                    maxRStart = min(i + 160 - 18, seqLen - 1)
                    
                    for rStart in range(minRStart, maxRStart + 1):
                        rObjList = getReverseValidVariants(rStart, sequences, False, memoRevPrimers)
                        if not rObjList: continue
                        
                        for rObj in rObjList[:2]:
                            rVariants, rLen = rObj['variants'], rObj['length']
                            ampliconSize = rStart + rLen - i
                            if ampliconSize < 70 or ampliconSize > 160: continue
                            
                            primerMaxTm = max(
                                max([float(calcTm(v)) for v in fVariants]),
                                max([float(calcTm(v)) for v in rVariants])
                            )

                            validProbes = []
                            # 寻找正向探针
                            for pStart in range(i + fLen + 1, i + fLen + 30):
                                if pStart >= rStart - 18: break
                                pObjList = getAllValidVariants(pStart, sequences, True, memoFwdProbes)
                                for pObj in pObjList:
                                    if pStart + pObj['length'] <= rStart - 1:
                                        minPTm = min([float(calcTm(v)) for v in pObj['variants']])
                                        if minPTm > primerMaxTm:
                                            p_copy = dict(pObj)
                                            p_copy.update({'start': pStart, 'isRev': False, 'gapToSameDir': pStart - (i + fLen)})
                                            validProbes.append(p_copy)
                            
                            # 寻找反向探针
                            for pEnd in range(rStart - 1, rStart - 31, -1):
                                if pEnd <= i + fLen + 18: break
                                for pLen in range(18, 37):
                                    pStart = pEnd - pLen
                                    if pStart >= i + fLen + 1 and (rStart - pEnd < 30):
                                        pObjList = getReverseValidVariants(pStart, sequences, True, memoRevProbes)
                                        for pObj in pObjList:
                                            if pObj['length'] == pLen:
                                                minPTm = min([float(calcTm(v)) for v in pObj['variants']])
                                                if minPTm > primerMaxTm:
                                                    p_copy = dict(pObj)
                                                    p_copy.update({'start': pStart, 'isRev': True, 'gapToSameDir': rStart - pEnd})
                                                    validProbes.append(p_copy)
                                                    
                            if not validProbes: continue
                            validProbes.sort(key=lambda x: x['totalScore'], reverse=True)
                            
                            for pObj in validProbes[:3]:
                                pVariants = pObj['variants']
                                primerMinTm = min(
                                    min([float(calcTm(v)) for v in fVariants]),
                                    min([float(calcTm(v)) for v in rVariants])
                                )
                                minPTm = min([float(calcTm(v)) for v in pVariants])
                                
                                crossSoftPenalty = 0
                                dimerRisk = False
                                for f in fVariants:
                                    for r in rVariants:
                                        if has3PrimeDimerRisk(f, r): dimerRisk = True
                                if dimerRisk: crossSoftPenalty += 20
                                
                                primerTmDiff = abs(primerMaxTm - primerMinTm)
                                if primerTmDiff > 2.0: crossSoftPenalty += (primerTmDiff - 2.0) * 5
                                if minPTm < primerMaxTm + 5.0: crossSoftPenalty += (primerMaxTm + 5.0 - minPTm) * 6
                                if pObj['gapToSameDir'] > 10: crossSoftPenalty += (pObj['gapToSameDir'] - 10) * 1.5
                                
                                mixF = -35 if len(fVariants) > 1 else 0
                                mixR = -35 if len(rVariants) > 1 else 0
                                
                                finalScore = fObj['totalScore'] + pObj['totalScore'] + rObj['totalScore'] + mixF + mixR - crossSoftPenalty
                                
                                baseP0Score = float(fObj['stats']['p0']) + (float(pObj['stats']['p0']) * 3) + float(rObj['stats']['p0'])
                                pP0 = float(pObj['stats']['p0'])
                                probeBonus = 50 if pP0 >= 99.0 else 0
                                probePenalty = (98.0 - pP0) * 10 if pP0 < 98.0 else 0
                                selfSoftPen = (fObj['score'] - (15 if fObj['length'] else 0)) + \
                                              (pObj['score'] - (15 if pObj['length'] else 0)) + \
                                              (rObj['score'] - (15 if rObj['length'] else 0))
                                displaySoftPen = -selfSoftPen + crossSoftPenalty
                                
                                allCandidates.append({
                                    'fwd': fVariants, 'rev': rVariants, 'probe': pVariants,
                                    'fStats': fObj['stats'], 'pStats': pObj['stats'], 'rStats': rObj['stats'],
                                    'size': ampliconSize, 'start': i, 'score': finalScore,
                                    'probeDir': "Reverse (-)" if pObj['isRev'] else "Forward (+)",
                                    'probeGap': pObj['gapToSameDir'],
                                    'details': {'base': baseP0Score, 'pBonus': probeBonus, 'pPenalty': -probePenalty, 'mixF': mixF, 'mixR': mixR, 'softPen': -displaySoftPen}
                                })
                
                # 【内存截断优化】：每当候选集庞大时执行瘦身，确保 Streamlit 内存安全
                if len(allCandidates) > 2000:
                    allCandidates.sort(key=lambda x: x['score'], reverse=True)
                    allCandidates = allCandidates[:1000]
            
            progress_bar.progress(96)
            status_text.markdown("**⏳ 正在整合排序并聚类黄金去重区...**")
            
            allCandidates.sort(key=lambda x: x['score'], reverse=True)
            globalLociGroups = []
            LOCUS_WINDOW = 50
            
            for cand in allCandidates:
                foundLocus = False
                for locus in globalLociGroups:
                    if abs(cand['start'] - locus['anchorStart']) <= LOCUS_WINDOW:
                        if len(locus['variants']) < 3: locus['variants'].append(cand)
                        foundLocus = True
                        break
                if not foundLocus:
                    globalLociGroups.append({'locusId': len(globalLociGroups) + 1, 'anchorStart': cand['start'], 'variants': [cand]})

            progress_bar.progress(100)
            status_text.markdown("**✅ 寻优完成！最终 DOE 清单已就绪。**")
            
            # HTML 报告渲染
            reportHTML = f"""
            <div class="report-box" style="display:block;">
                <div class="report-header">
                    <h2>多重 TaqMan 引探组合设计分析报告</h2>
                    <p>生成时间：{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
                    <p style="font-size: 13px; color: #95a5a6;">核心红线：探针 Tm 严格大于引物；探针 Tm 不足自动延展至 36bp；<b>上下游引物仅 1 个碱基差异且严禁落在 3' 端 5 个碱基内时，自动折叠为 IUPAC 兼并碱基；若落在 3' 端 5 bp 内或差异 ≥2 bp，严禁兼并，保留双条混合套</b>。</p>
                </div>
            """
            
            if not globalLociGroups:
                reportHTML += `<div class="candidate-card" style="border-left-color: #e74c3c;"><h4 style="color:#c0392b;">体系设计失败</h4><p>在该序列库中未能找到满足各项约束指标的有效靶区。</p></div>`
            else:
                def generateSeqRows(title, variantsArray):
                    html = ''
                    for idx, v in enumerate(variantsArray):
                        label = f"{title} {idx + 1}" if len(variantsArray) > 1 else title
                        html += f"""
                        <div class="seq-row">
                            <span class="seq-type">{label}:</span>
                            <span class="seq-string">5'- {v} -3'</span>
                            <span class="seq-stats">Len: {len(v)}bp | Tm: {calcTm(v)}°C | GC: {calcGC(v)}%</span>
                        </div>"""
                    return html
                    
                def getStatHTML(stats, isMix):
                    mix_badge = '<span class="mix-badge">混合套数扣分(-35)</span>' if isMix else ''
                    return f"""
                    <div class="mismatch-stats">
                        <span class="stat-badge bg-0">完全匹配(0): {stats['p0']}%, {stats['m0']}/{stats['total']}</span>
                        <span class="stat-badge bg-1">错配1碱基(1): {stats['p1']}%, {stats['m1']}/{stats['total']}</span>
                        <span class="stat-badge bg-2">错配2碱基(2): {stats['p2']}%, {stats['m2']}/{stats['total']}</span>
                        <span class="stat-badge bg-3">错配3碱基(≥3): {stats['p3']}%, {stats['m3p']}/{stats['total']}</span>
                        {mix_badge}
                    </div>"""

                for locus in globalLociGroups:
                    reportHTML += f"""
                    <div class="locus-group">
                        <h3 class="locus-title">
                            <span>🎯 独立黄金靶区 {locus['locusId']}</span>
                            <span style="font-size: 14px; color: #7f8c8d; font-weight: normal;">(参考起始坐标: {locus['anchorStart']})</span>
                        </h3>
                    """
                    for vIndex, cand in enumerate(locus['variants']):
                        isPrimary = (vIndex == 0)
                        cardClass = 'cand-primary' if isPrimary else 'cand-variant'
                        roleBadge = `<span class="role-badge role-main">主力优选</span>` if isPrimary else f'<span class="role-badge role-sub">微调备选 {vIndex}</span>'
                        
                        reportHTML += f"""
                        <div class="candidate-card {cardClass}">
                            <h4>
                                <span>{roleBadge} 综合得分: <span class="score-badge">{cand['score']:.1f} 分</span></span>
                                <span style="font-size: 13px; font-weight: normal; color: #95a5a6;">精确定位: {cand['start']} | 产物长度: {cand['size']} bp | 探针间距: {cand['probeGap']} bp</span>
                            </h4>
                            
                            <div class="seq-block">
                                {generateSeqRows('Forward', cand['fwd'])}
                                {getStatHTML(cand['fStats'], len(cand['fwd']) > 1)}
                            </div>
                            <div class="seq-block probe-block">
                                {generateSeqRows(f"Probe [{cand['probeDir']}]", cand['probe'])}
                                {getStatHTML(cand['pStats'], len(cand['probe']) > 1)}
                            </div>
                            <div class="seq-block">
                                {generateSeqRows('Reverse', cand['rev'])}
                                {getStatHTML(cand['rStats'], len(cand['rev']) > 1)}
                            </div>

                            <details>
                                <summary>🔍 展开查看当前候选体系全维度评估分</summary>
                                <div class="details-content">
                                    <div class="score-item"><span>基础匹配分 (探针 3 倍权重)</span><span class="score-plus">+{cand['details']['base']:.1f}</span></div>
                                    <div class="score-item"><span>探针匹配度奖励 (完美率 ≥99%)</span><span class="score-plus">+{cand['details']['pBonus']:.1f}</span></div>
                                    <div class="score-item"><span>探针错配惩罚 (低于98%十倍扣除)</span><span class="{'score-minus' if cand['details']['pPenalty'] != 0 else ''}">{cand['details']['pPenalty']:.1f}</span></div>
                                    <div class="score-item"><span>F/R 混合套数扣分 (-35分/次，极力优选单套)</span><span class="{'score-minus' if (cand['details']['mixF'] + cand['details']['mixR']) != 0 else ''}">{cand['details']['mixF'] + cand['details']['mixR']}</span></div>
                                    <div class="score-item" style="border-top: 1px solid #ddd; padding-top: 4px; font-weight: bold;"><span>偏离规则体系总扣分 (GC/温差/同向间距过长等)</span><span class="{'score-minus' if cand['details']['softPen'] != 0 else ''}">{cand['details']['softPen']:.1f}</span></div>
                                </div>
                            </details>
                        </div>"""
                    reportHTML += "</div>"
            reportHTML += "</div>"
            
            st.markdown(reportHTML, unsafe_allow_html=True)
            
            # CSV 导出生成
            if globalLociGroups:
                csv_content = "\ufeff靶区归属,变体角色,综合得分,寡核苷酸类型,序列 (5'->3'),长度 (bp),Tm (°C),GC (%),完美匹配(0),错配1碱基(1),错配2碱基(2),错配3碱基(≥3),预期产物长度 (bp),同向引物间距 (bp),精确起始坐标\n"
                for locus in globalLociGroups:
                    locusName = f"靶区_{locus['locusId']}"
                    for vIndex, cand in enumerate(locus['variants']):
                        role = "主力优选" if vIndex == 0 else f"备选_{vIndex}"
                        score = f"{cand['score']:.1f}"
                        
                        for i, seq in enumerate(cand['fwd']):
                            type_label = f"Forward_{i+1}" if len(cand['fwd']) > 1 else "Forward"
                            stats = cand['fStats']
                            csv_content += f"{locusName},{role},{score},{type_label},{seq},{len(seq)},{calcTm(seq)},{calcGC(seq)},{stats['p0']}% ({stats['m0']}/{stats['total']}),{stats['p1']}% ({stats['m1']}/{stats['total']}),{stats['p2']}% ({stats['m2']}/{stats['total']}),{stats['p3']}% ({stats['m3p']}/{stats['total']}),{cand['size']},-,{cand['start']}\n"
                        
                        for i, seq in enumerate(cand['probe']):
                            type_label = f"Probe({cand['probeDir']})_{i+1}" if len(cand['probe']) > 1 else f"Probe({cand['probeDir']})"
                            stats = cand['pStats']
                            csv_content += f"{locusName},{role},{score},{type_label},{seq},{len(seq)},{calcTm(seq)},{calcGC(seq)},{stats['p0']}% ({stats['m0']}/{stats['total']}),{stats['p1']}% ({stats['m1']}/{stats['total']}),{stats['p2']}% ({stats['m2']}/{stats['total']}),{stats['p3']}% ({stats['m3p']}/{stats['total']}),{cand['size']},{cand['probeGap']},{cand['start']}\n"
                        
                        for i, seq in enumerate(cand['rev']):
                            type_label = f"Reverse_{i+1}" if len(cand['rev']) > 1 else "Reverse"
                            stats = cand['rStats']
                            csv_content += f"{locusName},{role},{score},{type_label},{seq},{len(seq)},{calcTm(seq)},{calcGC(seq)},{stats['p0']}% ({stats['m0']}/{stats['total']}),{stats['p1']}% ({stats['m1']}/{stats['total']}),{stats['p2']}% ({stats['m2']}/{stats['total']}),{stats['p3']}% ({stats['m3p']}/{stats['total']}),{cand['size']},-,{cand['start']}\n"
                
                dateStr = datetime.datetime.now().strftime("%Y-%m-%d")
                st.download_button("📥 下载完整 DOE 清单 (Excel CSV)", data=csv_content.encode('utf-8-sig'), file_name=f"多重 TaqMan 引探组合设计分析报告{dateStr}.csv", mime="text/csv", use_container_width=True)
