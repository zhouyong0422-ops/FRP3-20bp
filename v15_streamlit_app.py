import streamlit as st

# 完整原HTML源码，替换境外chart.js为国内jsdelivr镜像，其余内容一字未改
full_html = """
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <title>V02 原料齐套 - 工业级柔性优先寻优系统 V15</title>
    <!-- 替换境外无法访问CDN为国内可访问镜像，图表功能完全不变 -->
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.8/dist/chart.umd.min.js"></script>
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background: #eef2f5; margin: 0; padding: 20px; color: #333; }
        .dashboard { max-width: 1200px; margin: auto; background: white; padding: 30px; border-radius: 12px; box-shadow: 0 6px 20px rgba(0,0,0,0.08); }
        h1 { color: #2c3e50; border-bottom: 3px solid #1abc9c; padding-bottom: 10px; margin-top: 0; display: flex; justify-content: space-between; align-items: flex-end;}
        .badge { font-size: 14px; background: #2c3e50; color: #f1c40f; padding: 5px 12px; border-radius: 15px; font-weight: bold; border: 1px solid #f1c40f; }
        
        .upload-section { display: flex; align-items: center; gap: 15px; margin-bottom: 15px; background: #f8f9fa; padding: 15px; border-radius: 8px; border: 2px dashed #bdc3c7; }
        .upload-btn { background: #27ae60; color: white; padding: 10px 20px; border-radius: 6px; cursor: pointer; font-weight: bold; transition: background 0.3s; }
        .upload-btn:hover { background: #219150; }
        #fileName { color: #7f8c8d; font-style: italic; }
        
        textarea { width: 100%; height: 100px; padding: 10px; border: 1px solid #ccc; border-radius: 6px; font-family: monospace; resize: vertical; box-sizing: border-box; }
        .run-btn { background: #2980b9; color: white; border: none; padding: 15px 24px; border-radius: 6px; cursor: pointer; font-size: 18px; font-weight: bold; transition: background 0.3s; margin-top: 15px; width: 100%; box-shadow: 0 4px 6px rgba(41, 128, 185, 0.3); }
        .run-btn:hover { background: #1f6391; }
        
        .chart-container { position: relative; height: 250px; width: 100%; margin-top: 20px; }
        
        .report-box { margin-top: 30px; background: #fff; padding: 30px; border: 1px solid #e0e0e0; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.05); display: none; position: relative;}
        .report-header { text-align: center; border-bottom: 2px dashed #ccc; padding-bottom: 20px; margin-bottom: 20px; position: relative; }
        .report-header h2 { margin: 0; color: #2c3e50; font-size: 24px;}
        .report-header p { margin: 5px 0 0 0; color: #7f8c8d; }
        
        .download-btn { position: absolute; right: 0; top: 0; background: #f39c12; color: white; border: none; padding: 10px 20px; border-radius: 6px; cursor: pointer; font-weight: bold; transition: 0.3s; box-shadow: 0 2px 4px rgba(243, 156, 18, 0.3);}
        .download-btn:hover { background: #d68910; }
        
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
        .seq-type { font-weight: bold; width: 85px; display: inline-block; color: #34495e;}
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
        
        #loadingMsg { display: none; margin-top: 15px; color: #2980b9; font-weight: bold; font-size: 16px; text-align: center; }
    </style>
</head>
<body>

<div class="dashboard">
    <h1>
        <span>🧬 自动化引物探针柔性优先寻优系统</span>
        <span class="badge">纯净单套优选 V15</span>
    </h1>
    
    <div class="upload-section">
        <label for="fileInput" class="upload-btn">📂 导入 FASTA 序列库</label>
        <input type="file" id="fileInput" accept=".fasta,.fas,.txt,.aln" style="display: none;" onchange="handleFileUpload(event)">
        <span id="fileName">死守4大绝对底线；取消短产物奖励（70-150bp严格一视同仁）；重罚混合引物（极力保障单套纯净输出）！</span>
    </div>

    <textarea id="fastaInput" placeholder="或在此处直接粘贴比对完成的 FASTA 序列..."></textarea>
    
    <button class="run-btn" onclick="runPipeline()">⚙️ 启动寻优：基于纯净单套准则出具靶区 DOE</button>
    <div id="loadingMsg">⏳ 正在进行全景矩阵与单套优先加权计算，请耐心等待...</div>

    <div id="analysisArea" style="display: none;">
        <h3 style="margin-top: 30px;">📊 靶标序列群变异强度扫描 (香农熵)</h3>
        <div class="chart-container"><canvas id="entropyChart"></canvas></div>
    </div>

    <div id="reportArea" class="report-box">
        <div class="report-header">
            <h2>V02 原料齐套：工业级多重引探 DOE 报告</h2>
            <p id="reportDate">生成时间：</p>
            <p style="font-size: 13px; color: #95a5a6;">规范内核：绝对死守 70–150 bp 产物（无偏向）、探针 5'端排 G 及长度边界；严重扣分遏制 F/R 双套混合，保证生产极简。</p>
            <button onclick="downloadCSV()" class="download-btn">📥 导出 Excel (含靶区分组)</button>
        </div>
        <div id="candidatesOutput"></div>
    </div>
</div>

<script>
    let chartInstance = null;
    let globalLociGroups = []; 

    function handleFileUpload(event) {
        const file = event.target.files[0];
        if (!file) return;
        document.getElementById('fileName').textContent = `已加载序列库: ${file.name}`;
        const reader = new FileReader();
        reader.onload = function(e) { document.getElementById('fastaInput').value = e.target.result; };
        reader.readAsText(file);
    }

    function calcTm(seq) {
        let g = (seq.match(/G/g) || []).length;
        let c = (seq.match(/C/g) || []).length;
        return (64.9 + 41 * (g + c - 16.4) / seq.length).toFixed(1);
    }

    function calcGC(seq) {
        let g = (seq.match(/G/g) || []).length;
        let c = (seq.match(/C/g) || []).length;
        return (((g + c) / seq.length) * 100).toFixed(1);
    }

    function reverseComplement(seq) {
        const map = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C' };
        return seq.split('').reverse().map(b => map[b] || b).join('');
    }

    function hasSecondaryStructureRisk(seq) {
        let len = seq.length;
        if (len < 12) return false;
        for (let i = 4; i <= 5; i++) {
            let head = seq.slice(0, i);
            let tailComp = reverseComplement(seq.slice(-i));
            if (head === tailComp) return true;
        }
        return false;
    }

    function has3PrimeDimerRisk(seq1, seq2) {
        let end1 = seq1.slice(-4);
        let end2Comp = reverseComplement(seq2.slice(-4));
        return end1 === end2Comp;
    }

    // =========================================================
    // 1. 绝对硬底线质控（触犯一票否决）
    // =========================================================
    function isHardValidOligo(seq, isProbe = false) {
        if (seq.includes('N') || /[^ATGC]/.test(seq)) return false;
        
        let gc = parseFloat(calcGC(seq));
        if (gc < 20 || gc > 80) return false;

        if (isProbe && seq.charAt(0) === 'G') return false;

        return true;
    }

    // =========================================================
    // 2. 柔性优先准则（触犯则执行降维惩罚打分）
    // =========================================================
    function calcSoftPenalties(seq, isProbe) {
        let penalty = 0;
        let gc = parseFloat(calcGC(seq));
        
        if (/([ATGC])\1{3,}/.test(seq)) penalty += 10;
        if (hasSecondaryStructureRisk(seq)) penalty += 15;

        if (isProbe) {
            if (gc < 30) penalty += (30 - gc) * 2;
            if (gc > 65) penalty += (gc - 65) * 2;
            
            let gCount = (seq.match(/G/g) || []).length;
            let cCount = (seq.match(/C/g) || []).length;
            if (gCount >= cCount) penalty += 8;
            
        } else {
            if (gc < 40) penalty += (40 - gc) * 2;
            if (gc > 60) penalty += (gc - 60) * 2;
            
            let end5 = seq.slice(-5);
            let endGC = (end5.match(/[GC]/g) || []).length;
            if (endGC < 1 || endGC > 2) penalty += 6;
            
            if (seq.slice(-1) === 'T') penalty += 12;
            if (/GG$|CC$/.test(seq)) penalty += 8;
        }
        return penalty;
    }

    function getTopVariants(startIndex, length, sequencesArray, maxVariants = 2) {
        let counts = {};
        let totalValid = 0;
        let totalSeq = sequencesArray.length;

        for(let i = 0; i < totalSeq; i++) {
            let seq = sequencesArray[i].substring(startIndex, startIndex + length);
            if(seq.includes('-') || seq.includes('N')) continue;
            counts[seq] = (counts[seq] || 0) + 1;
            totalValid++;
        }

        if (totalValid / totalSeq < 0.90) return [];
        if (Object.keys(counts).length === 0) return [];

        let sorted = Object.entries(counts).sort((a,b) => b[1] - a[1]);
        let variants = [sorted[0][0]];
        let coverage = sorted[0][1] / totalValid;

        // 仅当第一大序列覆盖率偏低，且第二大序列占比超 4% 时，才允许提取变体
        if (coverage < 0.97 && sorted.length > 1 && maxVariants > 1) {
            let secondCov = sorted[1][1] / totalValid;
            if (secondCov > 0.04) variants.push(sorted[1][0]);
        }
        return variants;
    }

    function getAllValidVariants(startIdx, sequences, isProbe) {
        let validList = [];
        const targetLengths = isProbe 
            ? [20, 22, 24, 18, 19, 21, 23, 25, 26, 27, 28, 29, 30]
            : [20, 21, 22, 19, 18, 23, 24, 25];
        
        for (let len of targetLengths) {
            if (startIdx + len > sequences[0].length) continue;
            let rawVariants = getTopVariants(startIdx, len, sequences, isProbe ? 1 : 2);
            if (rawVariants.length === 0) continue;

            let finalVariants = [];
            let allPassed = true;
            
            for (let v of rawVariants) {
                if (!isHardValidOligo(v, isProbe)) {
                    allPassed = false;
                    break;
                }
                finalVariants.push(v);
            }
            
            if (allPassed) {
                validList.push({ length: len, variants: finalVariants });
                break; 
            }
        }
        return validList;
    }

    function calculateMixMismatch(variantsArray, startIndex, sequencesArray) {
        let totalSeq = sequencesArray.length;
        let stats = { m0: 0, m1: 0, m2: 0, m3p: 0, total: totalSeq };
        const seqLen = variantsArray[0].length;

        for (let i = 0; i < totalSeq; i++) {
            let libSeqSnippet = sequencesArray[i].substring(startIndex, startIndex + seqLen);
            
            if (libSeqSnippet.includes('-') || libSeqSnippet.includes('N')) {
                stats.m3p++;
                continue; 
            }
            
            let bestMismatches = seqLen;
            for (let v = 0; v < variantsArray.length; v++) {
                let targetSeq = variantsArray[v];
                let mismatches = 0;
                for (let j = 0; j < seqLen; j++) {
                    if (targetSeq[j] !== libSeqSnippet[j]) mismatches++;
                }
                if (mismatches < bestMismatches) bestMismatches = mismatches;
            }
            
            if (bestMismatches === 0) stats.m0++;
            else if (bestMismatches === 1) stats.m1++;
            else if (bestMismatches === 2) stats.m2++;
            else stats.m3p++;
        }

        if(stats.total === 0) return { p0: '0.0', p1: '0.0', p2: '0.0', p3: '0.0', m0:0, m1:0, m2:0, m3p:0, total:0 };

        return {
            p0: parseFloat(((stats.m0 / stats.total) * 100).toFixed(1)),
            p1: ((stats.m1 / stats.total) * 100).toFixed(1),
            p2: ((stats.m2 / stats.total) * 100).toFixed(1),
            p3: ((stats.m3p / stats.total) * 100).toFixed(1),
            m0: stats.m0, m1: stats.m1, m2: stats.m2, m3p: stats.m3p, total: stats.total
        };
    }

    function runPipeline() {
        const text = document.getElementById('fastaInput').value.trim();
        if (!text) { alert("请导入包含变异库的文件！"); return; }
        
        document.getElementById('loadingMsg').style.display = 'block';
        document.getElementById('reportArea').style.display = 'none';
        
        setTimeout(() => { executeLogic(text); }, 100);
    }

    function executeLogic(text) {
        const lines = text.split('\n');
        let sequences = []; let currentSeq = "";
        for (let line of lines) {
            if (line.startsWith('>')) {
                if (currentSeq) sequences.push(currentSeq.toUpperCase());
                currentSeq = "";
            } else { currentSeq += line.trim(); }
        }
        if (currentSeq) sequences.push(currentSeq.toUpperCase());
        
        if (sequences.length < 2) {
            alert("文件格式有误或序列数不足2条。");
            document.getElementById('loadingMsg').style.display = 'none';
            return;
        }

        const seqLen = sequences[0].length;
        let entropies = []; 
        
        for (let i = 0; i < seqLen; i++) {
            let column = {}; let total = 0;
            for (let j = 0; j < sequences.length; j++) {
                const base = sequences[j][i];
                if (base && base !== '-') { column[base] = (column[base] || 0) + 1; total++; }
            }
            let entropy = 0;
            for (const base in column) {
                const p = column[base] / total;
                entropy -= p * Math.log2(p);
            }
            entropies.push(entropy);
        }

        renderChart(entropies);
        document.getElementById('analysisArea').style.display = 'block';

        let allCandidates = [];
        const minGap = 1; const maxGap = 15; 

        for (let i = 0; i < seqLen - 150; i++) {
            let fObjList = getAllValidVariants(i, sequences, false);
            if (fObjList.length === 0) continue;
            let fObj = fObjList[0];
            let fVariants = fObj.variants;
            let fLen = fObj.length;

            for (let gap1 = minGap; gap1 <= maxGap; gap1++) {
                let pStart = i + fLen + gap1;
                let pObjList = getAllValidVariants(pStart, sequences, true);
                if (pObjList.length === 0) continue;
                let pObj = pObjList[0];
                let pVariants = pObj.variants;
                let pLen = pObj.length;

                for (let gap2 = minGap; gap2 <= maxGap; gap2++) {
                    let rStart = pStart + pLen + gap2;
                    const targetLengths = [20, 21, 22, 19, 18, 23, 24, 25];
                    let rFound = false;
                    let rVariantsRaw = [];
                    let rVariants = [];
                    let rLen = 0;
                    
                    for (let testLen of targetLengths) {
                        if (rStart + testLen > seqLen) continue;
                        let rawVars = getTopVariants(rStart, testLen, sequences, 2);
                        if (rawVars.length === 0) continue;
                        
                        let allPassed = true;
                        let processedRev = [];
                        for(let rv of rawVars) {
                            let comp = reverseComplement(rv);
                            if (!isHardValidOligo(comp, false)) { allPassed = false; break; }
                            processedRev.push(comp);
                        }
                        
                        if (allPassed) {
                            rFound = true;
                            rVariantsRaw = rawVars;
                            rVariants = processedRev;
                            rLen = testLen;
                            break;
                        }
                    }
                    
                    if (!rFound) continue;
                    
                    // 【绝对硬底线】：产物长度严格 70–150 bp（只要在此范围内即可，无偏向加权）
                    let ampliconSize = rStart + rLen - i;
                    if (ampliconSize < 70 || ampliconSize > 150) continue;

                    let minFTm = Math.min(...fVariants.map(v => parseFloat(calcTm(v))));
                    let maxFTm = Math.max(...fVariants.map(v => parseFloat(calcTm(v))));
                    let minRTm = Math.min(...rVariants.map(v => parseFloat(calcTm(v))));
                    let maxRTm = Math.max(...rVariants.map(v => parseFloat(calcTm(v))));
                    let minPTm = Math.min(...pVariants.map(v => parseFloat(calcTm(v))));
                    
                    let primerMaxTm = Math.max(maxFTm, maxRTm);
                    let primerMinTm = Math.min(minFTm, minRTm);
                    
                    // ==========================================
                    // 柔性评分核算体系 (Soft Scoring Engine)
                    // ==========================================
                    let softPenalty = 0;
                    
                    for (let f of fVariants) softPenalty += calcSoftPenalties(f, false);
                    for (let p of pVariants) softPenalty += calcSoftPenalties(p, true);
                    for (let r of rVariants) softPenalty += calcSoftPenalties(r, false);

                    let dimerRisk = false;
                    for (let f of fVariants) {
                        for (let r of rVariants) { if (has3PrimeDimerRisk(f, r)) dimerRisk = true; }
                    }
                    if (dimerRisk) softPenalty += 20;

                    let primerTmDiff = Math.abs(primerMaxTm - primerMinTm);
                    if (primerTmDiff > 2.0) softPenalty += (primerTmDiff - 2.0) * 5;

                    if (minPTm < primerMaxTm + 5.0) {
                        softPenalty += (primerMaxTm + 5.0 - minPTm) * 6;
                    }

                    if (gap1 > 5) softPenalty += (gap1 - 5) * 2;
                    if (gap2 > 5) softPenalty += (gap2 - 5) * 2;

                    let fStats = calculateMixMismatch(fVariants, i, sequences);
                    let pStats = calculateMixMismatch(pVariants, pStart, sequences);
                    let rStats = calculateMixMismatch(rVariantsRaw, rStart, sequences);

                    let fP0 = parseFloat(fStats.p0);
                    let pP0 = parseFloat(pStats.p0);
                    let rP0 = parseFloat(rStats.p0);

                    // 基础匹配分 (探针 3 倍权重)
                    let baseScore = fP0 + (pP0 * 3) + rP0;
                    let probeBonus = (pP0 >= 99.0) ? 50 : 0;
                    let probePenalty = (pP0 < 98.0) ? (98.0 - pP0) * 10 : 0;

                    // 【重大调整】：对使用两套（混合引物）进行重罚 (-35分)，极大保障单套纯净优选
                    let mixF = fVariants.length > 1 ? -35 : 0;
                    let mixR = rVariants.length > 1 ? -35 : 0;
                    
                    // 最终计算汇总分 (取消了短扩增子偏向奖励 shortAmpliconBonus)
                    let totalScore = baseScore + probeBonus - probePenalty + mixF + mixR - softPenalty;

                    allCandidates.push({ 
                        fwd: fVariants, rev: rVariants, probe: pVariants, 
                        fStats: fStats, pStats: pStats, rStats: rStats,
                        size: ampliconSize, start: i, score: totalScore,
                        details: { base: baseScore, pBonus: probeBonus, pPenalty: -probePenalty, mixF: mixF, mixR: mixR, softPen: -softPenalty }
                    });
                }
            }
        }

        allCandidates.sort((a, b) => b.score - a.score);
        globalLociGroups = [];
        const LOCUS_WINDOW = 50; 

        for (let cand of allCandidates) {
            let foundLocus = false;
            for (let locus of globalLociGroups) {
                if (Math.abs(cand.start - locus.anchorStart) <= LOCUS_WINDOW) {
                    if (locus.variants.length < 3) { locus.variants.push(cand); }
                    foundLocus = true; break;
                }
            }
            if (!foundLocus) {
                globalLociGroups.push({ locusId: globalLociGroups.length + 1, anchorStart: cand.start, variants: [cand] });
            }
        }

        const now = new Date();
        document.getElementById('reportDate').textContent = `出具时间：${now.toLocaleString()}`;
        document.getElementById('loadingMsg').style.display = 'none';
        
        let reportHTML = "";
        if (globalLociGroups.length === 0) {
            reportHTML = `<div class="candidate-card" style="border-left-color: #e74c3c;"><h4 style="color:#c0392b;">⚠️ 体系设计失败</h4><p>在绝对硬底线限制下（产物70-150bp、探针5'不为G、长度18-30nt），该序列库未能找到完整无缺失(Gap<10%)的组合区域。</p></div>`;
        } else {
            const generateSeqRows = (title, variantsArray) => {
                let html = '';
                variantsArray.forEach((v, idx) => {
                    let label = variantsArray.length > 1 ? `${title} ${idx + 1}` : title;
                    html += `
                    <div class="seq-row">
                        <span class="seq-type">${label}:</span>
                        <span class="seq-string">5'- ${v} -3'</span>
                        <span class="seq-stats">Len: ${v.length}bp | Tm: ${calcTm(v)}°C | GC: ${calcGC(v)}%</span>
                    </div>`;
                });
                return html;
            };

            const getStatHTML = (stats, isMix) => `
                <div class="mismatch-stats">
                    <span class="stat-badge bg-0">完全匹配(0): ${stats.p0}%, ${stats.m0}/${stats.total}</span>
                    <span class="stat-badge bg-1">错配1碱基(1): ${stats.p1}%, ${stats.m1}/${stats.total}</span>
                    <span class="stat-badge bg-2">错配2碱基(2): ${stats.p2}%, ${stats.m2}/${stats.total}</span>
                    <span class="stat-badge bg-3">错配3碱基(≥3): ${stats.p3}%, ${stats.m3p}/${stats.total}</span>
                    ${isMix ? '<span class="mix-badge">混合套数扣分(-35)</span>' : ''}
                </div>
            `;

            globalLociGroups.forEach((locus) => {
                reportHTML += `
                <div class="locus-group">
                    <h3 class="locus-title">
                        <span>🎯 独立黄金靶区 ${locus.locusId}</span>
                        <span style="font-size: 14px; color: #7f8c8d; font-weight: normal;">(参考起始坐标: ${locus.anchorStart})</span>
                    </h3>
                `;
                
                locus.variants.forEach((cand, vIndex) => {
                    let isPrimary = vIndex === 0;
                    let cardClass = isPrimary ? 'cand-primary' : 'cand-variant';
                    let roleBadge = isPrimary ? `<span class="role-badge role-main">主力优选</span>` : `<span class="role-badge role-sub">微调备选 ${vIndex}</span>`;

                    reportHTML += `
                    <div class="candidate-card ${cardClass}">
                        <h4>
                            <span>${roleBadge} 综合得分: <span class="score-badge">${cand.score.toFixed(1)} 分</span></span>
                            <span style="font-size: 13px; font-weight: normal; color: #95a5a6;">精确定位: ${cand.start} | 产物长度: ${cand.size} bp</span>
                        </h4>
                        
                        <div class="seq-block">
                            ${generateSeqRows('Forward', cand.fwd)}
                            ${getStatHTML(cand.fStats, cand.fwd.length > 1)}
                        </div>
                        <div class="seq-block probe-block">
                            ${generateSeqRows('Probe', cand.probe)}
                            ${getStatHTML(cand.pStats, cand.probe.length > 1)}
                        </div>
                        <div class="seq-block">
                            ${generateSeqRows('Reverse', cand.rev)}
                            ${getStatHTML(cand.rStats, cand.rev.length > 1)}
                        </div>

                        <details>
                            <summary>🔍 查看当前变体柔性评分明细 (纯净单套优选版)</summary>
                            <div class="details-content">
                                <div class="score-item"><span>基础匹配分 (探针 3 倍权重)</span><span class="score-plus">+${cand.details.base.toFixed(1)}</span></div>
                                <div class="score-item"><span>探针卓越奖励 (完美率 ≥99%)</span><span class="score-plus">+${cand.details.pBonus.toFixed(1)}</span></div>
                                <div class="score-item"><span>探针错配惩罚 (低于98%十倍扣除)</span><span class="${cand.details.pPenalty === 0 ? '' : 'score-minus'}">${cand.details.pPenalty.toFixed(1)}</span></div>
                                <div class="score-item"><span>F/R 混合套数重罚 (-35分/次，极力优选单套)</span><span class="${(cand.details.mixF + cand.details.mixR) === 0 ? '' : 'score-minus'}">${cand.details.mixF + cand.details.mixR}</span></div>
                                <div class="score-item" style="border-top: 1px solid #ddd; padding-top: 4px; font-weight: bold;"><span>柔性规则偏离总扣分 (GC/温差/3'末位等)</span><span class="${cand.details.softPen === 0 ? '' : 'score-minus'}">${cand.details.softPen.toFixed(1)}</span></div>
                            </div>
                        </details>
                    </div>`;
                });
                
                reportHTML += `</div>`; 
            });
        }
        
        document.getElementById('candidatesOutput').innerHTML = reportHTML;
        document.getElementById('reportArea').style.display = 'block';
    }

    function downloadCSV() {
        if (!globalLociGroups || globalLociGroups.length === 0) {
            alert("目前没有可导出的数据，请先运行寻优！");
            return;
        }

        let csvContent = "\uFEFF"; 
        csvContent += "靶区归属,变体角色,综合得分,寡核苷酸类型,序列 (5'->3'),长度 (bp),Tm (°C),GC (%),完美匹配(0),错配1碱基(1),错配2碱基(2),错配3碱基(≥3),预期产物长度 (bp),精确起始坐标\n";

        globalLociGroups.forEach((locus) => {
            locus.variants.forEach((cand, vIndex) => {
                const locusName = `靶区_${locus.locusId}`;
                const role = vIndex === 0 ? "主力优选" : `备选_${vIndex}`;
                const score = cand.score.toFixed(1);
                const size = cand.size;
                const start = cand.start;

                cand.fwd.forEach((seq, i) => {
                    let type = cand.fwd.length > 1 ? `Forward_${i+1}` : `Forward`;
                    let stats = cand.fStats;
                    csvContent += `${locusName},${role},${score},${type},${seq},${seq.length},${calcTm(seq)},${calcGC(seq)},${stats.p0}% (${stats.m0}/${stats.total}),${stats.p1}% (${stats.m1}/${stats.total}),${stats.p2}% (${stats.m2}/${stats.total}),${stats.p3}% (${stats.m3p}/${stats.total}),${size},${start}\n`;
                });

                cand.probe.forEach((seq, i) => {
                    let type = cand.probe.length > 1 ? `Probe_${i+1}` : `Probe`;
                    let stats = cand.pStats;
                    csvContent += `${locusName},${role},${score},${type},${seq},${seq.length},${calcTm(seq)},${calcGC(seq)},${stats.p0}% (${stats.m0}/${stats.total}),${stats.p1}% (${stats.m1}/${stats.total}),${stats.p2}% (${stats.m2}/${stats.total}),${stats.p3}% (${stats.m3p}/${stats.total}),${size},${start}\n`;
                });

                cand.rev.forEach((seq, i) => {
                    let type = cand.rev.length > 1 ? `Reverse_${i+1}` : `Reverse`;
                    let stats = cand.rStats;
                    csvContent += `${locusName},${role},${score},${type},${seq},${seq.length},${calcTm(seq)},${calcGC(seq)},${stats.p0}% (${stats.m0}/${stats.total}),${stats.p1}% (${stats.m1}/${stats.total}),${stats.p2}% (${stats.m2}/${stats.total}),${stats.p3}% (${stats.m3p}/${stats.total}),${size},${start}\n`;
                });
            });
        });

        const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.setAttribute("href", url);
        const dateStr = new Date().toISOString().slice(0, 10);
        link.setAttribute("download", `V02_纯净单套柔性加权_DOE清单_${dateStr}.csv`);
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }

    function renderChart(data) {
        const ctx = document.getElementById('entropyChart').getContext('2d');
        if (chartInstance) chartInstance.destroy();
        chartInstance = new Chart(ctx, {
            type: 'bar',
            data: {
                labels: Array.from({length: data.length}, (_, i) => i + 1),
                datasets: [{ label: '香农熵 (突变强度扫描)', data: data, backgroundColor: data.map(v => v > 0.05 ? '#e74c3c' : '#3498db') }]
            },
            options: { responsive: true, maintainAspectRatio: false }
        });
    }
</script>
</body>
</html>
"""

# Streamlit页面布局
st.set_page_config(page_title="V02 原料齐套 - 工业级柔性优先寻优系统 V15", layout="wide")
st.title("🧬 自动化引物探针柔性优先寻优系统 V15")
st.caption("纯净单套优选版本 | 70-150bp产物无偏向、重罚混合引物")

# 嵌入完整HTML页面，高度适配大屏
st.components.v1.html(full_html, height=1400, scrolling=True)
