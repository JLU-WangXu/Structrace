# -*- coding: utf-8 -*-
import os
import requests
import numpy as np
import math
import random
import pandas as pd
from tqdm import tqdm
import logging

# --- 配置日志 ---
# 配置日志记录，将错误和信息输出到文件和控制台
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    handlers=[
                        logging.FileHandler("processing.log"),
                        logging.StreamHandler()
                    ])

# ==============================================================================
# 1. PDB 数据获取模块
# ==============================================================================

def get_all_pdb_ids():
    """从RCSB PDB通过新的搜索API获取所有PDB ID的列表。"""
    logging.info("正在从 RCSB PDB 的新版搜索API获取所有 PDB ID 列表...")
    
    # 新的API地址
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # 构建查询语句，意思是“找到所有存在的PDB条目”
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "operator": "exists",
                "attribute": "rcsb_entry_container_identifiers.entry_id"
            }
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True  # 这个选项告诉API返回所有结果，而不仅仅是第一页
        }
    }

    try:
        # 使用POST方法发送JSON查询
        response = requests.post(search_url, json=query, timeout=60)
        response.raise_for_status()  # 检查请求是否成功
        
        data = response.json()
        
        # 从返回的结果中提取PDB ID
        # 新API返回的每个条目是一个字典，ID在'identifier'字段里
        ids = [item['identifier'] for item in data['result_set']]
        
        logging.info(f"成功获取 {len(ids)} 个 PDB ID。")
        return ids
    except requests.exceptions.RequestException as e:
        logging.error(f"获取 PDB ID 列表失败: {e}")
        return []
    except KeyError:
        # 如果返回的JSON格式不对，捕获KeyError
        logging.error("解析PDB ID列表失败，返回的数据格式可能已更改。")
        return []


def download_pdb(pdb_id, directory):
    # """根据PDB ID下载PDB文件到指定目录。"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    save_path = os.path.join(directory, f"{pdb_id}.pdb")
    if os.path.exists(save_path):
        # logging.info(f"文件 {pdb_id}.pdb 已存在，跳过下载。")
        return save_path

    try:
        response = requests.get(url, timeout=20)
        response.raise_for_status()  # 如果下载失败（如404），则会抛出异常
        with open(save_path, 'w') as f:
            f.write(response.text)
        return save_path
    except requests.exceptions.RequestException as e:
        logging.warning(f"下载 {pdb_id}.pdb 失败: {e}")
        return None

# ==============================================================================
# 2. 核心水印算法模块 (基于您的代码重构)
# ==============================================================================

def extract_ca_from_lines(lines):
    # """从PDB文件行中提取所有CA原子信息。"""
    data = []
    for ln in lines:
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            try:
                b = float(ln[60:66])
                x, y, z = map(float, (ln[30:38], ln[38:46], ln[46:54]))
                data.append({"x": x, "y": y, "z": z, "b": b, "line": ln})
            except (ValueError, IndexError):
                # 忽略格式不正确的行
                continue
    return data

def embed_watermark(pdb_path, watermark_bits, output_path, start_freq=5, strength=0.02):
    """
    对一个PDB文件嵌入二进制水印。

    Args:
        pdb_path (str): 原始PDB文件路径。
        watermark_bits (str): 要嵌入的二进制字符串 (e.g., "1011")。
        output_path (str): 输出的带水印PDB文件路径。
        start_freq (int): FFT中开始修改的频率索引。
        strength (float): 嵌入强度。

    Returns:
        float: 计算出的全局CA RMSD，如果失败则返回 None。
    """
    try:
        with open(pdb_path, "r") as f:
            lines = f.readlines()

        all_ca = extract_ca_from_lines(lines)
        if not all_ca:
            raise ValueError("PDB文件中未找到CA原子。")

        bits_len = len(watermark_bits)
        needed_freqs = math.ceil(bits_len / 3)
        # N ≥ 2*(ceil(bits_len/3) + START_FREQ)
        min_top_n = 2 * (needed_freqs + start_freq)
        top_n = max(35, min_top_n)

        if len(all_ca) < top_n:
            raise ValueError(f"蛋白质残基数 ({len(all_ca)}) 过少，无法满足嵌入 {bits_len} 位所需的最少 {top_n} 个CA原子。")
        
        cas = sorted(all_ca, key=lambda x: x["b"], reverse=True)[:top_n]
        
        X = np.array([c["x"] for c in cas])
        Y = np.array([c["y"] for c in cas])
        Z = np.array([c["z"] for c in cas])
        
        FX, FY, FZ = np.fft.fft(X), np.fft.fft(Y), np.fft.fft(Z)

        # --- 嵌入过程 ---
        bits = [int(b) for b in watermark_bits]
        N = len(FX)
        FX1, FY1, FZ1 = FX.copy(), FY.copy(), FZ.copy()

        for i, bit in enumerate(bits):
            freq = start_freq + i // 3
            arr = [FX1, FY1, FZ1][i % 3]
            mag, phi = np.abs(arr[freq]), np.angle(arr[freq])
            mag_new = mag + (strength if bit else -strength)
            arr[freq] = mag_new * np.exp(1j * phi)
            arr[N - freq] = mag_new * np.exp(-1j * phi) # 保持共轭对称
        
        X1, Y1, Z1 = np.fft.ifft(FX1).real, np.fft.ifft(FY1).real, np.fft.ifft(FZ1).real

        # --- 写回PDB并计算RMSD ---
        new_lines = lines.copy()
        line_indices_map = {id(c["line"]): i for i, c in enumerate(cas)}
        
        # 建立一个从行内容到其在原始文件中的索引的映射，以处理重复行
        line_to_indices = {}
        for i, line in enumerate(lines):
            line_to_indices.setdefault(line, []).append(i)
        
        cas_line_indices_used = {}
        for i, c in enumerate(cas):
            ln = c["line"]
            # 获取该行所有出现的位置
            possible_indices = line_to_indices[ln]
            # 找到一个尚未被使用的位置索引
            idx_in_lines = -1
            for pi in possible_indices:
                if cas_line_indices_used.get(pi) is None:
                    idx_in_lines = pi
                    cas_line_indices_used[pi] = True
                    break
            
            if idx_in_lines == -1: continue # 以防万一找不到

            xs, ys, zs = f"{X1[i]:8.3f}", f"{Y1[i]:8.3f}", f"{Z1[i]:8.3f}"
            new_lines[idx_in_lines] = ln[:30] + xs + ys + zs + ln[54:]
        
        with open(output_path, "w") as f:
            f.writelines(new_lines)

        # --- 计算全局RMSD ---
        orig_coords = np.array([[c["x"], c["y"], c["z"]] for c in extract_ca_from_lines(lines)])
        encoded_coords = np.array([[c["x"], c["y"], c["z"]] for c in extract_ca_from_lines(new_lines)])
        
        if orig_coords.shape != encoded_coords.shape:
             raise ValueError("原始和编码后的CA原子数量不匹配，无法计算RMSD。")

        global_rmsd = np.sqrt(np.mean(np.sum((encoded_coords - orig_coords)**2, axis=1)))
        
        return global_rmsd

    except Exception as e:
        logging.error(f"在处理 {os.path.basename(pdb_path)} 时嵌入失败: {e}")
        return None

def decode_watermark(original_pdb_path, encoded_pdb_path, bit_length, start_freq=5):
    """
    从一个已编码的PDB文件中解码出水印。

    Args:
        original_pdb_path (str): 原始PDB文件路径。
        encoded_pdb_path (str): 已编码PDB文件路径。
        bit_length (int): 原始嵌入的水印比特数。

    Returns:
        str: 解码出的二进制字符串, 如果失败则返回 None。
    """
    try:
        with open(original_pdb_path, "r") as f0, open(encoded_pdb_path, "r") as f1:
            lines0, lines1 = f0.readlines(), f1.readlines()

        all_ca0 = extract_ca_from_lines(lines0)
        
        # 计算解码所需的top_n (必须与编码时完全一致)
        needed_freqs = math.ceil(bit_length / 3)
        min_top_n = 2 * (needed_freqs + start_freq)
        top_n = max(35, min_top_n)

        if len(all_ca0) < top_n:
             raise ValueError(f"蛋白质残基数 ({len(all_ca0)}) 过少，无法满足解码 {bit_length} 位所需的最少 {top_n} 个CA原子。")
        
        cas0 = sorted(all_ca0, key=lambda x: x["b"], reverse=True)[:top_n]
        
        # 关键：解码时必须使用原始PDB的B-factor顺序来提取已编码的坐标
        orig_lines_map = {c['line']: c for c in cas0}
        cas1_dict = {ln['line']: ln for ln in extract_ca_from_lines(lines1)}
        
        cas1 = []
        for c0 in cas0:
            if c0['line'] in cas1_dict:
                 cas1.append(cas1_dict[c0['line']])
            else: # 如果行被修改，通过坐标近似查找
                found = False
                for c1_key, c1_val in cas1_dict.items():
                    if c1_key[17:26] == c0['line'][17:26]: # 相同的残基序号和链ID
                        cas1.append(c1_val)
                        found = True
                        break
                if not found:
                    raise ValueError("无法在编码文件中找到与原始文件对应的CA原子。")


        X0 = np.array([c["x"] for c in cas0])
        Y0 = np.array([c["y"] for c in cas0])
        Z0 = np.array([c["z"] for c in cas0])
        X1 = np.array([c["x"] for c in cas1])
        Y1 = np.array([c["y"] for c in cas1])
        Z1 = np.array([c["z"] for c in cas1])

        FX0, FY0, FZ0 = np.fft.fft(X0), np.fft.fft(Y0), np.fft.fft(Z0)
        FX1, FY1, FZ1 = np.fft.fft(X1), np.fft.fft(Y1), np.fft.fft(Z1)

        # --- 解码过程 ---
        decoded_bits = []
        for i in range(bit_length):
            freq = start_freq + i // 3
            arr0 = [FX0, FY0, FZ0][i % 3]
            arr1 = [FX1, FY1, FZ1][i % 3]
            delta = np.abs(arr1[freq]) - np.abs(arr0[freq])
            decoded_bits.append("1" if delta > 0 else "0")
            
        return "".join(decoded_bits)

    except Exception as e:
        logging.error(f"在解码 {os.path.basename(encoded_pdb_path)} 时失败: {e}")
        return None

# ==============================================================================
# 3. 主流水线模块
# ==============================================================================

def generate_random_bit_string(length):
    """生成指定长度的随机二进制字符串。"""
    return ''.join(random.choice('01') for _ in range(length))

def main():
    # --- 参数配置 ---
    NUM_PDB_TO_PROCESS = 10000
    BIT_LENGTHS_TO_TEST = [4, 8, 12, 16]
    
    # --- 目录设置 ---
    OUTPUT_DIR = "watermarking_output"
    ORIGINALS_DIR = os.path.join(OUTPUT_DIR, "pdb_originals")
    WATERMARKED_DIR = os.path.join(OUTPUT_DIR, "pdb_watermarked")
    RESULTS_CSV_PATH = os.path.join(OUTPUT_DIR, "results.csv")

    os.makedirs(ORIGINALS_DIR, exist_ok=True)
    os.makedirs(WATERMARKED_DIR, exist_ok=True)
    
    # --- 获取PDB ID ---
    all_pdb_ids = get_all_pdb_ids()
    if not all_pdb_ids:
        logging.critical("无法获取PDB ID列表，程序退出。")
        return
        
    random.shuffle(all_pdb_ids) # 随机打乱，避免总是处理同一批蛋白
    pdb_ids_to_process = all_pdb_ids[:NUM_PDB_TO_PROCESS]

    # --- 初始化结果记录 ---
    results_data = []
    # 如果结果文件已存在，可以加载它以支持断点续传
    if os.path.exists(RESULTS_CSV_PATH):
        logging.info(f"发现已存在的结果文件 {RESULTS_CSV_PATH}，将在此基础上追加。")
        df_existing = pd.read_csv(RESULTS_CSV_PATH)
        processed_keys = set(df_existing[['PDB_ID', 'Bit_Length']].apply(lambda row: (row.PDB_ID, row.Bit_Length), axis=1))
        results_data = df_existing.to_dict('records')
    else:
        processed_keys = set()


    # --- 主处理循环 ---
    logging.info(f"开始处理 {len(pdb_ids_to_process)} 个PDB文件，每个测试 {len(BIT_LENGTHS_TO_TEST)} 种水印长度...")
    
    # 使用tqdm创建总进度条
    pbar = tqdm(total=len(pdb_ids_to_process) * len(BIT_LENGTHS_TO_TEST), desc="总体进度")

    for pdb_id in pdb_ids_to_process:
        # 下载原始PDB文件
        original_pdb_path = download_pdb(pdb_id, ORIGINALS_DIR)
        if not original_pdb_path:
            pbar.update(len(BIT_LENGTHS_TO_TEST)) # 更新进度条，跳过这个PDB
            continue
            
        for bits in BIT_LENGTHS_TO_TEST:
            pbar.set_description(f"处理 {pdb_id} ({bits}-bit)")

            # --- 断点续传检查 ---
            if (pdb_id, bits) in processed_keys:
                pbar.update(1)
                continue

            # 1. 生成随机水印
            watermark = generate_random_bit_string(bits)
            
            # 2. 定义输出文件名
            output_filename = f"{pdb_id}_{bits}bit_{watermark}.pdb"
            output_pdb_path = os.path.join(WATERMARKED_DIR, output_filename)
            
            # 3. 嵌入水印并获取RMSD
            global_rmsd = embed_watermark(original_pdb_path, watermark, output_pdb_path)
            
            if global_rmsd is None:
                # 嵌入失败，记录空结果
                verification_status = "Embed_Failed"
                decoded_watermark = "N/A"
            else:
                # 4. 解码验证
                decoded_watermark = decode_watermark(original_pdb_path, output_pdb_path, bits)
                if decoded_watermark is None:
                    verification_status = "Decode_Failed"
                elif decoded_watermark == watermark:
                    verification_status = "Success"
                else:
                    verification_status = "Mismatch"

            # 5. 记录结果
            result_entry = {
                "PDB_ID": pdb_id,
                "Bit_Length": bits,
                "Watermark_Embedded": watermark,
                "Watermark_Decoded": decoded_watermark,
                "Verification_Status": verification_status,
                "Global_CA_RMSD": round(global_rmsd, 5) if global_rmsd is not None else "N/A"
            }
            results_data.append(result_entry)
            
            # 6. 实时保存结果到CSV
            pd.DataFrame(results_data).to_csv(RESULTS_CSV_PATH, index=False)

            # 更新进度条
            pbar.update(1)

    pbar.close()
    logging.info("所有处理完成！")
    logging.info(f"原始PDB文件保存在: {ORIGINALS_DIR}")
    logging.info(f"水印PDB文件保存在: {WATERMARKED_DIR}")
    logging.info(f"详细分析结果保存在: {RESULTS_CSV_PATH}")


if __name__ == "__main__":
    main()