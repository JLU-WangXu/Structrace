import numpy as np
import math
import os
import requests
import traceback

print("所有库导入成功。")

# --- 步骤 0: 定义输入数据 ---

# --- 摘要文本 --- 
ABSTRACT_TO_EMBED = "Cleavage of amyloid precursor protein (APP) by the intramembrane protease g-secretase islinked to Alzheimer’s disease (AD). We report an atomic structure of human g-secretase incomplex with a transmembrane (TM) APP fragment at 2.6-angstrom resolution. The TMhelix of APP closely interacts with five surrounding TMs of PS1 (the catalytic subunit ofg-secretase). A hybrid b sheet, which is formed by a b strand from APP and two b strandsfrom PS1, guides g-secretase to the scissile peptide bond of APP between its TM and b strand.Residues at the interface between PS1 and APP are heavily targeted by recurring mutationsfrom AD patients. This structure, together with that of g-secretase bound to Notch, revealcontrasting features of substrate binding, which may be applied toward the design ofsubstrate-specific inhibitors"

# --- PDB 设置 ---
PDB_ID = "6iyc"
ORIGINAL_PDB_PATH = f"{PDB_ID}.pdb"
ENCODED_PDB_PATH = f"{PDB_ID}_abstract_watermarked.pdb"

print(f"摘要长度: {len(ABSTRACT_TO_EMBED)} 字符")
print(f"原始 PDB: {ORIGINAL_PDB_PATH}")
print(f"输出 PDB: {ENCODED_PDB_PATH}")


# --- 自动下载 PDB 文件 ---
if not os.path.exists(ORIGINAL_PDB_PATH):
    print(f"正在从 RCSB 下载 {ORIGINAL_PDB_PATH}...")
    url = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
    try:
        r = requests.get(url)
        r.raise_for_status() # 确保请求成功
        with open(ORIGINAL_PDB_PATH, 'w') as f:
            f.write(r.text)
        print(f"成功下载 {ORIGINAL_PDB_PATH}。")
    except Exception as e:
        print(f"错误: 下载 PDB 文件失败。 {e}")
else:
    print(f"文件 {ORIGINAL_PDB_PATH} 已存在，跳过下载。")


# --- 步骤 1: 定义核心算法函数 ---

# --- 常量 --- 
START_FREQ = 5      # 嵌入的起始频率 (避开低频)
STRENGTH   = 1.0    # 嵌入强度 (扰动幅度) -  原 0.02, 因四舍五入失败

# --- 函数 1: 提取所有重原子 ---
def extract_all_heavy_atoms(lines):
    """按文件顺序提取所有非氢 ATOM 原子及其 B 因子和坐标"""
    data = []
    for ln in lines:
        # 确保是 ATOM 记录，而不是 HETATM
        if ln.startswith("ATOM"):
            try:
                # 检查元素符号 (col 77-78)，排除氢 (H)
                element = ln[76:78].strip()
                if element == 'H':
                    continue
                
                atom_serial = ln[6:11] # 提取原子序列号 (e.g. " 123 ")
                b = float(ln[60:66])
                x,y,z = map(float,(ln[30:38],ln[38:46],ln[46:54]))
                data.append({"line":ln, "b":b, "x":x, "y":y, "z":z, "serial": atom_serial})
            except (ValueError, IndexError):
                # 跳过格式不正确的行或没有元素符号的行
                continue
    return data

# --- 函数 2: 智能选择嵌入原子 ---
def select_loop_atoms(ca_list, bits_len):
    """
    1) 自动计算嵌入 bits_len 所需的原子数 (N)
    2) 按 B 因子降序选取前 N 个 CA 原子
    """
    # 我们需要 bits_len 个比特。每个频率系数嵌入 3 个比特 (X, Y, Z)。
    needed_freqs = math.ceil(bits_len / 3)
    # 傅里叶变换的长度 N 必须满足: N >= 2 * (needed_freqs + START_FREQ)
    needed_N = max(35, 2 * (needed_freqs + START_FREQ))

    if len(ca_list) < needed_N:
        raise ValueError(f"符合条件的原子数太少 ({len(ca_list)})。 "
                         f"无法嵌入 {bits_len} bits (需要 {needed_N} 个原子)。")

    print(f"🔧 嵌入 {bits_len} bits (摘要 + 终止符)...")
    print(f"🔧 需要 {needed_freqs} 个频率系数, 自动选取 top_n = {needed_N} 个原子 (B-factor 最高)。")
    
    # 按 B 因子(b)降序排序，选取前 N 个
    loop = sorted(ca_list, key=lambda c: c["b"], reverse=True)[:needed_N]
    return loop

# --- 函数 3: FFT 嵌入消息 ---
def embed_message(FX, FY, FZ, message):
    """在三条 FFT 频谱振幅上嵌入 message + NULL 终止符 (0x00)"""
    # 组装比特流 (使用 utf-8 编码摘要, 忽略无法编码的字符, 最后 + 0x00 终结)
    msg_bytes = list(message.encode("utf-8", "ignore")) + [0]
    bits = [int(b) for byte in msg_bytes for b in f"{byte:08b}"]
    
    N = len(FX)
    max_bits = (math.floor(N/2) - START_FREQ) * 3
    if len(bits) > max_bits:
        raise ValueError(f"消息过长：需要 {len(bits)} bits, 但只有 {max_bits} bits 可用。")

    FX1, FY1, FZ1 = FX.copy(), FY.copy(), FZ.copy()
    for i, bit in enumerate(bits):
        freq = START_FREQ + i // 3
        arr = (FX1, FY1, FZ1)[i % 3]
        mag, phi = abs(arr[freq]), np.angle(arr[freq])
        
        # 1 -> 增加幅度; 0 -> 减小幅度
        mag_new = mag + (STRENGTH if bit else -STRENGTH)
        
        arr[freq]   = mag_new * np.exp(1j*phi)
        arr[N-freq] = mag_new * np.exp(-1j*phi) # 保持共轭对称以确保反变换为实数
        
    return FX1, FY1, FZ1

print("核心算法函数已定义。")


# --- 步骤 2: 执行编码 (嵌入摘要) ---
def run_encoding():
    print("\n--- 步骤 2: 开始编码（嵌入摘要）---")
    try:
        # 1. 读入 PDB
        with open(ORIGINAL_PDB_PATH) as f:
            lines = f.readlines()
        
        # 2. 提取所有重原子
        all_atoms = extract_all_heavy_atoms(lines)
        if not all_atoms:
            print(f"错误: 在 {ORIGINAL_PDB_PATH} 中未找到符合条件的 ATOM 原子。")
            return False
        else:
            print(f"在 {ORIGINAL_PDB_PATH} 中找到 {len(all_atoms)} 个重原子。")

        # 3. 根据消息长度和 B 因子选 loop 区域
        # 消息比特长度 = (字符数 * 8) + 8 (用于 0x00 终止符)
        bits_len = len(ABSTRACT_TO_EMBED.encode("utf-8", "ignore")) * 8 + 8
        loop_atoms = select_loop_atoms(all_atoms, bits_len)

        # 4. FFT → 嵌入 → IFFT
        X = np.array([c["x"] for c in loop_atoms])
        Y = np.array([c["y"] for c in loop_atoms])
        Z = np.array([c["z"] for c in loop_atoms])
        
        FX, FY, FZ = np.fft.fft(X), np.fft.fft(Y), np.fft.fft(Z)
        FX1, FY1, FZ1 = embed_message(FX, FY, FZ, ABSTRACT_TO_EMBED)
        
        X1 = np.fft.ifft(FX1).real
        Y1 = np.fft.ifft(FY1).real
        Z1 = np.fft.ifft(FZ1).real

        # 5. 写回新坐标 (使用基于序列号的稳健方法)
        # 创建一个 (序列号 -> 新坐标行) 的映射
        new_coords_map = {}
        for i, c in enumerate(loop_atoms): # loop_atoms 已按 B-factor 排序
            xs, ys, zs = f"{X1[i]:8.3f}", f"{Y1[i]:8.3f}", f"{Z1[i]:8.3f}"
            new_line = c["line"][:30] + xs + ys + zs + c["line"][54:]
            new_coords_map[c["serial"]] = new_line # 键是原子序列号 (e.g. " 123 ")

        # 遍历原始 PDB 的每一行
        new_lines = []
        for ln in lines:
            if ln.startswith("ATOM"):
                try:
                    # 提取这一行的序列号
                    current_serial = ln[6:11]
                    # 检查这个序列号是否在我们修改过的原子列表里
                    if current_serial in new_coords_map:
                        # 如果是，就用新行
                        new_lines.append(new_coords_map[current_serial])
                    else:
                        # 如果不是，就用原始行
                        new_lines.append(ln)
                except:
                    # 格式有问题的行，直接添加
                    new_lines.append(ln)
            else:
                # 非 ATOM 行 (如 REMARK, HETATM 等)，直接添加
                new_lines.append(ln)

        with open(ENCODED_PDB_PATH, "w") as f:
            f.writelines(new_lines)

        # 6. 报告扰动 (RMSD)
        disp = np.sqrt((X1-X)**2 + (Y1-Y)**2 + (Z1-Z)**2)
        rmsd_loop = np.sqrt(np.mean(disp**2))
        
        all_atoms_new = extract_all_heavy_atoms(new_lines)
        coords0 = np.array([[c["x"],c["y"],c["z"]] for c in all_atoms])
        coords1 = np.array([[c["x"],c["y"],c["z"]] for c in all_atoms_new])
        
        if coords0.shape != coords1.shape:
            print("警告: 编码前后原子数量不匹配，全局 RMSD 可能不准确。")
            min_len = min(len(coords0), len(coords1))
            coords0 = coords0[:min_len]
            coords1 = coords1[:min_len]
            
        global_rmsd = np.sqrt(np.mean(np.sum((coords1 - coords0)**2, axis=1)))
        
        print("="*50)
        print(f"✅ 编码完成，保存到：{ENCODED_PDB_PATH}")
        print(f"  扰动区域 RMSD (仅 {len(loop_atoms)} 个原子): {rmsd_loop:.6f} Å (Max {disp.max():.6f} Å)")
        print(f"  全局重原子 RMSD (所有 {len(all_atoms)} 个原子): {global_rmsd:.6f} Å")
        print("="*50)
        return True

    except Exception as e:
        print(f"编码过程中发生严重错误: {e}")
        traceback.print_exc()
        return False

# --- 步骤 3: 执行解码 (逆向工程) ---
def run_decoding():
    print("\n--- 步骤 3: 开始解码（逆向工程）---")
    try:
        # 1. 读两份 PDB
        with open(ORIGINAL_PDB_PATH) as f0, open(ENCODED_PDB_PATH) as f1:
            lines0 = f0.readlines()
            lines1 = f1.readlines()
        
        atoms0_all = extract_all_heavy_atoms(lines0)
        atoms1_all = extract_all_heavy_atoms(lines1)
        
        if not atoms0_all:
            raise FileNotFoundError(f"原始 PDB {ORIGINAL_PDB_PATH} 中没有 CA 原子或为空。")
        if not atoms1_all:
            raise FileNotFoundError(f"编码后 PDB {ENCODED_PDB_PATH} 中没有 CA 原子或为空。")

        # 2. 用与编码时 *完全相同* 的逻辑选出 loop0
        bits_len_to_decode = len(ABSTRACT_TO_EMBED.encode("utf-8", "ignore")) * 8 + 8
        loop0 = select_loop_atoms(atoms0_all, bits_len_to_decode)
        
        # --- 关键修复：重建 loop1 以匹配 loop0 的顺序 ---
        # 1. 创建一个编码后原子的查找表，以序列号为键
        atoms1_map_by_serial = {atom["serial"]: atom for atom in atoms1_all}
        
        # 2. 遍历 loop0 (已按 B-factor 排序)，并按此顺序构建 loop1
        loop1 = []
        for atom0 in loop0:
            serial_key = atom0["serial"]
            corresponding_atom1 = atoms1_map_by_serial.get(serial_key)
            
            if corresponding_atom1:
                loop1.append(corresponding_atom1)
            else:
                raise ValueError(f"解码失败：在编码文件中未找到序列号为 {serial_key} 的原子。")
        # --- 修复结束 ---

        if len(loop0) != len(loop1):
            raise ValueError("解码时原子数量不匹配，无法继续。")
        
        print(f"🔧 解码时同样选取 {len(loop0)} 个高 B-factor 原子")

        # --- 新增：打印前 5 个原子的 B-factor 和序列号以供验证 ---
        print("\n--- 验证 loop0 和 loop1 的原子是否匹配 (前 5 个) ---")
        print(" B-factor | loop0 (原始) 序列号 | loop1 (编码后) 序列号")
        print("---------------------------------------------------------")
        for i in range(5):
            b0 = loop0[i]['b']
            s0 = loop0[i]['serial']
            s1 = loop1[i]['serial']
            print(f" {b0:>8.3f} | {s0:>18} | {s1:>20}")
        print("---------------------------------------------------------")
        # --- 验证结束 ---


        # 3. FFT
        X0 = np.array([c["x"] for c in loop0])
        Y0 = np.array([c["y"] for c in loop0])
        Z0 = np.array([c["z"] for c in loop0])
        X1 = np.array([c["x"] for c in loop1])
        Y1 = np.array([c["y"] for c in loop1])
        Z1 = np.array([c["z"] for c in loop1])

        FX0, FY0, FZ0 = np.fft.fft(X0), np.fft.fft(Y0), np.fft.fft(Z0)
        FX1, FY1, FZ1 = np.fft.fft(X1), np.fft.fft(Y1), np.fft.fft(Z1)

        # 4. 差分读 bit 流，遇 NULL(0x00) 停
        bits = []
        N = len(FX0)
        max_bits_to_read = (N//2 - START_FREQ) * 3
        
        for i in range(max_bits_to_read):
            freq = START_FREQ + i // 3
            if freq >= N//2: # 停止，以防越界
                break
            
            # --- 关键修复：使用 FZ0/FZ1, 而不是 Z0/Z1 ---
            arr0 = (FX0, FY0, FZ0)[i % 3]
            arr1 = (FX1, FY1, FZ1)[i % 3]
            delta = abs(arr1[freq]) - abs(arr0[freq])
            bits.append(1 if delta > 0 else 0)

        # 5. 字节重组，0x00 终止
        
        # --- 新增：打印原始比特流 ---
        print(f"\n--- 原始比特流 (前 72 bits) ---")
        print("".join(map(str, bits[:72])))
        print("-----------------------------------")
        # --- 打印结束 ---

        msg_bytes = []
        for i in range(0, len(bits), 8):
            if i + 8 > len(bits):
                break # 字节不完整
            byte_str = "".join(str(b) for b in bits[i:i+8])
            byte = int(byte_str, 2)
            
            if byte == 0:
                break # 找到 NULL 终止符
            msg_bytes.append(byte)

        # --- 新增：打印重组的字节列表 ---
        print(f"\n--- 重组的字节列表 (前 10 字节) ---")
        print(msg_bytes[:10])
        print("-----------------------------------")
        # --- 打印结束 ---

        # 6. 解码
        decoded_message = bytes(msg_bytes).decode("utf-8", "ignore")
        
        print("\n" + "="*50)
        print("✅ 解码成功! 提取内容如下:")
        print("="*50)
        print(decoded_message)
        
        # 7. 验证
        if decoded_message == ABSTRACT_TO_EMBED:
            print("\n" + "="*50)
            print("🎉 验证成功: 解码后的摘要与原始摘要完全一致!")
            print("="*50)
        else:
            print("\n" + "="*50)
            print("⚠️ 验证失败: 解码后的摘要与原始摘要不匹配。")
            print("="*50)

    except FileNotFoundError as e:
        print(f"错误: 无法开始解码。文件未找到: {e}")
    except Exception as e:
        print(f"解码过程中发生错误: {e}")
        traceback.print_exc()

# --- 主程序入口 ---
if __name__ == "__main__":
    if os.path.exists(ORIGINAL_PDB_PATH):
        # 首先执行编码
        encoding_success = run_encoding()
        
        # 如果编码成功，则执行解码
        if encoding_success and os.path.exists(ENCODED_PDB_PATH):
            run_decoding()
        elif not encoding_success:
            print("编码失败，跳过解码。")
        else:
            print(f"编码输出文件 {ENCODED_PDB_PATH} 未找到，无法解码。")
    else:
        print(f"原始文件 {ORIGINAL_PDB_PATH} 未找到，请先下载。")