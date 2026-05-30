"""
CIF → PDB 转换模块

使用 Gemmi 读取 mmCIF 文件，输出 PDB 文件。
默认在 CIF 同目录生成同名 .pdb（不传 pdb_path 时）。

依赖: gemmi

用法:
    from cif2pdb import cif2pdb
    pdb = cif2pdb("structure.cif")             # → structure.pdb 在同目录
    pdb = cif2pdb("structure.cif", "a.pdb")    # → a.pdb
"""

import os
import gemmi


def cif2pdb(cif_path, pdb_path=None):
    """
    读取 mmCIF 文件，输出 PDB 文件。

    参数
    ----------
    cif_path : str
        .cif 文件路径。
    pdb_path : str, optional
        PDB 输出路径。不传时自动设为 cif_path 同目录同名 .pdb。

    返回
    -------
    str
        PDB 文件的绝对路径。
    """
    if not os.path.exists(cif_path):
        raise FileNotFoundError(f"文件不存在: {cif_path}")

    try:
        structure = gemmi.read_structure(cif_path)
    except Exception as e:
        raise ValueError(f"Gemmi 无法解析 CIF: {e}")

    if not structure:
        raise ValueError(f"CIF 中无模型数据: {cif_path}")

    if pdb_path is None:
        base, _ = os.path.splitext(cif_path)
        pdb_path = base + ".pdb"

    os.makedirs(os.path.dirname(pdb_path), exist_ok=True)
    structure.write_pdb(pdb_path)
    return os.path.abspath(pdb_path)
