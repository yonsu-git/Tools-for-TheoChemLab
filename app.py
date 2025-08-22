import streamlit as st
import math

# -------------------------
# 距離計算関数
# -------------------------
def distance(c1, c2):
    return math.dist(c1, c2)

# -------------------------
# fchkファイルを解析する関数（簡易版）
# -------------------------
def parse_fchk(file_content):
    numbers = []
    coordinates = []
    capture_coords = False
    for line in file_content.splitlines():
        if "Atomic numbers" in line:
            capture_numbers = True
            continue
        if capture_numbers:
            parts = line.split()
            if parts:
                numbers.extend(parts)
            if "Nuclear charges" in line:
                capture_numbers = False
        if "Cartesian Coordinates" in line:  # 適宜 fchk に合わせて修正
            capture_coords = True
            continue
        if capture_coords:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    coordinates.append(tuple(map(float, parts[:3])))
                except ValueError:
                    continue
    return numbers, coordinates

# -------------------------
# 最小距離計算関数
# -------------------------
def compute_min_distances(numbers, coordinates):
    atoms = [{"Z": z, "coord": c} for z, c in zip(numbers, coordinates)]
    target_atoms = atoms[-4:]
    other_atoms = atoms[:-4]
    
    min_distances = []
    for atom_i in target_atoms:
        min_dist = float("inf")
        for atom_j in other_atoms:
            if atom_i["Z"] == atom_j["Z"]:
                continue
            d = distance(atom_i["coord"], atom_j["coord"])
            if d < min_dist:
                min_dist = d
        min_distances.append(min_dist)
    return min_distances

# -------------------------
# Streamlit UI
# -------------------------
st.title("fchk 最小距離計算アプリ")

uploaded_file = st.file_uploader("fchkファイルをアップロードしてください", type=["fchk"])

if uploaded_file:
    content = uploaded_file.read().decode("utf-8")
    numbers, coordinates = parse_fchk(content)
    
    if len(numbers) < 4 or len(coordinates) < 4:
        st.warning("原子数が4個未満です。")
    else:
        min_distances = compute_min_distances(numbers, coordinates)
        st.write("最後の4個の原子それぞれの最小距離（同じ原子番号を除く）:")
        st.write(min_distances)
