import re
import os 
import math
from openpyxl import load_workbook
from pathlib import Path


file_path = r'c:\研究\test'
numbers = []
coordinates = []
atoms = []
capture = None

for path in Path(file_path).iterdir():

    with open(path, "r", encoding = "utf-8") as f:
        for line in f :
            # "Atomic Numbers"が出てきたら、次の行から取り出し
            if "Atomic numbers" in line:
                capture = "atomic"
                continue
            
            # "Nuclear charges"が出てきたら次のif文まで処理を飛ばす
            if "Nuclear charges" in line:
                capture = None
                continue
            
            # "Current cartesian coordinates"が出てきたら処理を開始
            if "Current cartesian coordinates" in line:
                capture = "coordinate"
                continue
            
            # "Number of symbols"が出てきたら処理を終了
            if "Number of symbols" in line:
                capture = None
                break
            
            if capture == "atomic":
                parts = line.strip().split() #文字列として原子番号を追加
                numbers.extend(parts)
            
            if capture == "coordinate":
                parts = line.strip().split()
                coordinates.extend(map(float, parts)) #数値として座標を追加
                
                # 得られた座標を(x,y,z)ごとにリストにする
                tri_coordinates = []
                for i in range(0, len(coordinates), 3):
                    tri_coordinates.append(tuple(coordinates[i:i+3]))

    # 源氏番号(Z)と座標(coord)をリストにして出力
    for z, coord in zip(numbers, tri_coordinates):
        atoms.append({"Z": z, "coord":coord})

    # 最小距離を求める    
    def distance(c1, c2):
        return math.dist(c1, c2)

    targets = 4 #最小距離を求めたい原子数
    target_atoms = atoms[-targets:] #最小距離を求めたい原子のみのリスト
    other_atoms = atoms[:-targets]
    results = []
            
    for i, atom_i in enumerate(target_atoms):
        min_dist = float("inf")
        closet_atom = None
        for j, atom_j in enumerate(other_atoms):
            # 同じ原子番号の計算は行わない
            if atom_i["Z"] == atom_j["Z"]:
                continue
            d = distance(atom_i["coord"], atom_j["coord"])
            if d < min_dist:
                min_dist = d
                closet_atom = atom_j
        results.append({
            "atoms" : atom_i,
            "closet_atom" : closet_atom,
            "min_distance" : min_dist
        })
            
            
    print(results)
    print(len(results))