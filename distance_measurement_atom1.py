import os
import sys
import fileinput
import math
import re
import numpy as np
from pathlib import Path

def main():
	# 引数１つ目のファイルを読み込む
	if Path(sys.argv[1]).exists():
		coordinates = []
		with open(sys.argv[1], "r") as f:
			file = f.read()

			target1 = 'Current cartesian coordinates'
			index1 = file.find(target1)
			target2 = 'Force Field'
			index2 = file.find(target2)
			coordinates_before = file[index1:index2]
			coordinates = re.findall(r"\S+", coordinates_before)
			del coordinates[0:5]
			# print(coordinates)
		# Current cartesian coordinates を読み取って、その後ろのNがいくつか読み取る
		coordinates_number = int(coordinates[0])
		coordinates_h2o_number = coordinates_number - 3
		atom_sum = int(coordinates_number / 3)
		atom_x = float(coordinates[coordinates_number - 2])
		atom_y = float(coordinates[coordinates_number - 1])
		atom_z = float(coordinates[coordinates_number])
		# print(f'atom_x = {atom_x}, atom_y= {atom_y}, atom_z = {atom_z}')

		# 原子番号の順番を元に、調べたい距離を調べる
		# atom = np.array([x,y,z])

		# 引数3個ずつ取って、atom x,y,zとの距離を出す
		# h2oのo原子の最短距離
		min_o = 100
		for i in range(1, coordinates_h2o_number, 9):
			# print(f'{i}番目')
			h2o_o_x = float(coordinates[i])
			h2o_o_y = float(coordinates[i + 1])
			h2o_o_z = float(coordinates[i + 2])
			# print(f'h2o_o_x = {h2o_o_x}, h20_o_y = {h2o_o_y}, h2o_o_z = {h2o_o_z}')
			distance = math.sqrt((atom_x - h2o_o_x)**2 + (atom_y - h2o_o_y)**2 + (atom_z - h2o_o_z)**2)
			if(min_o > distance):
				min_o = distance
		# h2oのh原子1つ目の最短距離
		min_h = 100	 
		for i in range(4, coordinates_h2o_number, 9):
			# print(f'{i}番目')
			h2o_h1_x = float(coordinates[i])
			h2o_h1_y = float(coordinates[i + 1])
			h2o_h1_z = float(coordinates[i + 2])
			# print(f'h2o_h1_x = {h2o_h1_x}, h2o_h1_y = {h2o_h1_y}, h2o_h1_z = {h2o_h1_z}')

			# distance = np.linalg.norm(atom) で距離出す方法もありそう
			distance = math.sqrt((atom_x - h2o_h1_x)**2 + (atom_y - h2o_h1_y)**2 + (atom_z - h2o_h1_z)**2)
			# print(distance)
			# 繰り返し、大きさ比べて小さい方を保存
			if(min_h > distance):
				min_h = distance
		# h2oのn原子1つ目の最短距離	
		for i in range(7, coordinates_h2o_number, 9):
			# print(f'{i}番目')
			h2o_h2_x = float(coordinates[i])
			h2o_h2_y = float(coordinates[i + 1])
			h2o_h2_z = float(coordinates[i + 2])
			# print(f'h2o_x = {h2o_x}, h2o_y = {h2o_y}, h2o_z = {h2o_z}')

			# distance = np.linalg.norm(atom) で距離出す方法もありそう
			distance = math.sqrt((atom_x - h2o_h2_x)**2 + (atom_y - h2o_h2_y)**2 + (atom_z - h2o_h2_z)**2)
			# print(distance)
			# 繰り返し、大きさ比べて小さい方を保存
			if(min_h > distance):
				min_h = distance

		# ボーア単位をオングストロームに直す
		min_o = min_o * 0.5291772083 
		min_h = min_h * 0.5291772083 
		print(f',min_o_distance = {min_o}, min_h_distance = {min_h}')
	else:
		print('Error')


if __name__ == "__main__":
    main()
