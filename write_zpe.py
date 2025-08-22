import os
import sys
import fileinput
import math
import re
import numpy as np
from pathlib import Path
import random

# avogadro constant
na = 6.02214199 * pow(10, 23)
		
# bohr radius
a0 = 5.291772083 * pow(10, -11)
# [bohr] ---> [angstrom]
b_to_a = a0 * pow(10, 10)

# boltzmann constant
kb = 1.380603 * pow(10, -23)	

# gas constant
R = 8.31446262

# [Hartree] ---> [J]
h_to_j = 4.3597482 * pow(10, -18)

# define a class named "Atom"
class Atom:
    def __init__(self, atomic_weight, coordinate1, coordinate2, coordinate3, Label):
        self.atom = []
        aw = atomic_weight
        self.atom.append(aw)
        x = coordinate1
        self.atom.append(x)
        y = coordinate2
        self.atom.append(y)
        z = coordinate3
        self.atom.append(z)
        label = Label
        self.atom.append(label)
        # atom = [aw, x, y, z, label]
    def get_atom(self):
        return self.atom

# get label from atomic_weight
def get_m_and_l(total_atoms, atomic_numbers):
    ATOMIC_WEIGHT = []
    LABEL = []
    for i in range(0, total_atoms):
        if atomic_numbers[i] == 1:
            ATOMIC_WEIGHT.append(float(1.0078250))
            LABEL.append('H')
        elif atomic_numbers[i] == 8:
            ATOMIC_WEIGHT.append(float(15.9949146))
            LABEL.append('O')
        elif atomic_numbers[i] == 6:
            ATOMIC_WEIGHT.append(float(12.0000000))
            LABEL.append('C')
    return ATOMIC_WEIGHT, LABEL

# get label from atomic_weight
def get_m_and_l_2(total_atoms, atomic_weight):
    ATOMIC_WEIGHT = []
    LABEL = []
    for i in range(0, total_atoms):
        if atomic_weight[i] == 1:
            ATOMIC_WEIGHT.append(float(1.0078250))
            LABEL.append('H')
        elif atomic_weight[i] == 16:
            ATOMIC_WEIGHT.append(float(15.9949146))
            LABEL.append('O')
        elif atomic_weight[i] == 12:
            ATOMIC_WEIGHT.append(float(12.0000000))
            LABEL.append('C')
    return ATOMIC_WEIGHT, LABEL

# get distance between ASW cluster and a radical
def measure_distance(N, info):
    Dist = []
    dist_o_min = 100
    dist_h_min = 100
    num_o = None
    num_h = None
    for i in range(0, N - 60):
        for j in range(0, 20):
            dist_o = math.sqrt((info[3*j].atom[1]-info[60+i].atom[1])**2 + (info[3*j].atom[2]-info[60+i].atom[2])**2 + (info[3*j].atom[3]-info[60+i].atom[3])**2) 
            if dist_o < dist_o_min:
                dist_o_min = dist_o
                num_o = 3*j+1
            dist_o = 100
            dist_h = math.sqrt((info[3*j+1].atom[1]-info[60+i].atom[1])**2 + (info[3*j+1].atom[2]-info[60+i].atom[2])**2 + (info[3*j+1].atom[3]-info[60+i].atom[3])**2) 
            if dist_h < dist_h_min:
                dist_h_min = dist_h
                num_h = 3*j+2
            dist_h = 100
            dist_h = math.sqrt((info[3*j+2].atom[1]-info[60+i].atom[1])**2 + (info[3*j+2].atom[2]-info[60+i].atom[2])**2 + (info[3*j+2].atom[3]-info[60+i].atom[3])**2) 
            if dist_h < dist_h_min:
                dist_h_min = dist_h
                num_h = 3*j+3
            dist_h = 100
        Dist.append([dist_o_min, num_o, dist_h_min, num_h])
        dist_o_min = 100
        dist_h_min = 100
        num_o = None
        num_h = None
    print('dist_o_min'.center(15), 'num_o'.center(5), 'dist_h_min'.center(15), 'num_h'.center(5))
    for i in range(0, N-60):
        print(str(round(Dist[i][0], 10)).center(15), str(Dist[i][1]).center(5), str(round(Dist[i][2], 10)).center(15), str(Dist[i][3]).center(5))
    return Dist

# define a fanction named "make_csvfile"
def make_csvfile(total_atoms, f_name, info, charge, multiplicity):
    f_name_input_number = int(f_name[-6])
    print(f_name_input_number)
    f_name_input_number += 1
    f_name_input = f_name[:-7] + '_' + str(f_name_input_number) + '.gjf'
    f_name_chk = f_name_input.replace('gjf', 'chk')
    with open(f_name_input, "w") as f:
        f.write('%nproc=8\n')
        f.write('%mem=24000mb\n')
        f.write(r'%chk=' + f_name_chk + '\n')
        f.write('#p wb97xd/6-311+(d,p) SCF(Xqc,MaxConv=60) Opt Freq Int(Grid=Ultra) GFinput NoSymm\n\n')
        f.write('RE:\n\n')
        f.write(str(charge) + ' ' + str(multiplicity) + '\n')
        for i in range(0, total_atoms):
            f.write(info[i].atom[4].ljust(8))
            f.write(str(round(info[i].atom[1], 14)).ljust(20))
            f.write(str(round(info[i].atom[2], 14)).ljust(20))
            f.write(str(round(info[i].atom[3], 14)).ljust(20))
            f.write('\n')
        f.write('\n\n\n')

def main():
    # read the file of the first argument
    if Path(sys.argv[1]).exists() and Path(sys.argv[2]).exists():
        f_name = None
        with open(sys.argv[1], "r") as f:
            f_name = (f.name).replace('Wat20.', '')
            f_name = (f_name).replace('_fin.out', '')

            # get the number of atoms
            target_n = ' NAtoms='
            found_n = False
            N = None
            line = f.readline()
            while found_n != True:
                if line.startswith(target_n):
                    # print(re.findall(r"\S+", line))
                    N = int(re.findall(r"\S+", line)[1])
                    # print('Number of Atoms = ', N)
                    found_n = True
                else:
                    line = f.readline()

            n_2 = 0
            if N % 10 == 0:
                n_2 = N/10
            else:
                n_2 = N // 10 + 1
            # print('n_2 : ', n_2)

	    # get the weight of atoms and labels
            AN = []
            target_an = ' IAtWgt='
	    # print('line = ', line)
            count_an = 0
            line = f.readline()
            while count_an < n_2:
                if line.startswith(target_an):
                    # print(line.rstrip("\n"))
                    AN_b = re.findall(r"\S+", line)
                    for i in range(1, len(AN_b)):
                        AN.append(int(AN_b[i]))
                    count_an += 1
                    # print(count_an)
                    line = f.readline()
                else:
                    line = f.readline()
            # print('The length of AN : ', len(AN))
            # print('AN')
            """
            for i in range(0, N):
                print(AN[i])
            """
            M = get_m_and_l_2(N, AN)[0]
            L = get_m_and_l_2(N, AN)[1]
            """
            for i in range(0, N):
                print('Atom', i+1, ' : AW=', M[i], ',L=', L[i])
            """

            target_ip = ' Initial Parameters'
            found_ip = 0
            line = f.readline()
            while found_ip != 2:
                if target_ip in line:
                    found_ip += 1
                    line = f.readline()
                    # print(line.rstrip("\n"))
                else:
                    line = f.readline()
            # get last cartesian coordinates(lcc)
            C = []
            # print(line)
            target_C = 'Coordinates'
            found_C = False
            while found_C != True:
                if target_C in line:
                    found_C = True
                    # print(line.rstrip("\n"))
                else:
                    line = f.readline()
            for i in range(0, 3):
                line = f.readline()
            for i in range(0, N):
                # print(line.rstrip("\n"))
                C_b = re.findall(r"\S+", line)
                for i in range(3, 6):
                    C.append(float(C_b[i]))
                line = f.readline()
            # print('The length of C : ', len(C))
            """
            for i in range(0, N):
                print('atom : ', i+1)
                print('x=', C[3*i], ' y= ', C[3*i+1], ' z= ', C[3*i+2])
            """

            SD1 = None
            target_SD1 = ' SCF Done:'
            found_SD1 = False
            while found_SD1 != True:
                if line.startswith(target_SD1):
                    found_SD1 = True
                else:
                    line = f.readline()
            SD1_b = re.findall(r"\S+", line)
            # print(SD1_b)
            SD1 = float(SD1_b[4])

            # get Mulliken charges
            MC1 = []
            target_mc1 = ' Mulliken charges and spin densities:'
            found_MC1 = False
            while found_MC1 != True:
                if line.startswith(target_mc1):
                    found_MC1 = True
                    line = f.readline()
                else:
                    line = f.readline()
            line = f.readline()
            for i in range(0, N):
                MC1_b = re.findall(r"\S+", line)
                MC1.append(float(MC1_b[2]))
                line = f.readline()
            """
            print('Mulliken Charges')
            for i in range(0, N):
                print(str(MC1[i]))
            """
            # get zero-point energies
            ZPE1 = 0
            line = f.readline()
            target_zpe1 = ' Sum of electronic and zero-point Energies='
            found_ZPE1 = False
            while found_ZPE1 != True:
                if line.startswith(target_zpe1):
                    found_ZPE1 = True
                else:
                    line = f.readline()
            ZPE1_b = re.findall(r"\S+", line)
            # print(line)
            # print('the length of line : ', len(ZPE1_b))
            # print(ZPE1_b)
            ZPE1 = float(ZPE1_b[6])
            line = f.readline()

            E_hartree1 = ZPE1
            E_kelvin1 = E_hartree1 * na * h_to_j / 1000 * 1000 / R
    
            """
            print(len(M))
            print(len(C))
            print(len(L))
            """
            info = []
            for i in range(0, N):
                info.append(Atom(M[i], C[3*i], C[3*i+1], C[3*i+2], L[i]))

        with open(sys.argv[2], "r") as f2:
            f_name_2 = (f2.name).replace('.out', '')

            csv_file_name = f2.name
            # get SCF energy
            SD2 = 0
            line2 = f2.readline()
            target_zpe2 = ' SCF Done'
            found_SD2 = False
            while found_SD2 != True:
                if line2.startswith(target_zpe2):
                    found_SD2 = True
                else:
                    line2 = f2.readline()
            SD2_b = re.findall(r"\S+", line2)
            SD2 = float(SD2_b[4])
            # print(SD2_b)

            # get Mulliken charges
            MC2 = []
            target_mc2 = ' Mulliken charges'
            found_MC2 = False
            while found_MC2 != True:
                if line2.startswith(target_mc2):
                    found_MC2 = True
                    line2 = f2.readline()
                else:
                    line2 = f2.readline()
            line2 = f2.readline()
            for i in range(0, 60):
                MC2_b = re.findall(r"\S+", line2)
                MC2.append(float(MC2_b[2]))
                line2 = f2.readline()
            """
            print('Mulliken Charges')
            for i in range(0, 60):
                print(str(MC2[i]))
            """

            # get zero-point energies
            ZPE2 = 0
            line2 = f2.readline()
            target_zpe2 = ' Sum of electronic and zero-point Energies='
            found_ZPE2 = False
            while found_ZPE2 != True:
                if line2.startswith(target_zpe2):
                    found_ZPE2 = True
                else:
                    line2 = f2.readline()
            ZPE2_b = re.findall(r"\S+", line2)
            # print(line2)
            # print('the length of line2 : ', len(ZPE2_b))
            # print(ZPE2_b)
            ZPE2 = float(ZPE2_b[6])
            line2 = f2.readline()

            E_hartree2 = ZPE2
            E_kelvin2 = E_hartree2 * na * h_to_j / 1000 * 1000 / R

        with open(sys.argv[3], "r") as f3:
            f_name_3 = (f3.name).replace('.out', '')

            # get SCF energy
            SD3 = 0
            line3 = f3.readline()
            target_zpe3 = ' SCF Done'
            found_SD3 = False
            while found_SD3 != True:
                if line3.startswith(target_zpe3):
                    found_SD3 = True
                else:
                    line3 = f3.readline()
            SD3_b = re.findall(r"\S+", line3)
            SD3 = float(SD3_b[4])
            # print(SD3_b)

            # get Mulliken charges
            MC3 = []
            target_mc3 = ' Mulliken charges and spin densities:'
            found_MC3 = False
            while found_MC3 != True:
                if line3.startswith(target_mc3):
                    found_MC3 = True
                    line3 = f3.readline()
                else:
                    line3 = f3.readline()
            line3 = f3.readline()
            for i in range(0, N-60):
                MC3_b = re.findall(r"\S+", line3)
                MC3.append(float(MC3_b[2]))
                line3 = f3.readline()
            """
            print('Mulliken Charges')
            for i in range(0, N-60):
                print(str(MC3[i]))
            """

            # get zero-point energies
            ZPE3 = 0
            line3 = f3.readline()
            target_zpe3 = ' Sum of electronic and zero-point Energies='
            found_ZPE3 = False
            while found_ZPE3 != True:
                if line3.startswith(target_zpe3):
                    found_ZPE3 = True
                else:
                    line3 = f3.readline()
            ZPE3_b = re.findall(r"\S+", line3)
            # print(line3)
            # print('the length of line3 : ', len(ZPE3_b))
            # print(ZPE3_b)
            ZPE3 = float(ZPE3_b[6])
            line3 = f3.readline()

            E_hartree3 = ZPE3
            E_kelvin3 = E_hartree3 * h_to_j * na / 1000 * 1000 / R

        # get adsorption energies
        E_ads_hartree = E_hartree2 + E_hartree3 - E_hartree1
        E_ads_kelvin = E_kelvin2 + E_kelvin3 - E_kelvin1
        print(E_ads_hartree)
        print(E_ads_kelvin)

        # get the distances of adsorption and decide the kind of adsorption (chemisorprion of physisorption)
        Dist = []
        Dist = measure_distance(N, info)

        # get the difference of Mulliken charges between before and after
        Dif = []
        for i in range(0, N-60):
            Dif.append(MC1[Dist[i][1]-1]-MC2[Dist[i][1]-1])
            Dif.append(MC1[Dist[i][3]-1]-MC2[Dist[i][3]-1])
            Dif.append(MC1[60+i]-MC3[i])

        with open("asw_zpe.csv", "a") as f4:
            f4.write(str(f_name).ljust(20))
            f4.write(str(round(E_ads_hartree, 10)).ljust(20))
            f4.write(str(round(E_ads_kelvin, 10)).ljust(20))
            for i in range(0, N-60):
                f4.write(str(round(Dist[i][0], 10)).ljust(20))
                f4.write(str(Dist[i][1]).ljust(20))
                f4.write(str(round(Dist[i][2], 10)).ljust(20))
                f4.write(str(Dist[i][3]).ljust(20))
                f4.write(str(round(Dif[3*i], 10)).ljust(20))
                f4.write(str(round(Dif[3*i+1], 10)).ljust(20))
                f4.write(str(round(Dif[3*i+2], 10)).ljust(20))
            f4.write('\n')   
    else:
            print('Error')


if __name__ == "__main__":
    main()
