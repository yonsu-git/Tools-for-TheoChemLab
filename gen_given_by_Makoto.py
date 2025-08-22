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
def get_m_and_l(total_atoms, nuclear_charges):
    ATOMIC_WEIGHT = []
    LABEL = []
    for i in range(0, total_atoms):
        if nuclear_charges[i] == 1:
            ATOMIC_WEIGHT.append(float(1.0078250))
            LABEL.append('H')
        elif nuclear_charges[i] == 8:
            ATOMIC_WEIGHT.append(float(15.9949146))
            LABEL.append('O')
        elif nuclear_charges[i] == 6:
            ATOMIC_WEIGHT.append(float(12.0000000))
            LABEL.append('C')
        elif nuclear_charges[i] == 7:
            ATOMIC_WEIGHT.append(float(14.00307))
            LABEL.append('N')
        elif nuclear_charges[i] == 15:
            ATOMIC_WEIGHT.append(float(30.973762))
            LABEL.append('P')

    return ATOMIC_WEIGHT, LABEL

# define a fanction named "gen_rad"
def gen_rad(N, N2, info, info_ten):
    # get a gravity center of ASW cluster (G_ASW)
    mt_w = 0; gx_w = 0; gy_w = 0; gz_w = 0
    for i in range(0, N):
        mt_w += info[i].atom[0]
        gx_w += info[i].atom[1]*info[i].atom[0]
        gy_w += info[i].atom[2]*info[i].atom[0]
        gz_w += info[i].atom[3]*info[i].atom[0]
    gx_w /= mt_w
    gy_w /= mt_w
    gz_w /= mt_w
    # print('G_ASW : x= ', gx_w, ', y= ', gy_w, ', z= ', gz_w)

    # get a gravity center of radical (G_rad)
    mt_r = 0; gx_r = 0; gy_r = 0; gz_r = 0
    for i in range(0, N2):
        mt_r += info_ten[i].atom[0]
        gx_r += info_ten[i].atom[1]*info_ten[i].atom[0]
        gy_r += info_ten[i].atom[2]*info_ten[i].atom[0]
        gz_r += info_ten[i].atom[3]*info_ten[i].atom[0]
    gx_r /= mt_r
    gy_r /= mt_r
    gz_r /= mt_r
    # print('G_rad : x= ', gx_r, ', y= ', gy_r, ', z= ', gz_r)

    # set the distance between G_ASW and G_rad [angstrom]
    r = 1
    r_cut = 3

    # the degrees between G_ASW and G_rad
    t1 = random.uniform(0, math.pi)
    p1 = random.uniform(0, math.pi*2)
    print('t1 : ', t1, ', p1 : ', p1)

    # the degrees of rotation of the radical
    t2 = random.uniform(0, math.pi)
    p2 = random.uniform(0, math.pi*2)

    theta_x = random.uniform(0, math.pi*2)
    theta_y = random.uniform(0, math.pi*2)
    theta_z = random.uniform(0, math.pi*2)
    print('theta_x= ', theta_x, ', theta_y= ', theta_y, ', theta_z= ', theta_z)

    comp_Grad = False
    while comp_Grad != True:
        # print('r : ', r)
        # the coordinate of G_rad
        x_rad = gx_w + r*math.sin(t1)*math.cos(p1)
        y_rad = gy_w + r*math.sin(t1)*math.sin(p1)
        z_rad = gz_w + r*math.cos(t1)
        # rotate the radical
        C_rot = []
        for i in range(0, N2):
            # move an atom just for the distance
            x_ten = info_ten[i].atom[1] - gx_r
            y_ten = info_ten[i].atom[2] - gy_r
            z_ten = info_ten[i].atom[3] - gz_r
            # print('Atom : x= ', x_ten, ', y= ', y_ten, ', z= ', z_ten)
            # print(r_ten)
            # rotate around the x axis
            x_1 = x_ten
            y_1 = y_ten*math.cos(theta_x) + z_ten*math.sin(theta_x)
            z_1 = y_ten*(-math.sin(theta_x)) + z_ten*math.cos(theta_x)
            # rotate around the y axis
            x_2 = x_1*math.cos(theta_y) + z_1*(-math.sin(theta_y))
            y_2 = y_1
            z_2 = x_1*math.sin(theta_y) + z_1*math.cos(theta_y)
            # rotate around the z axis
            x_rot = x_2*math.cos(theta_z) + y_2*math.sin(theta_z)
            y_rot = x_2*(-math.sin(theta_z)) + y_2*math.cos(theta_z)
            z_rot = z_2 
            x_rot += x_rad
            y_rot += y_rad
            z_rot += z_rad
            C_rot.append(x_rot)
            C_rot.append(y_rot)
            C_rot.append(z_rot)
            # print('Atom : x= ', x_rot, ', y= ', y_rot, ', z= ', z_rot)
            x_ten = 0; y_ten = 0; z_ten = 0
            x_rot = 0; y_rot = 0; z_rot = 0
        # print('distance : ', math.sqrt((C_rot[3] - C_rot[6])**2 + (C_rot[4] - C_rot[7])**2 + (C_rot[5] - C_rot[8])**2))
        dis = 0
        coordinate_rejection = False
        for i in range(0, N):
            for j in range(0, N2):
                dis = math.sqrt((C_rot[3*j] - info[i].atom[1])**2 + (C_rot[3*j + 1] - info[i].atom[2])**2 + (C_rot[3*j + 2] - info[i].atom[3])**2)
                if dis < r_cut:
                    coordinate_rejection = True
        if coordinate_rejection == True:
            r += 0.1
            C_rot = []
        elif coordinate_rejection == False:
            comp_Grad = True
            print('r : ', r)
            print('\n')
    # create a 'Atom' class named 'info_gen'
    info_gen = []
    for i in range(0, N2):
        info_gen.append(Atom(info_ten[i].atom[0], C_rot[3*i], C_rot[3*i+1], C_rot[3*i+2], info_ten[i].atom[4]))
    return info_gen

# define a fanction named "make_csvfile"
def make_gjffile(total_atoms, f_name1, f_name2, info, charge, spin, N_file):
    f_name_input = f_name1.replace('.fchk', '_' + f_name2.replace('.fchk', '_') + 'Try' + str(N_file).zfill(2) + '.gjf')
    f_name_chk = f_name_input.replace('gjf', 'chk')
    with open(f_name_input, "w") as f:
        f.write('%nproc=4\n')
        f.write('%mem=24000mb\n')
        f.write(r'%chk=' + f_name_chk + '\n')
        f.write('#p wb97xd/6-311+(d,p) pop=nbo SCF(Xqc,MaxConv=60) Opt Freq Int(Grid=Ultra) GFinput NoSymm\n\n')
        f.write('RE:\n\n')
        f.write(str(charge) + ' ' + str(spin) + '\n')
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
            f_name1 = f.name

            # set the total of asw cluster
            N = 60
            if N % 5 == 0:
                n = int(N/5)
                # print('n = ', n)
            else:
                n = math.floor(N/5 + 1)
                # print('n = ', n)
            N_2 = N * 3
            if N_2 % 5 == 0:
                n_2 = int(N_2/5)
            else:
                n_2 = math.floor(N_2/5 + 1)
            # print('n_2 = ', n_2)

            # get the weight of atoms
            NC = []
            target1 = 'Nuclear charges'
            line = f.readline()
            # print('line = ', line)
            found_NC = False
            while found_NC != True:
                line = f.readline()
                if line.startswith(target1):
                    found_NC = True
            line = f.readline()
            for i in range(0, n):
                # print(line.rstrip("\n"))
                NC_b = re.findall(r"\S+", line)
                for i in range(0, len(NC_b)):
                    NC.append(float(NC_b[i]))
                line = f.readline()
            # print('The length of NC : ', len(NC))
            M = get_m_and_l(N, NC)[0]
            L = get_m_and_l(N, NC)[1]
            """
            for i in range(0, N):
                print('Atom', i+1, ' : AW=', M[i], ',L=', L[i])
            """

            # get last cartesian coordinates (lcc) and last velocity coordinates (lvc)
            C = []
            # print(line)
            target2 = 'Current cartesian coordinates'
            found_C = False
            while found_C != True:
                if line.startswith(target2):
                    found_C = True
                    line = f.readline()
                else:
                    line = f.readline()
            for i in range(0, n_2):
                # print(line.rstrip("\n"))
                C_b = re.findall(r"\S+", line)
                for i in range(0, len(C_b)):
                    C.append(float(C_b[i])*b_to_a)
                line = f.readline()
            # print('The length of C : ', len(C))
            """
            for i in range(0, N):
                print('atom : ', i+1)
                print('x=', C[3*i], ' y= ', C[3*i+1], ' z= ', C[3*i+2])
            """
            info = []
            for i in range(0, N):
                info.append(Atom(M[i], C[3*i], C[3*i+1], C[3*i+2], L[i]))

        with open(sys.argv[2], "r") as f:
            f_name2 = f.name
            target0 = 'Number of atoms'
            found_N_mol = False
            while found_N_mol != True:
                line = f.readline()
                if line.startswith(target0):
                    found_N_mol = True
            N2_b = re.findall(r"\S+", line)
            N2 = int(N2_b[4])
            # print(N2)
            # set the total of asw cluster
            if N2 % 5 == 0:
                n2 = int(N2/5)
                # print('n = ', n)
            else:
                n2 = math.floor(N2/5 + 1)
                # print('n = ', n)
            N2_2 = N2 * 3
            if N2_2 % 5 == 0:
                n2_2 = int(N2_2/5)
            else:
                n2_2 = math.floor(N2_2/5 + 1)
            # print('n2_2 = ', n2_2)

            # get the weight of atoms
            NC2 = []
            target1 = 'Nuclear charges'
            line = f.readline()
            # print('line = ', line)
            found_NC2 = False
            while found_NC2 != True:
                line = f.readline()
                if line.startswith(target1):
                    found_NC2 = True
            line = f.readline()
            for i in range(0, n2):
                # print(line.rstrip("\n"))
                NC2_b = re.findall(r"\S+", line)
                for i in range(0, len(NC2_b)):
                    NC2.append(float(NC2_b[i]))
                line = f.readline()
            # print('The length of NC : ', len(NC2))
            M2 = get_m_and_l(N2, NC2)[0]
            L2 = get_m_and_l(N2, NC2)[1]
            """
            for i in range(0, N):
                print('Atom', i+1, ' : AW=', M2[i], ',L=', L2[i])
            """

            # get last cartesian coordinates (lcc) and last velocity coordinates (lvc)
            C2 = []
            # print(line)
            target2 = 'Current cartesian coordinates'
            found_C2 = False
            while found_C2 != True:
                if line.startswith(target2):
                    found_C2 = True
                    line = f.readline()
                else:
                    line = f.readline()
            for i in range(0, n2_2):
                # print(line.rstrip("\n"))
                C_b2 = re.findall(r"\S+", line)
                for i in range(0, len(C_b2)):
                    C2.append(float(C_b2[i])*b_to_a)
                line = f.readline()
            # print('The length of C : ', len(C2))
            """
            for i in range(0, N2):
                print('atom : ', i+1)
                print('x=', C2[3*i], ' y= ', C2[3*i+1], ' z= ', C2[3*i+2])
            """
            info_ten = []
            for i in range(0, N2):
                info_ten.append(Atom(M2[i], C2[3*i], C2[3*i+1], C2[3*i+2], L2[i]))

        # set the size of files
        size_file = 10

        # set charge and spin
        charge = 1
        spin = 1

        info_rad = []
        for i in range(0, size_file):
            info_rad = gen_rad(N, N2, info, info_ten)
            for j in range(0, N2):
                info.append(info_rad[j])
            make_gjffile(N + N2, f_name1, f_name2, info, charge, spin, i+1)
            info_rad = []
            del info[N:]

    else:
        print('Error')


if __name__ == "__main__":
    main()
