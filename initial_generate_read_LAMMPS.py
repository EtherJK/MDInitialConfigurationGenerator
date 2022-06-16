"""Generate the .in file to locate each atom for LAMMPS"""
import random
import math
import copy

#basic parameters
class Box():
    """A normal simulation box"""
    def __init__(self, box_value, bondary_type):
        self.x_length = box_value[0]
        self.y_length = box_value[1]
        self.z_length = box_value[2]
        self.bondary_type = bondary_type

class MoleculeStruct():
    """Molecule structure"""
    def __init__(self):
        self.atom_list = []
        self.bond_list = []
        self.angle_list = []
        self.dihedral_list = []
        self.improper_list = []
        self.bond_length = []

    def generate_bond_length(self, bond_length_list):
        for i in range(0, len(self.bond_list)):
            self.bond_length.append(bond_length_list[self.bond_list[i]])


class EachMolecule():
    """Actual molecule we are going to use"""
    def __init__(self, molecule_struct):
        self.struct = molecule_struct

    def generate_chain(self, box, charge_list, bond_length_list):
        """generate the atoms for a whole chain"""
        self.struct.generate_bond_length(bond_length_list)
        self.atoms = []
        for now_atom in range(0, len(self.struct.atom_list)):
            now_atom_type = self.struct.atom_list[now_atom]
            self.atoms.append(EachAtom(now_atom_type, charge_list[now_atom_type]))
            if now_atom == 0:
                self.atoms[now_atom].gen_origin_pos(box)
            else:
                self.atoms[now_atom].gen_process_pos(box, self.atoms[now_atom - 1], self.struct.bond_length[now_atom - 1])


class EachAtom():
    """Atoms in each molecule"""
    num_tot_atom = 0

    def __init__(self, atom_type, charge):
        self.atom_type = atom_type
        self.charge = charge
        EachAtom.num_tot_atom += 1
        self.num_atom = EachAtom.num_tot_atom

    def gen_origin_pos(self, box):
        """if this atom is the first atom for a chain, generate the position"""
        self.position = Position( random.uniform(-box.x_length/2, box.x_length/2),
                                  random.uniform(-box.y_length/2, box.y_length/2),
                                  random.uniform(-box.z_length/2, box.z_length/2),
                                  )

    def gen_process_pos(self, box, pre_atom, bond_length):
        """if this atom is NOT the first atom for a chain, generate the position"""
        _theta = math.pi * random.uniform(0, 1)
        _phi = 2 * math.pi * random.uniform(0, 1)
        _x_change = bond_length * math.sin(_theta) * math.cos(_phi)
        _y_change = bond_length * math.sin(_theta) * math.sin(_phi)
        _z_change = bond_length * math.cos(_theta)
        self.position = copy.copy(pre_atom.position)
        self.position.position_change(box, _x_change, _y_change, _z_change)

class Position():
    """position for each molecules"""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def position_change(self, box, x_change, y_change, z_change):
        """move position to a bond length away"""
        if box.bondary_type == "periodic":
            self.x = pos_change_in_box(box.x_length, self.x, x_change)
            self.y = pos_change_in_box(box.y_length, self.y, y_change)
            self.z = pos_change_in_box(box.z_length, self.z, z_change)
        else:
            print("Invalid Bondary Type ERROR")

def pos_change_in_box(box_length, x_origin, x_change):
    """change the position and keep it in box"""
    if (x_origin + x_change) > box_length/2:
        while x_origin + x_change > box_length/2:
            x_change -= box_length
    elif (x_origin + x_change) < -box_length/2:
        while x_origin + x_change < -box_length/2:
            x_change += box_length

    if ((x_origin + x_change) > -box_length/2) and ((x_origin + x_change) < box_length/2):
        #print (x_origin + x_change)
        return x_origin + x_change
    else:
        print("Change Position ERROR")
        return None

def read_value(lines, key_words, return_type = "defalt", number_type = "int"):
    """extract the value of key word from in file"""
    for line in lines:
        if line.startswith(key_words):
            words = line.split()
            if (len(words) == 2) and (not return_type == "list"):
                if number_type == "int":
                    return int(words[1])
                elif number_type == "float":
                    return float(words[1])
                elif number_type == "string":
                    return words[1]
            else:
                if number_type == "int":
                    return [float(word) for word in words[1:]]
                elif number_type == "float":
                    return [float(word) for word in words[1:]]
                elif number_type == "string":
                    return words[1:]

def split_mol_info(lines, mol_info, types):
    """there are infos for different molecules, take them apart into each mol_info[]"""
    lines = ["".join(line.split()) for line in lines]
    for type in range(0, types):
        try:
            now_line = lines.index("molecule" + str(type + 1))
        except ValueError:
            print("Molecule Index ERROR")
            return None
        else:
            now_line += 1
            while (now_line < len(lines)) and (not lines[now_line].startswith("molecule")):
                mol_info[type].append(lines[now_line])
                now_line += 1


def read_mol_info(mol_info, now_mol_type, num_mol, molecule_type):
    """use the info for each molecule to generate its parameter"""
    for line in mol_info:
        if line.startswith("number"):
            num_mol[now_mol_type] = int(line.strip().replace("number",""))
            break
        print("Not found 'number' ERROR")
        return None
    molecule_type[now_mol_type].atom_list = generate_list(mol_info, "atom")
    molecule_type[now_mol_type].bond_list = generate_list(mol_info, "bond")
    molecule_type[now_mol_type].angle_list = generate_list(mol_info, "angle")
    molecule_type[now_mol_type].dihedral_list = generate_list(mol_info, "dihedral")
    molecule_type[now_mol_type].improper_list = generate_list(mol_info, "improper")
    #check if the number of each is proper
    if ( len(molecule_type[now_mol_type].atom_list) - 1 != len(molecule_type[now_mol_type].bond_list) 
        and len(molecule_type[now_mol_type].bond_list) != 0) or\
        ( len(molecule_type[now_mol_type].bond_list) - 1 != len(molecule_type[now_mol_type].angle_list) 
        and len(molecule_type[now_mol_type].angle_list) != 0) or\
        ( len(molecule_type[now_mol_type].angle_list) - 1 != len(molecule_type[now_mol_type].dihedral_list) 
        and len(molecule_type[now_mol_type].dihedral_list) != 0):
        print("numbers of items not match ERROR in " + str(now_mol_type + 1) +" type")

def generate_list(lines, key_words):
    """to process the sentences like 'atom3*2,2,2*1'"""
    target_list = []
    for line in lines:
        if line.startswith(key_words):
            target_line = line.strip().replace(key_words,"")
            target_line = expand_brakets(target_line)


            parts = target_line.split(",")

            for part in parts:
                if part == '':
                    continue
                if "*" in part:
                    [num_repeat, type_repeat] = [int(i) for i in part.split("*")]
                    for _ in range (0, num_repeat):
                        target_list.append(type_repeat)
                else:
                    target_list.append(int(part))
    return target_list

def expand_brakets(line):
    res_line = ''
    now_pos = 0
    while now_pos < len(line):
        if line[now_pos] == ')':
            for prev_pos in range( -1, -len(res_line) - 1, -1):
                if res_line[prev_pos] =='(':
                    need_to_expand = res_line[prev_pos + 1:]
                    break
            
            for num_end in range(now_pos+1, len(line) + 1):
                if num_end == len(line) or not line[num_end].isdigit():
                    repeat_num = int(line[now_pos+1 : num_end])
                    break

            now_pos = num_end - 1

            for i in range(0, repeat_num - 1):
                res_line = res_line + ',' +need_to_expand
            res_line = rreplace(res_line, '(', '', 1)

        else:
            res_line = res_line + line[now_pos]

        now_pos += 1
    return res_line

def rreplace(input, old, new, *max):
    count = len(input)
    if max and str(max[0]).isdigit():
        count = max[0]
    while count:
        index = input.rfind(old)
        if index >= 0:
            chunk = input.rpartition(old)
            input = chunk[0] + new + chunk[2]
        count -= 1
    return input

        


def output_head(file, num_tot_atom, num_tot_bond, num_tot_angle, 
                num_tot_dihedral, num_tot_improper, molecule_type, box):
    """output the head part for LAMMPS in file"""
    file.write("This file generated polymers for LAMMPS\n")
    file.write("\n")
    file.write(str(num_tot_atom) + "\tatoms\n")
    file.write(str(num_tot_bond) + "\tbonds\n")
    file.write(str(num_tot_angle) + "\tangles\n")
    file.write(str(num_tot_dihedral) + "\tdihedrals\n")
    file.write(str(num_tot_improper) + "\timpropers\n")
    file.write("\n")
    tot_atom_list = []
    tot_bond_list = []
    tot_angle_list = []
    tot_dihedral_list = []
    tot_improper_list = []
    for i in range(0, len(molecule_type)):
        tot_atom_list += molecule_type[i].atom_list
        tot_bond_list += molecule_type[i].bond_list
        tot_angle_list += molecule_type[i].angle_list
        tot_dihedral_list += molecule_type[i].dihedral_list
        tot_improper_list += molecule_type[i].improper_list
    tot_atom_list = list(set(tot_atom_list))
    tot_bond_list = list(set(tot_bond_list))
    tot_angle_list = list(set(tot_angle_list))
    tot_dihedral_list = list(set(tot_dihedral_list))
    tot_improper_list = list(set(tot_improper_list))
    file.write(str(len(tot_atom_list)) + "\tatom types\n")
    file.write(str(len(tot_bond_list)) + "\tbond types\n")
    file.write(str(len(tot_angle_list)) + "\tangle types\n")
    file.write(str(len(tot_dihedral_list)) + "\tdihedral types\n")
    file.write(str(len(tot_improper_list)) + "\timproper types\n")
    file.write("\n")
    file.write(str(-box.x_length/2) + " " + str(box.x_length/2) + " xlo xhi\n")
    file.write(str(-box.y_length/2) + " " + str(box.y_length/2) + " ylo yhi\n")
    file.write(str(-box.z_length/2) + " " + str(box.z_length/2) + " zlo zhi\n")
    file.write("\n")
    file.write("Masses\n")
    file.write("\n")
    for i in range(0, len(tot_atom_list)):
        file.write(str(i + 1) + "\t 1\n")
    file.write("\n")

def output_atom(file, molecules):
    """output the Atoms part"""
    file.write("Atoms\n")
    file.write("\n")
    for now_mol in range(0, len(molecules)):
        for now_atom in range(0, len(molecules[now_mol].struct.atom_list)):
            file.write(str(molecules[now_mol].atoms[now_atom].num_atom) 
                        + "\t" + str(now_mol + 1) 
                        + "\t" + str(molecules[now_mol].struct.atom_list[now_atom] + 1)
                        + "\t" + str(molecules[now_mol].atoms[now_atom].charge) + "\t"
                        + "\t" + str(format(molecules[now_mol].atoms[now_atom].position.x, ".4f"))
                        + "\t" + str(format(molecules[now_mol].atoms[now_atom].position.y, ".4f"))
                        + "\t" + str(format(molecules[now_mol].atoms[now_atom].position.z, ".4f")) + "\n" )
    file.write("\n")

def output_bond(file, molecules):
    """output the bond part"""
    file.write("Bonds\n")
    file.write("\n")
    num_bond = 1
    for now_mol in range(0, len(molecules)):
        for now_atom in range(0, len(molecules[now_mol].struct.bond_list)):
            file.write(str(num_bond) 
                        + "\t" + str(molecules[now_mol].struct.bond_list[now_atom] + 1) 
                        + "\t" + str(molecules[now_mol].atoms[now_atom].num_atom)
                        + "\t" + str(molecules[now_mol].atoms[now_atom].num_atom + 1) + "\n")
            num_bond += 1
    file.write("\n")

def output_angle(file, molecules):
    """output the angle part"""
    file.write("Angles\n")
    file.write("\n")
    num_angle = 1
    for now_mol in range(0, len(molecules)):
        for now_atom in range(0, len(molecules[now_mol].struct.angle_list)):
            file.write(str(num_angle) 
                        + "\t" + str(molecules[now_mol].struct.angle_list[now_atom] + 1) 
                        + "\t" + str(molecules[now_mol].atoms[now_atom].num_atom)
                        + "\t" + str(molecules[now_mol].atoms[now_atom].num_atom + 1)
                        + "\t" + str(molecules[now_mol].atoms[now_atom].num_atom + 2) + "\n")
            num_angle += 1

#read file
filename_in = "polymer_info.txt"
with open(filename_in) as input_file:
    lines = input_file.readlines()
for line in lines:
    if (line.startswith("#")) or (not line.strip()):
        lines.remove(line)

box_value = read_value(lines, "box", "list", "float")
bondary_type = read_value(lines, "bondary_type", "defalt", "string")
box = Box(box_value, bondary_type)
num_molecule_type = read_value(lines, "molecule_type")
num_atom_type = read_value(lines, "atom_type")
atom_charge = read_value(lines, "atom_charge", "list")
num_bond_type = read_value(lines, "bond_type")
bond_length = read_value(lines, "bond_length", "list", "float")
angle_type = read_value(lines, "angle_type")

#generate the molecule structure
num_mol = [0 for _ in range(0,num_molecule_type)]
molecule_type = []
for _ in range(0, num_molecule_type):
    molecule_type.append(MoleculeStruct())

#read the molecule info and put it in molecule_type
mol_info = [[] for i in range(0, num_molecule_type)]
split_mol_info(lines, mol_info, num_molecule_type)
for now_mol_type in range(0, num_molecule_type):
    read_mol_info(mol_info[now_mol_type][:], now_mol_type, num_mol, molecule_type)

#count for the output_head
num_tot_mol = 0
num_tot_atom = 0
num_tot_bond = 0
num_tot_angle = 0
num_tot_dihedral = 0
num_tot_improper = 0
molecules = []
for now_type in range(0, num_molecule_type):
    num_tot_atom += num_mol[now_type] * len(molecule_type[now_type].atom_list)
    num_tot_bond += num_mol[now_type] * len(molecule_type[now_type].bond_list)
    num_tot_angle += num_mol[now_type] * len(molecule_type[now_type].angle_list)
    num_tot_dihedral += num_mol[now_type] * len(molecule_type[now_type].dihedral_list)
    num_tot_improper += num_mol[now_type] * len(molecule_type[now_type].improper_list)
    for now_mol in range(0, num_mol[now_type]):
        molecules.append( EachMolecule( molecule_type[now_type] ) )
        molecules[num_tot_mol].generate_chain(box, atom_charge, bond_length)
        num_tot_mol += 1

#begin output
filename_out = "polymer.in"
output_file = open(filename_out, 'w')
output_head(output_file, num_tot_atom, num_tot_bond, num_tot_angle, 
            num_tot_dihedral, num_tot_improper, molecule_type, box)
output_atom(output_file, molecules)
output_bond(output_file, molecules)
if num_tot_angle > 0:
    output_angle(output_file, molecules)
#output_dihedral(output_file, molecules)
#output_improper(output_file, molecules)
output_file.close()