"""
Generate the FIELD file to locate each atom for DL_Meso
Although many infos in polymer.txt is useless in generating DL_Meso FIELD file,
in order to keep consistance with LAMMPS .in generate file, still keeps them.
"""
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

    def generate_chain(self, charge_list, bond_length_list):
        """generate the atoms for a whole chain"""
        self.struct.generate_bond_length(bond_length_list)
        self.atoms = []
        for now_atom in range(0, len(self.struct.atom_list)):
            now_atom_type = self.struct.atom_list[now_atom]
            self.atoms.append(EachAtom(now_atom_type, charge_list[now_atom_type]))
            if now_atom == 0:
                self.atoms[now_atom].gen_origin_pos()
            else:
                self.atoms[now_atom].gen_process_pos(self.atoms[now_atom - 1], self.struct.bond_length[now_atom - 1])


class EachAtom():
    """Atoms in each molecule"""
    num_tot_atom = 0

    def __init__(self, atom_type, charge):
        self.atom_type = atom_type
        self.charge = charge
        EachAtom.num_tot_atom += 1
        self.num_atom = EachAtom.num_tot_atom

    def gen_origin_pos(self):
        """if this atom is the first atom for a chain, generate the position"""
        self.position = Position(0, 0, 0)

    def gen_process_pos(self, pre_atom, bond_length):
        """if this atom is NOT the first atom for a chain, generate the position"""
        _theta = math.pi * random.uniform(0, 1)
        _phi = 2 * math.pi * random.uniform(0, 1)
        _x_change = bond_length * math.sin(_theta) * math.cos(_phi)
        _y_change = bond_length * math.sin(_theta) * math.sin(_phi)
        _z_change = bond_length * math.cos(_theta)
        self.position = copy.copy(pre_atom.position)
        self.position.position_change(_x_change, _y_change, _z_change)

class Position():
    """position for each molecules"""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def position_change(self, x_change, y_change, z_change):
        """move position to a bond length away"""
        self.x = pos_change_in_box(self.x, x_change)
        self.y = pos_change_in_box(self.y, y_change)
        self.z = pos_change_in_box(self.z, z_change)


def pos_change_in_box(x_origin, x_change):
    """change the position and keep it in box"""
    return x_origin + x_change

def read_value(lines, key_words, return_type = "defalt", number_type = "int", return_line_number = False):
    """extract the value of key word from in file"""
    n = 0
    for line in lines:
        if line.startswith(key_words):
            words = line.split()
            if (len(words) == 2) and (not return_type == "list"):
                if number_type == "int":
                    result = int(words[1])
                elif number_type == "float":
                    result = float(words[1])
                elif number_type == "string":
                    result = words[1]
            else:
                if number_type == "int":
                    result = [float(word) for word in words[1:]]
                elif number_type == "float":
                    result = [float(word) for word in words[1:]]
                elif number_type == "string":
                    result = words[1:]
            if return_line_number == False:
                return result
            else:
                return [result, n]
        n += 1

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

def find_unbond(molecule_type, num_mol, num_tot_atom):
    unbond_atom = [0 for _ in range(0, num_tot_atom)]
    for mol_type in range(0, len(molecule_type)):
        if len(molecule_type[mol_type].atom_list) == 1:
            now_atom_type = molecule_type[mol_type].atom_list[0]
            unbond_atom[now_atom_type] += num_mol[mol_type]
    return unbond_atom

def read_interaction(line, num_atom_type):
    # useage: split the line into A_types, B_types and interaction_coef
    parts = line.split()
    A_types = parts[0]
    B_types = parts[1]
    def convert_words_into_types(types, num_atom_type):
        if types == "*":
            return range(0, num_atom_type)
        else:
            return [int(types)]
    A_types = convert_words_into_types(A_types, num_atom_type)
    B_types = convert_words_into_types(B_types, num_atom_type)
    coef = ""
    for n in range(2, len(parts)):
        coef += parts[n]
        coef += " "
    return [A_types, B_types, coef]

    





def output_head(file):
    file.write("Project Name\n\n")

def output_species(file, unbond_atom, molecule_type):
    tot_atom_list = []
    for i in range(0, len(molecule_type)):
        tot_atom_list += molecule_type[i].atom_list
    tot_atom_list = list(set(tot_atom_list))
    file.write("SPECIES " + str(len(tot_atom_list)) + "\n")
    for atom_type in range(0, len(tot_atom_list)):
        file.write(chr(ord("A") + atom_type) + " 1.0 " + str(atom_charge[atom_type]) + " " + str(unbond_atom[atom_type]) + "\n")
    file.write("\n")

def output_molecules(file, molecule_type, num_mol, molecules, bond_force, bond_coef, bond_length):
    n_molecule = 0
    for now_mol in range(0, len(molecule_type)):
        if len(molecule_type[now_mol].atom_list) > 1:
            n_molecule += 1
    if n_molecule > 0:
        file.write("MOLECULES " + str(n_molecule) + "\n")
    for now_mol in range(0, len(molecule_type)):
        if len(molecule_type[now_mol].atom_list) > 1:
            file.write("MOL_" + chr(ord("A") + now_mol) + "\n")
            file.write("nummols " + str(num_mol[now_mol]) + "\n")
            file.write("beads " + str(len(molecule_type[now_mol].atom_list)) + "\n")
        
            for now_atom in range(0, len(molecules[now_mol].struct.atom_list)):
                file.write(chr(molecules[now_mol].struct.atom_list[now_atom] + ord("A"))
                            + "\t" + str(format(molecules[now_mol].atoms[now_atom].position.x, ".4f"))
                            + "\t" + str(format(molecules[now_mol].atoms[now_atom].position.y, ".4f"))
                            + "\t" + str(format(molecules[now_mol].atoms[now_atom].position.z, ".4f")) + "\n" )
            file.write("bonds " + str(len(molecule_type[now_mol].bond_list)) + "\n")
            for now_atom in range(0, len(molecules[now_mol].struct.atom_list) - 1):
                file.write(bond_force + " " + str(now_atom + 1) + " " + str(now_atom + 2) + " " + str(bond_coef) + " " + str(bond_length[0]) + "\n")
            file.write("finish\n")
    file.write("\n")

def output_interaction(file, num_atom_type, inter_matrix, interactions):
    file.write("INTERACTIONS " + str(int((num_atom_type + 1) * num_atom_type / 2)) + "\n")
    for A_type in range(0, num_atom_type):
        for B_type in range(A_type, num_atom_type):
            file.write(chr(ord("A") + A_type) + " " + chr(ord("A") + B_type) + " " + interactions + " " + inter_matrix[A_type][B_type] + "\n")
            

#read file
filename_in = "polymer_info.txt"
with open(filename_in) as input_file:
    lines = input_file.readlines()
for line in lines:
    if (line.startswith("#")) or (not line.strip()):
        lines.remove(line)

num_molecule_type = read_value(lines, "molecule_type")
num_atom_type = read_value(lines, "atom_type")
atom_charge = read_value(lines, "atom_charge", "list")
num_bond_type = read_value(lines, "bond_type")
bond_length = read_value(lines, "bond_length", "list", "float")
angle_type = read_value(lines, "angle_type")
bond_force = read_value(lines, "bond_force", number_type="string")
bond_coef = read_value(lines, "bond_coef", number_type="float")

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
    for now_mol in range(0, 1):
        molecules.append( EachMolecule( molecule_type[now_type] ) )
        molecules[num_tot_mol].generate_chain(atom_charge, bond_length)
        num_tot_mol += 1

#generate the interaction matrix
[interactions, num_line] = read_value(lines, "interactions", number_type="string", return_line_number=True)
inter_matrix = [["" for i in range(0, num_atom_type)] for j in range(0, num_atom_type)]
while True:
    num_line += 1
    [A_types, B_types, interaction_coef] = read_interaction(lines[num_line], num_atom_type)
    for A_type in A_types:
        for B_type in B_types:
            inter_matrix[A_type][B_type] = interaction_coef
            inter_matrix[B_type][A_type] = interaction_coef
    next_line = lines[num_line + 1]
    if not (next_line.startswith("*") or next_line[0].isdigit()):
        break

#begin output
filename_out = "FIELD"
output_file = open(filename_out, 'w')
output_head(output_file)
unbond_atom = find_unbond(molecule_type, num_mol, num_tot_atom)
output_species(output_file, unbond_atom, molecule_type)
output_molecules(output_file, molecule_type, num_mol, molecules, bond_force, bond_coef, bond_length)
output_interaction(output_file, num_atom_type, inter_matrix, interactions)
output_file.write("\nCLOSE")


output_file.close()