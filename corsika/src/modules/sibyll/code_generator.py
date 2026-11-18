#!/usr/bin/env python3

# (c) Copyright 2018-2019 CORSIKA Project, corsika8@kit.edu
#
# See file AUTHORS for a list of contributors.
#
# This software is distributed under the terms of the 3-clause BSD license.
# See file LICENSE for a full version of the license.


import pickle, sys, itertools



def load_particledb(filename):
    '''
    loads the pickled particle_db (which is an OrderedDict)
    '''
    with open(filename, "rb") as f:
        particle_db = pickle.load(f)
    return particle_db



def read_sibyll_codes(filename, particle_db):
    '''
    reads to sibyll codes data file

    For particles known to sibyll, add 'sibyll_code' and 'sibyll_xsType' to particle_db
    '''
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if len(line)==0 or line[0] == '#':
                continue            
            identifier, sib_code, canInteractFlag, xsType = line.split()
            try:
                particle_db[identifier]["sibyll_code"] = int(sib_code)
                particle_db[identifier]["sibyll_xsType"] = xsType
                particle_db[identifier]["sibyll_canInteract"] = int(canInteractFlag)
            except KeyError as e:
                raise Exception("Identifier '{:s}' not found in particle_db".format(identifier))


            

def generate_sibyll_enum(particle_db):
    '''
     generates the enum to access sibyll particles by readable names
    '''
    output = "enum class SibyllCode : int8_t {\n"
    for identifier, pData in particle_db.items():
        if 'sibyll_code' in pData:
            output += "  {:s} = {:d},\n".format(identifier, pData['sibyll_code'])
    output += "};\n"
    return output

def generate_sibyll_can_interact(particle_db):
    '''
     generates the look-up table whether particles can interact in SIBYLL or not
    '''    
    string = "std::array<bool, {:d}> constexpr caninteract = {{\n".format(len(particle_db))    

    count = -1    
    for identifier, pData in particle_db.items():        
        count += 1
        if 'sibyll_canInteract' in pData:
            if pData['sibyll_canInteract'] != 0:
                string += "  true, // {:s} {:d}\n".format(identifier, count)
            else:
                string += "  false, // {:s} {:d}\n".format(identifier, count)
        else:
            string += "  false, // {:s} not implemented in SIBYLL \n".format(identifier)
    string += "};\n"
    return string


def generate_corsika2sibyll(particle_db):    
    '''
    generates the look-up table to convert corsika codes to sibyll codes
    '''
    string = "std::array<SibyllCode, {:d}> constexpr corsika2sibyll = {{\n".format(len(particle_db))
    for identifier, pData in particle_db.items():
        if pData['isNucleus']: continue
        if 'sibyll_code' in pData:
            string += "  SibyllCode::{:s}, \n".format(identifier)
        else:
            string += "  SibyllCode::Unknown, // {:s}\n".format(identifier + ' not implemented in SIBYLL')
    string += "};\n"
    return string
    


def generate_corsika2sibyll_xsType(particle_db):    
    '''
    generates the look-up table to convert corsika codes to sibyll codes
    '''
    string = "std::array<SibyllXSClass, {:d}> constexpr corsika2sibyllXStype = {{\n".format(len(particle_db))
    for identifier, pData in particle_db.items():
        if pData['isNucleus']: continue
        if 'sibyll_xsType' in pData:
            string += "  SibyllXSClass::{:s}, // {:s}\n".format(pData['sibyll_xsType'], identifier)
        else:
            string += "  SibyllXSClass::CrossSectionUnknown, // {:s}\n".format(identifier + ' not implemented in SIBYLL')
    string += "};\n"
    return string


def generate_sibyll2corsika(particle_db) :
    '''
    generates the look-up table to convert sibyll codes to corsika codes    
    '''
    string = ""
    
    minID = 0
    for identifier, pData in particle_db.items() :
        if 'sibyll_code' in pData:
            minID = min(minID, pData['sibyll_code'])

    string += "SibyllCodeIntType constexpr minSibyll = {:d};\n\n".format(minID)

    pDict = {}
    for identifier, pData in particle_db.items() :
        if 'sibyll_code' in pData:
            sib_code = pData['sibyll_code'] - minID
            pDict[sib_code] = identifier
    
    nPart = max(pDict.keys()) - min(pDict.keys()) + 1
    string += "std::array<corsika::Code, {:d}> constexpr sibyll2corsika = {{\n".format(nPart)
    
    for iPart in range(nPart) :
        if iPart in pDict:
            identifier = pDict[iPart]
        else:
            identifier = "Unknown"
        string += "  corsika::Code::{:s}, \n".format(identifier)
    
    string += "};\n"
    return string

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {:s} <particle_db.pkl> <sibyll_codes.dat>".format(sys.argv[0]), file=sys.stderr)
        sys.exit(1)
        
    print("code_generator.py for SIBYLL")
    
    particle_db = load_particledb(sys.argv[1])
    read_sibyll_codes(sys.argv[2], particle_db)
    
    with open("Generated.inc", "w") as f:
        print("// this file is automatically generated\n// edit at your own risk!\n", file=f)
        print(generate_sibyll_enum(particle_db), file=f)
        print(generate_sibyll_can_interact(particle_db), file=f)
        print(generate_corsika2sibyll(particle_db), file=f)
        print(generate_sibyll2corsika(particle_db), file=f)
        print(generate_corsika2sibyll_xsType(particle_db), file=f)
