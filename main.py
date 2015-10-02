# Classes for analysing dihedral angles from HOMSTRAD database entries

import Bio.PDB as pdb
import numpy as np
import re
import os
import itertools
from Bio import SeqIO

# TODO:
#   - Move the change of NoneType to np.NAN to getAlignedPhiPsi()
#   - Write parsing of pir/ali file in seperate function

class AlignStruc():
    '''Classes for each HOMSTRAD pdb'''
    def __init__(self, pdb_file_path, pir_file_path):
        self.parser = pdb.PDBParser(QUIET=True)
        self.ppb = pdb.PPBuilder()
        self.TriedFixed = False
        self.Working = False
        self.isExpanded = False

        # Load pir file
        # Needs own function
        handle = open(pir_file_path, "rU")
        strucParse = struc = SeqIO.parse(handle, "pir")
        self.pir_file_list = [(i.id, str(i.seq)) for i in struc]

        # Load pdb file
        self.family = self.parser.get_structure("Name", pdb_file_path)

    def __getitem__(self, key):
        if self.isExpanded == True:
            return self.subStrucList[self.subStrucDict[key]]
        else:
            return

    def printSeqs(self):
        for i in self.subStrucList:
            print i, "\n"

    def PirPdbMatch(self):
        ''' Checks if sequence of structures matches between pdb and pir file'''
        if (self.isMismatched() & self.isExpanded):
            print "Trying to fix"
            self.fixMismatches()

    def isMismatched(self):
        if (self.isExpanded):
            for i in self.subStrucList:
                if not (i.CheckConsistency()):
                    return True

        return False

    def isEqualLen(self):
        return len(self.pir_file_list) == len(self.family[0])

    def isWorking(self):
        return ((not self.isMismatched()) & self.isEqualLen())

    def fixMismatches(self):
        self.TriedFixed = True

        # Try to fix the order of the pir_list
        pirPermute = []
        for pep in self.subStrucList:
            for i, pir in enumerate(self.pir_file_list):
                pirSeq = pir[1].replace("\n", "").replace(" ", "").replace("-", "")
                if pirSeq == pep.getPdbSeq():
                    pirPermute.append(i)
        self.pir_file_list = [self.pir_file_list[i] for i in pirPermute]

        # Check if the fix worked
        if not self.isEqualLen():
            print "Couldn't fix"
            return

        self.GetContainedStructures()
        if not(self.isMismatched()):
            print "Success!"

    def GetContainedStructures(self):
        if not self.isEqualLen():
            print "Unequal lengths"
            return
        self.subStrucList = []
        self.subStrucDict = {}
        for strucId, structure in enumerate(self.family[0]):
            tempSubStruc = SingleStruc(structure, self.pir_file_list[strucId])
            self.subStrucList.append(tempSubStruc)
            self.subStrucDict[tempSubStruc.Name] = strucId
            self.isExpanded = True

    # NOT USED YET
    def UseIO(self, struc):
        self.io.set_structure(struc)
        self.io.save('temp.pdb')

    def parsePirFile(self):
        pass

    def getAlignedPhiPsi(self):
        phiPsiArray = []
        for i in self.subStrucList:
            phiPsiArray.append(i.getPhiPsi())
        return np.array(phiPsiArray)

    def getAlignedSeq(self):
        SeqArray = []
        for i in self.subStrucList:
            SeqArray.append(i.seq)
        return np.array(SeqArray)

    def getOutputCSV(self):
        PhiPsiAAArray = []
        for i in self.subStrucList:
            temp = np.hstack((i.Name, np.array(i.getPhiPsiAA()).flatten("F")))
            PhiPsiAAArray.append(temp)

        output = np.array(PhiPsiAAArray)
        header = ["%s_%s" % (j,i) for j in ("phi", "psi", "AA") for i in range(1, ((output.shape[1]-1)/3)+1)]
        header.insert(0, "PDB_ID")
        return np.vstack((header, PhiPsiAAArray))

class SingleStruc(AlignStruc):
    '''A single structure from an aligned HOMSTRAD pdb'''
    def __init__(self, structure, pir_seq):
        self.io = pdb.PDBIO()
        self.parser = pdb.PDBParser(QUIET=True)
        self.ppb = pdb.PPBuilder()

        # Figure out what kind of structure is send from AlignStruc
        self.struc = structure
        self.peptides = self.ppb.build_peptides(self.struc)

        # Name and seq from pir_seq
        self.Name = pir_seq[0]
        self.seq = pir_seq[1]

    def __str__(self):
        return "%s \n%s" % (self.seq, self.getPdbSeq())

    def CheckConsistency(self):
        '''Checks the alignment file against the pdb sequence'''
        pepSeq = self.getPdbSeq()
        isSame = pepSeq == self.seq.replace('-', '')
        return isSame

    def getPdbSeq(self):
        pepSeq = ''
        for pep in self.peptides:
            pepSeq += pep.get_sequence()
        return pepSeq

    def getSecStruc(self):
        '''Pipes the pdb structure to a perl script, that translates it to
        pdb'''
        # Make a temporary file to pipeline to the stride perl script
        # CHANGE TO USE ONLY ONE IO INSTANCE FROM PARENT CLASS
        self.io.set_structure(self.struc)
        self.io.save("temp.pdb")

        # Pipelining to perl
        os.system('stride temp.pdb -fTemp.str')
        os.system('perl stride2pdb temp.str > temp.ss')

        # FUNCTION NOT YET IMPLEMENTED
        pdb_header = self.listSecStruc('temp.ss')

        os.system('cat temp.ss temp.pdb > temp2.pdb')

        # Build peptide
        self.struc = self.parser.get_structure(self.Name, 'temp2.pdb')[0]
        self.peptides = self.ppb.build_peptides(self.struc)

        # Clean up temporary files
        os.remove('temp.ss')
        os.remove('temp.str')
        os.remove('temp.pdb')
        os.remove('temp2.pdb')

    def getPhiPsi(self):
        ''' Returns the aligned dihedral angles as a numpy array
        The angles are calculated by biopython'''
        alignedDihedral = np.empty([len(self.seq)], dtype=object)
        alignedDihedral[:] = np.NAN

        if (len(self.peptides) > 1):
            tempList = []
            offset = 0
            for i, pepchain in enumerate(self.peptides):

                # Compensates for the index offset in the pdb file
                if i > 0:
                #     # This loop needs to be checked for len(peptide) > 2 (more than one split in the protein)
                    Prev_endRes = self.peptides[i-1][-1].get_id()[1] + 1 - offset
                    offset = pepchain[0].get_id()[1] - Prev_endRes

                # flatten list
                phiPsi = pepchain.get_phi_psi_list()
                tempList.append([np.NAN]*offset)
                tempList.append(phiPsi)

            phiPsi = list(itertools.chain(*tempList))

        else:
            phiPsi = self.peptides[0].get_phi_psi_list()

        startPos = self.peptides[0][0].get_id()[1]-1
        count = 0
        for i, AA in enumerate(self.seq):
            if not (AA == '-'):
                if len(phiPsi) < count+1:
                    print count
                # Gave out of range error for "az" in Homstrad
                temp = phiPsi[count]
                alignedDihedral[i] = temp
                count += 1

        return alignedDihedral

    def listSecStruc(self, strucFile):
        '''Uses regex to find secondary structures in the file returned by
         stride2pdb'''
        # Import alignment
        with open(strucFile, "r") as myfile:
            align = myfile.read()

        temphelices = re.findall(r'(HELIX.*[^*]*?)', align)
        tempsheets = re.findall(r'(SHEET.*[^*]*?)', align)

        sheets = np.empty([0,3])
        helices = np.empty([0,3])

        if len(temphelices) > 0:
            helices = np.genfromtxt(temphelices, delimiter=None, usecols=(5, 8))
            # Id for expanding. 0 is coil, 1 is helix, 2 is sheet
            IDCol = np.empty((helices.shape[0], 1))
            IDCol[:] = 1
            helices = np.hstack((IDCol, helices))

        if len(tempsheets) > 0:
            sheets = np.genfromtxt(tempsheets, delimiter=None, usecols=(6, 9))
            # Id for expanding. 0 is coil, 1 is helix, 2 is sheet
            IDCol = np.empty((sheets.shape[0], 1))
            IDCol[:] = 2
            sheets = np.hstack((IDCol, sheets))

        return np.vstack((helices, sheets))

    def getPhiPsiAA(self):
        '''Returns a list of first dihedral angles, then amino acid for
        position'''
        tempList = []
        for i, dihedral in enumerate(self.getPhiPsi()):
            if type(dihedral) == tuple:
                if dihedral[0] is None:
                    tempList.append((("NA", dihedral[1]), self.seq[i]))
                elif dihedral[1] is None:
                    tempList.append(((dihedral[0], "NA"), self.seq[i]))
                else:
                    tempList.append((dihedral, self.seq[i]))
            else:
                tempList.append((("NA", "NA"), '-'))
        return [(a,b,c) for (a,b),c in tempList]
