# Comparative modeling by the automodel class
#
# Demonstrates how to build multi-chain models
#
import sys
sys.path.append('/vdb/User/wuxiaoqin/TCR-pMHC_Modelling/modeller-10.5/modlib/modeller')
from modeller import *
from modeller.automodel import *    # Load the automodel class
#from modeller import soap_protein_od

log.verbose()

class MyModel(automodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['A', 'B', 'C', 'D', 'E'],
                             renumber_residues=[1, 1, 1, 1, 1])
        # Another way to label individual chains:
        # self.chains[0].name = 'A'
        # self.chains[1].name = 'B'

env = Environ()
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

aln = Alignment(env)
code =[]
input_argvs = input("alignment sequence file(.ali):")
aln_file = input_argvs.split()[0]
model_num = input_argvs.split()[1]
#aln_file = input("alignment sequence file(.ali):")
read_aln = modfile.File(aln_file, 'r')
while aln.read_one(read_aln, alignment_format='PIR'):
	code.append(aln[0].code)
	
read_aln.close()
template_code = code[0]
target_code = code[1]

# Be sure to use 'MyModel' rather than 'automodel' here!
a = MyModel(env,
            alnfile  = aln_file ,     # alignment filename
            knowns   = template_code, # codes of the templates
            sequence = target_code,   # code of the target
            assess_methods = (assess.DOPE,
                              assess.GA341,
		              #soap_protein_od.Scorer()
			      )
	    )                       

a.starting_model= 1                # index of the first model
a.ending_model  = int(model_num)                # index of the last model
                                   # (determines how many models to calculate)
a.make()                           # do comparative modeling
