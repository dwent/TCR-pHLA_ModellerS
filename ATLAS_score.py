#!/usr/bin/env python

import argparse
import os
import numpy as np
import subprocess
import re

parser = argparse.ArgumentParser(description='score model TCR-pMHC complex using Rosetta\'s score app')
parser.add_argument('-i', help='PDB file of TCR-pMHC complex (ex. xxxx.pdb)', type=str, dest='i',required=True)
parser.add_argument('-ros', help='Rosetta3.10 bin path (ex. $YourRosettaPath/main/)', type=str, dest='ros',required=True)
parser.add_argument('-c', help='MHC alpha chain ID, peptide chain ID, TCR alpha chian ID and TCR beta chain ID of MHC-pMHC PDB file (ex. ACDE)', type=str, dest='c',required=True)
args = parser.parse_args()

def weight_file():
	'''
	Generate a feature weight file for Rosetta3 scoring function 
	'''
	features = ['dslf_ca_dih', 'dslf_cs_ang', 'dslf_ss_dih', 'dslf_ss_dst', 'fa_atr', 'fa_dun', 'fa_intra_rep', 'fa_pair', 'fa_rep', 'fa_sol', 'hbond_bb_sc', 'hbond_lr_bb', 'hbond_sc', 'hbond_sr_bb', 'omega', 'p_aa_pp', 'pro_close', 'rama', 'ref']
	w_1 = [(x+ ' 1.0\n') for x in features]
	w_2 = ''.join(w_1)
	out_fop = open('weights.txt', 'w')
	out_fop.write(w_2)
	out_fop.close()

def score(pdb, weights, label):
	'''
	Score PDB complex using Rosetta3 scoring function
	'''
	ros_path = args.ros
	score_cmd = [ros_path + 'source/bin/score.linuxgccrelease', '-database',
	ros_path + 'database/', '-s', pdb , '-out:file:scorefile', 
	label + '_score.sc', '-extrachi_cutoff', '1', '-ex1', '-ex2', '-ex3', '-score:weights',
	weights]
	process = subprocess.Popen(score_cmd)
	process.wait()

def GetPDBchain(pdb, chain_id):
	'''
	Extract protein chain from pdb file
	'''
	pdb_lines = ''    
	for line in open(pdb, 'r'):
		if ((line[0:4] == 'ATOM') and (line[21] == chain_id)):
			pdb_lines = pdb_lines + line
	pdb_lines = pdb_lines + 'TER\n'
	return(pdb_lines)

def dG_bind(COM_sc, TCR_sc, pMHC_sc):
	'''
	Calculate dG_bind scores = COM_score - (TCR_score + pMHC_score)
	'''
	COM = open(COM_sc, 'r')
	TCR = open(TCR_sc, 'r')
	pMHC= open(pMHC_sc, 'r')

	header = COM.readline().split()[2:21]
	TCR.readline()
	pMHC.readline()
	COM_scores = np.array(map(float, COM.readline().split()[2:21]))
	TCR_scores = np.array(map(float, TCR.readline().split()[2:21]))
	pMHC_scores = np.array(map(float, pMHC.readline().split()[2:21]))
	dG_scores = COM_scores - (TCR_scores + pMHC_scores)
	dG_hash = {}
	for i, energy in enumerate(header):
		dG_hash[energy] = dG_scores[i]
	COM.close()
	TCR.close()
	pMHC.close()
	return(dG_hash)

def ATLAS_dG_bind(dG_hash):
	'''
	Calculate ATLAS dG bind score                                                                         
	ATLAS_dG_bind = 2.1609*fa_atr + 0.1564*fa_pair + 0.0263*fa_rep + 1.6697*fa_sol - 0.3122*hbond_bb_sc - 0.5961*hbond_lr_bb + 0.4873*hbond_sc - 0.2401*hbond_sr_bb - 6.6377
	'''
	features = ['fa_atr', 'fa_pair', 'fa_rep', 'fa_sol', 'hbond_bb_sc', 'hbond_lr_bb', 'hbond_sc', 'hbond_sr_bb']
	ATLAS_dG_bind = 2.1609*dG_hash[features[0]] + 0.1564*dG_hash[features[1]] + 0.0263*dG_hash[features[2]] + 1.6697*dG_hash[features[3]] - 0.3122*dG_hash[features[4]] - 0.5961*dG_hash[features[5]] + 0.4873*dG_hash[features[6]] - 0.2401*dG_hash[features[7]] - 6.6377
	out_arr = [(x+ '\t' +str(dG_hash[x])+ '\n') for x in features]
	out_txt = ''.join(out_arr)
	out_txt = out_txt + 'ATLAS_bind_score' + '\t' +str(ATLAS_dG_bind)
	return(out_txt)

def main():
	pdb = args.i
	chains = args.c
	# Get protein chains
	chain_mhc_a = GetPDBchain(pdb, chains[0])
	chain_pep = GetPDBchain(pdb, chains[1])
	chain_tcr_a = GetPDBchain(pdb, chains[2])
	chain_tcr_b = GetPDBchain(pdb, chains[3])
	# TCR-pMHC complex
	com_fop = open('TCR-pMHC.pdb', 'w')
	com_fop.write(chain_mhc_a)
	com_fop.write(chain_pep)
	com_fop.write(chain_tcr_a)
	com_fop.write(chain_tcr_b)
	com_fop.close()
	# pMHC
	pmhc_fop = open('pMHC.pdb', 'w')
	pmhc_fop.write(chain_mhc_a)
	pmhc_fop.write(chain_pep)
	pmhc_fop.close()
	# TCR
	tcr_fop = open('TCR.pdb', 'w')
	tcr_fop.write(chain_tcr_a)
	tcr_fop.write(chain_tcr_b)
	tcr_fop.close()
	# Generate a feature weight file
	weight_file()
	# Score TCR-pMHC
	score('TCR-pMHC.pdb', 'weights.txt', 'COM')
	# Score TCR
	score('TCR.pdb', 'weights.txt', 'TCR')
	# Score pMHC
	score('pMHC.pdb', 'weights.txt', 'pMHC')			
	# Calculate dG bind
	dG_scores = dG_bind('COM_score.sc', 'TCR_score.sc', 'pMHC_score.sc')
	# ATLAS score
	ATLAS_score = ATLAS_dG_bind(dG_scores)
	out_f = pdb.replace('.pdb', '.sc')
	out_fop = open(out_f, 'w')
	out_fop.write(ATLAS_score)
	out_fop.close()
	# Remove temporary files
	os.remove('COM_score.sc')
	os.remove('TCR_score.sc')
	os.remove('pMHC_score.sc')
	os.remove('TCR-pMHC.pdb')
	os.remove('TCR.pdb')
	os.remove('pMHC.pdb')
	os.remove('weights.txt')
	print('###########Scoring for TCR-pMHC is finished!###########')
	print(ATLAS_score)
	print('#######################################################')

if __name__ == '__main__':
	main()
