#!/usr/bin/env python

import os

#####Function starts#####
ProtCode = {
    'GLY':'G',
    'ALA':'A',
    'VAL':'V',
    'LEU':'L',
    'ILE':'I',
    'SER':'S',
    'THR':'T',
    'CYS':'C',
    'MET':'M',
    'PRO':'P',
    'ASP':'D',
    'ASN':'N',
    'GLU':'E',
    'GLN':'Q',
    'LYS':'K',
    'ARG':'R',
    'HIS':'H',
    'PHE':'F',
    'TYR':'Y',
    'TRP':'W',
    'ASX':'B',
    'GLX':'Z',
    'UNK':'X'
    }

def pdb_seq(pdb, chain_id):
    id_arr = []
    id_aa = {}
    pro_seq = ''
    
    for line in open(pdb, 'r'):
        if((line[0:4] == 'ATOM') and (line[21] == chain_id)): 
            aa_id = line[23:26]
            id_arr.append(aa_id)
            id_aa[aa_id] = line[17:20]

    id_arr_1 = sorted(list(set(id_arr)))
    for aa_id in id_arr_1:
        aa_id_1 = aa_id.strip() 
        #aa = id_aa[aa_id]
        aa = ProtCode[id_aa[aa_id]]
        pro_seq = pro_seq + aa
    return(pro_seq)
#####Function stops######

model_dir = 'structure/'
seq_dir = 'sequence/'
os.system('mkdir ' +seq_dir)

pdb_arr = []

for root, dirs, files in os.walk(model_dir):
    pdb_arr = files

for pdb in pdb_arr:
    pdb_id = pdb.replace('.pdb', '')
    print(pdb_id)
    os.system('mkdir ' +seq_dir + pdb_id)
    ch_a = seq_dir + pdb_id + '/' + pdb_id+ '_A.fasta'
    ch_b = seq_dir + pdb_id + '/' + pdb_id+ '_B.fasta'
    ch_c = seq_dir + pdb_id + '/' + pdb_id+ '_C.fasta'
    ch_d = seq_dir + pdb_id + '/' + pdb_id+ '_D.fasta'
    ch_e = seq_dir + pdb_id + '/' + pdb_id+ '_E.fasta'
    
    fop_a = open(ch_a, 'w')
    fop_a.write('>' +pdb_id+ ' A' + '\n')
    fop_a.write(pdb_seq(model_dir+pdb, 'A') + '\n')
    fop_a.close()
    
    fop_b = open(ch_b, 'w')
    fop_b.write('>' +pdb_id+ ' B' + '\n')
    fop_b.write(pdb_seq(model_dir+pdb, 'B') + '\n')
    fop_b.close()

    fop_c = open(ch_c, 'w')
    fop_c.write('>' +pdb_id+ ' C' + '\n')
    fop_c.write(pdb_seq(model_dir+pdb, 'C') + '\n')
    fop_c.close()

    fop_d = open(ch_d, 'w')
    fop_d.write('>' +pdb_id+ ' D' + '\n')
    fop_d.write(pdb_seq(model_dir+pdb, 'D') + '\n')
    fop_d.close()

    fop_e = open(ch_e, 'w')
    fop_e.write('>' +pdb_id+ ' E' + '\n')
    fop_e.write(pdb_seq(model_dir+pdb, 'E') + '\n')
    fop_e.close()
    









