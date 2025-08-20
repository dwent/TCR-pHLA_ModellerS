#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################################
#Pipeline for modelling the protein structure of TCR-pMHC complex                       #
#Usage:TCR-pMHC_modelling.py [TCR-pMHC complex sequence] [Output dir] [Number of models]#
#########################################################################################

import sys
import os

#####Rosetta path#############################################################
ros_path = '/vdb/User/wuxiaoqin/source_Rosetta/rosetta_src_2018.33.60351_bundle/main/'  
##############################################################################

##TCR-pMHC complex modelling
#Get TCR-pMHC complex modelling script dir
scr_dir = ''

if (sys.argv[0] != 'TCR-pMHC_modelling.py'):
    scr_dir = sys.argv[0].replace('TCR-pMHC_modelling.py',  '')
else:
    scr_dir = './'

#Input
try:
    seq = sys.argv[1]  #TCR-pMHC complex sequence
    out_dir = sys.argv[2]  #Output dir
    model_num = sys.argv[3]  #Number of models
except:
    print('Input error! Some parameters are missing!')
    print('Usage:TCR-pMHC_modelling.py [TCR-pMHC complex sequence] [Output dir] [Number of models]')
    exit()
    

model_cmd ='perl ' +scr_dir+ 'Pipline_modeller-improve.pl ' +scr_dir+ 'Template ' +seq+ ' ' +out_dir+' ' +model_num+ ' ' +scr_dir
os.system(model_cmd)

print('ok')

##TCR-pMHC model scoring
model_score = {}
model_config = out_dir+ '/config.txt'  #Modelling config file
for li in open(model_config, 'r'):
    #Scoring for modelled pdb
    li_arr = li.split('     ')
    model = li_arr[0]    
    score_cmd = scr_dir+ 'ATLAS_score.py -i ' +model+ ' -ros ' +ros_path+ ' -c ACDE'
    os.system(score_cmd)

    #Get score of modelled pdb
    score = ''
    score_file = model.replace('.pdb', '.sc')
    for sc in open(score_file, 'r'):
        if (sc.find('ATLAS_bind_score') != -1):
            sc_arr = sc.split('\t')
            score = sc_arr[1].replace('\n', '')
            model_score[model] = round(float(score), 2)

##Output score results
#Sort modelled pdbs
model_score_sorted = sorted(model_score.items(), key=lambda x: x[1], reverse=False)
#Output
out_file = out_dir+ '/Score.txt'
out_fop = open(out_file, 'w')
out_fop.write('\t'.join(['Model', 'Score']) + '\n')
model_id = 0
for model,sc in model_score_sorted:
    model_id = model_id + 1
    #Change file name of modelled pdb
    model_new = out_dir + '/model' + str(model_id) + '.pdb'
    mv_cmd = 'mv ' + model + ' ' + model_new
    os.system(mv_cmd)
    #Change name of score file
    mv_cmd_sc = mv_cmd.replace('.pdb', '.sc')
    os.system(mv_cmd_sc)
    out_fop.write('\t'.join(['model'+str(model_id), str(sc)]) + '\n')
out_fop.close()
os.system('rm -f status.txt')
