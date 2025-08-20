from modeller import *
import sys, getopt

if __name__ == "__main__":
	# process command arguments
	template_pdb = ''
	chain_id_first = ''
	chain_id_last = ''
	target_seq = ''
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'ht:b:e:q:')
	except getopt.GetoptError:
		print('python3 salign.py -t <template_pdb> -b <chain_begin> -e <chain_end> -q <query sequence>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-t':
			template_pdb = arg
		elif opt == '-b':
			chain_id_first = arg
		elif opt == '-e':
			chain_id_last = arg
		elif opt == '-q':
			target_seq = arg
		elif opt == '-h':
			print("python3 salign.py -t <template_pdb> -b <chain_begin> -e <chain_end> -q <query sequence>")
			sys.exit(2)
	if template_pdb == '' or chain_id_first == '' or chain_id_last == '' or target_seq == '':
		print('python3 salign.py -t <template_pdb> -b <chain_begin> -e <chain_end> -q <query sequence>')
		sys.exit(2)

	env = environ()
	aln = alignment(env)
	#template_pdb = input("Template PDB file:")
	#chain_id_first = input("Template PDB chain ID first:")
	#chain_id_last = input("Template PDB chain ID last:")

	template_name = template_pdb.split(".")[-2]
	template_name = template_name.split("/")[-1]
	file_path = "/".join(template_pdb.split("/")[0:-1])
	print(template_name)
	mdl = model(env, file=template_pdb, model_segment=('FIRST:'+chain_id_first,'LAST:'+chain_id_last))
	if chain_id_first == chain_id_last:
		template_name = template_name+'_'+chain_id_first
	else:
		template_name = template_name+'_'+chain_id_first+'_'+chain_id_last

	aln.append_model(mdl, align_codes=template_name, atom_files=template_pdb)

	aln_target = alignment(env)
	#target_seq = input("target seq file(.ali):")
	target_name = ''
	read_target = modfile.File(target_seq, 'r')
	while aln_target.read_one(read_target, alignment_format='PIR'):
		target_name = aln_target[0].code
		break
	read_target.close()
	print(target_name+"--------------")

	aln.append(file=target_seq, align_codes=target_name) #align_codes=target_name
	#aln.align2d()
	aln.salign()
	aln.write(file=file_path+'/'+target_name+'-'+template_name+'.ali', alignment_format='PIR')
	aln.write(file=file_path+'/'+target_name+'-'+template_name+'.pap', alignment_format='PAP')
