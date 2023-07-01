#!/usr/bin/env python3

'''Output all the docking commands'''

import glob, sys, re, os

#refined-set is redocking
for fname in glob.glob('refined-set/*/*ligand_rdkit_250_etkdg_2dock*.sdf.gz')+glob.glob('refined-set/*/*ligand_dmcg_250_2dock*.sdf.gz'):  
    if 'docked' in fname:
        continue
    if 'tmp' in fname:
        continue
    m = re.search(r'refined-set/(\S+)/',fname)
    pdb = m.group(1)
    base = f'refined-set/{pdb}/{pdb}'
    print(f'gnina -r {base}_rec.pdb --autobox_ligand {base}_ligand.sdf -l {fname} --exhaustiveness 8 --cpu 1 --seed 0 -o {fname[:-7]}_docked.sdf.gz; obrms -f {base}_ligand.sdf {fname[:-7]}_docked.sdf.gz > {fname[:-7]}_docked.rms.txt')
    
    if '2dock_single_eranked_min.sdf.gz' in fname:
        #also try five times the sampling
        print(f'gnina -r {base}_rec.pdb --autobox_ligand {base}_ligand.sdf -l {fname} --exhaustiveness 40 --cpu 1 --seed 0 -o {fname[:-7]}_e5_docked.sdf.gz; obrms -f {base}_ligand.sdf {fname[:-7]}_e5_docked.sdf.gz > {fname[:-7]}_e5_docked.rms.txt')
        #alternative extra sampling
        cmd = ''
        for s in range(5):
            cmd += f'gnina -r {base}_rec.pdb --autobox_ligand {base}_ligand.sdf -l {fname} --exhaustiveness 8 --cpu 1 --seed {s} -o {fname[:-7]}_s{s}_tmp.sdf.gz; '
        cmd += f'zcat {fname[:-7]}_s?_tmp.sdf.gz | gzip > {fname[:-7]}_s5_docked.sdf.gz; rm  {fname[:-7]}_s?_tmp.sdf.gz ;'
        cmd += f'obrms -f {base}_ligand.sdf {fname[:-7]}_s5_docked.sdf.gz> {fname[:-7]}_s5_docked.rms.txt'
        print(cmd)

#crossdocking - the combinatorial number of pairs is downsampled in this file
for line in open('wierbowski_cd/ds_cd_input_pairs.txt'):
    rec,lig,reflig,_ = line.split()
    prefix = lig[:-4]
    rec_prefix = rec[:-4]
    if not os.path.exists(rec):
        sys.stderr.write(f'{rec} doesn\'t exist\n')
    if not os.path.exists(reflig):
        sys.stderr.write(f'{reflig} doesn\'t exist\n')
        
    for fname in glob.glob(f'{prefix}_dmcg_250_2dock*.sdf.gz')+glob.glob(f'{prefix}_rdkit_250_etkdg_2dock*.sdf.gz'):
        if 'docked' in fname:
            continue         
        if 'tmp' in fname:
            continue
        m = re.search(r'wierbowski_cd/\S+/(\S+)_LIG_aligned',fname)
        pdb = m.group(1)
        
        m = re.search(r'wierbowski_cd/\S+/(\S+)_PRO.pdb',rec)
        rec_prefix = m.group(1)
        dname = f'{fname[:-7]}_docked_{rec_prefix}'
        print(f'gnina -r {rec} --autobox_ligand {reflig} -l {fname} --exhaustiveness 8 --cpu 1 --seed 0 -o {dname}.sdf.gz; obrms -f {prefix}.sdf {dname}.sdf.gz > {dname}.rms.txt')
        if '2dock_single_eranked_min.sdf.gz' in fname:
            #also try five times the sampling
            dname = f'{fname[:-7]}_e5_docked_{rec_prefix}'

            print(f'gnina -r {rec} --autobox_ligand {reflig} -l {fname} --exhaustiveness 40 --cpu 1 --seed 0 -o {dname}.sdf.gz; obrms -f {prefix}.sdf {dname}.sdf.gz > {dname}.rms.txt')
            #alternative extra sampling
            cmd = ''
            for s in range(5):
                cmd += f'gnina -r {rec} --autobox_ligand {reflig} -l {fname} --exhaustiveness 8 --cpu 1 --seed {s} -o {fname[:-7]}_s{s}_tmp_{rec_prefix}.sdf.gz; '
            cmd += f'zcat {fname[:-7]}_s?_tmp_{rec_prefix}.sdf.gz | gzip > {fname[:-7]}_s5_docked_{rec_prefix}.sdf.gz; rm  {fname[:-7]}_s?_tmp_{rec_prefix}.sdf.gz ;'
            cmd += f'obrms -f {prefix}.sdf {fname[:-7]}_s5_docked_{rec_prefix}.sdf.gz > {fname[:-7]}_s5_docked_{rec_prefix}.rms.txt'
            print(cmd)

    
