import os
import time
import hail as hl
from pprint import pprint


from hail.plot import show
from gnomad_qc.v2.resources.sample_qc import *
from gnomad_qc.v2.resources import get_gnomad_data

import argparse
import logging
import csv
import sys
import pandas as pd
import numpy as np
import glob
import time
from tqdm import tqdm

from pprint import pprint

fn_meta=" _ input meta data _ "

fn_LCR="~/LCR-hs38.bed"

time.strftime('%c', time.localtime(time.time()))

hl.init(master='local[20]', spark_conf={'spark.driver.memory':'200g'},default_reference='GRCh38')
hl.plot.output_notebook()

### VQSR Filter

mt = hl.import_vcf(vcf_dir + file_vcf, 
                        force_bgz = True, reference_genome=vcf_version, 
                        contig_recoding = {f"{k}":f"chr{k}" for k in (list(range(1, 23)) + ['X', 'Y'])},
                        array_elements_required=False)

mt = mt.filter_rows(hl.len(mt.filters) == 0)

### Split Multi-allelic
mt = mt.select_entries(mt.AD, mt.DP, mt.GQ, mt.GT)
mt = hl.split_multi_hts(mt, permit_shuffle=True)
mt = mt.key_rows_by(**hl.min_rep(mt.locus, mt.alleles))

### LCR region
bed_LCR = hl.import_bed(fn_LCR)

mt = mt.annotate_rows(LCR_region = bed_LCR[mt.locus])
mt = mt.filter_rows(hl.is_defined(mt.LCR_region), keep=False)

mt=hl.sample_qc(mt, name='sample_qc')
    
## Genotype QC

mt = hl.read_matrix_table(path_input)

mt = mt.annotate_entries(AB = mt.AD[1] / hl.sum(mt.AD))
mt = mt.filter_entries((mt.GQ>=20))
mt = mt.filter_entries((mt.DP>=10) & (mt.DP<=200))

mt = mt.filter_entries((mt.GT.is_hom_ref()&(mt.AB <= 0.2))|
                           (mt.GT.is_het()&(mt.AB >= 0.2))&(mt.AB  <= 0.8)|
                           (mt.GT.is_hom_var()&(mt.AB >= 0.9))
                           , keep=True) 

mt=hl.sample_qc(mt, name='sample_qc')
mt=hl.variant_qc(mt, name='variant_qc')

mt = mt.annotate_rows(rsid = (hl.delimit([mt.locus.contig, hl.str(mt.locus.position), mt.alleles[0], mt.alleles[1]], ':')))

## DP 

mt = mt.filter_rows((hl.len(mt.alleles) == 2) & (hl.is_snp(mt.alleles[0], mt.alleles[1])) &
                            (mt.variant_qc.AF[1] > 0.001) &
                            (mt.variant_qc.call_rate > 0.99))
                            
mt=hl.sample_qc(mt, name='sample_qc')

mean_ht = mt.sample_qc.dp_stats.mean
call_ht = mt.sample_qc.n_called

mean_ht.export'~/' + prefix + '_sample_qc1_mean_dp_' + file_name)
call_ht.export('~/' + prefix  +'_sample_qc1_n_call_' + file_name)


## Sample quality

mt = mt.filter_rows((hl.len(mt.alleles) == 2) &
                    (mt.variant_qc.call_rate > 0.99))

mt=hl.sample_qc(mt, name='sample_qc')
        
mt_SNV = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
mt_Indel = mt.filter_rows(hl.is_indel(mt.alleles[0], mt.alleles[1]))
        
mt_SNV = hl.sample_qc(mt_SNV)
mt_Indel = hl.sample_qc(mt_Indel)
        
n_non_ref_SNV = mt_SNV.cols()[mt.s].sample_qc.n_non_ref
n_non_ref_Indel = mt_Indel.cols()[mt.s].sample_qc.n_non_ref
        
n_non_ref_SNV.export('~/' + prefix +'_sample_qc2_n_snp_' + file_name)
n_non_ref_Indel.export('~/' + prefix +'_sample_qc2_n_indel_' + file_name)

df_nrow=pd.DataFrame()
nrow = mt.count_rows()
tmp = pd.DataFrame({'POS':[file_name],
                                'nrow':[nrow]})
df_nrow = pd.concat([df_nrow, tmp], axis=0)
df_nrow.to_csv('~/' + prefix +'_sample_qc2_n_variant_' + file_name, sep='\t', index=False)

call_ht = mt.sample_qc.n_called
call_ht.export('~/' + prefix +'_sample_qc2_n_call_' + file_name)

n_insertion_ht = mt.sample_qc.n_insertion
n_insertion_ht.export('~/' + prefix +'_sample_qc2_n_insertion_' + file_name)

n_deletion_ht = mt.sample_qc.n_deletion
n_deletion_ht.export('~/' + prefix +'_sample_qc2_n_deletion_' + file_name)

n_het_ht = mt.sample_qc.n_het
n_het_ht.export('~/' + prefix +'_sample_qc2_n_het_' + file_name)

n_hom_var_ht = mt.sample_qc.n_hom_var
n_hom_var_ht.export('~/' + prefix +'_sample_qc2_n_hom_var_' + file_name)

n_transition_ht = mt.sample_qc.n_transition
n_transition_ht.export('~/' + prefix +'_sample_qc2_n_transition_' + file_name)

n_transversion_ht = mt.sample_qc.n_transversion
n_transversion_ht.export('~/' + prefix +'_sample_qc2_n_transversion_' + file_name)

n_singleton_ht = mt.sample_qc.n_singleton
n_singleton_ht.export('~/' + prefix +'_sample_qc2_n_singleton_' + file_name)

## Variant QC 

# +
remove_list ='~/' + prefix + "_sampleQC_final_list"
df_sample=pd.read_csv(remove_list, sep='\t')

samples_to_remove = df_sample.s.tolist()
samples_to_remove=str(samples_to_remove)
set_to_remove = hl.literal(samples_to_remove)

mt_variantQC1 = hl.read_matrix_table(path_input)
mt_variantQC1 = mt_variantQC1.filter_cols(~set_to_kinship.contains(mt_variantQC1['s']))
mt_variantQC1 = mt_variantQC1.filter_cols(~set_to_remove.contains(mt_variantQC1['s']))
mt_variantQC1 = mt_variantQC1.filter_cols(set_to_HC.contains(mt_variantQC1['s']))

mt_variantQC1= hl.variant_qc(mt_variantQC1, name='variant_qc')
mt_variantQC1_HWE = mt_variantQC1.filter_rows(mt_variantQC1.variant_qc.p_value_hwe < 1e-09)

row_key_hwe = mt_variantQC1_HWE.row_key
row_key_hwe.export('~/' + file_name + '.VQ1.GQ.UN.HC.' + pop_name + '.HWE.variant')

# +
## Reloading the matrix for removing variants without only QC failed samples

mt_variantQC1 = hl.read_matrix_table(path_input)
mt_variantQC1 = mt_variantQC1.filter_cols(~set_to_remove.contains(mt_variantQC1['s']))

# +
mt_variantQC1=hl.variant_qc(mt_variantQC1, name='variant_qc')
mt_variantQC1 = mt_variantQC1.filter_rows((mt_variantQC1.variant_qc.AC[1] > 0 ))
# +
mt_variantQC1=hl.variant_qc(mt_variantQC1, name='variant_qc')
mt_variantQC1 = mt_variantQC1.filter_rows((mt_variantQC1.variant_qc.call_rate >= 0.9))

mt_variantQC1 = mt_variantQC1.filter_rows(hl.is_defined(eur_ht[mt_variantQC1.row_key]), keep=False) 

mt_variantQC2_SNP = mt_variantQC2.filter_rows(hl.is_snp(mt_variantQC2.alleles[0], mt_variantQC2.alleles[1]))

HARD_CUTOFF = 
mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows(mt_variantQC2_SNP.info.QD >= HARD_CUTOFF)

HARD_CUTOFF = 
mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows(mt_variantQC2_SNP.info.SOR <= HARD_CUTOFF)

HARD_CUTOFF = 
mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows(mt_variantQC2_SNP.info.FS <= HARD_CUTOFF)

HARD_CUTOFF = 
mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows(mt_variantQC2_SNP.info.MQ >= HARD_CUTOFF)

HARD_CUTOFF = 
mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows((mt_variantQC2_SNP.info.MQRankSum > -HARD_CUTOFF)&
                                                 (mt_variantQC2_SNP.info.MQRankSum < HARD_CUTOFF))

HARD_CUTOFF = 
mt_variantQC2_SNP = mt_variantQC2_SNP.filter_rows((mt_variantQC2_SNP.info.ReadPosRankSum >= -HARD_CUTOFF) &
                                                  (mt_variantQC2_SNP.info.ReadPosRankSum <= HARD_CUTOFF))

mt_variantQC2_INDEL = mt_variantQC2.filter_rows(hl.is_indel(mt_variantQC2.alleles[0], mt_variantQC2.alleles[1]))

HARD_CUTOFF =
mt_variantQC2_INDEL = mt_variantQC2_INDEL.filter_rows(mt_variantQC2_INDEL.info.QD >= HARD_CUTOFF)

HARD_CUTOFF = 
mt_variantQC2_INDEL = mt_variantQC2_INDEL.filter_rows(mt_variantQC2_INDEL.info.FS <= HARD_CUTOFF)

HARD_CUTOFF = 
mt_variantQC2_INDEL = mt_variantQC2_INDEL.filter_rows((mt_variantQC2_INDEL.info.ReadPosRankSum > -HARD_CUTOFF)&
                                                      (mt_variantQC2_INDEL.info.ReadPosRankSum < HARD_CUTOFF))


mt_variantQC2_merged = mt_variantQC2_SNP.union_rows(mt_variantQC2_INDEL)
mt_variantQC2_merged.write(result_name, overwrite=True)

# Merge Export

files_joined = os.path.join(dir_output + 'VQ3/chr[0-9]*.VQ1.GQ.VQ2.VQ3.mt')
list_files = glob.glob(files_joined)

tables = [hl.read_matrix_table(list_files[j]) 
                for j in range(0, len(list_files))]
    
print(len(tables))
mt = hl.MatrixTable.union_rows(*tables)
hl.export_plink(mt, dir_output +'export/' + prefix,
                ind_id = mt.s, fam_id =mt.s)



