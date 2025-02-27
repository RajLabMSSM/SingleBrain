#Chip Arreay QC

## Split autosomal and X
input="Input file"
out_dir="Path output"
sample_name="Sample name"

ml plink

#vcf to bfile
plink --vcf ${input} --chr X --make-bed --out ${out_dir}/${sample_name}.chrX
plink --vcf ${input} --autosome --make-bed --out ${out_dir}/${sample_name}.auto


## Varint QC
input_dir="input directory"
output_dir="output directory"

ml plink

plink --bfile ${input_dir}/${sample_name}.auto \
	--maf 0.01 \
	--geno 0.02 \
	--keep-allele-order \
	--make-bed --out ${output_dir}/${sample_name}.auto.maf0.01.geno0.02

control="Control list"
plink --bfile ${output_dir}/${sample_name}.auto.maf0.01.geno0.02 \
	--hwe 0.000001 \
	--keep-fam $control \
	--keep-allele-order \
	--make-bed --out ${output_dir}/${sample_name}.auto.maf0.01.geno0.02.control

awk '{print $2}' ${output_dir}/${sample_name}.auto.maf0.01.geno0.02.control.bim > ${output_dir}/maf0.01.geno0.02.hwe10E_6.auto.bim_list



plink --bfile ${output_dir}/${sample_name}.auto.maf0.01.geno0.02 \
   --extract ${output_dir}/maf0.01.geno0.02.hwe10E_6.auto.bim_list \
   --keep-allele-order \
   --make-bed \
   --out ${output_dir}/${sample_name}.auto.maf0.01.geno0.02.hwe10E_6

## Sample Check

plink --bfile ${input_dir}/${sample_name}.chrX \
	--maf 0.01 \
	--geno 0.02 \
	--keep-allele-order \
	--make-bed --out ${output_dir}/${sample_name}.chrX.maf0.01.geno0.02


plink --bfile ${output_dir}/${sample_name}.auto.maf0.01.geno0.02.hwe10E_6 \
   --missing --het \
   --out ${output_dir}/${sample_name}.auto.maf0.01.geno0.02.hwe10E_6


plink --bfile ${output_dir}/${sample_name}.chrX.maf0.01.geno0.02 \
	--check-sex \
	--out ${output_dir}/${sample_name}.chrX

ml king
king -b ${output_dir}/${sample_name}.auto.maf0.01.geno0.02.hwe10E_6.bed \
         --kinship --ibs \
        --unrelated --degree 2 \
        --rpath ${output_dir}

plink --bfile ${output_dir}/${sample_name}.pca \
	--indep-pairwise 50 5 0.5 \
	--out ${output_dir}/${sample_name}.pca

plink --bfile ${output_dir}/${sample_name}.pca \
	--extract ${output_dir}/${sample_name}.pca.prune.in \
	--make-bed --out ${output_dir}/${sample_name}.pca.pruned


ml plink2
plink2 --bfile ${output_dir}/${sample_name}.pca.pruned \
	--freq --pca approx var-wts \
	--threads 20 \
	--out ${output_dir}/${sample_name}.pca


## Sample QC
cat "Low quality sample list" > ${path_dir}/sample_qc.list

plink --bfile ${output_dir}/VQ/${sample_name}.auto.maf0.01.geno0.02.hwe10E_6 \
   --remove ${path_dir}/sample_qc.list \
   --make-bed \
   --out ${output_dir}/SQ/${sample_name}.QC


## Split by chr and vcf
split_dir=''

for CHR in `seq 1 22`;
do

plink --bfile $file --allow-extra-chr --chr $CHR \
   --make-bed \
   --out ${split_dir}/sea_ad.hg19.chr$CHR

bgzip ${split_dir}/sea_ad.hg19.chr${CHR}.vcf

done

## imputation bot

${imputationbot} impute --refpanel topmed-r2 \
        --population mixed \
        --autoDownload \
        --files ${split_dir}/*vcf.gz



## Imputation filtering
imputation_dir=''

for i in {1..22}
do
echo chr${i}.info.gz
zcat ${imputation_dir}/chr${i}.info.gz| awk '{if($7>0.2) print $1}'| sed '1d' > ${imputation_dir}/info/chr${i}.info.Rsq0.2
done


## Filtering 
result_dir=''

for CHR in `seq 1 22`;
do
	plink --vcf ${imputation_dir}/chr${CHR}.dose.vcf.gz \
        	--make-bed \
        	--keep-allele-order \
        	--out ${imputation_dir}/chr$CHR

	plink2 --bfile ${imputation_dir}/chr${CHR} \
        	--make-bed \
        	--maf 0.005 \
        	--extract ${imputation_dir}/info/chr${CHR}.info.Rsq0.8 \
        	--import-dosage-certainty 0.9 \
        	--snps-only \
        	--rm-dup \
        	--out ${result_dir}/chr${CHR}

done


# Merge 

for chr in {1..22}
do
  imputed=${result_dir}/chr${chr}
  # make mergelist
  if [ $chr -ne 1 ];then
    if [ $chr -eq 2 ];then
      echo $imputed > ${result_dir}/merge${chr}.list
    fi
    if [ $chr -ne 2 ];then
      echo $imputed >> ${result_dir}/merge${chr}.list
    fi
  fi
done

plink --bfile ${result_dir}/chr1 \
	--merge-list ${result_dir}/merge.list \
	--keep-allele-order \
	--make-bed --out ${result_dir}/${sample_name}.hg38



#End


