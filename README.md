# ExomeDepthToVCF
----
This project is dedicated to the conversion of ExomeDepth output (CNV.calls) to VCF format and eventually the union with another VCF


## Usage
----
ExomeDepthToVCF requires ExomeDepth output file (in .csv format), BED file used in CNV calling process by the tool and reference genome fasta file. Optionally, some other parameters to determine CNV calls FILTER values could be specified. 

    RUNNING:
        
    python calls2vcf_conversion_v1.0.py 
    	--exomedepth CNV.calls 
    	--bed bed_file 
    	--fasta genome.fasta
    	--output output_file.vcf
    	--bedtools path_to_bedtools_bin
    	--sample	sample_name
    
    
    Optionally, some other parameters could be specified (default values are shown):
    
    python calls2vcf_conversion_v1.0.py 
    	--exomedepth CNV.calls 
    	--bed bed_file 
    	--fasta genome.fasta
    	--output output_file.vcf
    	--bedtools path_to_bedtools_bin
    	--sample	sample_name
    	--upper_threshold 1.2
    	--lower_threshold 0.8
    	--qual_threshold 6
    	

## Options

|   Option   |      Description      |  Required | Type | Default |
|-------|:------- |:---------:|:-------|:-------:|
| -e --exomedepth  |  ExomeDepth CNV calls output   | True | String | -  |
| -b --bed              |  BED file used as target for CNV calling                    | True      | String | -       |
| -f --fasta            |  Reference genome in fasta format                           | True      | String | -       |
| -o --output           |  Path to vcf output file                                    | True      | String | -       |
| -s --sample           |  Sample name                                                | True      | String | -       |
| -t --bedtools         |  Path to bedtools binary                                    | True      | String | -       |
| -ut --upper_threshold | Upper threshold for read-ratio based quality filter         | False     | Double | 1.2     |
| -lt --lower_threshold | Lower threshold for read-ratio based quality filter         | False     | Double | 0.8     |
| -q --qual_threshold   | Bayes factor threshold for ExomeDepth quality control value | False     | Double | 6       |


## Output
ExomeDepthToVCF generates VCF reporting CNV in a standardized format. 

FILTER field could reports 3 values:
- PASS: all filters pass
- reads_ratio_close_to_1: read ratio between case and reference is in range lower_threshold-upper_threshold
- low_bayes_factor: Bayes Factor computed by ExomeDepth as CNV quality value is less than qual_threshold

INFO field reports:
- SVTYPE: CNV type could be DEL or DUP
- SVLEN: CNV size
- READRATIO: read ration between case and controls for specific target
- BF: Bayes Factor as computed by ExomeDepth
- DP: Approximate depth of coverage (some reads may have been filtered)

NOTE:

In rare cases ExomeDepth reports two contiguous and partially overlapping CNV, one as deletion and one as duplication. ExomeDepthToVCF splits those CNV assigning type as follow:
- the overlapping target interval are assigned to the smaller interval type
- the remaining target intervals are asigned to the larger interval

For the refactor CNV, read ratio is calculated removing all reads belonging to the first interval. 
