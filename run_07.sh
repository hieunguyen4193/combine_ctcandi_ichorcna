output_version="20241229";
outdir="/media/hieunguyen/HNSD_mini/outdir";
for mode in all hypo hyper;do \
    for input_cancer_class in Liver Lung Breast CRC pan_cancer;do \
        for dataset_name LOD SPIKE_IN REPORT4 VALIDATION;do \
            echo -e "Working on mode: " $mode ", single/multi cancer tpe: " $input_cancer_class ", dataset " $dataset_name;
            python 07_count_reads_in_regions.py --output_version $output_version --dataset_name $dataset_name --input_cancer_class $input_cancer_class --outdir $outdir --mode $mode;\
                done;done;done