nextflow run main.nf -profile H1N1 \
  --input_dir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/H1N1 \
  --outdir results/H1N1

nextflow run main.nf -profile H3N2 \
  --input_dir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/H3N2 \
  --outdir results/H3N2

nextflow run main.nf -profile B_VIC \
  --input_dir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/B_VIC \
  --outdir results/B_VIC

nextflow run main.nf -profile B_YAM \
  --input_dir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/B_YAM \
  --outdir results/B_YAM

nextflow run main.nf -profile H5NX \
  --input_dir h5_example \
  --outdir results/H5NX



nextflow run main.nf -profile H1N1 \
  --input_dir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/H1N1 \
  --outdir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/results/H1N1

nextflow run main.nf -profile H3N2 \
  --input_dir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/H3N2 \
  --outdir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/results/H3N2

nextflow run main.nf -profile B_VIC \
  --input_dir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/B_VIC \
  --outdir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/results/B_VIC

nextflow run main.nf -profile B_YAM \
  --input_dir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/B_YAM \
  --outdir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/results/B_YAM

nextflow run main.nf -profile H5NX \
  -w /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/work \
  --input_dir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/H5NX \
  --outdir /Volumes/Lei_work/Projects/FLU/FLUComprehansiveData/results/H5NX
