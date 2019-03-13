# Betaduck
For data-handling on the PromethION beta device.

## Uses

### Tidying your data
* Bundle fast5 folds into gzipped tar files (much easier to then rsync via ssh or mount - see the commands below)
* Gzip fastq files and move to folder 'fastq'
* Move sequencing summary files to folder 'sequencing_summary'

![flowchart](images/tidying_flowchart.png)

### Plotting your sequencing data
* Plot datasets using the fastq and sequencing summary files. And provide a stats report. 
* Examples of some  plots produced are shown below.
  + Yield produced over time
![yield_plot](images/example.yield.png)
  + Histogram of read-length distribution
![hist_plot](images/example.weighted.hist.png)
  + Show the distribution in average read quality
![dist_plot](images/example.quality.hist.png)
  + Pair plot, showing the correlation between four variables provided for each read
![pair_plot](images/example.pair_plot.png)
  + Yield per channel as provided in the flowcell map.
![flowcell_map](images/example.flowcellmap.png)

## Installation through Docker.

### Benefits of running through Docker.
* Restrict memory and so that concurrent sequencing runs are not compromised.
* Version controlling for specific versions of MinKNOW
* Prerequisites are installed for you.
  * You will still need to install docker

### Downloading and installing Docker.
There is no need to repeat what the internet has already provided
[Click here](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-16-04) to see a nice straight forward guide to downloading docker.

You may need to log out and log back in again for the super user permissions to take place. If that is not feasible, you may need to precede the docker commands below with `sudo`

### Installing betaduck from docker
`docker pull alexiswl/betaduck:18.07.1-3`  
19.01.1-dev represents the current MinKNOW version.

## Tidying and compressing files

### betaduck config 
First we'll generate the config file. This is used as an input in the next command as to which files and folders to work on.  

For assistance with options use  
`docker run alexiswl/betaduck config --help`

**betaduck config parameters**
* --run-dir
  + Likely to be `/data/<run_id>/<sample_id>/<run_dir>`
* --output_yaml_file
  + Output config file generated. To be used in next script.  
  + Likely to be `/data/<run_id>/<sample_id>/<run_dir>/config.yaml`
* --sanitiser
  + Run the fastq sanitiser script before generating the config file.  
* --active
  + Use if the run is currently still going.
  + Removes the last fast5 folder it finds from the dataframe.

**docker parameters**  
You will also need to bind the /data volume to the container when executing the script.
Example:  
`docker run --volume /data:/data alexiswl/betaduck config ..parameters`


### betaduck tidy 
Most of the hard work has been done for us in the previous script.  

**betaduck tidy parameters**

* --config
  + Path to config file generated in the previous script
* --keep
  + Should we keep the original files
* --dry_run          
  + Don't actually run the subcommands. Just pretend.
* --overwrite        
  + Overwrite files if they already exist
* --threads
  + Number of folders to tar and zip simultaneously (default=1)
  
**docker parameters**  
Here is where docker shines, it can restrict the cpus and memory utilisations of a given container as to not blow up your system.  
Having said that, the PromethION beta device is pretty powerful machine. It has 96 cpus with 400 Gb of memory and all of this on SSD drives.
Use the `top` and `free` commands to view the current utilisations of your system.  
The `--cpus` parameter should be 1 higher than the `--threads` parameter used by betaduck. 
We can also trick the container into thinking the root directory is the run_dir. 
Example:  

`run_id=/data/<run_id>/<sample_id>/<run_dir>`  
`config_yaml=config.yaml`  
`docker run --memory=50g --cpus=7 --volume=${run_id}:/ alexiswl/betaduck tidy --config=${config.yaml} ..parameters`  

## Plotting datasets with betaduck
Now we get our rewards, some plots produced from the seaborn and matplotlib libraries

**betaduck plot parameters**

* --summary_dir 
  + Where the tidy command placed the summary files.
  + Likely `/data/basecalled/<sample>/<flowcell_port>/sequencing_summary/`
* --fastq_dir 
  + Where the tidy command placed the fastq files.
  + Likely `/data/basecalled/<sample>/<flowcell_port>/fastq`
* --name 
  + Name used in plots
  + Keep this to 10 characters or less
     + fig.tight_layout() is still in development.
* --plots_dir 
  + Where do you want to put these plots
  + Likely `/data/basecalled/<sample>/<flowcell_port>/plots`
* --threads
  + Number of threads to use when reading in fastq and summary datasets.
  + Number of threads used when generating the plots has been restricted to 1.

## Aligning data with betaduck
Using minimap2 v2.13 and samtools v.19, betaduck will align the fastq files to your reference genome.
Using the --filter_lambda option, the lambda genome is appended to the reference as an additional chromosome.
Reads that align to the lambda chromosome are then split from the rest of the rads that align to the genome.
## Rsyncing to an external location
Now that your nanopore data is tidy, you can rsync the data across using rsync.
You will need to specify which files to include and exclude in the rsync comamnd.
```
rsync --archive ' \ # Archive allows rsync to find files recursively
--include='*/' \ # Look for the following file endinds in subfolders
--include='*.fast5.tar.gz' \ # Find all the fast5 tars
--include='*.fastq.gz' \ # Find all the gzipped fastq files
--include='*.sequencing_summary.txt' \ # Find all the moved sequencing summary files
--include='*.png' \ # Carry over any plots that have been generated
--exclude='*' \ # Exclude all other files
--prune-empty-dirs \ # Don't download folders that won't have files in them (like reads/0 etc)
/path/to/data/basecalled/reads
/dest/directory
```

### Setting up your genome file
* Create the directory /data/genome on your PromethION.
* All genomes will sit under here with this structure
  * /data/genome/<genome_name>/genome.fa 
    * you may need to rename your reference file to genome.fa
  * minimap2 indexes will automatically be built the first time they're used under:
    * /data/genome/<genome_name>/genome.mmi
  * If the --filter_lambda option is used, two more files will be generated for each genome
    * /data/genome/<genome_name>/genome.w_lambda.fa
    * /data/genome/<genome_name>/genome.w_lambda.mmi 
* Betaduck will automatically download the lambda reference.
  * The lambda reference will be placed in:
  * /data/genome/lambda/genome.fa
  
**betaduck align parameters**
* --fastq_dir=/data/basecalled/<sample_name>/<flowcell_port>/fastq
* --genome_dir=/data/genome (where the reference directories are found)
* --genome=<genome_name> (Subdirectory of the genome_dir parameter)
* --output_dir=/data/basecalled/<sample_name>/<flowcell_port>/alignment
* --filter_lambda (Do you want the lambda sequences removed)
* --cs=long|short|none (Do you want the cs tag for each alignment)
* --md (Do you want the MD tag for each alignment)

This is also acceptable to do whilst the run is still being generated or the betaduck tidy script is running.
fast5.tar.gz files are written as fast5.tar.gz.tmp initially and then moved so they won't be listed by the rsync if they're still being written to.

### Additional options to rsync
`--remove-source-files` if you wish to remove the data from the PromethION device.

## Troubleshooting
* Info logs display in UTC time
  + Docker containers use the time of the container, not that of the outside.
  + You can fix this by mounting the /etc/localtime into the container
  + `docker run --volume /data:/data --volume /etc/localtime:/etc/localtime alexiswl/betaduck ..parameters`
