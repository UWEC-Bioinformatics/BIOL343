# Submitting Snakemake workflows with the Slurm job manager

1. Make sure you have the profile saved in `BIOL343/snakemake-profile`  
    a. This should include a config called `config.v8+.yaml`  
    b. This file includes the default configurations for Slurm submission - do not change it.  
        - ***Note***: the walltime is set to 24 hours. Please let Dr. Wheeler know if you think you'll require more time.  
2. Login to a head node of BOSE (i.e., don't start a VS Code job).  
    a. The easiest way to do this is to navigate to your BIOL343 folder in the OnDemand file explorer and click the button "Open in Terminal".  
3. Navigate to the folder that contains your `snakefile`.  
4. Run the command: `snakemake --profile ~/BIOL343/snakemake-profile/`
  >[!TIP]
  > Running a pipeline this way requires that you keep the SSH open for the entire analysis. To start a background process, run the following instead: `nohup snakemake --profile ~/GitHub/invision-tools/slurm-profile/ &`   
