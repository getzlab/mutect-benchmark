# GCloud commands

use `--zone` everytime or there will be additional confirmation.

* batch start/stop

  ```bash
  gcloud compute instances stop gce-worker{1..6} --zone us-central1-a
  ```

* clone by image

  ```bash
  gcloud compute instances create `eval echo gce-worker{1..$n_nodes}` --image $snap --preemptible --machine-type n1-standard-8 --zone us-east1-d --metadata-from-file startup-script=provision.sh || { echo "Error instantiating nodes!"; exit 1; 
  ```

* resize, only one instance at a time

  ```bash
  gcloud compute instances set-machine-type INSTANCE --machine-type n1-standard-1 --zone us-central1-a
  ```

  I am using `n1-stanrd-4` for Mutect2 and `n1-highmem-16` for GenomicsDBImport

# Slurm

* `slurmctld -c` hard reset, and restart by `slurmctld -f slurm.conf`
* `squeue | grep xxx`
* `scancel -n JOB_NAME`
* `scancel JOB_ID`

