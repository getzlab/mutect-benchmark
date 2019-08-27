# Google cloud commands

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

  ```
  gcloud compute instances set-machine-type INSTANCE --machine-type n1-standard-1 --zone us-central1-a
  ```

* 

