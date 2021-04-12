## Cross Layer QoS Control
This repository contains materials to reproduce the results for SC 2021 paper "It Takes Two to Tango: Reconciling Application and Storage to Manage I/O Interference on CGroups-based Storage".

## Prerequisite
1. We set up our experiment environment usings docker (v20.10.5). Please refer to [Docker documentation](https://docs.docker.com/engine/install/ubuntu/) for installation instructions. The image we use to create containers is "jrei/systemd-ubuntu:18.04";
2. We use ADIOS as the data management framework to refactor and write scientific data generated from simulations to storage systems. Please refer to the [ADIOS documentation](https://www.olcf.ornl.gov/center-projects/adios/) for installation instructions.

## How to reproduce the results
1. Create volume group and logical volume on HDD on host node which will be shared among all containers;
2. Create containers. Bind mount the HDD volume to each container;
3. We use three Applications for evaluation, you can find the these data sets under "/ADIOS/examples/C/global-array/larger\_data/current\_data/". And the decimation ratio can be set via editing "/ADIOS/examples/C/global-array/test\_xgc.xml" before running "/ADIOS/examples/C/global-array/test\_xgc" to refactor and write reduced base data to SSD and deltas to HDD;
4. Our cross layer qos control method can be run through "qos\_control.py". The parameters you need to set before run are listed below.

| Parameters | Description |
| ------ | ------ |
|application\_name| The name of application|
|deci\_ratio|Decimation ratio|
|time\_interval|The period of data analytics|
|fulldata\_len|The number of the elements of original data|
|reduced\_len|The number of the elements of reduced base data|
|error\_metric|The error control metrics (Can be "NRMSE" or "PSNR")|
|error\_target|The enforced error tolerance|
|priority|The priority of application|
|min\_priority|The minimum priority among all applciations|
|max\_priority|The maximum priority among all applciations|
|min\_accuracy|The lowest accuracy requirement|
|max\_accuracy|The highest accuracy requirement|
|error\_target\_list|All accuracy levels|
|docker\_path\_id|The container id|
|weight\_tag|Different weight function (Can be "size", "size\_priority", "size\_priority\_accuracy")|
5. After setting up all the parameters, you also need to set up the period and intensity of interference in each container. Then we are ready to run data analytics in one container and interference applications (noise.py) in the other containers;
6. The QoS run takes this form:
```sh
python qos_control.py [ERROR_CONTROL_BOOL] [WEIGHT_BOOL] [RUN_START_TIME]
```
ERROR\_CONTROL\_BOOL : The Boolean to specify whether the error control is on or not. (Can be "on" or "off")
WEIGHT\_BOOL: The Boolean to specify whether the dynamic weight adjustment is on or not. (Can be "on" or "off")
RUN\_START\_TIME: The start time of the code. (Can be the format "HOUR MINITE SECOND" or "now")
"

