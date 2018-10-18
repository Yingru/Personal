# Working_diary_2018-08

[TOC]

## 2018-07-30

#### 1.thermal boltzmann distribution
```python
def Boltzmann_E(E, temp, M):
    dfdE = np.exp(-E/temp) * np.sqrt(E**2 - M**2) * E
    normed = dfdE.sum()*(E[1] - E[0])
    return dfdE/normed
    
def Boltzmann_p(p, temp, M):
    E = np.sqrt(p**2 + M**2)
    dNdp = p/E * dNdE
    normed = dfdp.sum() * (p[1] - p[0])
    return dNdp/normed

# find the effective temperature of a system with equilibirum distribution
# being the boltzmann distribution
def effective_temperature(p0_list):
    from scipy.optmize import curve_fit
    dNdp0, bins = np.histogram(p0_list, bins=50, normed=True)
    p0_bins = 0.5*(bins[1:] + bins[:-1])
    def func(E, temp):
        return thermal_dfdE(E, temp, 1.3)
        
    popt, pcov = curve_fit(func, p0_bins, dNdp0)
    return popt  # array of parameters
```

## 2018-07-31
#### 1. comparing the diffusion induced radiation and absorption process (emitted/absorpted gluon energy): $k_0 < \mu \pi T$, where $\mu=0.1, 1$
    
#### 2. getting started with [Chameleon](https://www.chameleoncloud.org/), which is a configurable experimental enviroment for large-scale cloud research.
- [tutorial](https://chameleoncloud.readthedocs.io/en/latest/getting-started/index.html)
    - create an account -> create or join a project -> start using `Chameleon` 
    - log in through [UChicago](https://chi.uc.chameleoncloud.org/) 
    - or log in through [TACC](https://chi.tacc.chameleoncloud.org/)
    - reserve a node (when you create a reservation for one ot more nodes, only you and other users on your project will be able to use those nodes for the time specific), specify the name, and leasing time
    - launch bare-metal instance, specify name, source (OS images), flavor, key pair, lauch process including powering-up, loading the operating system, and booting up
    - associating an IP address (click a few button)
    - access your instance: 
    ```bash
    # change the permission on the file to user read/write only
    chmod 600 yx59chameleonkey.pem
    # add the key to your current SSH idenetiy
    ssh-add yx59chameleonkey.pem
    # log in to yout Chameleon instance via SSH using usrname@ip_adress
    ssh cc@ip_address
    ```
    - <font color='red'> `ssh ` failed, submitted a ticket 
      solved</font>: `ssh cc@ip_address`, the usrname by default is `cc`
     `Orchestration`: work with Chameleon's complex applications
- a few concept:
    - `instance`: virtual machine but on bare-metal nodes
    - `SU`: Service Unites, differs from traditional HPC or cloud service that charged in core-hours

- GUI interface [link](https://chameleoncloud.readthedocs.io/en/latest/technical/gui.html)
- command line interface [link](https://chameleoncloud.readthedocs.io/en/latest/technical/cli.html)
    - [ ] to be read
    
#### 3. write a container for HIC_HQ 
    - [x] trying `docker` container first
    - [ ] then `singularity` container
    - [ ] in Nersc, you can run docker container with `shifter` I believe


## 2018-08-01
#### 1. HQ meeting slides: [Diffusion induced radiation/absorption](./slides/LGV-vs-Lido_2018-08-01.pdf) 
> origin: /var/phy/project/nukeserv/yx59/summer2018/Langevin-vs-LBT/Scattering/build/src/Compare_between_LGV_Lido2/summary_0730-2018

<font color='red'> 

- [ ] figure out the large value of $\alpha_s$ in our *Lido* framework
- [ ] compre the momentum dependence of the two calibration, we should at least figuring out why are they different? how can they be different
- [ ] starting from comparing the prior range of the two parametrization. See if it is possible for both of them to have the similar prior.

</font>

#### 2. Keep testing `Chameleon`, updating the lease by one extra day
```bash
ssh cc@192.5.87.76

# check ubuntu version
lsb_release -a

# install docker and its dependencies
wget https://download.docker.com/linux/ubuntu/dists/trusty/pool/stable/amd64/docker-ce_17.03.2~ce-0~ubuntu-trusty_amd64.deb
sudo apt-get install libsystemd-journal0
sudo dpkg -i docker-ce_17.03.2~ce-0~ubuntu-trusty_amd64.deb

# testing the installation
sudo docker run hello-world

```

#### 3. update the [hic_HQ-osg](https://github.com/Yingru/hic_HQ-osg) git repository, now update all the url and clean the directory a bit

#### 4. write a container for HIC_HQ:
    - Before **Dockerfile**, I can first of all, write an **install_software.sh**, which could be much easier (testing the libraries requirements bit by bit)
    - write the **Dockerfile** in [github/hic_HQ](https://github.com/Yingru/hic_HQ)
    ```bash
    # build the container image
    sudo docker build -t hic_hq:v1 .
    ```
    - okay, find some drawback in the previous workflow, it is not easy to work with the container thing. (I should keep all the executables in one /bin folder, and all the parameters_file in the localhost)


## 2018-08-02
#### 1. update the previous workflow (modify the **CMakeLists** file), now updated to the git repository in [**container** branch](https://github.com/Yingru/hic_HQ-osg/tree/container)
#### 2. write the **Dockerfile**, and successfully run with `Docker`
```
# build the docker image
sudo docker build -t hic_hq:v1 .
cd workdir/

# to run the executable
sudo docker run -v `pwd`:/tmp/hic_HQ-osg/results hic_hq:v2 \
python3 run-events_cD.py args.conf 0
```
#### 3. Some useful lesson learned for `git submodule` (which is painful from the very beginning)
```
# of couse clone with submodule
git clone --recursive git_with_submodule_url

# if you made changes, remember to commit the changes in each of the submodule 
# and  push to remote
# and then go back to the parent repository, and add the commits again

# if you switched a different branch 
git submodule update
```
#### 4. reading the `Slurm` documenetion, goal for this week (use the container in NERSC)
    
#### 5. install `Singularity` on `Chameleon`
[tutorial](https://www.sylabs.io/guides/2.5.1/user-guide/quick_start.html)
```bash
# singularity dependencies
sudo apt-get install libarchive-dev python dh-autoreconf build-essential

# install the maste branch
git clone https://github.com/singularityware/singularity.git
cd singularity

git checkout vault/release-2.5

./autogen.sh
./configure --prefix=/usr/local
make
sudo make install

# remove an old version

```

## 2018-08-03
to do:
- [x] At least make a singularity image work 
- [ ] try `Slurm` in NERSC
- [ ] strentch goal, try `shifter` on NERSC


#### 1. `Singularity` practice
- the easiest way to do right now, is just convert singularity image from docker image <font color='red'> not working right now :( </font>
<font color='blue'> possible solution: do not trust /tmp folder, need to prove that </font>
```bash
sudo docker run \
-v /var/run/docker.sock:/var/run/docker.sock \
-v `pwd`:/output \
--privileged -t --rm \
singularityware/docker2singularity \
hic_hq:v1
```

- Or build singularity image from dockerhub <font color='red'> still cannot access /bin folder </font>
<font color='blue'> solutio: do not trust tmp/ folder! need to improve this part </font>
```bash
# as in previous days, installed singularity in Chameleon 
singularity --version # check version, currently 2.5.2
sudo docker images # list docker images

sudo apt-get update && sudo apt-get install squashfs-tools 
singularity pull docker://yingruxu/hic_hq
singularity 
```
    
- Or write a **Singularity** recipe from scratch. [Singularity recipe](https://github.com/Yingru/hic_HQ/blob/master/Singularity)
```bash
sudo singularity build hic_hq.simg Singularity
```


## 2018-08-04

#### 1. start a new lease at `Chameleon`, and update a new ip_address
- How are external file systems and paths handled in a Singularity Container?
> Because Singularity is based on container principals, when an application is run from within a Singularity container its default view of the file system is different from how it is on the host system. This is what allows the environment to be portable. This means that root (‘/’) inside the container is different from the host!
> 
> Singularity automatically tries to resolve directory mounts such that things will just work and be portable with whatever environment you are running on. This means that /tmp and /var/tmp are automatically shared into the container as is /home. Additionally, if you are in a current directory that is not a system directory, Singularity will also try to bind that to your container.
> 
>There is a caveat in that a directory must already exist within your container to serve as a mount point. If that directory does not exist, Singularity will not create it for you! You must do that. To create custom mounts at runtime, you should use the -B or --bind argument:
> `singularity run --bind /home/vanessa/Desktop:/data container.img`
> 

- Okay, caveat (not to use the system defined binding path), they will be overwritten.

- use `var/hic_HQ-osg` and finally able to run singularity container
    - `sudo singularity shell --writable -B $PWD:/var/hic_HQ-osg/results hic_hq_v1.img`
    - now the second issue is that this is a read-only image. But sure I can figure out a way to make it writable.



## 2018-08-06

#### 1. Loggin to NERSC, and test `Slurm`
```
ssh yx59@edison.nersc.gov   # or ssh yx59@cori.nersc.gov

# first of all, clean a bit of the conda env
conda env list  # simply delete the env folder

# nersc modify all the env in .bashrc.ext file

# install miniconda and and it to the path
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
export PATH="$HOME/miniconda3/bin:$PATH"
export LANG='en_US.UTF-8'

# create conda environment
# conda create --prefix  /global/common/software/<project>/<prefix> 
# (--name, name is not allowed with path)
conda create -p /global/common/software/m2730/hic_yx numpy scipy cython h5py
source activate /global/common/software/m2730/hic_yx 

```

#### 2. on NERSC, testing [hic-eventgen](https://github.com/morelandjs/hic-eventgen) locally
> one extra tips for nersc, by default, the color options are disabled. 
> but you can use `ls --color` to check the color option
> or alias command (for some reason, doesn't work for adding to .bashrc.ext file)
> `alias ls="ls --color"`

```
git clone --recursive https://github.com/morelandjs/hic-eventgen
cd hic-eventgen/nersc
./install  # make sure you are in conda env

# testing a job locally
run-events --nevents 10 --logfile $CSCRATCH/logs/local_pPb2.log --checkpoint $CSCRATCH/checkpoints/local_pPb2.pkl $CSCRATCH/events/local_pPb2.dat
```


- reading in the **results_file**: $CSCRATCH/events/local_pPb2.dat 
> It is stored in a binary format, and structured array (in my workflow, we stored as hdf5 format -- which I think is better)
```
import numpy as np
species = [
        ('pion',     211),
        ('kaon',     321),
        ('proton',  2212),
        ('Lambda',  3122),
        ('Sigma0',  3212),
        ('Xi',      3312),
        ('Omega',   3334),
]



# specify data type, to ensure consistenty across all the machines
float_t = '<f8'
int_t = '<i8'
complex_t = '<c16'


# special data types
nharmonic = 8
mean_pT_t = [('N', float_t), ('pT', float_t)]
flow_t = [('N', int_t), ('Qn', complex_t, nharmonic)]


dtype=[ ('trigger', (float_t, 2)),  \
        ('init_entropy', float_t),  \
        ('nsamples', int_t),        \
        ('dNch_deta', float_t),     \
        ('dET_deta', float_t),      \
        ('mean_pT', mean_pT_t),     \
        ('iden_dN_dy', [(s, float_t) for (s, _) in species]),   \
        ('iden_mean_pT', [(s, mean_pT_t) for (s, _) in species]),   \
        ('pT_fluct', [('N', int_t), ('sum_pT', float_t), ('sum_pTsq', float_t)]),   \
        ('flow', [('alice', flow_t), ('cms', flow_t)] ) \
        ]

events = np.fromfile('local_pPb2.dat', dtype=dtype)

nch = events['dNch_deta']
print('nch: ', nch)
```

#### 3.  Submit `Slrum` job, and check the status 

- [simple](https://github.com/Yingru/hic-eventgen/blob/master/nersc/examples/simple)

```
source activate /global/common/software/m2730/hic_yx

sbatch simple
squeue -u yx59  # check the status
sqs -u yx59    # also check the status and give you more information

```

- [design](https://github.com/Yingru/hic-eventgen/blob/master/nersc/examples/design), with **design-wrapper**
```
sbatch design

```


#### 4. On Nersc, testing [hic_HQ-osg](https://github.com/Yingru/hic_HQ-osg) locally 
- change the SHELL in cori from `zsh` back to `bash`
  Use NIM to change your default login shell. Login, then select **Change Shell** from the **Action** pull-down menu
  
- <font color='red'> previous I have some issue of compiling the source codes </font> <font color='blue'> now fixed! </font>
```
# create an virtual env
conda create -p /global/common/software/m2730/hic_HQ-yx numpy scipy h5py cython
source activate /global/common/software/m2730/hic_HQ-yx

# testing model by model
# bmsap
module swap PrgEnv-intel PrgEnv-gnu  # load gnu c, c++, fortran compiler
                                     # current version 7.3.0
  
module avail boost
module load boost/1.63 
cmake -DCMAKE_INSTALL_PREFIX="/" ..

# make sure you are in cori!
export PATH="$HOME/miniconda3/bin:$PATH"  
# for some reason, cori does not allow me to load path
# diffusion, which includes hdf5
# currently only work with cori (as hdf5 is not installed in Edison)
module load hdf5   # use cc, CC, ftn for compiler wrapper

module load gsl #(cori-2.1, edison-2.3)
cmake -DCMAKE_INSTALL_PREFIX="/" ..
# 
```
- for `diffusion, vishnew`, which need hdf5-fortran library. Only cori is working right now. As edison only has installed `cray-hdf5`. you need to load `hdf5` and use `gfortran` compiler explicity

- [compiler wrapper](https://software.ecmwf.int/wiki/download/attachments/68163173/Compilers%20and%20libraries.pdf?api=v2
), [cross plat form compiler](https://www.hpc.kaust.edu.sa/sites/default/files/files/public//KSL/150520-User_Workshop/KSL_ProgEnv.pdf)
    - application that will run in parallel on Cray XC should be compiled with standard language wrappers
    - `cc` wrapper for C compiler, `CC` wapper for C++ compiler, `ftn` wrapper around the Fortran compiler
    - do not use `crayftn, craycc, ifort, icpc, gcc, g++ ...` unless you want an executable for the login node
    - with `GNU automake` or `cmake`, add specifier `-host=x86_64-unknown-linux-gnu` for the configure tool
    - with `cmake`, provide the `CMAKE_SYSTEM_NAME` and the used compilers in   toolchain file or when incovling cmake
    `cmake -DCMAKE_SYSTEM_NAME=linux \`
    `-DCMAKE_C_COMPILER=cc \ `
    `-DCMAKE_CC_COMPILIER=CC \ `
    `-DCMAKE_FORTRAN_COMPILER=ftn `


## 2018-08-07
<font color='red'> 
todos:

- [ ] [linux comand line](http://linuxcommand.org/index.php)
- [x] need to wrap up with the NERSC stuff, it takes too loog
- [ ] finnaly start with something physical
- [ ] [slurm totorial](https://support.ceci-hpc.be/doc/_contents/QuickStart/SubmittingJobs/SlurmTutorial.html)
- [ ] improve container by specifying the library version

</font>

#### 1. Keep testing the hic-eventgen framework with [taskfarmer](https://github.com/Yingru/hic-eventgen/blob/master/nersc/slurm/design-taskfarmer)
- `export THREADS=32`, don't really know what's that for
- just keep in mind, [`taskfarmer`](http://www.nersc.gov/users/data-analytics/workflow-tools/taskfarmer/) requires a full compute node, so you will need to request minimum 2 nodes in your batch script. The total amount of time you request should be equal to (numer of taks*task time)/$THREADS, not simply the time required to run one task.
- design-test: edison-10142208, taskfarmer-test: edison-10142255

#### 2. wrap-up the compiler for **hic_HQ-osg**
- choose `gnu` compiler instead (which works for all the old fortran codes)
`ftn` for some reason cannot properly compile hdf5_fortran

- a clean start:
(for some reason, cmake find_hdf5 cannot find the HDF5 library in my conda env in diffusion (maybe need to re-write the diffusion CMakeLists.txt))
the actual work is first find hdf5 in my conda env `/hic_HQ-yx/lib/hdf5` and then the hdf5 in /usr/common/software/hdf5`

    - for trento => diffusion
    ```
    ssh yx59@cori.nersc.gov
    alias ls="ls --color"
    source activate /global/common/software/m2730/hic_HQ-yx

    module swap PrgEnv-intel PrgEnv-gnu
    export CC=gcc CXX=g++ FC=gfortran
    module load boost hdf5

    cd trento
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX="/" ..
    make
    
    cd diffusion
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX="/" ..
    make
    ```

- <font color='blue'> Finally, it compiles! </font>
- Now need to modify the workflow ...

#### 3. Let's do some physics!
- make some comparison between the Langevin and Lido prior parameters range
- one more thing I need to check is, whether what I did for the $\alpha_s$ interpolation is correct?


## 2018-08-08

#### 1. Reading Duke-PHSD-Nantes paper
    Updated the latest version on dropbox
#### 2. Keep doing some physics? Which is better? Interpolate $\alpha_S$ first or interpolate $D2\pi T$ first? Let's try it out?
- compare two types of interpolation: 
    - first: interpolate $D2\pi T$ first and interpolate $\alpha_S$ later
    - second: interpolate $\alpha_S$ first, and calculate $D2\pi T$
    - even though they do not differ much, the second way should be the more correct one
    - <font color='red'> testing `curve_fit`? (ValueError?) </font>
    - <font color='blue'> fixed, dimension converting problem! use `func` to calculate `ydata` instead of directly use `qhat`. </font>



## 2018-08-09
Todos:
- [x] <font color='red'> starting from comparing the prior range of the two parametrization. See if it is possible for both of them to have the similar prior. </font> </font color='blue'> checked, it is not possible </font>
- [ ] at two another calibration. One is with the similar format of the Lido model, and at least with the pre-equlibirum stage, as well as the updated medium evolution.
- [x] First thing's first, let's update the medium background!

#### 1. debugging `curve_fit`?
- finally, it is the dimension conversion problem!
[analysis_0807](/var/phy/project/nukeserv/yx59/summer2018/Langevin-vs-LBT/analysis_0807_prior)

- calculate $\hat{q}$ in **Lido** and **LGV**
**LGV**: 
\begin{align*}
(D2\pi T)^{\rm linear} &= \alpha * (1 + \beta * (T / T_c -1)) \\
(D2\pi T)^{\rm pQCD} &= \frac{8\pi}{\hat{q}^{\rm pQCD}/T^3} \\
(D2\pi T) &= \frac{1}{1 + (\gamma^2 p)^2} (D2\pi T)^{\rm linear} + \frac{(\gamma^2 p)^2}{1 + (\gamma^2 p)^2} (D2\pi T)^{\rm pQCD}
\end{align*}

**Lido**:
\begin{align*}
\frac{\hat{q}}{T^3} = \kappa_D(x_D + \frac{1 - x_D}{ET}) + \frac{\hat{q}_{\rm el}}{T^3}
\end{align*}

```python
# LGV
self.D2piT_pQCD = 8*np.pi/self.qT3_pQCD
self.D2piT_linear = qhatMin * (1. + qhatSlope * (self.temp/0.154 - 1.))
factor = (qhatPower**2 * DL.momentum) ** 2
self.D2piT_tuned = 1./(1.+factor) * self.D2piT_linear \
                          + factor/(1.+factor) * self.D2piT_pQCD
self.qT3 = 8*np.pi/self.D2piT_tuned


# Lido
a = f['Boltzmann/{}/rate/scalar'.format(process[0])].attrs
E = np.linspace(a['low-0'], a['high-0'], a['shape-0'])
T = np.linspace(a['low-1'], a['high-1'], a['shape-1'])
kappa_zT3 = np.zeros([len(E), len(T)])
kappa_TT3 = np.zeros([len(E), len(T)])
for pro in process:
    dNdt = f['Boltzmann/{}/rate/scalar/0'.format(pro)].value
    dpzdt = f['Boltzmann/{}/rate/vector/3'.format(pro)].value
    dpx2dt = f['Boltzmann/{}/rate/tensor/5'.format(pro)].value
    dpz2dt = f['Boltzmann/{}/rate/tensor/15'.format(pro)].value
    kappa_zT3 += (dpz2dt - dpzdt**2 / dNdt)/T**3
    kappa_TT3 += dpx2dt/T**3

kappa_diffusion = A + B/np.outer(E, T)
qT3 = 2 * (kappa_TT3 + kappa_diffusion)


```
- Question? For some reason, the **Lido** table is generating something different from what I observed previously? <font color='red'> need to confirm with Weiyao </font> <font color='blue'> confirmed, two modifications on *matrix_element.cpp* and *initialize_mD* </font>
- What does it mean by comparablilty? Right now the format and shape of **LGV** qhat is different from Lido one, but how is that related to the fact we are suppose to compare to calibrated results from different formulism? Is it necessary to have exactly same format to verify the validity of comparison. What is the prior the two calibration is different, how are suppose to justify that?

- update the medium background of the Langevin model.
- But my ultimate question is, how can we distinguish between two types of energy loss if they both describes similar Raa and v2, is there anyway for the machine learning techniques to distinguish between:
    - collision vs. collisional + radiative
    - improved Langevin vs. Lido


#### 2. update the scattering submodule. (now I basically understand, creating another branch to track the scattering case of hic_HQ-osg)
```
git clone --recursive git@github.com:Yingru/hic_HQ-osg.git
cd hic_HQ-osg
git branch
git submodule update

git branch Scattering_2018summer
git checkout Scattering_2018summer

git submodule add git@github.com:Yingru/Duke-Lido.git models/Duke-Lido
mv Duke-Lido Scattering
git checkout run_cD_2018summer

git submoduke add https://github.com/keweiyao/fortranformat.git models/fortranformat

git add *changes*
git commit -m 'for the scattering but with the old framework'

git push origin Scattering_2018summer
```


## 2018-08-10
todos:
- [ ] new calibration with new medium background as well as new parameterization?
- [ ] question? What is the difference between a Langevin dynamics and a LBT dynamics?

1. OSG: submit calibrated LGV and LBT comparison for LBT-collision only case. Goal, use any kind of classifier to distinguish those two models.
2. write a `binder` to do the prior, posterior fitting
3. update the `hic_HQ-osg`, to include the freestreaming version
```
git checkout container
git submodule update 
git branch freestreaming
git checkout freestreaming
git submodule update


# update the freestreaming submodule
git submodule add https://github.com/Duke-QCD/freestream.git models/freestreaming
git add .gitmodules
git add models/frestreaming
git commit -m 'add frzout submodule'
git push origin freestreaming   
# note, there are some problem in OSG that https is forbidden, you have to update the url every time to push some. pull is okay

# add the frzout submodule
git submodule add https://github.com/Duke-QCD/frzout.git model/frzout



# update trento
git remote add Duke-QCD https://github.com/Duke-QCD/trento.git
git fetch
git checkout remotes/Duke-QCD/master
git branch master_fs
git checkout master_fs 
git push origin master_fs

# to build trento requires
# cmake>-=3.4, boost>=1.50, hdf5
module load cmake/3.8.0 boost/1.62.0-cxx11  hdf5/1.8.20-cxx11
module list

# load library path for cmake
export CMAKE_INCLUDE_PATH=$CPATH
export CMAKE_LIBRARY_PATH=$LIBRARY_PATH

```


## 2018-08-13
todos: 
- [ ] updating the freestreaming freework, including a new trento, a new HQ production, a new sampler
- [ ] strength goal, updating the HQ hadronization process
- [ ] further testing `slurm` on NERSC, which part I did is wrong, and it is not working?

1. updating the `trento` in the `freestreaming` branch
2. debugging why the `slurm` submitting job is not working?
```
# run interactively, takes sometime for waiting for resources
salloc -N 1 -q debug -C haswell -t 00:30:00 -L SCRATCH
srun run-events # works

# okay, let me change the design file instead (design-wrapper is not a good idea)
```
    - some notes on the nersc job: 
    - 1 node, 1 CPU per task => 64 cpus per node, 64 task per srun

3. Lesson learned, giving up following Jonah and Scott's script, instead of creating your own! this is much easier!


## 2018-08-15
#### 1. Finally, a working version of `hic-eventgen` on NERSC!
- 1.1 [`taskfarmer tutorial`](http://www.nersc.gov/users/data-analytics/workflow-tools/taskfarmer/)
```
#====wrapper.sh================
cd $CSCRATCH/taskfarmer   # important! taskfarmer landed in /tmp
                          # need to cd to directory
python calBinary.py $1 $2 $3

#====tasks======================
bash wrapper.sh 0 1 2
bash wrapper.sh 1 2 3
bash warpper.sh 2 3 4

#====Batch script: batch.sl======
#!/bin/bash
#SBATCH -N 2 -c 64
#SBATCH -p debug
#SBATCH -t 00:05:00
#SBATCH -C haswell

cd $CSCRATCH/taskfamer
export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=32
runcommands.sh tasks.txt
#==================================

# to submit taskfarmer job
module load taskfarmer
sbatch batch.sl

```

- 1.2 to submit tasks
    - files: [*design-wrapper*](https://github.com/Yingru/hic-eventgen/blob/master/nersc/pPb_yx/design-wrapper), [*generate_tasks*](https://github.com/Yingru/hic-eventgen/blob/master/nersc/pPb_yx/generate_tasks), [*batch.sl*](https://github.com/Yingru/hic-eventgen/blob/master/nersc/pPb_yx/batch.sl), [*slurm_design*](https://github.com/Yingru/hic-eventgen/blob/master/nersc/pPb_yx/slurm_design)
    - to run a single task locally
    `bash design-wrapper $inputfile $taskID`
    - to submit a slurm single job
    `sbatch slurm_design`
    - to submit taskfarmer job
    `bash generate_tasks`
    `module load taskfarmer`
    `sbatch batch.sl`
    `sqs -u yx59`
    
    

#### 2. Continue modifying the `hic_HQ`, adding the pre-equlibrium stages in the module
- updated: freestreaming (added freestreaming process)
- updated: diffusion evolution (now able to read different HQ inputs, and different medium background, such that it can do evolution during pre-equlibirum stages)
- updated: UrQMD (now able to read in new frzout format, and run lightOnly case)

#### 3. Testing the docker and singularity container on *Chameleon*, I think mine works for now (just get refresh of what happened)
following the github [hic_HQ](https://github.com/Yingru/hic_HQ)
```
# a few more notes on python install
 sudo apt-get install python-dev python-pip libhdf5-dev
 sudo apt-get install numpy scipy h5py
```


## 2018-08-16
#### 1. Keep testing on `singularity`
[tutorial-2.5.2](https://www.sylabs.io/guides/2.5.2/user-guide/singularity_and_docker.html#tldr-too-long-didnt-read)
```
# a few more notes
# convert a singularity image to writable singularity image
singularity shell -B $PWD:/var/hic_HQ-osg/results hic_hq_write.img
# invoke to current directory?

sudo singularity shell -B $PWD:/var/hic_HQ-osg/results hic_hq_write.img
# invoke to /root?

sudo singularity exec --writable -B $PWD:/var/hic_HQ-osg/results hic_hq.img /var/hic_HQ-osg/bin/trento Pb Pb 10 --output initial.hdf5

sudo singularity exec --writable -B $PWD:/var/hic_HQ-osg/results hic_hq.img  python3 /var/hic_HQ-osg/results/run-events_cD.py /var/hic_HQ-osg/results/args.conf 0
```

> Question, what is the difference between `docker` and `singularity` and `shifer`?
> a very good explanation: http://geekyap.blogspot.com/2016/11/docker-vs-singularity-vs-shifter-in-hpc.html
> 

#### 2. Add args.parse on the new workflow
- testing the new workflow, the first version seems work
```
# work locally
git clone --recursive 

# python3.6 is only available in some PPA
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt-get update
sudo apt-get install python3.6
```

#### 3. Meething with Chameleon people
1:00 - 2:00, meeting with chameleon Paul and Cong, and they will help me to set up with the `condor` on chameleon stuff

## 2018-08-17
#### 1. testing new workflow on ubuntu14.04 instance? 
after rebuild the instance (with the same ip-address), I got the error as:
```
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@    WARNING: REMOTE HOST IDENTIFICATION HAS CHANGED!     @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
IT IS POSSIBLE THAT SOMEONE IS DOING SOMETHING NASTY!
Someone could be eavesdropping on you right now (man-in-the-middle attack)!
It is also possible that a host key has just been changed.
The fingerprint for the RSA key sent by the remote host is
51:82:00:1c:7e:6f:ac:ac:de:f1:53:08:1c:7d:55:68.
Please contact your system administrator.
Add correct host key in /Users/isaacalves/.ssh/known_hosts to get rid of this message.
Offending RSA key in /Users/isaacalves/.ssh/known_hosts:12
RSA host key for 104.131.16.158 has changed and you have requested strict checking.
Host key verification failed.
```
The reason is the finger-print for the host has changed,
you can actually delete the associated ip-address, and regen the know_host again
```
ssh-keygen -R ip_address_associated_with_host/host_name
ssh cc@ip_address # it will ask you to re-add the known-host
```

> problem: Ubuntu14.04 only have python3.4 in the repository, but frzout requires 3.5+ python?
```
# https://askubuntu.com/questions/865554/how-do-i-install-python-3-6-using-apt-get
# ubuntu14.04 repostory does not have python3.5+ registered
# use J Fernhough's PPA
# CAUTION: do not under any circumstances remove python3.4 by `sudo apt remove python3.4`
# python is fundamentally baked into ubuntu, and you could break your ubuntu install
# if you want python3 to map to python3.6, create a symlink instead!
sudo add-apt-repository ppa:jonathonf/python-3.6
sudo apt-get update
sudo apt-get install python3.6
# link python3 to python3.6
rm /usr/bin/python3
ln -s /usr/bin/python3.6 /usr/bin/python3

```


> Err, use `conda` instead
```
# install conda in silent mode
wget https://repo.continuum.io/miniconda/Miniconda3-4.5.4-Linux-x86_64.sh -O ~/Miniconda3.sh
bash ~/Miniconda3.sh -b -p $HOME/miniconda3
export PATH="$HOME/miniconda3/bin:$PATH"

# in principle, you can create a virtual env, but I am not sure how virtual environment is functional in container (let's leave that to be fixed?)
```
> install gnu compiler 6.4+ in ubuntu 14
> ubuntu 14 is using gnu 4,8 (cause issue while compiling)
> https://gist.github.com/application2000/73fd6f4bf1be6600a2cf9f56315a2d91
> 

> cmake install 3.4 on Ubuntu14:
> https://askubuntu.com/questions/610291/how-to-install-cmake-3-2-on-ubuntu?rq=1
```
sudo apt-get install build-essential
wget http://www.cmake.org/files/v3.4/cmake-3.4.0.tar.gz
tar xf cmake-3.4.0.tar.gz
cd cmake-3.2.2
./configure
make 
export PATH="`pwd`/cmake-3.4.1-Linux-x86_64/bin:$PATH" 

```


### 3. try on Ubuntu 16 <font color='red'> Worked! </font>
```
# gnu compiler is 5.4
# 1. Install cmake 3.4
sudo apt-get install build-essential
wget http://www.cmake.org/files/v3.4/cmake-3.4.0.tar.gz
tar xf cmake-3.4.0.tar.gz
cd cmake-3.4.0
./configure
make 
export PATH="`pwd`/bin:$PATH" 

# 2. install libhdf5 (1.8.16)
sudo apt-get install libhdf5-dev libhdf5-serial-dev

# 3. install boost (1.58.0)
sudo apt-get install libboost-all-dev


# 4. install gsl (2.1)
sudo apt-get install gsl-bin libgsl-dbg libgsl-dev

# 5. install python3.5 
sudo apt-get install python3-dev liblapack-dev libatlas-base-dev
sudo apt-get install python3-numpy python3-scipy python3-h5py
sudo apt-get install cython3
sudo apt-get install -y python3-setuptools 

# 6. upgrade scipy
sudo apt-get install python3-pip
sudo pip3 install --upgrade scipy
```


## 2018-08-20 : 2018-08-21
#### 1. Try new framework on Ubuntu 14? 
<font color='red'> forget about it, unless I can upgrade python3 (which is linked to python3.4) to python3.5, there is barely anything I can do.

Now have a walk around with `activate conda-env` while compiling python packages, and `deactivate conda-env` while done...
</font>

- reserve node -> launch instances -> (question: what is the different between maximum host =3 and maximum host = 1?)
- loggin into instance
```
ssh cc@192.5.87.178
git clone --recursive https://github.com/Yingru/hic_HQ-osg.git
cd hic_HQ-osg/
git checkout freestreaming
git submodule init
git submodule update

## update and install necessary library
sudo apt-get install update

## install cmake 3.4
sudo apt-get install build-essential
wget http://www.cmake.org/files/v3.4/cmake-3.4.0.tar.gz
tar xf cmake-3.4.0.tar.gz
cd cmake-3.4.0
./configure
make 
export PATH="$PWD/bin:$PATH"

# install libhdf5 (.8.11)
sudo apt-get install libhdf5-dev libhdf5-serial-dev

# install boost (1.54.0)
sudo apt-get install libboost-all-dev

# 4. install gsl (1.16)
sudo apt-get install gsl-bin libgsl0-dbg libgsl0-dev libgsl0ldbl

# 5. python3.6
# there are two options
# 5.1 if you want to use conda environment
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"

conda create --name py36 python=3.6
source activate py36
conda install numpy scipy h5py cython

# for some reason, the CMake cannot find the correct HDF5 library in a conda environment, therefore rearange the PATH
# export PATH="/usr/bin:$PATH"
# no, it does not work!! because in that way, python3 goes back to the system default python3
# wait, is it because I installed hdf5 instead of h5py?

# or a walking around is:
# activate the environment only when compiling the python setup.py
# which needs some hard-coded script (not sure if it will eventually work)
# yes it actually worked!

# 5.2 If you want to upgrade python3.4 => python3.5 (tried but failed, python3 is more baked in ubuntu than you may think, there are a lot of other libraries such as setup-tools not compatable with python3.5)

```

### write a new `Dockerfile` (with Ubuntu 16)
- best practice in write `Dockerfile` (https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#add-or-copy)

Because of image size matters, using `ADD` to fetch packages from URLs is strongly discouraged, you should use `curl` or `wgrt` instead. 
```
RUN mkdir -p /var/cmake \
    && curl -SL http://www.cmake.org/files/v3.4/cmake-3.4.0.tar.gz \
    | tar xf cmake-3.4.0.tar.gz /var/cmake
```

### try `shifter` on NERSC
- `shifter` image name: `source:image_name:image_version`, typically the source is docker and image name and version are defined by user.

- `root` is squashed on shifter images, so the software should be installed in a way that is executable to someone with user-level permissions. Additionally, iamges are mounted read-only NERSC, so software should be configured to output to NERSC filesystem, like $SCRATCH or project. 
```
# pull image down on Cori from dockerhub
shifterimg -v pull docker:yingruxu/hic_hq:latest
shifterimg images | grep 'yingruxu'
salloc -N 1 --image yingruxu/hic_hq -C haswell
```

**the /etc and /var directories are reserved for use by the system and wil be overwritten when the image is mounted

project must be accessed in a shifter image by using its full path /global/project instead of just /project**


## 2018-08-22:2018-08-24
Todos:
- [x] finalize the workflow of container, make it not writable
- [x] try the container in NERSC?
- [ ] try the slurm coupled with shifter
- [ ] PHSD project (need to run more events), and modify the paper
- [ ] I think now it is also not a bad time to finalize what I have learned from the different comparison between two types of energy loss 

### 1. processing the OSG data
- first, compare if the container works the same as the previous results, `run95_reproduce_run90`
- second, compare the LGV-LBT results
`run94_LGV-vs-LBT`, `run93_LGV-vs-LBT`

### 2. finalize the workflow
```
# now all the binary installed in $pkgname/bin, all the non-change data installed in $pkgname/share, all the configure files in $pkgname/results (which is expected to be changed for user case)
export PATH=$pkgname/bin:$PATH
export PYTHONPATH=$pkgname/lib/python*/site*....
export XDG_DATA_HOME=$pkgname/share:$XDG_DATA_HOME

# python can access the XDG standards by
os.environ.get('XDG_DATA_HOME')

# dockerfile mount to /usr/local/hic_HQ-osg/results
sudo docker run -v `pwd`:/usr/local/hic_HQ-osg/results hic_hq:v1 python3 run-events_cD.py args.conf 0
```

### 3. test on NERSC
```
# pull the image deom docker
shifterimg -v pull docker:yingruxu/hic_hq:freestreaming

# check images
shifterimg images | grep 'yingruxu'

# entering the shifter image
shifter --image=yingruxu/hic_hq:freestreaming
echo $PATH
python3 run-events_cD.py args.conf 0

# run the image (just as docker) 
shifter --image=yingruxu/hic_hq:freestreaming python3 run-events_cD.py args.conf 0

# run slurm interactive job
salloc -N 1 -C haswell -p debug --image=yingruxu/hic_hq:freestreaming -t 00:30:00
shifter python3 run-events_cD.py args.conf 4

# run slurm single job
#========== singleBatch.sl ========================
#!/bin/bash
#SBATCH --image=yingruxu/hic_hq:freestreaming
#SBATCH --nodes=1
#SBATCH --partition=regular
#SBATCH --constaint haswell
#SBATCH --cpus-per-task=1
#SBATCH --time 00:30:00
#SBATCH --license=cscratch1

job=$SLURM_JOB_ID
srun shifter python3 run-events_cD.py args.conf $job
```

- note, in `slurm`, specify the filesystem by `-L cscratch1` is able to make the job start only is the filesystem is available, in which way we can prevent some failures or descreased performance
- for submit `shifter` job, you need to specify the licence, as the shifter images are stores in `Cori` global scratch filesystem, and in `Edison` scaatch3 filesytem.

## 2018-08-24
- submitted a job on OSG, only if I can get all the soft and hard results done, and do a new calibration 
    - [x] local shifter image done
    - [x] local design-wrapper done
    - [x] slrum interactively done
    - [ ] slrum a single node, one cpu per task


- how can I get the task ID in the output or error file name for a batch job?
> if you want seperate the output by task, you will need to build a script containing this specification.
> ```
> ## ====== test file =====
> #!/bin/sh
> echo begin_test
> srun -o out_%j_%t hostname
> ## =======================
>
> sbatch -n7 -o out_%j test  # submitted batch job 65541
> 
> ls -l out*
> out_65541, out_65541_0, out_65541_1, out_65541_2 ...
> 
> cat out_65541  # begin_test
> 
> cat out_65542  # tdev2
> 
>
> ```


- [x] step 0: interactively salloc
- [x] step 1: partition=debug full events, shifter no srun
- [x] step 2: partition=regular full events, shifter no srun
- [ ] step 3: partition=regular full events, srun -n 32 shifter ? (not passing?)


- two ideas, use taskfarmer, which does not have srun
- but then hw to proceed??

```
# interactively but with more cores??
#Wait for it to complete
salloc -N 2 -C haswell -p regular -A ntrain --reservation=sc17_shifter --image <mydockerid>/hellompi:latest
# Wait for prepare_compilation_report
# Cori has 32 physical cores per node with 2 hyper-threads per core.  
# So you can run up to 64 tasks per node.
srun -N 2 -n 128 shifter /app/hello

```

## 2018-08-26
1. test locally on NERSC
- `nersc/install` have some problem loading frzout (which I don't understand), should figure it out later
- `install_frzout` do a second run up
- `design-wrapper` load `$PATH, $LD_LIBRATY_PATH, $XDG_DATA_HOME`
- `result_local` run a test run
- `taskfarmer` tested as well, failed even with original workflow (?!) why, is there any problem with taskfarmer?????! I am really pissed.


## 2018-08-27
1. PHSD project, debugging it??
2. processing OSG data?
3. provide the Indian group the information for nucleus?
- `trento/src/collision.cxx`
- `for (auto a in nucleus_A), output a; for (auto b in nucleus_B), output b`
- add random seed generator also for the nucleus position 

```cpp
// ============== nucleus.h ========================
// stores an ensemble of nucleons and randomly samples their positions
// interator interface through: begin(), end()
class Nucleus {
  public: 
      static NucleusPtr create(const std::string& species, double nucleon_dmin=0);
      virtual double radius() const=0;
      void sample_nucleons(double offset);
      size_type size() const noexecept;
      iterator begin() noexcept;
      iterator end() noexecept;
      const_iterator cbegin();
      const_iterator cend();    
}

class Nucleon{
  public: 
      double x() const;
      doubel y() const;
      double z() const;
      bool is_participant() const;
}
```
- generate initial nucleus position A and B
```
./trento Pb Pb 10 --cross-section 7.0 --normalization 18.5 -p 0.0 
--fluctuation 0.9 --nucleon-width 0.96  --nucleon-min-dist 1.280 
--grid-step 0.1 --grid-max 15.05 --output initial_minBias.hdf5


# b=6 (1234), b=2 (12345), b=15 (123456)
./trento Pb Pb 100 --cross-section 7.0 --normalization 18.5 -p 0.0 \
--fluctuation 0.9 --nucleon-width 0.96  --nucleon-min-dist 1.280 \
--grid-step 0.1 --grid-max 15.05 --b-min 8.0 --b-max 8.0 \
--output initial_b8.hdf5 --random-seed 1234
```

## 2018-08-28
1. possible problem for NERSC 
`srun -n 20` 

2. calculate the path towards equilibrium for both Boltzmann and Langevin dynamics
3. `python setup.py build_ext -i` ?? why should i spend so much time to debug the old code???


## 2018-08-30
1. verify the workflow for taskfarmer on NERSC, it seems working right now!
    - create a dockerfile, upload to dockerhub and build docker image
    - pull shifter image from docker to NERSC (note, NERSC cannot delete docker image, therefore, you have to rename the docker image everytime you pull from dockerup)
    ```
    # pull image down on Cori from dockerhub
    shifterimg -v pull docker:yingruxu/hic_hq:latest
    shifterimg images | grep 'yingruxu'
    ```
    - create job-wrapper file, taskfarmer submit file
    ```
    # ======= run-wrapper ===========
    #!/bin/bash

    workdir=$CSCRATCH/PbPb5020_taskfarmer/results_shifter_verify_20180830
    cd $workdir
    shifter python3 $workdir/run-events_cD.py $workdir/args.conf $1

    # ======= generate task ==========
    #!/bin/bash

    ntasks=10

    for (( n = 0; n < $ntasks; ++n )); do
      echo bash run-wrapper $n >> tasks
    done
    
    
    # ====== taskfarmer.sl =============
    #!/bin/bash
    #SBATCH --image=yingruxu/hic_hq:freestreaming
    #SBATCH --nodes=2
    #SBATCH --partition=regular
    #SBATCH --constraint haswell 
    #SBATCH --cpus-per-task=64
    #SBATCH --time 00:30:00
    #SBATCH --license=cscratch1 
    #SBATCH --mail-type ALL
    #SBATCH --mail-user yx59@duke.edu

    export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
    export THREADS=32

    cd $CSCRATCH/PbPb5020_taskfarmer/results_shifter_verify_20180830

    runcommands.sh tasks

    
    ```
    - submit jobs
    ```
    module load taskfarmer
    sbtach taskfarmer.sl
    ```

2. debug the `srun` error? Is is because I have the sameName intermediate files?
    if I can get the Freestream.h5 files back, which would be this is the true
    <font color='red'> I think I could just give up 
    does not fix it at all
    and also, I probably need to use the MICH version??    
    </font>
    
    
3. processing the OSG data (find a bug, now fixed)
4. need to run one set of PbPb 2.76 TeV data and check how good the medium is??