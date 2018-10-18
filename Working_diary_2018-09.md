# Working_diary_2018-09

## 2018-09-04
1. debugging HQ initialization position??
[freestreaming-HQ](/var/phy/project/nukeserv/yx59/summer2018/run96_PbPb5020_fs/debug_HQPosition)
    - one of the caveat is that, when transforming from `freestreaming` to `hydro` evolution, the system is phrase transforming from `ideal gas` to `QGP`, therefore, although the energy density is continous, the actual temperature has changed quite a lot
    - Okay, figured it out, need to transform `Ncoll` density while sample HQ positions.
    - <font color='red'> correcting again! 
    The Fortran column/row is reversed from python/c++ column/row.
    Therefore you don't need to transform Ncoll in run event
    </font>
    - Okay, HQ evolution is still wrong! I think the biggest problem is for initialize the HQs and evolve HQs inside the medium??
    - option 1: evolve HQs by t: p_r0(i,j) + deltat
2. With `freestreaming`, check the soft hadron observables
3. Now back to my comparison between 'Langevin-vs-Linearized-Boltzmann' project


## 2018-09-05
1. Question: what is the reason that new calculation $v2\{2\}$ is so much larger than $v2(\rm EP)$?
    - debugging? (compare pp_angle participant  vs. pp_angle event-plane)
    - light hadron calculated EP plane? 
    - light hadron calculated EP plane with light hadron cut?
    - how about Dmeson calculated EP plane? what do you think??
    - <font color='red'> finally find out the error! $p_x, p_y$ inversed! </font>

2. Also, it really matters that we should account for all the events that trento generates...
   
## 2018-09-06
1. [x] keep debug $v_2(EP)$, $v_2\{2\}$
    - `condor_run90-posterior0_PbPb5020-debug2`: this is for with Freestream and correct HQ initialization
    - `condor_run90-posterior0_PbPb5020_noFs`: without Freestream and maybe not correct HQ initialization?
    - No light freestream? How about ignore the light quark flow initialization?
    - now `condor_submit` a new job `condor_run90-posterior0_PbPb5020-debug3`, hopefully this time it will make sense
    
    
2. meeting with HF?
    - heavy quark develop flow through initial-state effect or final-state effect?

3. try on NERSC:
```
# clone source file
git clone --recursive https://github.com/Yingru/hic_HQ-osg.git
cd hic_HQ-osg/
git branch
git checkout freestreaming
git submodule init
git submoduke update

# create conda virtual env to install the software (need to specify python3.5 or 3.6) , there are some problem for python3.7
echo $PATH
echo $LANG  # make sure have the correct path and LANG, see last previous notes for install miniconda and setting LANG

conda create -p /global/common/software/m2730/hic_HQ-yx59 numpy scipy cython h5py python=3.6
source activate /global/common/software/m2730/hic_HQ-yx59

# install the source files
cd nersc/
./intall
./intall_frzout   # there is some problem for install frzout, which I should fix it later


# after intallation, you can deactivate the virtual env, as long as you specify the $PATH
source deactivate /global/common/software/m2730/hic_HQ-yx59

# ==== run-wrapper =========
export CONDA_PREFIX=/global/common/software/m2730/hic_HQ-yx59
export PATH=$CONDA_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export XDG_DATA_HOME=$CONDA_PREFIX/share

python3 run-events_cD.py args.conf $1

```

## 2018-09-07
1. update the hic_HQ/add the condor folder and test chameleon
```
# the yx59chameleonkey.pem is expired? I have to re-add it
ssh-add yx59chameleonkey.pem
```

2. try the container in NERSC
```
# pull image down on Cori from dockerhub
shifterimg -v pull docker:yingruxu/hic_hq:freestreaming

# 1. locally
# cd $CSCRATCH/PbPb5020_2018-0907/results_test1_localShifter
bash run-wrapper test1

# 2. taskfarmer submit (1 sample)
# cd $CSCRATCH/PbPb5020_2018-0907/results_test2_taskFarmer
module load taskfarmer
sbatch taskfarmer.sl
```

## 2018-09-10
1. R-install lhs packages 
```
sudo apt-get install r-base
sudo apt-get install r-cran-lhs
```

2. generate design-points for new calibration
`python3 design.py --inputfile-dir design_run97`

3. finally, submit jobs through OSG (let me do that first, and keep testing NERSC) -- I really should do that earlier (just for productivity)
also tested on NERSC (now let's see if that works??)

4. machine learning aspect (distinguish LGV and LBT??)


## 2018-09-17
- Calibrate to new parametrization of HQ evolution model:
    - `emcee` need to use the v2.2.x version, in order to be consistent with my current code
    - the posterior for ($\alpha_s, \alpha, \beta, \gamma, \xi$) is very similar to what I get for previous calculations


## 2018-09-19
1. How much does experimental error matters in the analysis?
    - case 1: if there is no Experimental error?
    - case 2: add experimental error back one-by-one ?
    
> kind of solved: but I need to figure out what is the effect of the extra_std to the ?? Because it should not be like that?
> Is it because my emulator Z is not correctly transformed?

2. Send out the projected estimation of the diffusion coefficients to CMS group.


## 2018-09-24
1. Finally, hopefully this is the last iteration of the PHSD paper.
2. The effect of correlations on final observables??


## 2018-09-27
1. Finally, submit the paper to arxiv today.
2. <font color='red'> todoList 
 - [ ] uncertainty invest: from a very simple example, to a more realistic investigate on Langevin dynamics but calibrate on different formular, Langevin calibrate on Boltzmann. 
 - [ ] Just more observables to calculate. It might be boring, but the experimentalist like to see those: like different Meson species, event-engine observables
 - [ ] LGV vs. LBT comparison, the second point makes sense
 - [ ] pPb collisions, I think it is a good idea to make some calculation on pPb collisions now.
</font>

3. submitted 4 validation cases, to validate the calibration, also, test if I can calculate the event-engine observables.
