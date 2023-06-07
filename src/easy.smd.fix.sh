#!/bin/bash

# file: easy_smd.slurm

#SBATCH --job-name=gromaces_smd
#SBATCH --partition=gpu
#SBATCH --nodelist=c003
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --gres=gpu:0
#SBATCH --mem=400G
#SBATCH --output=gromacs_smd_%j.log
#SBATCH --error=gromacs_smd_%j.err
#SBATCH --distribution=cyclic:block

export OMP_NUM_THREADS=36

echo TIME subset start $(date)
source /share/apps/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate md_env
## gromacs ver 2022.04

#------------------
box="2.0"
ion1="NA"
ion2="CL"
gpu="0"
lig="atp"
workdir="./"
GMX="/usr/local/gromacs/bin/gmx"

# SET NUMBER OF MPI TASKS 
#mpi="1"
mpi=$(($SLURM_NTASKS_PER_NODE*$SLURM_JOB_NUM_NODES))
# SET NUMBER OF OPENMP THREADS
#omp="36"
omp=$SLURM_CPUS_PER_TASK

sed -n '/MODEL        1/,/TER/p' clusters.pdb | head -n -1 | tail -n+2 > ref-complex.pdb
grep ' A ' ref-complex.pdb >| protein.pdb
echo TER >> protein.pdb
grep ' B ' ref-complex.pdb >| ${lig}.pdb
echo TER >> ${lig}.pdb

## Prepare Ligand Structure
antechamber -i ${lig}.pdb -fi pdb -o ${lig}_antechamber.mol2 -fo mol2 -c bcc -rn ${lig^^} -s 2 -at gaff2 -pf y -m 1
parmchk2 -i ${lig}_antechamber.mol2 -f mol2 -o ${lig}_antechamber.frcmod

cat <<EOF > tl.in
source oldff/leaprc.ff99SB
source leaprc.gaff2
loadamberparams ${lig}_antechamber.frcmod
lig=loadmol2 ${lig}_antechamber.mol2
check lig
saveamberparm lig ${lig}.prmtop ${lig}.inpcrd
quit
EOF

tleap -f tl.in

acpype -p ${lig}.prmtop -x ${lig}.inpcrd -d
cp ./${lig^^}.amb2gmx/${lig^^}_GMX.gro ./${lig}.gro
cp ./${lig^^}.amb2gmx/${lig^^}_GMX.top ./${lig}.top
#mv ./${lig^^}.amb2gmx/posre_${lig^^}.itp ./posre_${lig}.itp
sed -n '/atomtypes/,/moleculetype/{/moleculetype/b;p}' ${lig}.top > ${lig}.prm
sed -n '/moleculetype/,/system/{/system/b;p}' ${lig}.top > ${lig}.itp
rm -rf ${lig^^}.amb2gmx
rm -rf QOUT punch tl.in
rm -rf ${lig}.out ${lig}.gesp ${lig}.chk ${lig}.gjf
rm ${lig}.prmtop ${lig}.inpcrd ${lig}_antechamber.mol2 ${lig}_antechamber.frcmod leap.log

## Prepare Protein Structure
$GMX pdb2gmx -f protein.pdb -o protein_processed.gro -ff amber99sb -water spc -ignh -nochargegrp -p topol.top -i posre_protein.itp
mv topol.top protein.top
sed -n '/moleculetype/,/restraint/{/restraint/b;p}' protein.top > protein.itp
sed -i 's/Protein_chain_A/Protein /' protein.itp

echo "--------- Prepared Complex.gro Files---------------"

## Prepare Complex Structure [not position restraint for ligand]

(echo "0 & ! a H*"; echo "name 3 ${lig^^}-H"; echo "q") | $GMX make_ndx -f ${lig}.gro -o index_${lig}.ndx
echo "3" | $GMX genrestr -f ${lig}.gro -n index_${lig}.ndx -o posre_${lig}.itp -fc 1000 1000 1000

cp -r /share/home/dw_user0001/auto_smd/src/topol.top ./$workdir/
cp -r /share/home/dw_user0001/auto_smd/src/complex_topol.py ./$workdir/
python3 complex_topol.py protein_processed.gro ${lig}.gro
sed -i '/; Include water topology/i\; Include Ligand position restraints\n#ifdef POSRES_LIGAND\n#include "'"./posre_${lig}.itp"'"\n#endif\n' topol.top

## Steered MD

cp /share/home/dw_user0001/auto_smd/src/mdp/* ./$workdir/
sed /JZ4/s//${lig^^}/ md_pull.mdp -i
sed /JZ4/s//${lig^^}/ npt_umbrella.mdp -i
sed /JZ4/s//${lig^^}/ md_umbrella.mdp -i
sed /JZ4/s//${lig^^}/ npt.mdp -i 

$GMX editconf -f complex.gro -o newbox.gro -bt cubic -d 12
$GMX solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
$GMX grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo 15 | $GMX genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

(echo "1 | 13"; echo "q") | $GMX make_ndx -f solv_ions.gro -o index.ndx

$GMX grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
$GMX mdrun -v -deffnm em -gpu_id 0 

$GMX grompp -f npt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 5
$GMX mdrun -v -deffnm npt -gpu_id 0

$GMX grompp -f md_pull.mdp -c npt.gro -p topol.top -r npt.gro -n index.ndx -t npt.cpt -o pull.tpr -maxwarn 5
$GMX mdrun -deffnm pull -gpu_id 1 -pf pullf.xvg -px pullx.xvg 

grep -v [@#] pullf.xvg | awk '{print $2}' > forces
grep -v [@#] pullx.xvg | awk '{print $2}' > x
paste x forces > pull_force_vs_displacement.xvg
rm x forces

## get distance

mkdir -p pull_out
echo "0" | $GMX trjconv -s pull.tpr -f pull.xtc -o pull_out/conf.gro -sep
cd pull_out
count=$(ls | wc -l)

for (( i=0; i<$count; i++ ));do
$GMX distance -s ../pull.tpr -f conf${i}.gro -select 'com of group "Protein" plus com of group "'"${lig^^}"'"' -n ../index.ndx -oall dist${i}.xvg
d=`tail -n 1 dist${i}.xvg | awk '{print $2}'`
echo "${i} ${d}" >> distances.dat
rm dist${ps}.xvg
done

#Selecting frames for umbrella sampling input

space=0.2
range=`cat distances.dat| awk '{if(min==""){min=max=$2}; if($2>max) {max=$2}; if($2<min) {min=$2}; total+=$2; count+=1} END {print min" "sp" "max}' sp="$space"`
for xx in `seq -f "%f" ${range}`;do 
aa=`awk -v c=2 -v t=${xx} '{a[NR]=$c}END{
        asort(a);d=a[NR]-t;d=d<0?-d:d;v = a[NR]
        for(i=NR-1;i>=1;i--){
                m=a[i]-t;m=m<0?-m:m
                if(m<d){
                    d=m;v=a[i]
                }
        }
        print v
}' distances.dat`
grep $aa distances.dat |head -n1 >> list
done

#sed -i '1d' list

#renaming input frames and deleting unnecessary frames
aa=1
awk '!a[$0]++' list|while read i;do num=`echo $i|awk '{print $1}'`;cp conf${num}.gro umb${aa}.gro;aa=`expr $aa + 1`;done

## umbrella sampling
cp /share/home/dw_user0001/auto_smd/src/mdp/*_umbrella.mdp ./
cp /share/home/dw_user0001/auto_smd/index.ndx ./
cp /share/home/dw_user0001/auto_smd/atp.prm ./
cp /share/home/dw_user0001/auto_smd/posre_protein.itp ./
cp /share/home/dw_user0001/auto_smd/posre_atp.itp ./
cp /share/home/dw_user0001/auto_smd/atp.itp ./
cp /share/home/dw_user0001/auto_smd/protein.itp ./
cp /share/home/dw_user0001/auto_smd/topol.top ./
sed /JZ4/s//${lig^^}/ md_umbrella.mdp -i
sed /JZ4/s//${lig^^}/ npt_umbrella.mdp -i

for ii in $(seq 1 `awk '!a[$0]++' list |wc -l`);do
$GMX grompp -f npt_umbrella.mdp -c umb${ii}.gro -r umb${ii}.gro -p topol.top -n index.ndx -o npt${ii}.tpr -maxwarn 5
$GMX mdrun -deffnm npt${ii} -gpu_id 1
$GMX grompp -f md_umbrella.mdp -c npt${ii}.gro -t npt${ii}.cpt -p topol.top -r npt${ii}.gro -n index.ndx -o umbrella${ii}.tpr -maxwarn 5
$GMX mdrun -deffnm umbrella${ii} -gpu_id 1 
echo "umbrella${ii}.tpr" >> tpr-files.dat
echo "umbrella${ii}_pullf.xvg" >> pullf-files.dat
done

wait
$GMX wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal
#$GMX wham -it tpr-files.dat -if pullf-files.dat -o profile_shift.xvg -hist -unit kCal -zprof0 0.5
xmgrace -nxy histo.xvg

