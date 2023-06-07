#!/bin/bash

grep 'ATOM  ' cluster.pdb >| Protein.pdb
grep 'HETATM' cluster.pdb >| Ligand.pdb

#grep JZ4 3 HTB_clean .pdb > jz4.pdb
#gmx pdb2gmx -f 3 HTB_clean .pdb -o 3 HTB_processed .gro

echo "--------- Prepared Complex.gro Files---------------"

## Prepare Ligand Structure

perl sort_mol2_bonds.pl ${lig}.mol2 ${lig}_fix.mol2
antechamber -i ${lig}_fix.mol2 -fi mol2 -o ${lig}_antechamber.mol2 -fo mol2 -at gaff2 -c bcc -rn ${lig^^} -pf y -s 2 -m 1
parmchk2 -i ${lig}_antechamber.mol2 -f mol2 -o ${lig}_antechamber.frcmod

cat <<EOF > tl.in
source leaprc.gaff
loadamberparams ${lig}_antechamber.frcmod
lig=loadmol2 ${lig}_antechamber.mol2
check lig
saveamberparm lig ${lig}.prmtop ${lig}.inpcrd
quit
EOF

tleap -f tl.in
#amb2gro_top_gro.py -p ${lig}.prmtop -c ${lig}.inpcrd -t ${lig}_GMX.top -g ${lig}_GMX.gro -b ${lig}_GMX.pdb
acpype -p ${lig}.prmtop -x ${lig}.inpcrd -d
rm sqm.in sqm.out sqm.pdb tl.in
mv ./${lig^^}.amb2gmx/${lig^^}_GMX.gro ./${lig}.gro
mv ./${lig^^}.amb2gmx/${lig^^}_GMX.top ./
mv ./${lig^^}.amb2gmx/posre_${lig^^}.itp ./posre_${lig}.itp
mv ${lig^^}_GMX.top ${lig}.top
sed -n '/moleculetype/,/system/{/system/b;p}' ${lig}.top > ${lig}.itp
sed -n '/atomtypes/,/moleculetype/{/moleculetype/b;p}' ${lig}.top > ${lig}.prm
rm -rf ${lig^^}.amb2gmx
rm ${lig}.prmtop ${lig}.inpcrd ${lig}_antechamber.mol2 ${lig}_antechamber.frcmod leap.log

## Prepare Receptor Structure

pdbfixer receptor.pdb --keep-heterogens=none --replace-nonstandard --add-residues --output protein.pdb
gmx pdb2gmx -f protein.pdb -o protein_processed.gro  -ff amber99sb-ildn -water spc -ignh -p topol.top -i posre_protein.itp
mv topol.top protein.top
sed -n '/moleculetype/,/Include/{/Include/b;p}' protein.top > protein.itp
sed -i 's/Protein_chain_A/Protein /' protein.itp
cp -r /share/home/dw_user0001/auto_md/src/topol.top ./$workdir/

## Prepare Complex Structure

python3 complex_topol.py protein_processed.gro ${lig}.gro
#cat ${lig}.top | grep -A 50 atomtypes | sed -e '/^$/,$d' >> topol.top

## Steered MD

gmx editconf -f complex.gro -o newbox.gro -center 3.280 2.181 2.4775 -box 6.560 4.362 12
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1

gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -nt 10 -v -deffnm em

gmx grompp -f npt.mdp -c em.gro -r em.gro -p topol.top -o npt.tpr
gmx mdrun -nt 10 -v -deffnm npt

gmx make_ndx -f npt.gro

    > r 1-121
    > name 19 Chain_B
    > r 145-288
    > name 20 Chain_A
    > q

gmx grompp -f md_pull.mdp -c npt.gro -p topol.top -r npt.gro -n index.ndx -t npt.cpt -o pull.tpr
gmx mdrun - deffnm pull

grep -v [@#] pullf.xvg | awk '{print $2}' > forces
grep -v [@#] pullx.xvg | awk '{print $2}' > x
paste x forces > pull_force_vs_displacement.xvg
rm x forces

## get distance
echo "0" | gmx trjconv -s pull.tpr -f pull.xtc -o conf.gro -sep

for (( i=0; i<1001; i++ ));do
gmx distance -s pull.tpr -f conf${i}.gro -select 'com of group "Chain_A" plus com of group "Chain_B"' -n index.ndx -oall dist${i}.xvg
done

for (( i=0; i<1001; i++ ));do
d=`tail -n 1 dist${i}.xvg | awk '{print $2}'`
echo "${i} ${d}" >> summary_distances.dat
rm dist${i}.xvg
done

## umbrella sampling

gmx grompp -f npt_umbrella.mdp -c conf6.gro -r conf6.gro -p topol.top -n index.ndx -o npt0.tpr
gmx mdrun - deffnm npt0

gmx grompp -f md_umbrella.mdp -c nptXXX.gro -t nptXXX.cpt -p topol.top -r nptXXX.gro -n index.ndx -o umbrellaXXX.tpr  -maxwarn 2
gmx mdrun - deffnm umbrellaXXX

gmx wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal
xmgrace -nxy histo.xvg

gmx wham -it tpr-files.dat -if pullf-files.dat -o profile_shift.xvg -hist -unit kCal -zprof0 0.5
