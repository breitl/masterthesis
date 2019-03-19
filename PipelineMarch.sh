#!/bin/bash
#SBATCH -o ThirdRun_%j.out
#SBATCH -N 1
#SBATCH -B 2:8:1
#SBATCH -t 24:00:00
 
#source /etc/profile.d/modules.sh
module load gompi/2017b
source /home/breitlpa/shortcuts.sh

#Preprocessing: Create PHYLIP, matrix and partition file from NEXUS, sanity check with RAxML
../nexconvert Rabosky.scincids.nex
../reducer Rabosky.scincids.prematrix Rabosky.scincids.prepartition Rabosky.scincids.prephy
../raxml-ng --msa Rabosky.scincids.phy --model Rabosky.scincids.partition --prefix Rabosky.scincids.parse --threads 16 --force

#Step 1: Create linked, unlinked and scaled trees with RAxML
../raxml-ng --msa Rabosky.scincids.phy --model Rabosky.scincids.partition --prefix Rabosky.scincids.linked --brlen linked --threads 16 --force
../raxml-ng --msa Rabosky.scincids.phy --model Rabosky.scincids.partition --prefix Rabosky.scincids.unlinked --brlen unlinked --threads 2 --force
../raxml-ng --msa Rabosky.scincids.phy --model Rabosky.scincids.partition --prefix Rabosky.scincids.scaled --brlen scaled --threads 16 --force

#Step 2: Calculate trees on terrace with terraphast
../terraphast Rabosky.scincids.linked.raxml.bestTree Rabosky.scincids.matrix > Rabosky.scincids.linked.terraphastoutput
../terraphast Rabosky.scincids.unlinked.raxml.bestTree Rabosky.scincids.matrix > Rabosky.scincids.unlinked.terraphastoutput
../terraphast Rabosky.scincids.scaled.raxml.bestTree Rabosky.scincids.matrix > Rabosky.scincids.scaled.terraphastoutput

#Step 3: Crop terraphast output to only trees (delete first 3 lines)
tail -n +4 Rabosky.scincids.linked.terraphastoutput > Rabosky.scincids.linked.treesOnTerrace
tail -n +4 Rabosky.scincids.unlinked.terraphastoutput > Rabosky.scincids.unlinked.treesOnTerrace
tail -n +4 Rabosky.scincids.scaled.terraphastoutput > Rabosky.scincids.scaled.treesOnTerrace

#Step 4: Calculate likelihoodscores for all trees on terrace (for later plotting)
../raxml-ng --evaluate --msa Rabosky.scincids.phy --tree Rabosky.scincids.linked.treesOnTerrace --model Rabosky.scincids.partition --prefix Rabosky.scincids.linked.loglhs --brlen linked --threads 16 --force > Rabosky.scincids.linked.loglhscoresOutput
../raxml-ng --evaluate --msa Rabosky.scincids.phy --tree Rabosky.scincids.unlinked.treesOnTerrace --model Rabosky.scincids.partition --prefix Rabosky.scincids.unlinked.loglhs --brlen unlinked --threads 2 --force > Rabosky.scincids.unlinked.loglhscoresOutput
../raxml-ng --evaluate --msa Rabosky.scincids.phy --tree Rabosky.scincids.scaled.treesOnTerrace --model Rabosky.scincids.partition --prefix Rabosky.scincids.scaled.loglhs --brlen scaled --threads 16 --force > Rabosky.scincids.scaled.loglhscoresOutput

#Step 4b: Calculate likelihoodscores under unlinked tree
../raxml-ng --evaluate --msa Rabosky.scincids.phy --tree Rabosky.scincids.unlinked.treesOnTerrace --model Rabosky.scincids.partition --prefix Rabosky.scincids.unlinked-scaled.loglhs --brlen scaled --threads 16 --force > Rabosky.scincids.unlinked-scaled.loglhscoresOutput
../raxml-ng --evaluate --msa Rabosky.scincids.phy --tree Rabosky.scincids.unlinked.treesOnTerrace --model Rabosky.scincids.partition --prefix Rabosky.scincids.unlinked-linked.loglhs --brlen linked --threads 16 --force > Rabosky.scincids.unlinked-linked.loglhscoresOutput

#Step 5: Calculate RF-distance
cat Rabosky.scincids.linked.raxml.bestTree Rabosky.scincids.unlinked.raxml.bestTree Rabosky.scincids.scaled.raxml.bestTree> Rabosky.scincids.allBestTrees
../raxml-ng --rfdist --tree Rabosky.scincids.allBestTrees --prefix Rabosky.scincids.RF --threads 16 --force
