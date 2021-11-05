# Iterating with RifDock

This page describes the various ways you can use rifdock to achieve your
goals. 

This was originally an internal page at the Baker Lab, but we posted it here. With time this will clean up.

The example flags files can be found here:

https://www.github.com/rifdock/rifdock/tree/master/data/example_flags/

This site provides background info on RifDock and the Hierarchical
search:

[<https://rif.readthedocs.io/en/devel/legacy.html>](https://rif.readthedocs.io/en/devel/legacy.html "wikilink")


The most complete list of rifdock flags can be found by looking at the code where they are defined. 
If you scroll down far enough, most of the flags have descriptions. Don't fear the code, you can do it!

https://www.github.com/rifdock/rifdock/tree/master/apps/rosetta/rif_dock_test.hh

## Hierarchical search overview

The Hierarchical search is invoked whenever you do not pass the flag:
`-xform_pos`

The Hierarchical search was Will\'s grand plan for RifDock. The basic
idea is that one may perform a branch-and-bound search of the docking
space by zooming in on promising regions and ignoring regions that do
not score well. A critical parameter here is the **beam_size**. This is
the total number of conformations that RifDock stores and is adjusted
with **-beam_size_M**. The larger this is set, the less likely RifDock
is to lose designs. But be warned, the runtime of the Hierarchical
search is often linear with this number.

During the search stage, the worst scoring designs may be eliminated
with the following flag:

` -rif_dock:global_score_cut -10.0          # After HSearch and after HackPack, anything worse than this gets thrown out`

There are two primary ways that the Hierarchical search is used

### Hierarchical search total de novo

This is where you use the hierarchical search without **requirements**.

This is the original RifDock and the total score calculated by RifDock
is based on the sum of the best rotamer energies at each position. While
this should work in practice, several problems have been discovered:

* At low resolutions, non-physical conformations may overwhelm good, real solutions. One can imagine situations where rotamers cross each other in their position on target to their position on the scaffold. Poly-TRP designs are abundant in the low resolutions.
* The correct solution can slip through the cracks. The bins for the side-chains and the bins for the conformational sampling are not perfectly aligned. It may be the case that a great interaction at one resolution is not available at a higher resolution due to a rotamer being on the edge of a bin

#### Ways to improve total de novo

There are ways to combat these problems however:

* **Use a larger beam size** - This solution makes your runs take longer and doesn't completely solve the problem, but you'd likely get a lot better results with a beam size of 100M
* **Constraints** - If you know some key characteristics of what you want your interface to look like, the constraints may help you. See the constraints section
* **-tether_to_input_position** - This is totally cheating, but if you know where you want your scaffold to dock, you can put it there before-hand and set this to the RMSD away it's allowed to be. One could likely get very good results by first doing a run without this flag and then taking every output and feeding it back in with this flag.
    * However, if you're going to do that, you can also use **-seed_with_these_pdbs**. See the Seeding Position section

### Hierarchical search with requirements

The hierarchical search performs much better with highly constrained
searches. This section applies if you use the flags
**-requirements** or **-require_satisfaction** (and maybe also
if you use constraints).

The reason for the better performance is that the issue with
non-physical solutions outcompeting good solutions no longer is an
issue. This is because typically, these runs never exceed the beam size.
One issue does remain however:

* The correct solution can slip through the cracks. The bins for the side-chains and the bins for the conformational sampling are not perfectly aligned. It may be the case that a great interaction at one resolution is not available at a higher resolution due to a rotamer being on the edge of a bin

Brian C accidentally found a fix for this. One can run many hierarchical
searches simultaneously each offset a little from the other. In this
way, if a rotamer slips through the cracks in one search, it will likely
be found in another.

Adding these flags will fix the issue:

```
# WARNING, these flags may greatly increase runtime! Only use these if you aren't getting enough output.
-apply_seeding_xform_after_centering
-seeding_by_patchdock false
-seeding_pos 1000_random_16A_xforms.x
```
https://github.com/rifdock/rifdock/tree/master/data/1000_random_16A_xforms.x


## Seeding position overview

Seeding position searches are invoked whenever you pass the flag
**-xform_pos**.

In this mode, you tell RifDock all the places it should start a search,
and then it searches exhaustively around each of the positions. The flag
**-xform_pos** is the list of positions around the input position that
RifDock will look. A set of default **-xform_pos** exists here:

https://github.com/rifdock/rifdock/tree/master/data/example_xform_pos/

The filename indicates the maximum angular sampling in degrees. These
have a cartesian step of 0.7A and an angular step of 2.3 degrees. Brian
thinks that remaking these at finer resolution would provide better
results.

People typically use the 10 degree file.

During the search stage, the worst scoring designs may be eliminated
with the following flag:

` -rif_dock:global_score_cut -10.0          # After HSearch and after HackPack, anything worse than this gets thrown out`

There are many places to obtain seeding positions from:

### Seeding positions from patchdock

See the full tutorial for this at (find the SI for Cao 2021 Nature. You want the .tar.gz SI then look for cao_2021_protocol/)

In short, as long as one runs patchdock using the centered target pdb
from the rifgen output folder, the patchdock output files may be used as
seeding positions.

The input to rifdock would look like this:

```
-xform_pos my_xform_file.x
-seeding_pos my_patchdock_output_file.out
```

Two extra options exists to limit the patchdock outputs:

```
-rif_dock:patchdock_min_sasa  1000        # Only take patchdock outputs with more than this sasa
-rif_dock:patchdock_top_ranks 2000        # Only take the first N patchdock outputs
```

#### Seeding positions from pdbs

This is a rare scenario, but if you have a bunch of docked structures,
you can use the docks as seeding positions. A critical requirement here
is that the docked scaffold must have exactly the same backbone as the
pdb you are passing to rifdock.

When passing pdbs, ensure that they are all aligned to the target in the
rifgen output folder and that you have removed the target.

```
-rif_dock:seed_with_these_pdbs *.pdb       # List of scaffolds floating in space above the target that you
                                           #   would like to use instead of numeric seeding positions. The
                                           #   target shouldn't be present and the scaffold must match exactly
                                           #   Use this instead of -seeding_pose
-rif_dock:seed_include_input true          # Include the input pdb as one of the pdbs for -seed_with_these_pdbs
```

### Seeding positions from xforms

This is an advanced topic, but it may be the case that you know where
you want to put the scaffold and you aren\'t using patchdock.

In this situation, write your transforms to a file in 12 number format
(Rotation row1, rotation row2, rotation row3, translation) that move the
scaffold from it\'s input pdb to the location above the target in the
rifgen output folder. Then pass these flags:

```
-seeding_pos my_xforms_file.x
-seeding_by_patchdock false
```

## The standard RifDock Protocol

RifDock begins with a large scale docking search which is coordinated
either with the Hierarchical Search or the Seeding Position search. But
this is only the first step. There are then two more steps that you
scaffolds pass through: HackPack and Rosetta Score and Min.

### HackPack

This stage involves packing the rotamers found during the search stage.
All rotamers are rescored with the scaffold in its current position and
a quick monte-carlo is performed in order to determine which rotamers
fit together.

This stage is critical and should not be skipped, the following
determines what percent of search results go on to HackPack:

```
-rif_dock:hack_pack_frac  0.1              # What fraction of your HSearch results (that passed global_score_cut)
                                           # do you want to HackPack?
```

After HackPack, the following flag is used to filter results: (This is
the same as during the Hierarchical Search)

```
-rif_dock:global_score_cut -10.0          # After HSearch and after HackPack, anything worse than this gets thrown out
```

However, if the seeding positions are used, two filters are applied
immediately before **-rif_dock:global_score_cut**:

```
-rif_dock:cluster_score_cut -6.0          # After HackPack, what results should be thrown out before applying
                                         # -keep_top_clusters_frac
-rif_dock:keep_top_clusters_frac 0.5      # After applying the cluster_score_cut, what fraction of remaining 
                                           # seeding positions should survive?
```

### Rosetta Score and Min

After a set of side-chains has been decided. RifDock will load your
outputs as Rosetta poses and score them, and then score them after
minimizing. This step almost always leads to better downstream results,
however, it is rather slow. For this reason, Brian recommends disabling
it during debugging and enabling it for production runs.

There are three primary flags that affect this step:

```
-rif_dock:rosetta_score_fraction 0.006               # These two flags greaty affect runtime!!!!!
-rif_dock:rosetta_min_fraction 0.07                  # Choose wisely, higher fractions give more, better
                                                   # output at the cost of runtime
-rif_dock:rosetta_score_cut -10.0                    # After RosettaScore, anything with a score worse 
                                                   # than this gets thrown out
```

Setting fractions can be hard though. The following flags exist to help
you decide how many to score. Often you would set the previous two
fractions to 1.0 and then use these:

```
-rif_dock:rosetta_min_at_least 30                    # Make sure at least this many survive the rosetta_min_fraction
-rif_dock:rosetta_min_at_most 300                    # Make sure no more than this get minned
-rif_dock:rosetta_score_at_most  3000                # Make sure that no more than this many go to rosetta score
```

## I want to do X

This section is a guide on how to do specific things with rifdock.

### There\'s a hydrophobic pocket I want to target

* One method to do this is with the rifdock requirements. Here you would use the hydrophobic option.
* Another method is to force the scaffold to be near the pocket. This can be done with constraints.
* Finally, the Hydrophobic Filters were purpose built for this task. At the time of writing, this is probably the best solution.

### I have a specific hotspot I want to target

* If your hotspot is just a pocket or a specific hydrogen bond you want to form, the requirements may work for you.
* Alternatively, if you have your hotspot in a pdb and you want rifdock to match this, you can use the User Hotspots framework.

## Hydrophobic Filters Information

These filters are applied during HackPack and are used to guide RifDock
towards designs that make good hydrophobic contacts. These flags are
primarily designed to work with Seeding Position searches.

```
# These are rather experimental flags. You'll have to play with the values.
# Hydrophobic ddG is roughly fa_atr + fa_rep + fa_sol for hydrophobic residues.

#-hydrophobic_ddg_cut -12                  # All outputs must have hydrophobic ddG at least this value
#-require_hydrophobic_residue_contacts 5   # All outputs must make contact with at least this many target hydrophobics

#-one_hydrophobic_better_than -3           # Require that at least one rifres have a hydrophobic ddG better than this
#-two_hydrophobics_better_than -3          # Require that at least two rifres have a hydrophobic ddG better than this
#-three_hydrophobics_better_than -1        # Require that at least three rifres have a hydrophobic ddG better than this

# This next flag affects the *_hydrophobics_better_than flags. A rifres can only be counted towards 
#  those flags if it passes this one.
#-hydrophobic_ddg_per_atom_cut -0.4        # Require that hydrophobics for the *_hydrophobics_better_than 
                                           #  flags have at least this much ddG per side-chain heavy atoms.

#-hydrophobic_target_res 1,15,29,35        # If you want your selection of hydrophobic residues to 
                                         #  include only a subset of the ones you selected for the target_res, 
                                           #  place that selection here with commas.
```

There are a lot of ways you can use these to success.

* Setting **-require_hydrophobic_residue_contacts** may alone be enough for you to target your hydrophobic region, especially if you use **-hydrophobic_target_res**
* The **-*_hydrophobics_better_than** flags can also be used with great success. If you find your hotspots aren't making enough contacts, try increasing **-hydrophobic_ddg_per_atom_cut**

### Simple Heuristic to set these

If you\'re not sure what to do, the following heuristic might work.

* Set **-global_score_cut** -6, **-cluster_score_cut** -6, and **-rosetta_score_cut** -10. (Don't forget to disable Rosetta Score and Min for debugging speed)
* Increase **-hydrophobic_ddg_cut** until you stop getting outputs. Note this number and set it back to 0.
* Increase **-require_hydrophobic_residue_contacts** until you stop getting outputs. Note this number and set it back to 0.
* Increase one of the **-*_hydrophobics_better_than** flags until you stop getting outputs. Note this number and set it back to 0.
* Set the previous 3 flags to 75% of the value that caused you to stop getting outputs.

## User Hotspot Information

If you have \"hotspots\" in pdb format that you want to incorporate into
your designs, RifDock has tools for you to do this. The bad news is that
you need to re-run rifgen.

RifGen Flags:

```
###################### Hotspot configuration #################################
#   (use this either with or without apores, donres, and accres)

# Pick one of the two following methods for hotspot input:

# Hotspot input multiple distinct groups
# -hotspot_groups group0.pdb group1.pdb group2.pdb group3.pdb

# Hotspot input every hotspot is a group
# -hotspot_groups all_my_hotspots.pdb
# -single_file_hotspots_insertion

# -hotspot_sample_cart_bound 1.5   # How much do you want your hotspots to move left/right/up/down
# -hotspot_sample_angle_bound 15   # What angular deviation from your hotspot will you accept

# -hotspot_nsamples 100000  # How many times should the random sampling be done. 100000 - 1000000 is good

# -hotspot_score_thresh -0.5 # What score must a hotspot produce in order to be added to the RIF
# -hotspot_score_bonus -4    # Be careful, rifdock has a maximum score of -9
                             #  do not exceed this (this gets added to the hotspot score)
```

When you give a hotspot to RifGen, it does a few things:

* Add the exact rotamer into the RifDock rotamers (chi angles) so that that exact rotamer will appear inside the RIF
* Build inverse rotamers for each of your hotspots keeping the last 3 atoms of the hotspot fixed
* Randomly perturb the hotspot so as to completely fill in the RIF. The amount of sampling is controlled with **-hotspot_sample_cart_bound**, **-hotspot_sample_angle_bound**, and **-hotspot_nsamples**

The reason for the hotspot groups is that RifDock will allow you to
require individual hotspots to be placed. So for instance, if you have
three locations on your target where you have multiple hotspots, you
would make 3 pdb files, one for each location. If however, each of your
input residues is a unique hotspot, then you may pass them all in one
file and use **-single_file_hotspots_insertion**.

If you plan to require individual hotspots (or X hotspots) in RifDock,
you must use the following flag:

```
-rifgen::rif_type RotScoreSat_1x16  # you already have this flag, just change the value
```

Something terrible has happened during the development of RifDock
however. The system that allows you to require individual hotspots was
originally designed to be used to require hbonds. If you build hbonding
residues AND build hotspots, these systems will overlap. This means that
that in the system, item 0 is both your first hotspot and hbond 0. This
seriously needs fixed, but at the moment, you can limit the damage by
adding this flag:

```
-min_hb_quality_for_satisfaction -0.999999 # If using require_satisfaction (or buried unsats). How good does a
                                           #   hydrogen bond need to be to "count"? The scale is from -1.0 to 0 where -1.0
                                             #   is a perfect hydrogen bond.
```

When it comes time to use RifDock, there are a few flags you can use to
control the behavior of your hotspots:

```
#-require_satisfaction 4                   # Require at least this many hbonds, hotspots, or "requirements"

#-user_rotamer_bonus_constant 0            # Anything that makes a hydrogen-bond, is a hotspot, or 
                                           # is a "requirement" gets this bonus
#-user_rotamer_bonus_per_chi 0             # Anything that makes a hydrogen-bond, is a hotspot, or
                                           # is a "requirement" gets this bonus * number of chis

#-requirements 0,1,2,8                     # Require that certain satisfactions be required in all outputs
                                           # If one runs a standard RifDock, these will be individual hydrogen bonds
                                           #       to specific atoms
                                           # If one uses hotspots during rifgen, these will correspond the the 
                                           #       hotspots groups
                                           #   However, due to some conflicts, these will also overlap with hydrogen
                                           #       bonds to specific atoms
                                           # Finally, if one uses a -tuning_file, these will correspond to the
                                           #       "requirements" set there
```

## Requirement/tuning_file Information

Longxing adding a really cool system that allows one to specify specific
interactions during RifGen that you would like to later require in
RifDock. One begins by re-running RifGen with a **-tuning_file**.

Here is an example tuning file:

```
DONOR_DEFINITION
309 SER TYR
366 ASN GLN
369 ASN GLN
370 SER TYR
END_DONOR_DEFINITION

ACCEPTOR_DEFINITION
309 SER TYR
366 ASN GLN
369 ASN GLN
370 SER TYR
END_ACCEPTOR_DEFINITION

HBOND_DEFINITION
O 311 TYR TRP THR SER
HE1 335 SER THR TYR GLU GLN ASP ASN
END_HBOND_DEFINITION
BIDENTATE_DEFINITION
O 311 HG1 311
END_BIDENTATE_DEFINITION
REQUIREMENT_DEFINITION
1 APOLAR LEU PHE VAL ILE TRP MET
O 284 9.0
CG2 370 10.0
O5 31 9.0
O 34 9.0
END_APOLAR
2 HBOND HE1 335
3 HBOND O 311
4 BIDENTATE O 311 HG1 311
END_REQUIREMENT_DEFINITION
```

The tuning files doe not have the same issue with overlapping with
hydrogen bonds because all of the hydrogen bond info is cleared when the
tuning file is loaded.

The tuning file consists of multiple sections. All sections are optional.

### DONOR_DEFINITION

For a given target hbond donor residue, specify which aa types may be used to make hbonds.
This is mainly used to speedup the rifgen process. Hbonds with non-listed aa types aren't
excluded per-se, they just aren't sampled explicitly. Specify no aa-types to effectively disable
hbond sampling to that residue number. (But note that hbonds may still appear from bidentates interactions with other atoms)

Format:

```
DONOR_DEFINITION
  <RESNUM> <NAME3> <NAME3> <NAME3> ...
END_DONOR_DEFINITION
```

### ACCEPTOR_DEFINITION

For a given target hbond acceptor residue, specify which aa types may be used to make hbonds.
This is mainly used to speedup the rifgen process. Hbonds with non-listed aa types aren't
excluded per-se, they just aren't sampled explicitly. Specify no aa-types to effectively disable
hbond sampling to that residue number. (But note that hbonds may still appear from bidentates interactions with other atoms)

Format:

```
ACCEPTOR_DEFINITION
  <RESNUM> <NAME3> <NAME3> <NAME3> ...
END_ACCEPTOR_DEFINITION
```

### HBOND_DEFINITION

Define which residue types are allowed to hydrogen bond to a specific
atom. If you specify an atom here, only aa-types you list will
be stored in the rif if an hbond is detected.
    
(It's still possible to end up with accidental hbonds though at the very end of rifdock
 when rotamers are drawn from the rif and placed on the scaffold (because they move a little)).

Format:
    
```
HBOND_DEFINITION
  <ATOM_NAME> <RESNUM> <NAME3> <NAME3> <NAME3> ...
END_HBOND_DEFINITION
```

### BIDENTATE_DEFINITION

Limit target atoms to only making a very specific bidentate. All hbonding rotamers to these
atoms will be discarded unless the given bidentate is formed.

Format:

```
BIDENTATE_DEFINITION
    <ATOM_NAME1> <RESNUM1> <ATOM_NAME2> <RESNUM2>
BIDENTATE_DEFINITION
```

### HOTSPOT_DEFINITION

There's a way to do this. But only Longxing knows how...

### REQUIREMENT_DEFINITION

Here you specify the individual numbered requirements you want to use.
The number you select here will be used later in RifDock. There are
several things you can specify:

#### APOLAR

No guarentees this is correct. See if it works for you (we need Longxing to verify).

Specify that you want a hydrophobic residue of specific type to make an
interaction with a hydrophobic pocket. Interacting with the pocket is
defined as when the hydrophobic residue\'s CA is within X
Angstroms of all the atoms specified. (I mean, probably. Whow knows. Seems like a good fit for the input format though)

Format:

```
<REQ_NUM> APOLAR <NAME3> <NAME3> <NAME3>
  <ATOM_NAME1> <RESNUM1> <DISTANCE1>
  <ATOM_NAME2> <RESNUM2> <DISTANCE2>
  <ATOM_NAME3> <RESNUM3> <DISTANCE3>
END_APOLAR
```

#### HBOND

Specify that you want a certain hbond to exist. You do not have to
declare this in the HBOND definition.

Format:

```
<REQ_NUM> HBOND <ATOM> <RESNUM>
```

#### BIDENTATE

Specify that you want a certain bidentate hbond to exist.

Format:

```
<REQ_NUM> BIDENTATE <ATOM_NAME1> <RESNUM1> <ATOM_NAME2> <RESNUM2>
 ```

#### HOTSPOT
  
There's a way to do this. But only Longxing knows how...

Specify that you want a certain hotspot to exist.

```
<REQ_NUM> HOTSPOT ???????????????????
```

### Using the Requirements

After you\'ve built your rifgen, the same options that control hotspots
may be used to control the requirements. Longxing originally intended
for **-requirements** to be used, but a few alternati Brian wants to
add something to allow a per-requirement bonus.

```
#-requirements 0,1,2,8                     # Require that that all of these requirements be present in the output
#-sat_score_bonus 0:-2,1:-1.5              # Give a bonus score for each of the requirments.
#-sat_score_override 0:-2,1:-1.5           # Totally override the original score of the requirement rif residue
#-requirement_groups ......                # A way to group requirements and require a few of them. Space separated.
                                              # I want at least 3 of these requirements: 3:1,5,8,10,23. 
                                              # I want less than 2 of these: -2:4,6,7.
                                              # Negative requirements means not this requirement.
                                              # Make sure these 4 things don't happen simultanesouly: 4:-4,-6,-7,-8
                                                   # (You could also do this) -4:4,6,7,8
                                           
#-require_satisfaction 4                   # Require at least this many hbonds, hotspots, or "requirements" (repurposed flag)
```
  
This last one's complicated. You can specify that certain scaffold positions must be the ones to satisfy a requirement.

The first step is to label each position with a PDBInfoLabel. If you're not using Rosetta, inject this into your scaffold 
PDB for each position you want to label. Feel free to repeat labels depending on what you're trying to accomplish.

```
REMARK PDBinfo-LABEL:    1 LABEL_OF_YOUR_CHOICE ANOTHERLABEL THIRDLABEL
REMARK PDBinfo-LABEL:    43 ANOTHERLABEL THIRDLABEL
```
  
Then, these two flags will work their magic.
```
-pdbinfo_requirements ANOTHERLABEL:6    # Pairs of pdbinfo_label:req1,req2,req3 that specify that a residue with this pdbinfo_label must satisfy these requirements/sats.
                                        # Space separated for multiple
-num_pdbinfo_requirements_required 3    # Minimum number of pdbinfo_requirements to satisfy. -1 for all.
```

## Constraint Information

The constraints are a useful system also added by Longxing before the
tuning file. These are mostly atom-pair constraints that throw out docks
that fail the filter.

Put these into a file and pass it to the flag: -cst_files

For all cases, positive numbers throw out things that are too far away
and negative numbers throw out things that are too close together.
Descriptions below are given for the positive number.

### AtomPair

Force a scaffold atom to be near a target atom.

Format:

```
AtomPair <TARGET_ATOM_NAME> <TARGET_RESNUM> <SCAFF_ATOM> <SCAFF_RESNUM> <DISTANCE> 
```

### AtomToScaffold

Force an atom of the target to be within X Angstroms of any of the
scaffold backbone atoms.

Format:

```
AtomToScaffold <TARGET_ATOM> <TARGET_RESNUM> <DISTANCE>
```

### AtomToTarget

Force an atom of the scaffold to be within X Angstroms of any of the
target backbone atoms.

Format:

```
AtomToTarget <SCAFF_ATOM> <SCAFF_RESNUM> <DISTANCE>
```

### RayToScaffold

Draw a ray from the target out to a certain length. Ensure that at least
one scaffold backbone atoms is within X Anstroms of this ray. The ray is
defined by two target atoms and the length is defined starting from the
FROM_ATOM.

Format:

```
RayToScaffold <TO_ATOM> <TO_RESNUM> <FROM_ATOM> <FROM_RESNUM> <LENGTH_OF_RAY> <DISTANCE>
```

### RayToTarget

Draw a ray from the scaffold out to a certain length. Ensure that at
least one target backbone atoms is within X Anstroms of this ray. The
ray is defined by two scaffold atoms and the length is defined starting
from the FROM_ATOM.

Format:

```
RayToTarget <TO_ATOM> <TO_RESNUM> <FROM_ATOM> <FROM_RESNUM> <LENGTH_OF_RAY> <DISTANCE>
```
  
  
