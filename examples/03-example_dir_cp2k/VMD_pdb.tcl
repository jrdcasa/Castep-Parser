proc newRep { sel type color rep imol } {

    mol selection $sel
    mol representation $type 
    mol addrep $imol
    mol showrep $imol $rep on
    mol modcolor $rep $imol $color

}

display projection orthographic
axes location off
color Display Background white
display depthcue off
color Labels Bonds black

set listFiles [lsort [glob ./01-PDBs_OUT/*.pdb]]
foreach ifile $listFiles {


    mol addfile $ifile  type pdb waitfor -1
}
    set imol1 [molinfo top]
set nframes [molinfo top get numframes]
mol delrep 0 $imol1
set rep 0
newRep "all" "DynamicBonds 1.900000 0.100000 12.000000" "name" $rep $imol1 
set rep 1
newRep "all" "VDW 0.200000 12.000000" "name" $rep $imol1 

animate goto start
pbc box

