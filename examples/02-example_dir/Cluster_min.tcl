proc newRep { sel type color rep imol size} {

    mol selection $sel
    mol representation $type $size
    mol addrep $imol
    mol showrep $imol $rep on
    mol modcolor $rep $imol $color

}

display projection orthographic
axes location off
color Display Background white
display depthcue off
color Labels Bonds black

mol new  Cluster_min.xyz type xyz waitfor -1
set imol1 [molinfo top]
set nframes [molinfo top get numframes]
mol delrep 0 $imol1
set rep 0
newRep "all" "cpk" "name" $rep $imol1 0.7
set rep 1
newRep "all" "cpk" "name" $rep $imol1 0.7
mol drawframes $imol1 $rep {0:1:1000}
mol showrep $imol1 $rep 0
mol modcolor $rep $imol1 Timestep

animate goto start

