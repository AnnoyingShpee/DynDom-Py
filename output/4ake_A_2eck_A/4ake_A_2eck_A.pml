reinitialize
load 4ake_A_2eck_A.pdb
bg_color white
color grey
select region0, resi 0-26
select region0, region0 + resi 61-113
select region0, region0 + resi 157-209
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 114-154
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 27-59
set_color colour2 = [255,255,0  ]
color colour2, region2
select region3, resi 113-114
set_color colour3 = [0  ,255,0  ]
color colour3, region3
select region4, resi 154-170
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 59-63
set_color colour5 = [0  ,255,0  ]
color colour5, region5
set dash_gap, 0
set dash_radius, 0.2
